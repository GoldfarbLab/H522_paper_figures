library(here)
library(limma)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(ComplexHeatmap)
library(ConsensusClusterPlus)
library(org.Hs.eg.db)
library(annotate)
library(biogridr) # devtools::install_github("npjc/biogridr")
library(CellSurfaceTools) # install.packages('devtools') # devtools::install_github("GoldfarbLab/CellSurfaceTools")

source(here("common.R"))
select <- get(x = "select", pos = "package:dplyr") # deal with function masking

log.stats.threshold = -log10(0.05)
log.fc.threshold = log2(1.3)

################################################################################
# Read data
################################################################################
data <- read_tsv(here("data/MS/Phospho (STY)Sites_05.txt"), guess_max=20000)
design <- read_csv(here("data/MS/Experimental Design H522 Paper.csv"))
correction_factors <- read_tsv("data_processed/protein_loading_correction_factors.tsv")

#data <- data %>% select(-grep("Rep 1", colnames(data), value = T)) 
#correction_factors <- correction_factors[10:27,]
#design <- design %>% filter(Replicate != 1)

data <- data %>% select(-grep("Reporter intensity corrected 10.*Rep", colnames(data), value = T)) 
# get quant columns
reporter.intensities <- grep("Reporter intensity corrected.*Rep", colnames(data), value = T)

#FILTERING DATA FOR NAs, CONTAMINANTS, NONZERO ROWS
data$`Gene names`[is.na(data$`Gene names`)] <- data$`Proteins`[is.na(data$`Gene names`)] #replace NA gene names with Protein ID column (Majority protein IDs)
filtered.data <- filter(data, is.na(data$"Reverse"), is.na(data$"Potential contaminant")) %>%
  filter(`Localization prob` >= 0.75) %>%
  select("Leading proteins", "Gene names", "Positions within proteins", "id", all_of(reporter.intensities))

# Reshape so that the separate cases for 1 phosphosite, 2 phosphosites, and 3+ phosphosites are in separate rows
phosphosites.long <- as_tibble(reshape(as.data.frame(filtered.data), 
                                       timevar="multiplicity",
                                       idvar="id", 
                                       varying=reporter.intensities,
                                       sep="___",
                                       direction='long'))

reporter.intensities <- grep("Reporter intensity corrected.*Rep", colnames(phosphosites.long), value = T)
phosphosites.long <- phosphosites.long %>% group_by(`Leading proteins`, `Positions within proteins`, `Gene names`) %>%
  summarise_at(reporter.intensities, sum)

# filter for valid values
rep1.valid <- phosphosites.long$`Reporter intensity corrected 9 Rep 1` > 0 # bridge
rep2.valid <- phosphosites.long$`Reporter intensity corrected 9 Rep 2` > 0
rep3.valid <- phosphosites.long$`Reporter intensity corrected 9 Rep 3` > 0

phosphosites.long <- as_tibble(phosphosites.long) %>% filter( (rep1.valid + rep2.valid + rep3.valid) >= 1)

intensities.rep1 <- grep("Reporter intensity corrected.*Rep 1", colnames(phosphosites.long), value = T) #filter data into Reps 
intensities.rep2 <- grep("Reporter intensity corrected.*Rep 2", colnames(phosphosites.long), value = T)
intensities.rep3 <- grep("Reporter intensity corrected.*Rep 3", colnames(phosphosites.long), value = T)

rep1.valid <- rowSums(select(phosphosites.long, all_of(intensities.rep1)) > 0) >= 3 #finding at least 1 nonzero value in each Rep 
rep2.valid <- rowSums(select(phosphosites.long, all_of(intensities.rep2)) > 0) >= 3
rep3.valid <- rowSums(select(phosphosites.long, all_of(intensities.rep3)) > 0) >= 3

valid.data <- filter(phosphosites.long, (rep1.valid + rep2.valid + rep3.valid) >= 1)

reporter.intensities <- grep("Reporter intensity corrected.*Rep", colnames(valid.data), value = T)

#Replace 0s with min of dataset 
gnames <- valid.data$"Gene names"
normalization.data <- select(valid.data, all_of(reporter.intensities))
normalization.data[normalization.data == 0] <- min(normalization.data[normalization.data > 0])
normalization.data <- sweep(normalization.data, 2, t(as.matrix(correction_factors)), "/")

#REP1
normalization.data.rep1 <- select(normalization.data, all_of(intensities.rep1))
normalization.data.rep1 <- sweep(normalization.data.rep1, 1, normalization.data.rep1[, 9], "/") #dividing all rows by reference channel (TMT9)

#REP2
normalization.data.rep2 <- select(normalization.data, all_of(intensities.rep2))
normalization.data.rep2 <- sweep(normalization.data.rep2, 1, normalization.data.rep2[, 9], "/")

#REP3
normalization.data.rep3 <- select(normalization.data, all_of(intensities.rep3))
normalization.data.rep3 <- sweep(normalization.data.rep3, 1, normalization.data.rep3[, 9], "/")

normalized.data <- cbind(normalization.data.rep1, normalization.data.rep2, normalization.data.rep3)

#NaN and Inf Values? 
normalized.data <- log2(normalized.data)

#remove Ref Channel (TMT9) after normalization (after normalization: all 0s) 
ref.channel <- grep("Reporter intensity corrected 9", colnames(normalized.data), value = T)
normalized.data <- select(normalized.data, -all_of(ref.channel))

# batch correction
batch <- c(rep(1,8), rep(2,8), rep(3,8))

for (i in 1:nrow(normalized.data)) {
  normalized.data[i,] <- removeBatchEffect(normalized.data[i,], batch)
}

design <- filter(design, Condition !="Bridge") #have to remove Ref Channel in design file too 

################################################################################
# Annotations
################################################################################
# replace column names
quant.colnames <- paste(design$Condition, "Rep", design$Replicate)
colnames(normalized.data) <- quant.colnames

# add back names
phosphosites <- cbind(select(valid.data, "Leading proteins", "Gene names", "Positions within proteins"), normalized.data)

write_tsv(valid.data, '~/Downloads/phosphosites_valid.txt')

################################################################################
# Stats
################################################################################
fit <- lmFit(select(phosphosites, 
                    grep("Mock - 96 hr", quant.colnames, value = T), 
                    grep("96 hr", quant.colnames, value = T)), method="robust", 
             model.matrix(~ c(rep("Mock - 96 hr", 3), rep("96 hr", 3))))

fit <- eBayes(fit)
phosphosites.stats <- topTable(fit, n=nrow(phosphosites))
phosphosites <- cbind(phosphosites.stats, phosphosites[as.numeric(rownames(phosphosites.stats)), ])

diff.phosphosites <- phosphosites %>% filter(abs(`logFC`) >= log.fc.threshold & adj.P.Val <= 0.05)

################################################################################
# Clustering
################################################################################

phosphosites.z.scored <- t(scale(t(as.matrix(select(diff.phosphosites, all_of(quant.colnames)))), center=T, scale=T))
#rownames(phosphosites.z.scored) <- phosphosites$`Gene names`

clustering.results = ConsensusClusterPlus(t(phosphosites.z.scored),
                                          maxK=3,
                                          reps=10,
                                          pItem=0.8,
                                          pFeature=1,
                                          #clusterAlg="km",
                                          #distance="euclidean",
                                          clusterAlg="hc",
                                          innerLinkage="complete",
                                          finalLinkage="ward.D2",
                                          distance="pearson",
                                          plot="pdf", 
                                          title="figures/phos_complete_ward2_pearson_1000", 
                                          seed=1262118388.71279)

num.clusters = 3
row.clusters <- clustering.results[[num.clusters]][["consensusClass"]]


diff.phosphosites$cluster <- row.clusters



################################################################################
# Generate figures
################################################################################
figPCA <- plotPCA(normalized.data, design)
print(figPCA)

exit()


figHeatmap <- plotHeatmap(phosphosites.z.scored, phosphosites, design)

colnames(phosphosites.z.scored) <- paste(design$Condition, design$Replicate)
colnames(phosphosites.z.scored)[which(str_detect(colnames(phosphosites.z.scored), "1$"))] <- "1"
colnames(phosphosites.z.scored)[which(str_detect(colnames(phosphosites.z.scored), "3$"))] <- "3"
# reorder columns
phosphosites.z.scored <- phosphosites.z.scored[, c(1,9,17, 2,10,18, 3,11,19, 4,12,20, 5,13,21, 6,14,22, 7,15,23, 8,16,24)]

heatmap.phosphosites <- Heatmap(phosphosites.z.scored,
                            # labels
                            column_title = paste("Differentially Expressed Proteins (n = ", nrow(phosphosites.z.scored), ")", sep=""),
                            column_title_gp = gpar(fontsize=7),
                            show_row_names = F,
                            row_names_gp = gpar(fontsize = 5),
                            column_names_gp = gpar(fontsize = 6),
                            name = "z-score",
                            # colors #-max(abs(quant)), max(abs(quant)
                            col = colorRamp2(seq(-3, 3, length=256), rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))),
                            # legends
                            show_heatmap_legend = F,
                            heatmap_legend_param = list(color_bar = "continuous",
                                                        title_gp = gpar(fontsize = 6),
                                                        labels_gp = gpar(fontsize = 6),
                                                        grid_width = unit(2,"mm")),
                            # clustering
                            cluster_columns = F,
                            clustering_distance_rows = "pearson",
                            cluster_row_slices = F,
                            show_row_dend = F,
                            #cluster_rows = row.clusters,
                            split = row.clusters,
                            # labels
                            row_title_rot = 0,
                            row_title_gp = gpar(fontsize=7),
                            #row_title = NULL,
)


draw(heatmap.phosphosites)
