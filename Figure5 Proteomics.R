library(here)
library(limma)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(ComplexHeatmap)
library(ConsensusClusterPlus)
library(org.Hs.eg.db)
library(annotate)
library(CellSurfaceTools) # install.packages('devtools') # devtools::install_github("GoldfarbLab/CellSurfaceTools")

source(here("common.R"))
select <- get(x = "select", pos = "package:dplyr") # deal with function masking

#source(here("RequantifyProteomics.R")) # only necessary to run once

log.stats.threshold = -log10(0.05)
log.fc.threshold = log2(1.3)

################################################################################
# PCA plot
################################################################################
plotPCA <- function(data, design)
{
  
  data.t <- as.data.frame(t(data))
  pca <- prcomp(data.t, center=T, scale=F)
  eigen <- pca$sdev^2
  variance <- (eigen/sum(eigen)) * 100
  variance <- format(round(variance, 2), nsmall=2) # show 2 digits after decimal
  
  plotting.data <- cbind(data.frame(pca$x), design)
  y.range <- max(plotting.data$PC2) - min(plotting.data$PC2)
  min.x <- min(plotting.data$PC1)
  max.x <- max(plotting.data$PC1)
  x.range <- max.x - min.x
  
  plotting.data$Label <- plotting.data$Condition
  plotting.data$Label[which(design$Replicate != 1)] <- ""
  
  p <- (ggplot(plotting.data, aes(x=PC1, y=PC2, shape=as.factor(Infected), color=Time, label=Label)) 
        + geom_point(size=2, alpha=0.8)
        + geom_text(size=2, nudge_y = y.range/20)
        
        + scale_x_continuous(limits=c(min.x - x.range/10, max.x + x.range/10))
        + scale_shape_discrete(name = "Condition")
        + scale_color_discrete(guide='none')
        #+ scale_color_brewer(guide='none', type="qual", palette="Set1")
        
        + xlab(paste("PC1 (", variance[1], "% explained variance)", sep=""))
        + ylab(paste("PC2 (", variance[2], "% explained variance)", sep=""))
        + ggtitle(label = "Principal Component Analysis (PCA) - Proteome",
                  subtitle = paste("n =", formatC(nrow(normalized.data), big.mark=","), "proteins"))
        
        + theme.basic
        + theme(legend.position = c(0.9, 0.2), 
                plot.subtitle = element_text(size=7, hjust = 0.5, face = "italic"))
  )
}

################################################################################
# Heatmap
################################################################################
plotHeatmap <- function(proteins.z.scored, diff.proteins, design)
{
  colnames(proteins.z.scored) <- paste(design$Condition, design$Replicate)
  colnames(proteins.z.scored)[which(str_detect(colnames(proteins.z.scored), "1$"))] <- "1"
  colnames(proteins.z.scored)[which(str_detect(colnames(proteins.z.scored), "3$"))] <- "3"
  # reorder columns
  proteins.z.scored <- proteins.z.scored[, c(1,9,17, 2,10,18, 3,11,19, 4,12,20, 5,13,21, 6,14,22, 7,15,23, 8,16,24)]
  
  label.sars2.pos <-  which(diff.proteins$SARS_CoV_2 == T)
  label.int.pos <-  which(!is.na(diff.proteins$Bait))
  label.pos <- c(label.sars2.pos, label.int.pos)
  
  labels.sars2 <- diff.proteins$`Gene names`[label.sars2.pos]
  labels.int <- paste(diff.proteins$`Gene names`[label.int.pos], " (", diff.proteins$Bait[label.int.pos], ")", sep="")
  labels <- c(labels.sars2, labels.int)
  
  heatmap.proteins <- Heatmap(proteins.z.scored,
                           # labels
                           column_title = paste("Differentially Expressed Proteins\nn = ", nrow(proteins.z.scored), sep=""),
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
                           #split=7,
                           split = row.clusters,
                           # labels
                           row_title_rot = 0,
                           row_title_gp = gpar(fontsize=7),
                           row_title = NULL,
                           # size
                           width=ncol(proteins.z.scored)*0.2
  )
  
  cluster.rowAnnot <- rowAnnotation(block = anno_block(gp = gpar(fill = brewer.pal(num.clusters, "Set3"), col = NA)),
                                    width = unit(2, "mm"))
  
  gene.rowAnnot = rowAnnotation(gene.name = anno_mark(at = label.pos, 
                                                      labels = labels,
                                                      labels_gp = gpar(fontsize = 5), 
                                                      link_gp = gpar(lwd = 0.5),
                                                      link_width = unit(4, "mm"),
                                                      padding = unit(0.5, "mm")))
  
  ht_list = heatmap.proteins + cluster.rowAnnot + gene.rowAnnot
  
  gb_heatmap = grid.grabExpr(draw(ht_list), height=5, width=2.5)
}

################################################################################
# Profile plots of clusters
################################################################################
plotClusterProfiles <- function(diff.proteins.averaged.condition.fc.mock)
{
  num.in.cluster <- (diff.proteins.averaged.condition.fc.mock %>%
                       group_by(cluster) %>% 
                       summarise(n=n() / (nlevels(Condition)-1)))
  
  cluster.labels <- paste("Cluster ", num.in.cluster$cluster, " (n=", num.in.cluster$n, ")", sep="")
  names(cluster.labels) <- seq(1:num.clusters)
  
  representative.profiles <- (diff.proteins.averaged.condition.fc.mock %>% 
                                group_by(Condition, cluster) %>% 
                                summarise_at(c("FC"),  mean))
  
  p <- (ggplot(diff.proteins.averaged.condition.fc.mock, aes(x=Condition, 
                                         y=FC, 
                                         group=interaction(`Gene names`, cluster),
                                         color=as.factor(cluster),
                                         facet=cluster,
                                         linetype=as.factor(SARS_CoV_2)))
        
        + geom_line(size=0.5, alpha=0.75)
        + geom_line(size=0.5, data=representative.profiles, aes(x=Condition, y=FC, group=as.factor(cluster)), linetype="solid", color=grey)
        
        + scale_y_continuous(name = "log2(fold-change) over 4h mock", breaks=seq(-5,5,1))
        + scale_x_discrete(name = "Hours post-infection")
        
        + facet_grid(cluster ~ ., scales="free", 
                     labeller = labeller(cluster = cluster.labels))
        
        + scale_color_brewer(type="qual",palette="Set3")
        
        + theme.basic
        + theme(legend.position = "none",
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  )
}

################################################################################
# Profile plots of SARS-CoV-2 proteins
################################################################################
plotCoVProfiles <- function(proteins)
{
  SARS2.prots <- (proteins %>%
                    filter(SARS_CoV_2 == T & !str_detect(Condition, "Mock")) %>%
                    mutate(Condition = as.numeric(str_replace(Condition, " hr", ""))) %>%
                    mutate(Label = `Gene names`))
  SARS2.prots$Label[which(SARS2.prots$Condition != "96")] <- ""
  
  p <- (ggplot(SARS2.prots, aes(x=Condition, y=FC, group=`Gene names`, label=Label))
        
        + geom_line(size=0.5, alpha=0.75, color=medium.grey)
        + geom_point(size=1, color=medium.grey)
        + geom_text(size=2)
        
        + ggtitle("SARS-CoV-2 Proteins")
        
        + scale_y_continuous(name = "log2(Intensity / 4h Mock)", breaks = seq(0,5,1))
        + scale_x_continuous(name = "Hours post-infection", 
                             breaks = c(4, 12, 24, 48, 72, 96),
                             expand = c(0.05, 0))
        
        + theme.basic
        + theme(legend.position = "none")
  )
}


################################################################################
# Volcano plot 96 hrs
################################################################################
plotVolcano <- function(proteins)
{
  proteins$logFC <- -proteins$logFC
  min.x <- min(proteins$logFC)
  max.x <- max(proteins$logFC)
  x.range <- max.x - min.x
  x.limit <- max(abs(min.x), abs(max.x))
  
  proteins$log.adj.P.Val <- -log10(proteins$adj.P.Val)
  y.range <- max(proteins$log.adj.P.Val) - min(proteins$log.adj.P.Val)
  
  proteins$significance <- "none"
  
  levels(proteins$significance) <- c("none", "sig fold-change", "sig statistic", "sig both", "SARS-CoV-2")
  
  proteins$significance[abs(proteins$logFC) >= log.fc.threshold] = "sig fold-change"
  proteins$significance[proteins$log.adj.P.Val >= -log10(0.05)] = "sig statistic"
  proteins$significance[abs(proteins$logFC) >= log.fc.threshold & proteins$log.adj.P.Val >= -log10(0.05)] = "sig both"
  proteins$significance[proteins$SARS_CoV_2] = "SARS-CoV-2"
  
  proteins$Label <- proteins$`Gene names`
  proteins$Label <- str_replace(proteins$`Gene names`, "SARS_CoV_2", "")
  proteins$Label[!proteins$SARS_CoV_2] <- ""
  
  p <- (ggplot(data=proteins, aes(x=logFC, y=log.adj.P.Val, label=Label, color=significance))
        + geom_point(size=0.5) 
        + geom_text(nudge_y = 0.05 * y.range, nudge_x = 0.05 * x.range, size=2, color=COV2.color)
        
        + geom_vline(xintercept = -log.fc.threshold, linetype="dashed", color=grey)
        + geom_vline(xintercept = log.fc.threshold, linetype="dashed", color=grey)
        + geom_hline(yintercept = -log10(0.05), linetype="dashed", color=grey)
        
        + scale_color_manual(values=c("none"="#AAAAAA55", "sig fold-change"="#AAAAAA55", "sig statistic"="#AAAAAA55", "sig both"="#1f78b455", "SARS-CoV-2"="#fb8072"))
        
        + scale_x_continuous(limits=c(min.x - x.range/20, x.limit + x.range/20),
                             breaks = seq(-3,6,1))
        
        + ylab("-log10(q-value)")
        + xlab("log2(fold-change)")
        
        + ggtitle("96 hpi vs 96 hr mock")
        
        + theme.basic
        + theme(legend.position = "none")
  )
}


################################################################################
# Read data
################################################################################
data <- read_tsv(here("data_processed/requantifiedProteins.txt"), guess_max=10000)
design <- read_csv(here("data/MS/Experimental Design H522 Paper.csv"))
SARS.interactors <- select(read_csv(here("annotations/SARS2_interactome.csv")), c("Bait", "PreyGeneName"))
H522.mutations <- read_tsv(here("annotations/H522_mutations.tsv"))
#TMT9: Reference Channel 
#TMT10: remove

#FILTERING DATA FOR NAs, CONTAMINANTS, NONZERO ROWS
data$`Gene names`[is.na(data$`Gene names`)] <- data$`Majority protein IDs`[is.na(data$`Gene names`)] #replace NA gene names with Protein ID column (Majority protein IDs)
filtered.data <- filter(data, is.na(data$"Only identified by site"), is.na(data$"Reverse"), is.na(data$"Potential contaminant"), data$"Razor + unique peptides" > 1) #filter data for contaminants, cols ES (Only ID'd by site) ET (Reverse)  EU (Potential Contaminant) and K (at least 1 unique peptide)

#condense data to just the useful columns, ie Gene Names and Reporter Intensities (AND REMOOVE TMT 10)
reporter.intensities <- grep("Reporter intensity corrected", colnames(filtered.data), value = T)
reporter.counts <- grep("Reporter intensity count", colnames(filtered.data), value = T)

#At least 1 value for each Rep (aka not all 0s)
counts.rep1 <- grep("Reporter intensity count.*Rep 1", colnames(filtered.data), value = T) #filter data into Reps 
counts.rep2 <- grep("Reporter intensity count.*Rep 2", colnames(filtered.data), value = T)
counts.rep3 <- grep("Reporter intensity count.*Rep 3", colnames(filtered.data), value = T)

intensities.rep1 <- grep("Reporter intensity corrected.*Rep 1", colnames(filtered.data), value = T) #filter data into Reps 
intensities.rep2 <- grep("Reporter intensity corrected.*Rep 2", colnames(filtered.data), value = T)
intensities.rep3 <- grep("Reporter intensity corrected.*Rep 3", colnames(filtered.data), value = T)

rep1.valid <- rowSums(select(filtered.data, all_of(counts.rep1)) > 1) > 0 #finding at least 1 nonzero value in each Rep 
rep2.valid <- rowSums(select(filtered.data, all_of(counts.rep2)) > 1) > 0
rep3.valid <- rowSums(select(filtered.data, all_of(counts.rep3)) > 1) > 0

filtered.data <- select(filtered.data, "Gene names", "Protein IDs", "Leading razor protein", all_of(reporter.intensities))

valid.data <- filter(filtered.data, (rep1.valid + rep2.valid + rep3.valid) >= 2)
#valid.data <- filter(filtered.data, (rep1.valid + rep2.valid + rep3.valid) >= 2 | str_detect(`Protein IDs`, "SARS_CoV_2"))

#NORMALIZE DATA 
#1- take column sums of quantifications from each REP 
#2- divide sums by max of those sums 
#3- correct the data by dividing each column by its percent
#4- Then divide all rows by the ref channel
#5- Combine back into one table
#5- Then take the log2 of the dataset 
#6- Combine back into one table with gene names and protein names 

#Replace 0s with min of dataset 
gnames <- valid.data$"Gene names"
normalization.data <- select(valid.data, all_of(reporter.intensities))
normalization.data[normalization.data == 0] <- min(normalization.data[normalization.data > 0])

#REP1
normalization.data.rep1 <- select(normalization.data, all_of(intensities.rep1))
normalization.factor.rep1 <- colSums(normalization.data.rep1)
normalization.factor.rep1 <- normalization.factor.rep1/max(normalization.factor.rep1)
normalization.data.rep1 <- sweep(normalization.data.rep1, 2, normalization.factor.rep1, `/`)
normalization.data.rep1 <- sweep(normalization.data.rep1, 1, normalization.data.rep1[, 9], "/") #dividing all rows by reference channel (TMT9)

#REP2
normalization.data.rep2 <- select(normalization.data, all_of(intensities.rep2))
normalization.factor.rep2 <- colSums(normalization.data.rep2)
normalization.factor.rep2 <- normalization.factor.rep2/max(normalization.factor.rep2)
normalization.data.rep2 <- sweep(normalization.data.rep2, 2, normalization.factor.rep2, `/`)
normalization.data.rep2 <- sweep(normalization.data.rep2, 1, normalization.data.rep2[, 9], "/")

#REP3
normalization.data.rep3 <- select(normalization.data, all_of(intensities.rep3))
normalization.factor.rep3 <- colSums(normalization.data.rep3)
normalization.factor.rep3 <- normalization.factor.rep3/max(normalization.factor.rep3)
normalization.data.rep3 <- sweep(normalization.data.rep3, 2, normalization.factor.rep3, `/`)
normalization.data.rep3 <- sweep(normalization.data.rep3, 1, normalization.data.rep3[, 9], "/")

normalized.data <- cbind(normalization.data.rep1, normalization.data.rep2, normalization.data.rep3)

#NaN and Inf Values? 
normalized.data <- log2(normalized.data)


#remove Ref Channel (TMT9) after normalization (after normalization: all 0s) 
ref.channel <- grep("Reporter intensity corrected 9", colnames(normalized.data), value = T)
normalized.data <- select(normalized.data, -all_of(ref.channel))
design <- filter(design, Condition !="Bridge") #have to remove Ref Channel in design file too 


################################################################################
# Annotations
################################################################################
# replace column names
quant.colnames <- paste(design$Condition, "Rep", design$Replicate)
colnames(normalized.data) <- quant.colnames
# add back names
proteins <- cbind(select(valid.data, "Gene names", "Protein IDs", "Leading razor protein"), normalized.data)
# is it a SARS_CoV_2 protein?
proteins$"SARS_CoV_2" <- ifelse(str_detect(proteins$`Protein IDs`, "SARS_CoV_2"), T, F)
# is it an interactor of a SARS_CoV_2 protein?
proteins <- left_join(proteins, SARS.interactors, by=c("Gene names" = "PreyGeneName"))
# is it mutated?
proteins <- left_join(proteins, H522.mutations, by=c("Gene names" = "GeneName"))
# is it a cell surface or plasma membrane protein?
proteins <- left_join(proteins, human_surface_and_plasma_membrane_protLevel, by=c("Leading razor protein" = "UniProt"))

# update out of date gene names
proteins$`Gene names`[proteins$`Gene names` == "FAM64A"] = "PIMREG"
proteins$`Gene names`[proteins$`Gene names` == "FAM92A1"] = "FAM92A"
proteins$`Gene names`[proteins$`Gene names` == "KIAA0101"] = "PCLAF"
proteins$`Gene names`[proteins$`Gene names` == "C10orf35"] = "FAM241B"
proteins$`Gene names`[proteins$`Gene names` == "MB21D1"] = "CGAS"
proteins$`Gene names`[proteins$`Gene names` == "RNF219"] = "OBI1"
proteins$`Gene names`[proteins$`Gene names` == "C11orf84"] = "SPINDOC"
proteins$`Gene names`[proteins$`Gene names` == "WHSC1"] = "NSD3"
proteins$`Gene names`[proteins$`Gene names` == "C19orf66"] = "SHFL"
proteins$`Gene names`[proteins$`Gene names` == "MFI2"] = "MELTF"
proteins$`Gene names`[proteins$`Gene names` == "Q6ZT62-2;Q6ZT62"] = "BARGIN"
proteins$`Gene names`[proteins$`Gene names` == "C17orf104"] = "MEIOC"
proteins$`Gene names`[proteins$`Gene names` == "LEPREL2"] = "P3H3"
proteins$`Gene names`[proteins$`Gene names` == "FAM129A"] = "NIBAN1"
proteins$`Gene names`[proteins$`Gene names` == "SELH"] = "SELENOH"
proteins$`Gene names`[proteins$`Gene names` == "FAM127A;FAM127C;FAM127B"] = "RTL8C"
proteins$`Gene names`[proteins$`Gene names` == "MUT"] = "MMUT"
proteins$`Gene names`[proteins$`Gene names` == "LOH12CR1"] = "BORCS5"
proteins$`Gene names`[proteins$`Gene names` == "CTGF"] = "CCN2"
proteins$`Gene names`[proteins$`Gene names` == "CASC5"] = "KNL1"


################################################################################
# Stats
################################################################################
fit <- lmFit(select(proteins, 
                    grep("Mock - 96 hr", quant.colnames, value = T), 
                    grep("96 hr", quant.colnames, value = T)),
             model.matrix(~ c(rep("Mock - 96 hr", 3), rep("96 hr", 3))))

fit <- eBayes(fit)
proteins.stats <- topTable(fit, n=nrow(proteins))
proteins <- cbind(proteins.stats, proteins[as.numeric(rownames(proteins.stats)), ])

diff.proteins <- proteins %>% filter(abs(`logFC`) >= log.fc.threshold & adj.P.Val <= 0.05)

# calc enrichment of interactors
num.proteins <- sum(proteins$SARS_CoV_2 == F)
num.interactors <- sum(!is.na(proteins$Bait))
num.diff <- sum(diff.proteins$SARS_CoV_2 == F)
num.diff.interactors <- sum(diff.proteins$SARS_CoV_2 == F & !is.na(diff.proteins$Bait))
phyper(num.diff.interactors, num.interactors, num.proteins-num.interactors, num.diff, lower.tail = T, log.p = FALSE)
# 0.0924 not significant

################################################################################
# Clustering
################################################################################
diff.proteins.avg <- select(diff.proteins, all_of(quant.colnames)) %>% 
  mutate("Mock - 4 hr" = rowMeans(select(., grep("Mock - 4 hr", quant.colnames, value=T)))) %>%
  mutate("4 hr" = rowMeans(select(., grep("^4 hr", quant.colnames, value=T)))) %>%
  mutate("12 hr" = rowMeans(select(., grep("12 hr", quant.colnames, value=T)))) %>%
  mutate("24 hr" = rowMeans(select(., grep("24 hr", quant.colnames, value=T)))) %>%
  mutate("48 hr" = rowMeans(select(., grep("48 hr", quant.colnames, value=T)))) %>%
  mutate("72 hr" = rowMeans(select(., grep("72 hr", quant.colnames, value=T)))) %>%
  mutate("96 hr" = rowMeans(select(., grep("96 hr", quant.colnames, value=T)))) %>%
  mutate("Mock - 96 hr" = rowMeans(select(., grep("Mock - 96 hr", quant.colnames, value=T)))) %>%
  select(-all_of(quant.colnames))
  


proteins.z.scored.avg <- t(scale(t(as.matrix(diff.proteins.avg)), center=T, scale=T))
proteins.z.scored <- t(scale(t(as.matrix(select(diff.proteins, all_of(quant.colnames)))), center=T, scale=T))
rownames(proteins.z.scored) <- diff.proteins$`Gene names`

clustering.results = ConsensusClusterPlus(t(proteins.z.scored.avg),
                                          maxK=15,
                                          reps=1000,
                                          pItem=0.8,
                                          pFeature=1,
                                          #clusterAlg="km",
                                          #distance="euclidean",
                                          clusterAlg="hc",
                                          innerLinkage="complete",
                                          finalLinkage="ward.D2",
                                          distance="pearson",
                                          plot="pdf", 
                                          title="figures/complete_ward2_pearson_1000", 
                                          seed=1262118388.71279)

num.clusters = 8
row.clusters <- clustering.results[[num.clusters]][["consensusClass"]]

# reorder clusters
pos.1 <- which(row.clusters == 1)
pos.2 <- which(row.clusters == 2)
pos.3 <- which(row.clusters == 3)
pos.4 <- which(row.clusters == 4)
pos.5 <- which(row.clusters == 5)
pos.6 <- which(row.clusters == 6)
pos.7 <- which(row.clusters == 7)
pos.8 <- which(row.clusters == 8)
row.clusters[pos.6] <- 1
row.clusters[pos.3] <- 2
row.clusters[pos.2] <- 3
row.clusters[pos.1] <- 4
row.clusters[pos.8] <- 5
row.clusters[pos.4] <- 6
row.clusters[pos.5] <- 7
row.clusters[pos.7] <- 8

diff.proteins$cluster <- row.clusters
proteins <- left_join(proteins, select(diff.proteins, "Protein IDs", "cluster"), by="Protein IDs")



################################################################################
# Profiles
################################################################################
# average replicates
proteins.averaged.condition <- (proteins %>%
                                  pivot_longer(all_of(quant.colnames), names_to="Condition", values_to="FC") %>%
                                  mutate(Condition = str_extract(Condition, ".*(?=( Rep ))")) %>%
                                  group_by(`Gene names`, Condition) %>%
                                  summarise_at(c("FC","cluster","SARS_CoV_2"), mean, na.rm=T))

proteins.averaged.condition.fc.mock <- (proteins.averaged.condition %>%
                                          left_join(filter(., Condition == "Mock - 4 hr"), by="Gene names", suffix=c("",".mock")) %>%
                                          mutate(FC = FC - FC.mock) %>%
                                          mutate(Condition = factor(Condition, levels=c("Mock - 4 hr", "4 hr", "12 hr", "24 hr", "48 hr", "72 hr", "96 hr", "Mock - 96 hr"))) %>%
                                          filter(Condition != "Mock - 4 hr"))

diff.proteins.averaged.condition <- (diff.proteins %>%
                                       pivot_longer(all_of(quant.colnames), names_to="Condition", values_to="FC") %>%
                                       mutate(Condition = str_extract(Condition, ".*(?=( Rep ))")) %>%
                                       group_by(`Gene names`, Condition) %>%
                                       summarise_at(c("FC","cluster","SARS_CoV_2"), mean, na.rm=T))

diff.proteins.averaged.condition.fc.mock <- (diff.proteins.averaged.condition %>%
                                               left_join(filter(., Condition == "Mock - 4 hr"), by="Gene names", suffix=c("",".mock")) %>%
                                               mutate(FC = FC - FC.mock) %>%
                                               mutate(Condition = factor(Condition, levels=c("Mock - 4 hr", "4 hr", "12 hr", "24 hr", "48 hr", "72 hr", "96 hr", "Mock - 96 hr"))) %>%
                                               filter(Condition != "Mock - 4 hr"))




################################################################################
# Generate figures
################################################################################
figPCA <- plotPCA(normalized.data, design)
figHeatmap <- plotHeatmap(proteins.z.scored, diff.proteins, design)
figClusterProfiles <- plotClusterProfiles(diff.proteins.averaged.condition.fc.mock)
figCOVProfiles <- plotCoVProfiles(proteins.averaged.condition.fc.mock)
figVolcano <- plotVolcano(proteins)


F5.top <- arrangeGrob(figPCA, figCOVProfiles, figVolcano,
                  nrow = 1,
                  ncol = 3)
F5.bottom <- arrangeGrob(figHeatmap, figClusterProfiles,
                         widths = unit(c(2.2, 0.95, 1.35, 1.35), "in"),
                         nrow = 1,
                         ncol = 4)

#grid.draw(F5.top)  # to view the plot
saveFig(F5.top, "Figure5_top", 9, 6.85)
saveFig(F5.bottom, "Figure5_bottom", 5, 6.85)

################################################################################
# Work in progress
################################################################################

first.name.universe <- separate(tibble(name=proteins$`Gene names`), name, "first.name", sep=";", remove=F, extra="drop")$first.name
geneIds.universe <- mapIds(org.Hs.eg.db, first.name.universe, 'ENTREZID', 'SYMBOL', multiVals = "first")
geneIds.universe <- geneIds.universe[!is.na(geneIds.universe)]

msig <- msigdbr(species = "Homo sapiens") %>% filter(gs_cat == "H" | (gs_cat == "C2" & gs_subcat == "CP:REACTOME") | (gs_cat == "C5" & gs_subcat %in% c("BP", "CC")))  %>% select(gs_name, entrez_gene)
enrichment.results <- tibble()

for (cluster.id in 1:num.clusters) {
  first.name.cluster <- separate(tibble(name=filter(diff.proteins, cluster == cluster.id)$`Gene names`), name, "first.name", sep=";", remove=F, extra="drop")$first.name
  geneIds.cluster <- mapIds(org.Hs.eg.db, first.name.cluster, 'ENTREZID', 'SYMBOL', multiVals = "first")
  
  em <- enricher(gene=geneIds.cluster,
                 universe=geneIds.universe,
                 TERM2GENE=msig,
                 minGSSize = 10,
                 qvalueCutoff=0.1)
  
  if (!is.null(em)) {
    result <- em@result
    result$cluster <- cluster.id
    result <- filter(result, qvalue <= 0.1)
    if (nrow(result) > 0) {
      result$genes <- apply(result, 1, function(x) {str_c(sort(getSYMBOL(unlist(str_split(x[["geneID"]],"/")), data='org.Hs.eg')),collapse=",")} )
      enrichment.results <- rbind(enrichment.results, result)
    }
  }
}

# convert to matrix format
enrichment.results = pivot_wider(enrichment.results, id_cols = Description, names_from = cluster, values_from=c(pvalue, p.adjust, qvalue, Count, GeneRatio, BgRatio, genes))

enrichment.results$category <- ""
enrichment.results$category[str_detect(enrichment.results$Description, "^HALLMARK")] <- "Hallmark"
enrichment.results$category[str_detect(enrichment.results$Description, "^REACTOME")] <- "Reactome"
enrichment.results$category[str_detect(enrichment.results$Description, "^GO")] <- "GO"

write_tsv(filter(enrichment.results, category == "Hallmark"), str_c("~/Downloads/enrichment_hallmark.tsv"), na = "")
write_tsv(filter(enrichment.results, category == "Reactome"), str_c("~/Downloads/enrichment_reactome.tsv"), na = "")
write_tsv(filter(enrichment.results, category == "GO"), str_c("~/Downloads/enrichment_GO.tsv"), na = "")

enrichment.results.pruned <- filter(enrichment.results, Description %in% c("REACTOME_CELL_CYCLE_CHECKPOINTS", 
                                                                           "REACTOME_SIGNALING_BY_RHO_GTPASES", 
                                                                           "REACTOME_DNA_REPAIR", 
                                                                           "REACTOME_DNA_REPLICATION", 
                                                                           "REACTOME_REGULATION_OF_TP53_ACTIVITY", 
                                                                           "REACTOME_SUMOYLATION",
                                                                           
                                                                           "GO_INTRINSIC_COMPONENT_OF_PLASMA_MEMBRANE", 
                                                                           "GO_SPINDLE", 
                                                                           "GO_ANAPHASE_PROMOTING_COMPLEX", 
                                                                           "GO_CULLIN_RING_UBIQUITIN_LIGASE_COMPLEX", 
                                                                           "GO_UBIQUITIN_LIGASE_COMPLEX", 
                                                                           "GO_REPLICATION_FORK", 
                                                                           "GO_PHAGOCYTIC_VESICLE_MEMBRANE", 
                                                                           "GO_ENDOCYTIC_VESICLE_MEMBRANE", 
                                                                           "GO_ORGANELLE_ENVELOPE_LUMEN",
                                                                           
                                                                           "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
                                                                           "HALLMARK_INTERFERON_GAMMA_RESPONSE"))

enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "intrinsic component of plasma membrane")] <- "Plasma Membrane"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_INTERFERON_ALPHA_RESPONSE")] <- "Interferon alpha signaling"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_INTERFERON_GAMMA_RESPONSE")] <- "Interferon gamma signaling"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_TNFA_SIGNALING_VIA_NFKB")] <- "TNFA signaling via NFKB"

num.in.cluster <- (diff.proteins.averaged.condition.fc.mock %>%
                     group_by(cluster) %>% 
                     summarise(n=n() / (nlevels(Condition)-1)))

enrichment.results.pruned$ratio <- enrichment.results.pruned$Count / num.in.cluster$n[enrichment.results.pruned$cluster]

print(ggplot(data = enrichment.results.pruned, aes(x=as.factor(cluster), y=Description, size=Count, color=-log10(qvalue)))
      + geom_point()
      + xlab("Cluster")
      + ylab("")
      + theme.basic
      + theme(aspect.ratio = NULL, axis.title.y=element_blank()))




#colnames(proteins.annot)[8:31] <- paste(design$Condition, "Rep", design$Replicate)
#write_tsv(proteins.annot, "~/Box/H522-paperoutline/Figures/Figure5/proteinGroups_log2_relative_to_bridge.tsv")
#write_tsv(diff.proteins.melted, "~/Downloads/clusters.tsv")

################################################################################
# Write processed data
################################################################################
write_tsv(proteins, here("data_processed/proteinsNormedToBridge.txt"))






################################################################################
# Temp exploration
################################################################################
# proteins <- select(proteins, -c('logFC', 'adj.P.Val', 'AveExpr', 't', 'P.Value', 'B'))
# 
# threshold <- log(2)
# 
# # 4hr mock vs 4hr
# fit <- lmFit(select(proteins, 
#                     grep("^Mock - 4 hr", quant.colnames, value = T), 
#                     grep("^4 hr", quant.colnames, value = T)),
#              model.matrix(~ c(rep("Mock - 4 hr", 3), rep("4 hr", 3))))
# 
# fit <- eBayes(fit)
# proteins.stats.4hr <- topTable(fit, n=nrow(proteins))
# proteins.4hr <- cbind(proteins.stats.4hr, proteins[as.numeric(rownames(proteins.stats.4hr)), ])
# 
# diff.proteins.4hr <- proteins.4hr %>% filter(abs(`logFC`) >= threshold & adj.P.Val <= 0.05)
# 
# # 4hr vs 12 hr
# fit <- lmFit(select(proteins, 
#                     grep("^4 hr", quant.colnames, value = T), 
#                     grep("^12 hr", quant.colnames, value = T)),
#              model.matrix(~ c(rep("4 hr", 3), rep("12 hr", 3))))
# 
# fit <- eBayes(fit)
# proteins.stats.12hr <- topTable(fit, n=nrow(proteins))
# proteins.12hr <- cbind(proteins.stats.12hr, proteins[as.numeric(rownames(proteins.stats.12hr)), ])
# 
# diff.proteins.12hr <- proteins.12hr %>% filter(abs(`logFC`) >= threshold & adj.P.Val <= 0.05)
# 
# # 12hr vs 24 hr
# fit <- lmFit(select(proteins, 
#                     grep("^12 hr", quant.colnames, value = T), 
#                     grep("^24 hr", quant.colnames, value = T)),
#              model.matrix(~ c(rep("12 hr", 3), rep("24 hr", 3))))
# 
# fit <- eBayes(fit)
# proteins.stats.24hr <- topTable(fit, n=nrow(proteins))
# proteins.24hr <- cbind(proteins.stats.24hr, proteins[as.numeric(rownames(proteins.stats.24hr)), ])
# 
# diff.proteins.24hr <- proteins.24hr %>% filter(abs(`logFC`) >= threshold & adj.P.Val <= 0.05)
# 
# # 24hr vs 48 hr
# fit <- lmFit(select(proteins, 
#                     grep("^24 hr", quant.colnames, value = T), 
#                     grep("^48 hr", quant.colnames, value = T)),
#              model.matrix(~ c(rep("24 hr", 3), rep("48 hr", 3))))
# 
# fit <- eBayes(fit)
# proteins.stats.48hr <- topTable(fit, n=nrow(proteins))
# proteins.48hr <- cbind(proteins.stats.48hr, proteins[as.numeric(rownames(proteins.stats.48hr)), ])
# 
# diff.proteins.48hr <- proteins.48hr %>% filter(abs(`logFC`) >= threshold & adj.P.Val <= 0.05)
# 
# # 48hr vs 72 hr
# fit <- lmFit(select(proteins, 
#                     grep("^48 hr", quant.colnames, value = T), 
#                     grep("^72 hr", quant.colnames, value = T)),
#              model.matrix(~ c(rep("48 hr", 3), rep("72 hr", 3))))
# 
# fit <- eBayes(fit)
# proteins.stats.72hr <- topTable(fit, n=nrow(proteins))
# proteins.72hr <- cbind(proteins.stats.72hr, proteins[as.numeric(rownames(proteins.stats.72hr)), ])
# 
# diff.proteins.72hr <- proteins.72hr %>% filter(abs(`logFC`) >= threshold & adj.P.Val <= 0.05)
# 
# # 72hr vs 96 hr
# fit <- lmFit(select(proteins, 
#                     grep("^72 hr", quant.colnames, value = T), 
#                     grep("^96 hr", quant.colnames, value = T)),
#              model.matrix(~ c(rep("72 hr", 3), rep("96 hr", 3))))
# 
# fit <- eBayes(fit)
# proteins.stats.96hr <- topTable(fit, n=nrow(proteins))
# proteins.96hr <- cbind(proteins.stats.96hr, proteins[as.numeric(rownames(proteins.stats.96hr)), ])
# 
# diff.proteins.96hr <- proteins.96hr %>% filter(abs(`logFC`) >= threshold & adj.P.Val <= 0.05)
# 
# 
# diff.union <- rbind(diff.proteins.4hr, diff.proteins.12hr, diff.proteins.24hr, diff.proteins.48hr, diff.proteins.72hr, diff.proteins.96hr)
# 
# diff.venn <- full_join(diff.union, diff.proteins, by="Gene names")
# both <- sum(!is.na(diff.venn$`Protein IDs.x`) & !is.na(diff.venn$`Protein IDs.y`))
# print(sum(!is.na(diff.venn$`Protein IDs.x`)) - both) # pairwise only
# print(sum(!is.na(diff.venn$`Protein IDs.y`)) - both) # 96 only
# print(both)
# pairwise.only <- filter(diff.venn, is.na(`Protein IDs.y`))
# 
# library(timecourse)
# 
# timecourse.infected.data <- select(normalized.data, -grep("Mock", colnames(normalized.data), value=T))
# assay.1 <- rep(c("1", "2", "3"), each = 6)
# time.grp.1 <- rep(1:6, 3)
# reps <- rep(3, nrow(timecourse.infected.data))
# timecourse.infected.results <- mb.long(timecourse.infected.data, method="1D", times=6, reps=reps, rep.grp = assay.1, time.grp = time.grp.1) 
# # get top 265
# plotProfile(timecourse.infected.results, type="b", gnames=gnames, legloc=c(2,15), pch=c("1","2","3"), xlab="Time Point", rank = 1, col = c("pink", "black", "green")) 
# #use gid= to look for a specific protein's profile 
# sig.timecourse <- tibble("Gene names" = gnames[order(-timecourse.infected.results[["HotellingT2"]])[1:265]])
# sig.timecourse <- left_join(sig.timecourse, diff.proteins)
# print(sig.timecourse$`Gene names`[which(is.na(sig.timecourse$`B`))])

