library(here)
library(limma)
library(topGO)
library(timecourse)
library(ComplexHeatmap)
library(cluster)
library(ConsensusClusterPlus)

source(here("common.R"))
select <- get(x = "select", pos = "package:dplyr") # deal with function masking

source(here("RequantifyProteomics.R")) # only necessary to run once

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
                           column_title_gp = gpar(fontface = "bold", fontsize=7),
                           show_row_names = F,
                           row_names_gp = gpar(fontsize = 5),
                           column_names_gp = gpar(fontsize = 6),
                           name = "z-score",
                           # colors #-max(abs(quant)), max(abs(quant)
                           col = colorRamp2(seq(-3, 3, length=256), rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))),
                           # legends
                           show_heatmap_legend = F,
                           heatmap_legend_param = list(color_bar = "continuous",
                                                       title_gp = gpar(fontsize = 8),
                                                       labels_gp = gpar(fontsize = 6),
                                                       grid_width = unit(2,"mm")),
                           # clustering
                           cluster_columns = F,
                           clustering_distance_rows = "euclidean",
                           cluster_row_slices = F,
                           show_row_dend = F,
                           #cluster_rows = row.clusters,
                           #split=7,
                           split = row.clusters,
                           # labels
                           row_title_rot = 0,
                           row_title_gp = gpar(fontsize=8),
                           #row_title = NULL,
                           # size
                           width=ncol(proteins.z.scored)*0.2
  )
  
  cluster.rowAnnot <- rowAnnotation(block = anno_block(gp = gpar(fill = brewer.pal(num.clusters, "Set3"), col = NA)),
                                    width = unit(2, "mm"))
  
  gene.rowAnnot = rowAnnotation(gene.name = anno_mark(at = label.pos, 
                                                      labels = labels,
                                                      labels_gp = gpar(fontsize = 6), 
                                                      link_width = unit(4, "mm"),
                                                      padding = unit(0.5, "mm")))
  
  ht_list = heatmap.proteins + cluster.rowAnnot + gene.rowAnnot
  
  gb_heatmap = grid.grabExpr(draw(ht_list), height=6, width=4)
}

################################################################################
# Profile plots of clusters
################################################################################
plotClusterProfiles <- function(diff.proteins)
{
  num.in.cluster <- (diff.proteins.averaged.condition.fc.mock %>%
                       group_by(cluster) %>% 
                       summarise(n=n() / (nlevels(Condition)-1)))
  
  cluster.labels <- paste("Cluster ", num.in.cluster$cluster, " (n=", num.in.cluster$n, ")", sep="")
  names(cluster.labels) <- seq(1:num.clusters)
  
  representative.profiles <- (diff.proteins.averaged.condition.fc.mock %>% 
                                group_by(Condition, cluster) %>% 
                                summarise_at(c("FC"),  mean))
  
  p <- (ggplot(diff.proteins, aes(x=Condition, 
                                         y=FC, 
                                         group=interaction(`Gene names`, cluster),
                                         color=as.factor(cluster),
                                         facet=cluster,
                                         linetype=as.factor(SARS_CoV_2)))
        
        + geom_line(size=0.5, alpha=0.75)
        + geom_line(size=0.5, data=representative.profiles, aes(x=Condition, y=FC, group=as.factor(cluster), linetype="solid"), color="black")
        
        + scale_y_continuous(name = "log2(Intensity / 4h Mock)", breaks=seq(-5,5,0.5))
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
        
        + geom_line(size=0.5, alpha=0.75, color=COV2.color)
        + geom_point(size=1, color=COV2.color)
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
SARS.interactors.krogan <- select(read_csv(here("annotations/SARS2_interactome.csv")), c("Bait.Krogan", "PreyGeneName"))
SARS.interactors.mann <- select(read_csv("annotations/Mann_Interactors_Caco2.csv"), c("Bait.Mann", "gene_name"))
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

filtered.data <- select(filtered.data, "Gene names", "Protein IDs", all_of(reporter.intensities))

valid.data <- filter(filtered.data, (rep1.valid + rep2.valid + rep3.valid) >= 2 | str_detect(`Protein IDs`, "SARS_CoV_2"))

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
proteins <- cbind(select(valid.data, "Gene names", "Protein IDs"), normalized.data)
# is it a SARS_CoV_2 protein?
proteins$"SARS_CoV_2" <- ifelse(str_detect(proteins$`Protein IDs`, "SARS_CoV_2"), T, F)
# is it an interactor of a SARS_CoV_2 protein?
proteins <- left_join(proteins, SARS.interactors.krogan, by=c("Gene names" = "PreyGeneName"))
proteins <- left_join(proteins, SARS.interactors.mann, by=c("Gene names" = "gene_name"))
proteins$Bait <- str_c(str_replace_na(proteins$Bait.Krogan, replacement = ""), 
                       str_replace_na(proteins$Bait.Mann, replacement = ""), sep = ";")
proteins$Bait <- apply(proteins, 1, function(x) {
  y <- c(x["Bait.Krogan"], x["Bait.Mann"])
  str_c(base::unique(y[!is.na(y)]), collapse = ";")})
proteins <- proteins %>% mutate_all(na_if,"") #replace "" with NAs in Bait Column 

# interferon response 
interferon.response.alpha <- read_csv(here("annotations/interferon_response_alpha.txt"))
interferon.response.beta <- read_csv(here("annotations/interferon_response_beta.txt"))
interferon.response.gamma <-  read_csv(here("annotations/interferon_response_gamma.txt"))
interferon.regulation.type1 <-  read_csv(here("annotations/regulation_of_type_I_interferon_mediated_signaling_pathway.txt"))
interferon.regulation.type2.immune.response <-  read_csv(here("annotations/regulation_type_2_immune_response.txt"))
antigen.processing.presentation <- read_csv(here("annotations/antigen_processing_and_presentation.txt"))

summarized.gene.names <- aggregate(proteins["Bait"], proteins["Gene names"], 
               FUN = function(X) paste(unique(X), collapse=", "))


proteins <- left_join(proteins, interferon.response.alpha, by=c("Gene names" = "GO_CELLULAR_RESPONSE_TO_INTERFERON_ALPHA"))
proteins <- left_join(proteins, interferon.response.beta, by=c("Gene names" = "GO_CELLULAR_RESPONSE_TO_INTERFERON_BETA"))
proteins <- left_join(proteins, interferon.response.gamma, by=c("Gene names" = "GO_RESPONSE_TO_INTERFERON_GAMMA"))
proteins <- left_join(proteins, interferon.regulation.type1, by=c("Gene names" = "GO_REGULATION_OF_TYPE_I_INTERFERON_MEDIATED_SIGNALING_PATHWAY"))
proteins <- 
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
                                          clusterAlg="km",
                                          distance="euclidean",
                                          plot="pdf", 
                                          title="figures/km_1000", 
                                          seed=1262118388.71279)

num.clusters = 9
row.clusters <- clustering.results[[num.clusters]][["consensusClass"]]

# reorder clusters
# 7->1, 2->2, 3->3, 1->4, 4->5, 9->6, 5->7, 6->8, 8->9
pos.1 <- which(row.clusters == 1)
pos.2 <- which(row.clusters == 2)
pos.3 <- which(row.clusters == 3)
pos.4 <- which(row.clusters == 4)
pos.5 <- which(row.clusters == 5)
pos.6 <- which(row.clusters == 6)
pos.7 <- which(row.clusters == 7)
pos.8 <- which(row.clusters == 8)
pos.9 <- which(row.clusters == 9)
row.clusters[pos.7] <- 1
row.clusters[pos.2] <- 2
row.clusters[pos.3] <- 3
row.clusters[pos.1] <- 4
row.clusters[pos.4] <- 5
row.clusters[pos.9] <- 6
row.clusters[pos.5] <- 7
row.clusters[pos.6] <- 8
row.clusters[pos.8] <- 9

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
                         widths = unit(c(4, 1.5, 2), "in"),
                         nrow = 1,
                         ncol = 3)

#grid.draw(F5.top)  # to view the plot
saveFig(F5.top, "Figure5_top", 9, 7.5)
saveFig(F5.bottom, "Figure5_bottom", 6, 7.5)

################################################################################
# Work in progress
################################################################################

#colnames(proteins.annot)[8:31] <- paste(design$Condition, "Rep", design$Replicate)
#write_tsv(proteins.annot, "~/Box/H522-paperoutline/Figures/Figure5/proteinGroups_log2_relative_to_bridge.tsv")


#sampleGOdata.1 <- new("topGOdata",
#                    description = "Simple session", ontology = "BP",
#                    allGenes = proteins.annot$Gene.names, geneSel = diff.proteins$Gene.names,
#                    nodeSize = 10,
#                    annot = annFUN.org)

#write_tsv(diff.proteins.melted, "~/Downloads/clusters.tsv")

################################################################################
# Write processed data
################################################################################
write_tsv(proteins, here("data_processed/proteinsNormedToBridge.txt"))




