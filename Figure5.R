library(plyr)
library(here)
library(edgeR)
source(here("common.R"))


#library(stringr)
#library(plotly)

# these are raw counts per gene, no normalization or transforamtion
load(here("data/Figure 5/rep1_hum_covid_counts.RData"))
load(here("data/Figure 5/rep2_hum_covid_counts.RData"))
#RNA is in log10 viral/human
load(here("data/Figure 5/exp_set_up_covid_rna_Rep1.RData"))
load(here("data/Figure 5/exp_set_up_covid_rna_Rep2.RData"))

load(here("data/Figure 5/limma_voom_deg_data.RData"))
load(here("data/Figure 5/norm_lcpm.RData")) 

# fix sample swap
exp_set_up_covid_rna_Rep1[which(exp_set_up_covid_rna_Rep1$dig_one==1 &
                                  exp_set_up_covid_rna_Rep1$time==48), ]$time = NA
exp_set_up_covid_rna_Rep1[which(exp_set_up_covid_rna_Rep1$dig_one==1 &
                                  exp_set_up_covid_rna_Rep1$time==72 ), ]$time = 48
exp_set_up_covid_rna_Rep1[which(exp_set_up_covid_rna_Rep1$dig_one==1 &
                                  is.na(exp_set_up_covid_rna_Rep1$time)), ]$time = 72

pal <- colors.MOI

data <- rbind(exp_set_up_covid_rna_Rep1, exp_set_up_covid_rna_Rep2)

data$time <- as.numeric(as.character(data$time))

lim_data <- data[, names(data) %in% c("my_moi", "time", "logR_hum_cov")]

cdata <- ddply(lim_data, c("time", "my_moi"), summarise,
               N    = length(logR_hum_cov),
               mean = mean(logR_hum_cov),
               sd   = sd(logR_hum_cov),
               se   = sd / sqrt(N)
)

min.x <- min(cdata$time)
max.x <- max(cdata$time)
x.range <- max.x - min.x
x.limit <- max(abs(min.x), abs(max.x))

SARS2.profile <- ggplot(cdata, aes(x = time, y = mean, color = my_moi))+
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 2.5))+
  geom_line()+
  geom_point()+
  scale_color_manual(values = pal)+
  xlab("Time (hours post infection)")+
  ylab("log10(SARS-CoV-2/H.Sapiens)")+
  ggtitle("SARS-CoV-2 Viral mRNA Expression")+
  labs(color = "MOI")+
  scale_x_continuous(limits=c(min.x - x.range/20, x.limit + x.range/20), breaks = c(0, 4, 24, 48, 72, 96)) +
  theme.basic+
  theme(legend.position = "none", aspect.ratio = 0.75)





###########################################################################################################################
# PCA
###########################################################################################################################


#bring replicates together in one data structure
raw_counts <- cbind(rep1_hum_covid_counts, rep2_hum_covid_counts)
meta_data <- rbind(exp_set_up_covid_rna_Rep1, exp_set_up_covid_rna_Rep2)
meta_data$replicate <- c(rep("Rep1", time = 24), rep("Rep2", time = 24))
meta_data$sample <- paste(meta_data$replicate, meta_data$my_moi, meta_data$time, sep = "_")
names(raw_counts) <- as.character(meta_data$sample)

meta_data <- subset(meta_data, meta_data$my_moi %in% c("0.25",  "1.0"))
raw_counts <- raw_counts[, names(raw_counts) %in% c(meta_data$sample)]

### Set up for Limma Package Data Structure Use
groups <- c()
for(i in meta_data$time){
  if(i == 0){
    groups <- c(groups, 0)
  }else{
    groups <- c(groups, 1)
  }
}

group_name_vector <- c()
for(i in groups){
  if(i == 0){
    group_name_vector <- c(group_name_vector, "mock")
  }else{
    group_name_vector <- c(group_name_vector, "sars")
  }
}


group_name_vector <- as.factor(group_name_vector)
genes <- as.character(row.names(raw_counts))
samples <- as.character(names(raw_counts))
raw_count_matrix <- as.matrix(raw_counts)

#cols_to_filter <- which(meta_data$my_moi < 0.1)# | meta_data$replicate == "Rep1")
#group_name_vector <- group_name_vector[-cols_to_filter]
#samples <- samples[-cols_to_filter]
#raw_count_matrix <- raw_count_matrix[, -cols_to_filter]

##### Create LIMMA DEGList
my_DGEList <- DGEList(raw_count_matrix, genes = genes, samples = samples, group = group_name_vector , lib.size = colSums(raw_count_matrix))
cpm <- cpm(my_DGEList)
lcpm <- cpm(my_DGEList, log=TRUE)

###### I adjusted this value below to where the peak near -8 in the susequent plot goes away. 
keep.exprs <- rowSums(cpm>=5)>=8 #35 and 24
filt_my_DGEList <- my_DGEList[keep.exprs, , keep.lib.sizes = FALSE]
filt_lcpm <- cpm(filt_my_DGEList, log=TRUE)


######## now will calculate normalization factors
norm_filt_my_DGEList <- calcNormFactors(filt_my_DGEList, method = "TMM")
norm_lcpm <- cpm(norm_filt_my_DGEList, log=TRUE)
#### AT THIS POINT norm_lcpm is both filtered of low expresion transcripts and normalized. 


meta_data_high <- subset(meta_data, meta_data$my_moi > 0.1)
#meta_data_high$time <- as.numeric(meta_data_high$time)

norm_lcpm_df <- data.frame(norm_lcpm)
norm_lcpm_high <- norm_lcpm_df[, names(norm_lcpm_df) %in% as.character(meta_data_high$sample)]


data <- as.matrix(norm_lcpm_high)
data.t <- as.data.frame(t(data))
pca <- prcomp(data.t, center=T, scale=F)
eigen <- pca$sdev^2
variance <- (eigen/sum(eigen)) * 100
variance <- format(round(variance, 2), nsmall=2) # show 2 digits after decimal

design <- meta_data_high
row.names(design) <- design$sample

plotting.data <- cbind(data.frame(pca$x), design)
y.range <- max(plotting.data$PC2) - min(plotting.data$PC2)
min.x <- min(plotting.data$PC1)
max.x <- max(plotting.data$PC1)
x.range <- max.x - min.x

plotting.data$Label <- plotting.data$sample

pca_plot <- (ggplot(plotting.data, aes(x=PC1, y=PC2, shape=as.factor(my_moi), color=time, label=Label)) 
             + geom_point(size=2, alpha=0.8)
             #+ geom_text(size=4, nudge_y = y.range/20)
             
             + scale_x_continuous(limits=c(min.x - x.range/10, max.x + x.range/10))
             + scale_shape_discrete(name = "MOI")
             + scale_color_discrete(name = "time")
             #+ scale_color_brewer(guide='none', type="qual", palette="Set1")
             
             + xlab(paste("PC1 (", variance[1], "% explained variance)", sep=""))
             + ylab(paste("PC2 (", variance[2], "% explained variance)", sep=""))
             + ggtitle(label = "Principal Component Analysis",
                       subtitle = paste("n =", formatC(nrow(norm_lcpm_high), big.mark=","), "Genes"))
             
             + theme.basic
             + theme(legend.position = "none", 
                     aspect.ratio = 0.75,
                     plot.subtitle = element_text(size=6, hjust = 0.5, face = "italic"))
)



###########################################################################################################################
# Volcano
###########################################################################################################################


proteins <- subset(limma_voom_deg_data, limma_voom_deg_data$time == 96)

log.stats.threshold = -log10(0.005)
log.fc.threshold = 2

IFN_genes <- c("ISG15", "IRF9", "MX1", "IFI35", "OAS3")

proteins$logFC <- -proteins$logFC
min.x <- min(proteins$logFC)
max.x <- max(proteins$logFC)
x.range <- max.x - min.x
x.limit <- max(abs(min.x), abs(max.x))

proteins$log.adj.P.Val <- -log10(proteins$adj.P.Val)
y.range <- max(proteins$log.adj.P.Val) - min(proteins$log.adj.P.Val)

proteins$significance <- "none"

levels(proteins$significance) <- c("none", "sig fold-change", "sig statistic", "sig both", "SARS-CoV-2", "IFN")

proteins$significance[abs(proteins$logFC) >= log.fc.threshold] = "sig fold-change"
proteins$significance[proteins$log.adj.P.Val >= -log10(0.05)] = "sig statistic"
proteins$significance[abs(proteins$logFC) >= log.fc.threshold & proteins$log.adj.P.Val >= -log10(0.05)] = "sig both"
proteins$significance[proteins$genes %in% IFN_genes] = "IFN"



sars2_genes <- read.table(here("data/Figure 5/sars_cov2_genes_brief.txt"))

covid_genes <- as.character(sars2_genes$V1)


my_dict <- as.character(sars2_genes$V2)
names(my_dict) <- as.character(sars2_genes$V1)

my_labels <- c()
for(i in as.character(proteins$genes)){
  if(i %in%  covid_genes){
    my_labels <- c(my_labels, my_dict[[i]])
    proteins$significance[length(my_labels)] = "SARS-CoV-2"
  }else if (i %in% IFN_genes) {
    my_labels <- c(my_labels, i)
  }else{
    my_labels <- c(my_labels, "")
  }
}



proteins$Label <- my_labels

#proteins$Label <- str_replace(proteins$`Gene names`, "SARS_CoV_2", "")
#proteins$Label[!proteins$SARS_CoV_2] <- ""

volcano <- (ggplot(data=proteins, aes(x=logFC, y=log.adj.P.Val, label = Label, color=significance))
      + geom_point(size=0.5) 
      + geom_text(nudge_y = 0.05 * y.range, nudge_x = 0.05 * x.range, size=2)
      
      + geom_vline(xintercept = -log.fc.threshold, linetype="dashed", color=grey)
      + geom_vline(xintercept = log.fc.threshold, linetype="dashed", color=grey)
      + geom_hline(yintercept = -log10(0.05), linetype="dashed", color=grey)
      
      + scale_color_manual(values=c("none"="#AAAAAA55", "sig fold-change"="#AAAAAA55", "sig statistic"="#AAAAAA55", 
                                    "sig both"="#1f78b455", "SARS-CoV-2"="#fb8072", "IFN"="#39b54a"))
      
      + scale_x_continuous(limits=c(min.x - x.range/20, x.limit + x.range/20), breaks = seq(-9,18,3))
      
      + ylab("-log10(q-value)")
      + xlab("log2(fold-change)")
      
      + ggtitle("96 hours post infection vs mock")
      
      + theme.basic
      + theme(legend.position = "none",
              aspect.ratio = 0.75)
)






arranged.QC <- arrangeGrob(SARS2.profile, pca_plot, volcano,
                           nrow = 1,
                           ncol = 3)



saveFig(arranged.QC, "Figure5_B", 9, 5.7)















####################################################################################################################################################
# HEATMAP
####################################################################################################################################################

DE_Genes <- limma_voom_deg_data
##### 
##### padjusted for all genes to be considered "significant"
#####
p_adjusted_threshold <- 0.005
logFC_threshold <- 2

#DE_Genes <- subset(DE_Genes, DE_Genes$time < 96)
Sig_Genes <- subset(DE_Genes, DE_Genes$multi_sample_qvalue < p_adjusted_threshold )
Sig_Genes <- subset(Sig_Genes, abs(Sig_Genes$logFC) > logFC_threshold  )

Sig_Genes <- names(table(Sig_Genes$genes))

##### 
##### padjusted for all genes to be PLOTTED IN MASS HEAT MAP
#####

time_sort <- test_annotation_data[order(test_annotation_data$Time), ]

norm_lcpm <- data.frame(norm_lcpm)
high_moi_norm_lcpm <- norm_lcpm[, names(norm_lcpm) %in% c(time_sort$Sample)]
high_moi_norm_lcpm <- high_moi_norm_lcpm[, as.character(time_sort$Sample)]
time_sort_plot_matrix <- high_moi_norm_lcpm   

#[row.names(high_moi_norm_lcpm) %in% subset(out_data_impulseDE, out_data_impulseDE$padj < p_adjusted_threshold)$Gene, ]

#pheatmap(log(time_sort_plot_matrix, base = 2), cluster_cols = FALSE, scale = "row", show_rownames = FALSE)

plot_genes <- c(Sig_Genes, tail(row.names(norm_lcpm), n = 12))

COVID_genes <- tail(row.names(norm_lcpm), n = 12)
#plot_genes <- subset(plot_genes, plot_genes$isTransient == FALSE)

time_sort_plot_matrix_all <- time_sort_plot_matrix[row.names(time_sort_plot_matrix) %in% plot_genes, ]

### average the two replicates for heat map plotting, and oder by time post infection
plot_matrix_all <- time_sort_plot_matrix_all[, c(1,2,  5,6,  9,10,  13,14,  17,18, 21,22 ) ] + time_sort_plot_matrix_all[, (c(1,2,  5,6,  9,10,  13,14,  17,18, 21,22  ) + 2)]
plot_matrix_all <- plot_matrix_all/2

scaled_plot_matrix_all <- scale(t(plot_matrix_all), scale = TRUE, center = TRUE)

scaled_plot_matrix_all <- t(scaled_plot_matrix_all)

#png(file = "all_sig_genes_heatmap.png", height = 900, width = 900)
#print(pheatmap(plot_matrix_all, cluster_cols = FALSE, scale = "row", border_color = NA))
#dev.off()


annotation <- subset(time_sort, time_sort$Batch %in% c("B1"))

annotation$MOI <- rep(c("MOI0.25", "MOI1.0"), times = 6)
################################################################################
# Clustering
################################################################################


diff.proteins.avg <- plot_matrix_all

proteins.z.scored.avg <- t(scale(t(as.matrix(diff.proteins.avg)), center=T, scale=T))


clustering.results = ConsensusClusterPlus(t(proteins.z.scored.avg),
                                          maxK=13,
                                          reps=50,
                                          pItem=0.8,
                                          pFeature=1,
                                          #clusterAlg="km",
                                          #distance="euclidean",
                                          clusterAlg="hc",
                                          innerLinkage="complete",
                                          finalLinkage="ward.D2",
                                          distance="pearson",
                                          plot="pdf", 
                                          title=here(), 
                                          seed=1262118388.71279)

num.clusters = 12
row.clusters_orig <- clustering.results[[num.clusters]][["consensusClass"]]

row.clusters <- row.clusters_orig
# reorder clusters
pos.1 <- which(row.clusters == 1)
pos.2 <- which(row.clusters == 2)
pos.3 <- which(row.clusters == 3)
pos.4 <- which(row.clusters == 4)
pos.5 <- which(row.clusters == 5)
pos.6 <- which(row.clusters == 6)
pos.7 <- which(row.clusters == 7)
pos.8 <- which(row.clusters == 8)
pos.9 <- which(row.clusters == 9)
pos.10 <- which(row.clusters == 10)
pos.11 <- which(row.clusters == 11)
pos.12 <- which(row.clusters == 12)
row.clusters[pos.1] <- 7# 1
row.clusters[pos.2] <- 6# 2
row.clusters[pos.3] <- 7# 1
row.clusters[pos.4] <- 3 #4
row.clusters[pos.5] <- 7# 1
row.clusters[pos.6] <- 7# 1
row.clusters[pos.7] <- 4# 5
row.clusters[pos.8] <- 6# 2
row.clusters[pos.9] <- 4# 5
row.clusters[pos.10] <- 1# 7
row.clusters[pos.11] <- 5# 3
row.clusters[pos.12] <- 2# 6

diff.proteins <- plot_matrix_all
diff.proteins$cluster <- row.clusters


#proteins <- left_join(proteins, select(diff.proteins, "Protein IDs", "cluster"), by="Protein IDs")

################################################################################
# SET UP SARS COV2 LABELS
################################################################################

sars2_genes <- read.table("sars_cov2_genes_brief.txt")

my_dict <- as.character(sars2_genes$V2)
names(my_dict) <- as.character(sars2_genes$V1)





################################################################################
# Heatmap
################################################################################

colnames(proteins.z.scored.avg) <- paste(annotation$Time, annotation$MOI)

# reorder columns

#label.sars2.pos <-  which(diff.proteins$SARS_CoV_2 == T)
#label.int.pos <-  which(!is.na(diff.proteins$Bait))
#label.pos <- c(label.sars2.pos, label.int.pos)
label.pos <- which(row.names(diff.proteins) %in% COVID_genes)

#labels.sars2 <- diff.proteins$`Gene names`[label.sars2.pos]
#labels.int <- paste(diff.proteins$`Gene names`[label.int.pos], " (", diff.proteins$Bait[label.int.pos], ")", sep="")
#labels <- c(labels.sars2, labels.int)
holder <- row.names(diff.proteins)[label.pos]
labels <- c()
for(i in holder){
  labels <- c(labels, my_dict[[i]])
}

heatmap.proteins <- Heatmap(proteins.z.scored.avg,
                            # labels
                            column_title = paste("Differentially Expressed Genes\nn = ", nrow(proteins.z.scored.avg), sep=""),
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
                            width=ncol(proteins.z.scored.avg)*0.2
)

cluster.rowAnnot <- rowAnnotation(block = anno_block(gp = gpar(fill = brewer.pal(num.clusters, "Set3"), col = NA)),
                                  width = unit(3, "mm"))

gene.rowAnnot = rowAnnotation(gene.name = anno_mark(at = label.pos, 
                                                    labels = labels,
                                                    labels_gp = gpar(fontsize = 7), 
                                                    link_gp = gpar(lwd = 0.5),
                                                    link_width = unit(4, "mm"),
                                                    padding = unit(0.5, "mm")))

ht_list = heatmap.proteins + cluster.rowAnnot + gene.rowAnnot

gb_heatmap = grid.grabExpr(draw(ht_list), height=5, width=2.5)

annotated.degs.lcpm.clusters <- diff.proteins
colnames(annotated.degs.lcpm.clusters) <- c(paste(annotation$Time, annotation$MOI), "cluster")
annotated.degs.lcpm.clusters$gene <- row.names(annotated.degs.lcpm.clusters)
write.table(annotated.degs.lcpm.clusters, file = "Supplemental_Table_normLCPM_DEGs_by_Cluster_Group.csv", row.names = FALSE, sep = ",")

################################################################################
# Enrichment heatmap
################################################################################

diff.proteins$First_GeneID <- as.character(row.names(diff.proteins))
# plotEnrichment <- function(proteins, diff.proteins, num.clusters){

geneIds.universe <- as.character(row.names(norm_lcpm))

msig <- msigdbr(species = "Homo sapiens") %>% filter(gs_cat == "H" | (gs_cat == "C2" & gs_subcat == "CP:REACTOME") | (gs_cat == "C5" & gs_subcat %in% c("BP", "CC")))  %>% select(gs_name,human_gene_symbol)


enrichment.results <- tibble()

for (cluster.id in c(1, 2, 3,4,5,6, 7)) {
  geneIds.cluster <- filter(diff.proteins, cluster == cluster.id)$First_GeneID
  geneIds.cluster <- geneIds.cluster[which(!is.na(geneIds.cluster))]
  
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


enrichment.results$category <- ""
enrichment.results$category[str_detect(enrichment.results$Description, "^HALLMARK")] <- "Hallmark"
enrichment.results$category[str_detect(enrichment.results$Description, "^REACTOME")] <- "Reactome"
enrichment.results$category[str_detect(enrichment.results$Description, "^GO")] <- "Gene Ontology"

num.in.cluster <- data.frame(table(diff.proteins$cluster))
names(num.in.cluster) <- c("cluster", "n")

enrichment.results$ratio <- enrichment.results$Count / num.in.cluster$n[enrichment.results$cluster]

# convert to matrix format
enrichment.results_wide = pivot_wider(enrichment.results, id_cols = c(Description, category), names_from = cluster, values_from=c(pvalue, p.adjust, qvalue, Count, GeneRatio, BgRatio, ratio, genes))

createDir(here("tables"))
write_tsv(filter(enrichment.results_wide, category == "Hallmark"), here("Supplemental_Table_Cluster_Gene_Enrichment_Hallmark_sigs.tsv"), na = "")
write_tsv(filter(enrichment.results_wide, category == "Reactome"), here("Supplemental_Table_Cluster_Gene_Enrichment_Reactome_sigs.tsv"), na = "")
#write_tsv(filter(enrichment.results_wide, category == "Gene Ontology"), here("tables/enrichment_GO.tsv"), na = "")

include_sets <- c(
  
  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",                            
  "REACTOME_G1_S_SPECIFIC_TRANSCRIPTION",                                
  "HALLMARK_E2F_TARGETS",                                                
  #"REACTOME_NUCLEAR_EVENTS_KINASE_AND_TRANSCRIPTION_FACTOR_ACTIVATION_", 
  #"REACTOME_RAF_INDEPENDENT_MAPK1_3_ACTIVATION",                         
  #"REACTOME_G0_AND_EARLY_G1",                                            
  "REACTOME_MAPK_TARGETS_NUCLEAR_EVENTS_MEDIATED_BY_MAP_KINASES",        
  "REACTOME_ASPARAGINE_N_LINKED_GLYCOSYLATION",                          
  #"REACTOME_CARGO_CONCENTRATION_IN_THE_ER",                              
  "REACTOME_TRANSPORT_TO_THE_GOLGI_AND_SUBSEQUENT_MODIFICATION",         
  #"REACTOME_ER_TO_GOLGI_ANTEROGRADE_TRANSPORT,",                          
  #"REACTOME_COPII_MEDIATED_VESICLE_TRANSPORT",                           
  "HALLMARK_ANDROGEN_RESPONSE",                                          
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",                          
  #"REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",          
  #"REACTOME_NONSENSE_MEDIATED_DECAY_NMD_",                               
  #"REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",                          
  #"REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS",                
  #"REACTOME_SELENOAMINO_ACID_METABOLISM",                                
  "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
  #"REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE",
  
  
  #from cluster 6
  #"REACTOME_NGF_STIMULATED_TRANSCRIPTION"                                                
  "HALLMARK_KRAS_SIGNALING_UP",                                                          
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",                                                 
  #"REACTOME_NUCLEAR_EVENTS_KINASE_AND_TRANSCRIPTION_FACTOR_ACTIVATION_1"                 
  #"REACTOME_G1_S_SPECIFIC_TRANSCRIPTION1"                                                
  #"REACTOME_MITOTIC_G1_PHASE_AND_G1_S_TRANSITION"                                        
  #"HALLMARK_INFLAMMATORY_RESPONSE",                                                     
  #"REACTOME_REMOVAL_OF_THE_FLAP_INTERMEDIATE_FROM_THE_C_STRAND"                          
  "HALLMARK_TGF_BETA_SIGNALING",                                                      
  "REACTOME_SIGNALING_BY_NTRKS",                                                      
  "REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING",                                
  #"REACTOME_PROCESSIVE_SYNTHESIS_ON_THE_LAGGING_STRAND"                                  
  #"REACTOME_EXTRA_NUCLEAR_ESTROGEN_SIGNALING"                                            
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",                                      
  #"HALLMARK_UV_RESPONSE_UP"                                                              
  "REACTOME_BINDING_AND_UPTAKE_OF_LIGANDS_BY_SCAVENGER_RECEPTORS",                      
  #"REACTOME_PROCESSIVE_SYNTHESIS_ON_THE_C_STRAND_OF_THE_TELOMERE"                        
  #"REACTOME_TRANSCRIPTION_OF_E2F_TARGETS_UNDER_NEGATIVE_CONTROL_BY_DREAM_COMPLEX",    
  #"HALLMARK_ESTROGEN_RESPONSE_LATE"                                                      
  #"REACTOME_ESTROGEN_DEPENDENT_NUCLEAR_EVENTS_DOWNSTREAM_OF_ESR_MEMBRANE_SIGNALING"      
  #"REACTOME_LAGGING_STRAND_SYNTHESIS"                                                    
  #"REACTOME_PCNA_DEPENDENT_LONG_PATCH_BASE_EXCISION_REPAIR"                              
  "REACTOME_SIGNALING_BY_EGFR_IN_CANCER",                                                 
  #"REACTOME_UPTAKE_AND_ACTIONS_OF_BACTERIAL_TOXINS"                                      
  #"REACTOME_ATTENUATION_PHASE"                                                           
  #"REACTOME_SIGNALING_BY_ERBB2_IN_CANCER"                                                
  #"REACTOME_RESOLUTION_OF_AP_SITES_VIA_THE_MULTIPLE_NUCLEOTIDE_PATCH_REPLACEMENT_PATHWAY"
  #"REACTOME_DOWNREGULATION_OF_ERBB2_SIGNALING"                                           
  #"REACTOME_G0_AND_EARLY_G11"                                                            
  #"REACTOME_HSF1_ACTIVATION",                                                    
  "HALLMARK_P53_PATHWAY",                                                         
  #"REACTOME_S_PHASE"                                                                     
  "REACTOME_CELL_CYCLE_MITOTIC",                                                    
  #"REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES"                                      
  #"REACTOME_HSF1_DEPENDENT_TRANSACTIVATION"                                              
  #"REACTOME_DNA_STRAND_ELONGATION"                                                       
  #"REACTOME_TELOMERE_C_STRAND_LAGGING_STRAND_SYNTHESIS"                                  
  #"REACTOME_SIGNALING_BY_NUCLEAR_RECEPTORS"                                              
  "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS"                      
  #"REACTOME_TRANSPORT_OF_BILE_SALTS_AND_ORGANIC_ACIDS_METAL_IONS_AND_AMINE_COMPOUNDS"    
  #"REACTOME_RESOLUTION_OF_ABASIC_SITES_AP_SITES_"                                        
  #"HALLMARK_E2F_TARGETS1"                                                                
  #"REACTOME_BASE_EXCISION_REPAIR"                                                        
  #"REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP_"                             
  #"REACTOME_SIGNALING_BY_ERBB2"                                                          
  #"REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_CELL_CYCLE_GENES"                            
  #"REACTOME_DNA_REPLICATION"           
  
  
)

enricher_gene_sets <- subset(msig, msig$gs_name %in% include_sets)
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING")] <- "Interferon Alpha Beta Signaling"
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "REACTOME_G1_S_SPECIFIC_TRANSCRIPTION")] <- "G1/S Specific Transcription"
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "HALLMARK_E2F_TARGETS")] <- "E2F Targets"
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "REACTOME_MAPK_TARGETS_NUCLEAR_EVENTS_MEDIATED_BY_MAP_KINASES")] <- "MAPK Targets"
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "REACTOME_ASPARAGINE_N_LINKED_GLYCOSYLATION")] <- "Asparagine Linked Glycosylation"
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "REACTOME_TRANSPORT_TO_THE_GOLGI_AND_SUBSEQUENT_MODIFICATION")] <- "Transport to Golgi"
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "HALLMARK_ANDROGEN_RESPONSE")] <- "Androgen Response"    
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION")] <- "Translation Elongation"
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE")] <- "Protein Targeting to Membrane"    
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "HALLMARK_KRAS_SIGNALING_UP")] <- "KRAS Signaling"    
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB")] <- "TNFA Signaling via NFkB"    
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "HALLMARK_TGF_BETA_SIGNALING")] <- "TGFB Signaling"    
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "REACTOME_SIGNALING_BY_NTRKS")] <- "NTRK Signaling"

enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING")] <- "IL4 and IL13 Signaling"
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")] <- "Epithelial to Mesencymal Transition"
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "REACTOME_BINDING_AND_UPTAKE_OF_LIGANDS_BY_SCAVENGER_RECEPTORS")] <- "Scavenger Receptor Binding"
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "REACTOME_SIGNALING_BY_EGFR_IN_CANCER")] <- "EGFR Signaling"
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "HALLMARK_P53_PATHWAY")] <- "TP53 Pathway"
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "REACTOME_CELL_CYCLE_MITOTIC")] <- "Mitosis"
enricher_gene_sets$gs_name[which(enricher_gene_sets$gs_name == "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS")] <- "Replication Stress / ATR Activation"



save(enricher_gene_sets, file = "enricher_gene_sets.RData")


enrichment.results.pruned <- filter(enrichment.results_wide, Description %in% c(include_sets))#, subset(enrichment.results, enrichment.results$cluster %in% c(4,5,6,7))$ID)) ## enrichment.results.wide #
#enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_CELL_CYCLE_CHECKPOINTS")] <- "Cell cycle checkpoints"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING")] <- "Interferon Alpha Beta Signaling"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_G1_S_SPECIFIC_TRANSCRIPTION")] <- "G1/S Specific Transcription"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_E2F_TARGETS")] <- "E2F Targets"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_MAPK_TARGETS_NUCLEAR_EVENTS_MEDIATED_BY_MAP_KINASES")] <- "MAPK Targets"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_ASPARAGINE_N_LINKED_GLYCOSYLATION")] <- "Asparagine Linked Glycosylation"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_TRANSPORT_TO_THE_GOLGI_AND_SUBSEQUENT_MODIFICATION")] <- "Transport to Golgi"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_ANDROGEN_RESPONSE")] <- "Androgen Response"    
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION")] <- "Translation Elongation"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE")] <- "Protein Targeting to Membrane"    
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_KRAS_SIGNALING_UP")] <- "KRAS Signaling"    
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_TNFA_SIGNALING_VIA_NFKB")] <- "TNFA Signaling via NFkB"    
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_TGF_BETA_SIGNALING")] <- "TGFB Signaling"    
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_SIGNALING_BY_NTRKS")] <- "NTRK Signaling"

enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING")] <- "IL4 and IL13 Signaling"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")] <- "Epithelial to Mesencymal Transition"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_BINDING_AND_UPTAKE_OF_LIGANDS_BY_SCAVENGER_RECEPTORS")] <- "Scavenger Receptor Binding"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_SIGNALING_BY_EGFR_IN_CANCER")] <- "EGFR Signaling"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_P53_PATHWAY")] <- "TP53 Pathway"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_CELL_CYCLE_MITOTIC")] <- "Mitosis"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS")] <- "Replication Stress / ATR Activation"

enrichment.results.pruned <- arrange(enrichment.results.pruned, category)

em.q.values <- select(enrichment.results.pruned, grep("qvalue", colnames(enrichment.results.pruned), value=T))


em.q.values$qvalue_6 <- rep(NA, times = length(em.q.values$qvalue_2))
em.q.values$qvalue_1 <- rep(NA, times = length(em.q.values$qvalue_2))
#em.q.values$qvalue_7 <- rep(NA, times = length(em.q.values$qvalue_1))

em.q.values <- em.q.values[, c("qvalue_1", "qvalue_2", "qvalue_3", "qvalue_4", "qvalue_5", "qvalue_6", "qvalue_7")]

rownames(em.q.values) <- as.character(enrichment.results.pruned$Description)
colnames(em.q.values) <- c(1,2,3,4,5,6,7)

#rownames(em.q.values) <- as.character(enrichment.results.pruned$Description)
#blank_row <- rep(1, times = 7)
#names(blank_row) <- names(em.q.values)
#row.names(blank_row) <- c("blank")
#em.q.values <- rbind(em.q.values, blank_row)

ratio.values <- select(enrichment.results.pruned, grep("ratio", colnames(enrichment.results.pruned), value=T))
ratio.values$ratio_6 <- rep(NA, times = length(ratio.values$ratio_2))
#ratio.values$ratio_6 <- rep(NA, times = length(ratio.values$ratio_1))
ratio.values$ratio_1 <- rep(NA, times = length(ratio.values$ratio_2))

ratio.values <- ratio.values[, c("ratio_1", "ratio_2", "ratio_3", "ratio_4", "ratio_5", "ratio_6", "ratio_7")]

#ratio.values <- ratio.values[ ,c(1,2,3,5,4,6,7) ]
colnames(ratio.values) <- c(1,2,3,4,5,6,7)
ratio.values <- data.frame(ratio.values)



em.q.values <- -log10(as.matrix(em.q.values))
max.q = max(em.q.values, na.rm=T)
min.q = min(em.q.values, na.rm=T)
col_fun = colorRamp2(c(min.q, max.q-25, max.q), c("#300101", "#e41a1c", "#e41a1c"))

col.annot <- HeatmapAnnotation(cluster = anno_simple(colnames(em.q.values),
                                                     pt_gp = gpar(fontsize = 6),
                                                     height = unit(1, "mm"),
                                                     col = structure(brewer.pal(num.clusters, "Set3"), 
                                                                     names = colnames(em.q.values))),
                               show_annotation_name = F,
                               show_legend = F,
                               annotation_name_gp = gpar(fontsize = 6))


h <- Heatmap(em.q.values, 
             # labels
             column_title = "Enriched gene sets",
             column_title_gp = gpar(fontsize=7),
             row_names_gp = gpar(fontsize = 7),
             column_names_gp = gpar(fontsize = 6),
             column_names_rot = 0,
             row_title_gp = gpar(fontsize=7),
             #name = "-log10(q-value)",
             # legends
             bottom_annotation = col.annot,
             show_heatmap_legend = T,
             col = col_fun,
             heatmap_legend_param = list(color_bar = "continuous",
                                         title_gp = gpar(fontsize = 6),
                                         labels_gp = gpar(fontsize = 6),
                                         grid_width = unit(2,"mm")),
             
             cluster_rows = F,
             cluster_columns = F,
             show_row_names = T,
             show_column_names = T,
             #split = enrichment.results.pruned$category,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.rect(x = x, y = y, width = width, height = height, 
                         gp = gpar(fill = "white", col = "#EEEEEE"))
               grid.circle(x = x, y = y, r = sqrt(ratio.values[i, j]) * 0.1,
                           gp = gpar(fill = col_fun(em.q.values[i, j]), col = NA))
             }
)

gb_heatmap = grid.grabExpr(draw(h), height=5, width=2.5)
#  }






################################################################################
# Profile plots of clusters
################################################################################
#plotClusterProfiles <- function(diff.proteins.averaged.condition.fc.mock)
#{
load("limma_voom_deg_data.RData")
fc_data <- limma_voom_deg_data
fc_data <- subset(fc_data, !(fc_data$genes %in% sars2_genes$V1))


cluster_dictionary <- diff.proteins$cluster
names(cluster_dictionary) <- row.names(diff.proteins)

add_cluster <- c()
for(i in fc_data$genes){
  add_cluster <- c(add_cluster, cluster_dictionary[i])
}

fc_data$cluster <- add_cluster
fc_data_cluster <- fc_data[complete.cases(fc_data), ]
fc_data_cluster$Condition <- fc_data_cluster$time

fc_data_cluster_mock <- subset(fc_data_cluster, fc_data_cluster$time == 96)
fc_data_cluster_mock$time <- rep(0, times = length(fc_data_cluster_mock$time))
fc_data_cluster_mock$Condition <- fc_data_cluster_mock$time
fc_data_cluster_mock$logFC <- rep(0, times = length(fc_data_cluster_mock$time))



### add mock to plottind DF
fc_data_cluster <- rbind(fc_data_cluster, fc_data_cluster_mock)
fc_data_cluster$logFC <- fc_data_cluster$logFC * -1


num.in.cluster <- data.frame(table(diff.proteins$cluster))
names(num.in.cluster) <- c("cluster", "n")

cluster.labels <- paste("Cluster ", num.in.cluster$cluster, " (n=", num.in.cluster$n, ")", sep="")
names(cluster.labels) <- seq(1:7)

representative.profiles <- (fc_data_cluster %>% 
                              group_by(time, cluster) %>% 
                              summarise_at(c("logFC"),  mean))

#diff.proteins.averaged.condition.fc.noMock <- filter(diff.proteins.averaged.condition.fc.mock, Mock==F)
#diff.proteins.averaged.condition.fc.MockOnly <- filter(diff.proteins.averaged.condition.fc.mock, Mock==T)


profile_plots <- (ggplot(fc_data_cluster, aes(x=as.numeric(time), 
                                              y=logFC, 
                                              group=interaction(`genes`, cluster),
                                              color=as.factor(cluster),
                                              #linetype=as.factor(Mock),
                                              facet=cluster))
                  
                  #    + geom_line(size=0.25, data = diff.proteins.averaged.condition.fc.MockOnly, color="#CCCCCC")
                  + geom_line(size=0.5, alpha=0.5)
                  + geom_line(size=0.5, data=representative.profiles, 
                              aes(x=as.numeric(time), y=logFC, group=interaction(as.factor(cluster))),
                              color=grey)
                  + geom_point(size=0.5, data=representative.profiles, 
                               aes(x=as.numeric(time), y=logFC, group=NULL), 
                               color=grey)
                  
                  
                  + scale_y_continuous(name = "log2(fold-change over mock)", breaks=seq(-6,8,2))
                  + scale_x_continuous(name = "Hours post-infection", 
                                       breaks = c(4, 12, 24, 48, 72, 96),
                                       expand = c(0.05, 0))
                  
                  + facet_grid(cluster ~ ., scales="free", 
                               labeller = labeller(cluster = cluster.labels, size = 6))
                  
                  + scale_color_brewer(type="qual",palette="Set3")
                  
                  + theme.basic
                  + theme(legend.position = "none",
                          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                          axis.text=element_text(size=8))
)
#}


