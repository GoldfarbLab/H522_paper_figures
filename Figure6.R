library(here)
library(msigdbr)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(annotate)
source(here("common.R"))
select <- get(x = "select", pos = "package:dplyr") # deal with function masking

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
                aspect.ratio = 0.75,
                plot.subtitle = element_text(size=7, hjust = 0.5, face = "italic"))
  )
}

################################################################################
# Heatmap
################################################################################
plotHeatmap <- function(proteins.z.scored, diff.proteins, design, num.clusters)
{
  colnames(proteins.z.scored) <- paste(design$Condition, design$Replicate)
  colnames(proteins.z.scored)[which(str_detect(colnames(proteins.z.scored), "1$"))] <- "1"
  colnames(proteins.z.scored)[which(str_detect(colnames(proteins.z.scored), "3$"))] <- "3"
  # reorder columns
  proteins.z.scored <- proteins.z.scored[, c(1,9,17, 2,10,18, 3,11,19, 4,12,20, 5,13,21, 6,14,22, 7,15,23, 8,16,24)]
  
  #label.sars2.pos <-  which(diff.proteins$SARS_CoV_2 == T)
  #label.int.pos <-  which(!is.na(diff.proteins$Bait))
  #label.pos <- c(label.sars2.pos, label.int.pos)
  label.pos <- which(diff.proteins$SARS_CoV_2 == T)
  
  #labels.sars2 <- diff.proteins$`Gene names`[label.sars2.pos]
  #labels.int <- paste(diff.proteins$`Gene names`[label.int.pos], " (", diff.proteins$Bait[label.int.pos], ")", sep="")
  #labels <- c(labels.sars2, labels.int)
  labels <- diff.proteins$`Gene names`[label.pos]
  
  heatmap.proteins <- Heatmap(as.matrix(proteins.z.scored),
                           # labels
                           column_title = paste("Differentially Expressed Proteins (n = ", nrow(proteins.z.scored), ")", sep=""),
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
                           split = diff.proteins$cluster,
                           # labels
                           row_title_rot = 0,
                           row_title_gp = gpar(fontsize=7),
                           #row_title = NULL,
                           # size
                           width=ncol(proteins.z.scored)*0.2
  )
  
  cluster.rowAnnot <- rowAnnotation(block = anno_block(gp = gpar(fill = brewer.pal(num.clusters, "Set3"), col = NA)),
                                    width = unit(1, "mm"))
  
  gene.rowAnnot = rowAnnotation(gene.name = anno_mark(at = label.pos, 
                                                      labels = labels,
                                                      labels_gp = gpar(fontsize = 5), 
                                                      link_gp = gpar(lwd = 0.5),
                                                      link_width = unit(4, "mm"),
                                                      padding = unit(0.5, "mm")))
  
  ht_list = heatmap.proteins + cluster.rowAnnot + gene.rowAnnot
  
  gb_heatmap = grid.grabExpr(draw(ht_list), height=3.5, width=2.5)
}

################################################################################
# Profile plots of clusters
################################################################################
plotClusterProfiles <- function(diff.proteins.averaged.condition.fc.mock, num.clusters)
{
  num.in.cluster <- (diff.proteins.averaged.condition.fc.mock %>%
                       group_by(cluster) %>% 
                       summarise(n=n() / (nlevels(Condition))))
  
  cluster.labels <- paste("Cluster ", num.in.cluster$cluster, " (n=", num.in.cluster$n, ")", sep="")
  names(cluster.labels) <- seq(1:num.clusters)
  
  representative.profiles <- (diff.proteins.averaged.condition.fc.mock %>% 
                                group_by(Time, cluster, Mock) %>% 
                                summarise_at(c("FC"),  mean))
  
  diff.proteins.averaged.condition.fc.noMock <- filter(diff.proteins.averaged.condition.fc.mock, Mock==F)
  diff.proteins.averaged.condition.fc.MockOnly <- filter(diff.proteins.averaged.condition.fc.mock, Mock==T)
  
  
  p <- (ggplot(diff.proteins.averaged.condition.fc.noMock, aes(x=as.numeric(Time), 
                                         y=FC, 
                                         group=interaction(`Gene names`, cluster, Mock),
                                         color=as.factor(cluster),
                                         linetype=as.factor(Mock),
                                         facet=cluster))
        
        #+ geom_line(size=0.25, data = diff.proteins.averaged.condition.fc.MockOnly, color="#CCCCCC")
        + geom_line(size=0.25, alpha=0.5)
        + geom_line(size=0.4, data=representative.profiles, 
                    aes(x=as.numeric(Time), y=FC, group=interaction(as.factor(cluster), Mock)),
                    color=grey)
        + geom_point(size=0.4, data=representative.profiles, 
                     aes(x=as.numeric(Time), y=FC, group=NULL), 
                     color=grey)
        
        
        + scale_y_continuous(name = "log2(fold-change over 4h mock)", breaks=seq(-5,5,1))
        + scale_x_continuous(name = "Hours post-infection", 
                             breaks = c(4, 12, 24, 48, 72, 96),
                             expand = c(0.05, 0))
        
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
plotCoVProfiles <- function(proteins, proteins.avg)
{
  SARS2.prots <- (proteins %>%
                    filter(SARS_CoV_2 == T & !str_detect(Condition, "Mock")) %>%
                    mutate(Condition = as.numeric(str_replace(Condition, " hr", ""))) %>%
                    mutate(Label = `Gene names`) %>%
                    mutate(FC.new = 2^(FC)))

  SARS2.prots.avg <- (proteins.avg %>%
                        filter(SARS_CoV_2 == T & !str_detect(Condition, "Mock")) %>%
                        mutate(Condition = as.numeric(str_replace(Condition, " hr", ""))) %>%
                        mutate(Label = `Gene names`))
  
  dfwc_between <- summarySE(data=SARS2.prots, measurevar="FC.new", groupvars=c("Condition","Label"), na.rm=FALSE, conf.interval=.95)
  SARS2.prots.avg <- SARS2.prots.avg %>% 
    left_join(dfwc_between, by = c("Label", "Condition"), suffix=c("","y")) %>%
    filter(!(Label %in% c("Nsp8")))
  
  SARS2.prots.avg$Label[which(SARS2.prots.avg$Condition != "96")] <- ""

  p <- (ggplot(SARS2.prots.avg, aes(x=Condition, y=FC, group=`Gene names`, label=Label))
        
        + geom_ribbon(aes(ymin = log2(FC.new-se), ymax = log2(FC.new+se)), fill = "black", alpha=0.075)
        #+ geom_errorbar(width=3, color=light.grey,
        #                aes(ymin=log2(FC.new-se),
        #                    ymax=log2(FC.new+se)))

        + geom_line(aes(color=`Gene names`), size=0.35) #medium.grey
        + geom_point(aes(color=`Gene names`), size=0.75) #medium.grey
        + geom_text(aes(color=`Gene names`, x=Condition+10), size=2)
        
        + ggtitle("SARS-CoV-2 Proteins")
        
        + scale_color_manual(values=colors.SARS.proteins)
        
        + scale_y_continuous(name = "log2(Intensity / 4hr Mock)", breaks = seq(0,6,1), limits=c(-0.4,6.3))

        + scale_x_continuous(name = "Hours post-infection", 
                             breaks = c(4, 12, 24, 48, 72, 96),
                             expand = c(0.05, 0))
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 0.75)
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
  
  proteins <- proteins %>% mutate(significance = case_when(SARS_CoV_2 == T ~ "SARS-CoV-2",
                                                           abs(logFC) >= proteomics.log.fc.threshold & log.adj.P.Val >= proteomics.log.stats.threshold ~ "sig both",
                                                           abs(logFC) >= proteomics.log.fc.threshold ~ "sig fold-change",
                                                           log.adj.P.Val >= proteomics.log.stats.threshold ~ "sig statistic",
                                                           TRUE ~ "none"), 
                                  Label = case_when(SARS_CoV_2 == T ~ str_replace(`Gene names`, "SARS_CoV_2", ""),
                                                    TRUE ~ ""))
  
  p <- (ggplot(data=proteins, aes(x=logFC, y=log.adj.P.Val, label=Label, color=significance))
        + geom_point(size=0.5) 
        
        + geom_vline(xintercept = -proteomics.log.fc.threshold, linetype="dashed", color=grey)
        + geom_vline(xintercept = proteomics.log.fc.threshold, linetype="dashed", color=grey)
        + geom_hline(yintercept = proteomics.log.stats.threshold, linetype="dashed", color=grey)
        
        + geom_text(nudge_y = 0.05 * y.range, nudge_x = 0.05 * x.range, size=2, color=COV2.color)
        
        + scale_color_manual(values=c("none"="#AAAAAA55", "sig fold-change"="#AAAAAA55", "sig statistic"="#AAAAAA55", "sig both"="#1f78b455", "SARS-CoV-2"="#fb8072"))
        
        + scale_x_continuous(limits=c(min.x - x.range/20, x.limit + x.range/20),
                             breaks = seq(-3,6,1))
        
        + ylab("-log10(q-value)")
        + xlab("log2(fold-change)")
        
        + ggtitle("96 hpi vs 96 hr mock")
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 0.75)
  )
}

################################################################################
# Enrichment heatmap
################################################################################
plotEnrichment <- function(proteins, diff.proteins, num.clusters)
{
  geneIds.universe <- as.character(proteins$First_GeneID[which(!is.na(proteins$First_GeneID))])
  
  msig <- msigdbr(species = "Homo sapiens") %>% filter(gs_cat == "H" | (gs_cat == "C2" & gs_subcat == "CP:REACTOME") | (gs_cat == "C5" & gs_subcat %in% c("BP", "CC")))  %>% select(gs_name, entrez_gene)
  enrichment.results <- tibble()
  
  for (cluster.id in 1:num.clusters) {
    geneIds.cluster <- filter(diff.proteins, cluster == cluster.id)$First_GeneID
    geneIds.cluster <- geneIds.cluster[which(!is.na(geneIds.cluster))]
    
    em <- enricher(gene=geneIds.cluster,
                   universe=geneIds.universe,
                   TERM2GENE=msig,
                   minGSSize = 10,
                   qvalueCutoff=0.05)
    
    if (!is.null(em)) {
      result <- em@result
      result$cluster <- cluster.id
      result <- filter(result, qvalue <= 0.05)
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
  
  num.in.cluster <- (diff.proteins %>%
                       group_by(cluster) %>% 
                       summarise(n=n()))
  
  enrichment.results$ratio <- enrichment.results$Count / num.in.cluster$n[enrichment.results$cluster]
  
  # convert to matrix format
  enrichment.results_wide = pivot_wider(enrichment.results, id_cols = c(Description, category), names_from = cluster, values_from=c(pvalue, p.adjust, qvalue, Count, GeneRatio, BgRatio, ratio, genes))
  
  createDir(here("tables"))
  write_tsv(filter(enrichment.results_wide, category == "Hallmark"), here("tables/enrichment_hallmark.tsv"), na = "")
  write_tsv(filter(enrichment.results_wide, category == "Reactome"), here("tables/enrichment_reactome.tsv"), na = "")
  write_tsv(filter(enrichment.results_wide, category == "Gene Ontology"), here("tables/enrichment_GO.tsv"), na = "")
  
  enrichment.results.pruned <- filter(enrichment.results_wide, Description %in% c("REACTOME_CELL_CYCLE_CHECKPOINTS", 
                                                                                  "REACTOME_SIGNALING_BY_RHO_GTPASES", 
                                                                                  "REACTOME_DNA_REPAIR", 
                                                                                  "REACTOME_DNA_REPLICATION", 
                                                                                  "REACTOME_REGULATION_OF_TP53_ACTIVITY", 
                                                                                  "REACTOME_SUMOYLATION",
                                                                                  "REACTOME_NUCLEAR_SIGNALING_BY_ERBB4",
                                                                                  "REACTOME_DNA_STRAND_ELONGATION",
                                                                                  
                                                                                  "GO_INTRINSIC_COMPONENT_OF_PLASMA_MEMBRANE", 
                                                                                  "GO_SPINDLE", 
                                                                                  "GO_ANAPHASE_PROMOTING_COMPLEX", 
                                                                                  "GO_UBIQUITIN_LIGASE_COMPLEX", 
                                                                                  "GO_REPLICATION_FORK", 
                                                                                  "GO_PHAGOCYTIC_VESICLE_MEMBRANE", 
                                                                                  "GO_ENDOCYTIC_VESICLE_MEMBRANE", 
                                                                                  "GO_ORGANELLE_ENVELOPE_LUMEN",
                                                                                  "GO_MICROTUBULE_CYTOSKELETON_ORGANIZATION",
                                                                                  "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_ANTIGEN",
                                                                                  "GO_INTRINSIC_COMPONENT_OF_ENDOPLASMIC_RETICULUM_MEMBRANE",
                                                                                  "GO_SEMAPHORIN_PLEXIN_SIGNALING_PATHWAY",
                                                                                  "GO_REGULATION_OF_B_CELL_PROLIFERATION",
                                                                                  "GO_RESPONSE_TO_WOUNDING",
                                                                                  
                                                                                  "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
                                                                                  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                                                                                  "HALLMARK_KRAS_SIGNALING_DN",
                                                                                  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                                                                  "HALLMARK_E2F_TARGETS",
                                                                                  "HALLMARK_G2M_CHECKPOINT"))
  
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_CELL_CYCLE_CHECKPOINTS")] <- "Cell cycle checkpoints"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_SIGNALING_BY_RHO_GTPASES")] <- "Signaling by Rho GTPases"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_DNA_REPAIR")] <- "DNA repair"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_DNA_REPLICATION")] <- "DNA replication"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_REGULATION_OF_TP53_ACTIVITY")] <- "Regulation of TP53 activity"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_SUMOYLATION")] <- "Sumoylation"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_NUCLEAR_SIGNALING_BY_ERBB4")] <- "Nuclear signaling by ERBB4"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_DNA_STRAND_ELONGATION")] <- "DNA strand elongation"
  
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "GO_INTRINSIC_COMPONENT_OF_PLASMA_MEMBRANE")] <- "Plasma membrane"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "GO_SPINDLE")] <- "Spindle"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "GO_ANAPHASE_PROMOTING_COMPLEX")] <- "Anaphase promoting complex"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "GO_UBIQUITIN_LIGASE_COMPLEX")] <- "Ubiquitin ligase complex"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "GO_REPLICATION_FORK")] <- "Replication Fork"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "GO_PHAGOCYTIC_VESICLE_MEMBRANE")] <- "Phagocytic vesicle membrane"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "GO_ENDOCYTIC_VESICLE_MEMBRANE")] <- "Endocytic vesicle membrane"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "GO_ORGANELLE_ENVELOPE_LUMEN")] <- "Organelle envelope lumen"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "GO_MICROTUBULE_CYTOSKELETON_ORGANIZATION")] <- "Microtubule cytoskeleton organization"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_ANTIGEN")] <- "Antigen processing and presentation"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "GO_INTRINSIC_COMPONENT_OF_ENDOPLASMIC_RETICULUM_MEMBRANE")] <- "Endoplasmic reticulum membrane"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "GO_SEMAPHORIN_PLEXIN_SIGNALING_PATHWAY")] <- "Semaphorin-plexin signaling pathway"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "GO_REGULATION_OF_B_CELL_PROLIFERATION")] <- "Regulation of B-cell proliferation"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "GO_RESPONSE_TO_WOUNDING")] <- "Response to wounding"
  
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_INTERFERON_ALPHA_RESPONSE")] <- "Interferon alpha response"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_INTERFERON_GAMMA_RESPONSE")] <- "Interferon gamma response"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_KRAS_SIGNALING_DN")] <- "Down-regulated by KRAS signaling"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")] <- "Epithelial-mesenchymal transition"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_E2F_TARGETS")] <- "E2F targets"
  enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_G2M_CHECKPOINT")] <- "G2M checkpoints"
  
  enrichment.results.pruned <- arrange(enrichment.results.pruned, category)
  em.q.values <- select(enrichment.results.pruned, grep("qvalue", colnames(enrichment.results.pruned), value=T))
  colnames(em.q.values) <- seq(1:num.clusters)
  ratio.values <- select(enrichment.results.pruned, grep("ratio", colnames(enrichment.results.pruned), value=T))
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
          row_names_gp = gpar(fontsize = 5),
          column_names_gp = gpar(fontsize = 6),
          column_names_rot = 0,
          row_title_gp = gpar(fontsize=7),
          row_labels = enrichment.results.pruned$Description,
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
}

plotCorrelation <- function(proteins.averaged.condition.fc.mock, limma_voom_deg_data, proteins)
{
  proteins.clust <- proteins %>% select(`Gene names`, cluster) %>% filter(!is.na(cluster)) %>% unique()
  # format proteomics data
  protein.abundance <- proteins.averaged.condition.fc.mock %>% 
    filter(Condition %in% c("4 hr", "24 hr", "48 hr", "72 hr", "96 hr")) %>%
    select(`Gene names`, Time, FC) %>%
    spread(Time, FC)
  # format proteomics data
  rna <- limma_voom_deg_data %>% 
    filter(time %in% c(4, 24, 48, 72, 96)) %>%
    select(`Gene names` = genes, time, logFC) %>%
    spread(time, logFC)
  
  protein.rna <- inner_join(protein.abundance, rna, by="Gene names", suffix=c(".protein", ".rna"))
  protein.rna <- protein.rna %>% left_join(proteins.clust, by="Gene names")
  
  for (i in 1:nrow(protein.rna)){
    protein.rna$pcc[i] <- cor(t(protein.rna[i, 2:6]), t(protein.rna[i, 7:11]))
  }
  
  hist.all <- hist(protein.rna$pcc, breaks=seq(-1,1,0.1), plot=F)
  hist.diff <- hist((protein.rna %>% filter(!is.na(cluster)))$pcc, breaks=seq(-1,1,0.1), plot=F)
  
  hist.data <- tibble("count" = c(hist.all$counts, hist.diff$counts), 
                      "pcc" = c(hist.all$mids, hist.diff$mids), 
                      "class"=c(rep("All proteins", length(hist.all$mids)), 
                                rep("Differentially expressed proteins", length(hist.diff$mids))
                                )
  )
  
  p <- (ggplot(hist.data, aes(x=pcc, y=count))
        + geom_col(fill="#EEEEEE", color=grey, width=0.1, size=0.5)
        
        + facet_wrap(~ class, scales="free_y", ncol=1)
        
        #+ scale_x_continuous(limits=c(min.x - x.range/20, x.limit + x.range/20),
        #                     breaks = seq(-3,6,1))
        
        + ylab("Count")
        + xlab("Pearson's correlation coefficient")
        
        + ggtitle("RNA vs protein correlations")
        
        + theme.basic
        + theme(aspect.ratio = 0.5,
                plot.title = element_text(size = 7,
                                          margin = margin(0.1, 0, 0.1, 0, "cm")),
                strip.text.x = element_text(size = 6,
                                            margin = margin(0.1, 0, 0.1, 0, "cm")))
  )
}

plotCorrelationEnrichment <- function(proteins.averaged.condition.fc.mock, limma_voom_deg_data, proteins)
{
  proteins.clust <- proteins %>% select(`Gene names`, cluster, First_GeneID) %>% filter(!is.na(cluster)) %>% unique()
  # format proteomics data
  protein.abundance <- proteins.averaged.condition.fc.mock %>% 
    filter(Condition %in% c("4 hr", "24 hr", "48 hr", "72 hr", "96 hr")) %>%
    select(`Gene names`, Time, FC) %>%
    spread(Time, FC)
  # format proteomics data
  rna <- limma_voom_deg_data %>% 
    filter(time %in% c(4, 24, 48, 72, 96)) %>%
    select(`Gene names` = genes, time, logFC) %>%
    spread(time, logFC)
  
  protein.rna <- inner_join(protein.abundance, rna, by="Gene names", suffix=c(".protein", ".rna"))
   
  protein.rna <- protein.rna %>% left_join(proteins.clust, by="Gene names") %>%
    filter(!is.na(cluster))
  
  for (i in 1:nrow(protein.rna)){
    protein.rna$pcc[i] <- cor(t(protein.rna[i, 2:6]), t(protein.rna[i, 7:11]))
  }
  
  protein.rna <- arrange(protein.rna, -pcc)
  
  geneList <- protein.rna$pcc
  names(geneList) <- as.character(protein.rna$First_GeneID)
  
  msig <- msigdbr(species = "Homo sapiens") %>% filter(gs_cat == "H" | (gs_cat == "C2" & gs_subcat == "CP:REACTOME") | (gs_cat == "C5" & gs_subcat %in% c("BP", "CC")))  %>% select(gs_name, entrez_gene)
  
  enriched <- GSEA(geneList, TERM2GENE=msig, verbose=FALSE, pvalueCutoff = 0.5)
  
  p1 <- gseaplot(enriched, geneSetID = 1, by = "runningScore", title = enriched$Description[1], base_size=6)
  p2 <- gseaplot(enriched, geneSetID = 2, by = "runningScore", title = enriched$Description[2], base_size=6)
  
  arranged.enrich <- arrangeGrob(p1, p2,
                             nrow = 2,
                             ncol = 1)
}

################################################################################



################################################################################
# Read data
################################################################################
design <- read_csv(here("data/Figure 6/MS/Experimental Design H522 Paper.csv")) %>% filter(Condition !="Bridge") # remove Ref Channel
normalized.data <- read_tsv(here("data_processed/proteinsQuantNormalizedToBridge.txt"))
proteins <- read_tsv(here("data_processed/proteins.txt"), guess_max = 10000)
proteins.z.scored <- read_tsv(here("data_processed/proteinsZscored.txt"))
diff.proteins <- read_tsv(here("data_processed/diffProteins.txt"))
proteins.averaged.condition.fc.mock <- read_tsv(here("data_processed/proteinsFCvsMockAveraged.txt"))
diff.proteins.averaged.condition.fc.mock <- read_tsv(here("data_processed/diffProteinsFCvsMockAveraged.txt"))
proteins.fc.mock <- read_tsv(here("data_processed/proteinsFCvsMock.txt"), guess_max = 70000)
#load(here("data/Figure 6/limma_voom_deg_data.RData"))
limma_voom_deg_data <- read_csv(here("data/Figure 6/RNA_DEG_table_1_1_26.csv"))
################################################################################



################################################################################
# Generate figures
################################################################################
panel.PCA <- plotPCA(normalized.data, design)
panel.heatmap <- plotHeatmap(proteins.z.scored, diff.proteins, design, proteomics.num.clusters)
panel.clusterProfiles <- plotClusterProfiles(diff.proteins.averaged.condition.fc.mock, proteomics.num.clusters)
panel.COV2Profiles <- plotCoVProfiles(proteins.fc.mock, proteins.averaged.condition.fc.mock)
panel.volcano <- plotVolcano(proteins)
panel.enrichment <- plotEnrichment(proteins, diff.proteins, proteomics.num.clusters)
panel.correlation <- plotCorrelation(proteins.averaged.condition.fc.mock, limma_voom_deg_data, proteins)
panel.correlation.enrichment <- plotCorrelationEnrichment(proteins.averaged.condition.fc.mock, limma_voom_deg_data, proteins)

arranged.QC <- arrangeGrob(panel.PCA, panel.COV2Profiles, panel.volcano,
                  nrow = 1,
                  ncol = 3)

#grid.draw(F5.top)  # to view the plot
saveFig(arranged.QC, "Figure6_BCD", 9, 5.7)
saveFig(panel.heatmap, "Figure6_E", 3.5, 2.3)
saveFig(panel.clusterProfiles, "Figure6_F", 6, 1.3)
saveFig(panel.enrichment, "Figure6_G", 3.1, 2.5)
saveFig(panel.correlation, "Figure6_H", 2.1, 2)
saveFig(panel.correlation.enrichment, "Figure6_H2", 6.3, 5)
################################################################################
