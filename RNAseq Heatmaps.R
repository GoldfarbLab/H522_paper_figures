library(matrixStats)
library(methods)
library(cluster)
interferon.type1.rna <- read_csv(here("type_1_interferon.csv"))
normalized.immune.heatmap.data.rna <- read_csv(here("RNANormalizedImmuneHeatmapData.csv"))
rnadata <- read_tsv(here("uq_normalized_expr copy.txt"))

#Subset RNAseq DATA to have: WASHU_VERO_AGM, SRA_CACO-2, SRA_CALU-3, and UNC ones and matching Gene Names with Interferon Type I 
rnadata <- select(rnadata,"SYMBOL", "ENTREZID", "SRA_CACO2", "SRA_CALU3", starts_with("UNC"))
cell.lines.type1 <- filter(rnadata, rnadata$ENTREZID %in% interferon.type1.rna$GeneID)
cell.lines.type1 <- cell.lines.type1 %>% select(-ENTREZID) #%>% column_to_rownames(var="SYMBOL")

#filter out ones that don't change
variance <- rowVars(as.matrix(cell.lines.type1))
cell.lines.type1 <- cbind(cell.lines.type1, variance) %>% filter(variance >=.5)  #%>% select(-variance) 

#z-score cell lines
cell.lines.type1.zscored <- t(scale(t(as.matrix(cell.lines.type1)), center=T, scale=T))
cell.lines.type1.zscored[is.nan(cell.lines.type1.zscored)] <- NA
cell.lines.type1.zscored <- na.omit(cell.lines.type1.zscored)

#ANNOTATIONS select ones for annotations for ACE2, TMPRSS2, FURIN, CTSB, CTSL, NRP1 - unfinished
#annotations.cell.lines.type1 <- rnadata[c("ACE2", "TMPRSS2", "FURIN", "CTSB", "CTSL", "NRP1"),]
covid.annotations <- c("ACE2", "TMPRSS2", "FURIN", "CTSB", "CTSL", "NRP1")
annotations.cell.lines.type1<- filter(rnadata, SYMBOL %in% covid.annotations) %>%
  column_to_rownames(var="SYMBOL") %>% select(-ENTREZID)
#turn it into a heatmap, not z scored tho 
column.annotation <- HeatmapAnnotation(receptors = anno_simple(t(as.matrix(annotations.cell.lines.type1))),
                                       #name = column.annotation,
                                       which = 'col',
                                       col = list(receptors = colorRamp2(seq(-3, 3, length=256), 
                                                                         rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)))),
                                       annotation_width = unit(c(1, 4), 'cm'),
                                       gap = unit(1, 'mm'))

#Subset Normalized Heatmap Data to be just Type1 as well 
protein.data.type1 <- filter(normalized.immune.heatmap.data.rna, normalized.immune.heatmap.data.rna$First_GeneID %in% interferon.type1.rna$GeneID)
protein.data.type1 <- protein.data.type1 %>% column_to_rownames(var="GeneNames") %>% select(-First_GeneID) %>% 
  select("4 hr", "12 hr", "24 hr", "48 hr", "72 hr", "96 hr", "significant") 
#heatmap.data.type1[is.na(heatmap.data.type1)] <- 0 
protein.data.type1 <- as.matrix(protein.data.type1)

#join them 
joined.heatmap.data <- full_join(as_tibble(cell.lines.type1.zscored, rownames = "GeneNames"), 
                                 as_tibble(protein.data.type1, rownames = "GeneNames"), 
                                 by = "GeneNames") %>% column_to_rownames("GeneNames")

rnavariance <- as.integer(!is.na(joined.heatmap.data$variance)) 
joined.heatmap.data$rnavariance <- rnavariance
joined.heatmap.data <- select(joined.heatmap.data, -variance) %>% 
  filter(joined.heatmap.data$significant == 1 | joined.heatmap.data$rnavariance == 1) %>%
  select(-rnavariance, -significant)

#So now that they're the same size, split them back into sep tables to put in heatmap 
heatmap.cell.line.data <- select(joined.heatmap.data, "SRA_CACO2", "SRA_CALU3", starts_with("UNC")) %>% as.matrix(.)
heatmap.data.type1 <- select(joined.heatmap.data, "4 hr", "12 hr", "24 hr", "48 hr", "72 hr", "96 hr") %>% as.matrix(.)

notna.cell.lines <- heatmap.cell.line.data
  notna.cell.lines[is.na(notna.cell.lines)] <- 10

  #Cluster 
row.clusters <- hclust(dist(notna.cell.lines)) #row clustering 

#MAKE HEATAMP
#should be in a function but let's work outside one for a bit 
#plotRNAHeatmap <- function(heatmap.data, cell.lines.data)
#in this case, heatmap.data = heatmap.data.type1 and cell.lines.data = cell.lines.type1
#Make a heatmap of heatmap.data.type1 
rna.heatmap.of.proteins <- Heatmap(heatmap.data.type1,
                               # labels
                               column_title = paste("Type 1 Interferon Response Proteinsin H522\nn = ", nrow(heatmap.data.type1), sep=""),
                               column_title_gp = gpar(fontsize=7),
                               show_row_names = T,
                               row_names_gp = gpar(fontsize = 7),
                               row_names_side = "right",
                               column_names_gp = gpar(fontsize = 7),
                               name = "z-score",
                               col = colorRamp2(seq(-3, 3, length=256), rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))),
                               # legends
                               show_heatmap_legend = F,
                               heatmap_legend_param = list(color_bar = "continuous",
                                                           title_gp = gpar(fontsize = 6),
                                                           labels_gp = gpar(fontsize = 6),
                                                           grid_width = unit(2,"mm")),   
                               #REST UNEDITED 
                               # clustering
                               cluster_columns = F,
                               clustering_distance_rows = "euclidean",
                               cluster_rows = F,
                               show_row_dend = F,
                               #cluster_rows = row.clusters,
                               #split=7, #probs dont wanna split
                               #split = row.clusters,
                               # labels
                               row_title_rot = 0,
                               row_title_gp = gpar(fontsize=7),
                               
                               # size
                               width=ncol(heatmap.data.type1)*0.1)

#Make a heatmap of rnadata.type1 -- what should this look like? 
rna.cell.line.heatmap<- Heatmap(heatmap.cell.line.data,
                                   # labels
                                   column_title = paste("Type 1 Inferon Response Proteins\nn = ", nrow(cell.lines.type1), sep=""),
                                   column_title_gp = gpar(fontsize=7),
                                   show_row_names = T,
                                   row_names_gp = gpar(fontsize = 6),
                                   row_names_side = "right",
                                   column_names_gp = gpar(fontsize = 7),
                                   name = "z-score",
                                   #na_col = grey,
                                   #col = colorRamp2(seq(-3, 3, length=256), rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))),
                                   # legends
                                   show_heatmap_legend = F,
                                   heatmap_legend_param = list(color_bar = "continuous",
                                                               title_gp = gpar(fontsize = 6),
                                                               labels_gp = gpar(fontsize = 6),
                                                               grid_width = unit(2,"mm")),   
                                   #REST UNEDITED 
                                   # clustering
                                   cluster_columns = T,
                                   clustering_distance_rows = "euclidean",
                                   show_row_dend = F,
                                   cluster_rows = row.clusters,
                                   #split=7, #probs dont wanna split
                                   #split = row.clusters,
                                   # labels
                                   row_title_rot = 0,
                                   row_title_gp = gpar(fontsize=7),
                                   
                                   # size
                                   width=ncol(cell.lines.type1)*0.2,
                                   #annotations 
                                   #bottom_annotation_height=unit(1.0,"cm"), 
                                   bottom_annotation=column.annotation)
#add annotation row for ACE2, TMPRSS2, FURIN, CTSB, CTSL, NRP1 
rna.heatmap.list <-  rna.cell.line.heatmap + rna.heatmap.of.proteins 
#interferon_heatmap <- draw(heatmap.list)
type1_interferon_response_heatmap <- grid.grabExpr(draw(rna.heatmap.list), height=6.5, width=2)
