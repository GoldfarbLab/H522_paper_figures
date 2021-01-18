library(here)
library(limma)
library(ConsensusClusterPlus)
library(org.Hs.eg.db)
library(annotate)
library(biogridr) # devtools::install_github("npjc/biogridr")
library(CellSurfaceTools) # install.packages('devtools') # devtools::install_github("GoldfarbLab/CellSurfaceTools")
library(tidygraph) #devtools::install_github('thomasp85/tidygraph')

source(here("common.R"))
select <- get(x = "select", pos = "package:dplyr") # deal with function masking

#source(here("RequantifyProteomics.R")) # only necessary to run once


################################################################################
# Read data
################################################################################
data <- read_tsv(here("data_processed/requantifiedProteins.txt"), guess_max=10000)
design <- read_csv(here("data/Figure 6/MS/Experimental Design H522 Paper.csv"))

BioGRID.SARS2 <- read_tsv("annotations/BIOGRID-PROJECT-covid19_coronavirus_project-GENES-4.2.193.projectindex.txt", na="-") %>%
  filter(`VIRUS VALUES` == 'SARS-CoV-2')

BioGRID.SARS2.interactions <- read_tsv("annotations/BIOGRID-PROJECT-covid19_coronavirus_project-INTERACTIONS-4.2.193.tab3.txt", na="-", guess_max=100000) %>%
  mutate("minGene" = pmin(`Entrez Gene Interactor A`, `Entrez Gene Interactor B`), 
         "maxGene" = pmax(`Entrez Gene Interactor A`, `Entrez Gene Interactor B`),
         "humanGene" = ifelse(`Organism ID Interactor A` == 9606, `Entrez Gene Interactor A`, `Entrez Gene Interactor B`),
         "SARS2Gene" = ifelse(`Organism ID Interactor A` == 9606, `Entrez Gene Interactor B`, `Entrez Gene Interactor A`),
         "SARS2 Interactor" = ifelse(`Organism ID Interactor A` == 9606, `Official Symbol Interactor B`, `Official Symbol Interactor A`)) %>%
  filter(!is.na(minGene)) %>%
  filter(minGene != maxGene) %>%
  filter(SARS2Gene %in% BioGRID.SARS2$`ENTREZ GENE ID`) %>%
  filter(`Experimental System Type` == "physical") %>%
  group_by(humanGene, SARS2Gene, `SARS2 Interactor`) %>%
  summarize(num_unique_systems = n_distinct(`Experimental System`), 
            num_unique_publications = n_distinct(`Publication Source`), 
            is_mv = (n_distinct(`Experimental System`) > 1 | n_distinct(`Publication Source`) > 1)) %>%
  filter(is_mv == T) %>%
  select(humanGene, SARS2Gene, `SARS2 Interactor`) %>%
  mutate(`SARS2 Interactor` = str_to_title(`SARS2 Interactor`))

BioGRID.SARS2.interactions.per.human.gene <- BioGRID.SARS2.interactions %>%
  group_by(humanGene) %>%
  mutate("SARS2 Interactor" = paste0(`SARS2 Interactor`, collapse = ";")) %>%
  unique()

SARS.interactors.krogan <- select(read_csv(here("annotations/SARS2_interactome.csv")), c("Bait.Krogan", "PreyGeneName"))

SARS.interactors.mann <- read_csv("annotations/Mann_Interactors_Caco2.csv") %>% 
  filter(str_detect(bait_full_name, "SARS_CoV2")) %>% 
  select(c("Bait.Mann", "gene_name")) %>% 
  group_by(gene_name) %>% 
  summarise(Bait.Mann = str_c(base::unique(Bait.Mann), collapse = ";"))

SARS.interactors.liang <- read_tsv(here("annotations/Liang_interactors.txt")) %>% 
  mutate(Bait.Liang = Bait) %>% 
  select("Bait.Liang", "Gene name") %>%
  group_by(`Gene name`) %>% 
  summarise(Bait.Liang = str_c(base::unique(Bait.Liang), collapse = ";"))


H522.mutations <- read_tsv(here("annotations/H522_mutations.tsv"))
uniprot.mapping <- read_tsv(here("annotations/uniprot_mapping.tsv.zip"))
corum <- read_tsv(here("annotations/allComplexes.txt"))
biogrid <- read_tsv(here("annotations/BIOGRID-MV-Physical-HUMAN-4.2.193.tab3.txt")) %>%
  filter(`Entrez Gene Interactor A` != `Entrez Gene Interactor B`)
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

#valid.data <- filter(filtered.data, (rep1.valid + rep2.valid + rep3.valid) >= 2)
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
proteins <- cbind(select(valid.data, "Gene names", "Protein IDs", "Leading razor protein"), normalized.data)
# remove isoform suffix
proteins <- separate(proteins, "Leading razor protein", c("Leading canonical"), sep="-", remove=F, extra="drop")
# fill Entrez GeneID
proteins <- left_join(proteins, uniprot.mapping, by=c("Leading canonical" = "UniProt"))
# Fill missing GeneIDs for diff proteins: HLA-A, HLA-B, BARGIN, SLC5A3
proteins$First_GeneID[which(proteins$`Gene names` == "HLA-A")] = 3105
proteins$First_GeneID[which(proteins$`Gene names` == "HLA-B")] = 3106
proteins$First_GeneID[which(proteins$`Gene names` == "SLC5A3")] = 6526


# is it a SARS_CoV_2 protein?
proteins$"SARS_CoV_2" <- ifelse(str_detect(proteins$`Protein IDs`, "SARS_CoV_2"), T, F)
# is it an interactor of a SARS_CoV_2 protein?

proteins <- left_join(proteins, SARS.interactors.krogan, by=c("Gene names" = "PreyGeneName"))
proteins <- left_join(proteins, SARS.interactors.mann, by=c("Gene names" = "gene_name"))
proteins <- left_join(proteins, SARS.interactors.liang, by=c("Gene names" = "Gene name"))
proteins <- left_join(proteins, BioGRID.SARS2.interactions.per.human.gene, by=c("First_GeneID" = "humanGene"))
proteins$Bait <- str_c(str_replace_na(proteins$Bait.Krogan, replacement = ""), 
                       str_replace_na(proteins$Bait.Mann, replacement = ""),
                       str_replace_na(proteins$Bait.Liang, replacement = ""),
                       sep = ";")
proteins$Bait <- apply(proteins, 1, function(x) {
  baits <- unlist(str_split(x["Bait"], ";"))
  baits <- baits[baits != ""]
  str_c(base::unique(baits), collapse = ";")})
proteins <- proteins %>% mutate_all(na_if,"") #replace "" with NAs in Bait Column 

# interferon response 
interferon.response.alpha <- read_csv(here("annotations/interferon_response_alpha.txt")) %>% mutate(Interferon.Alpha = TRUE)
interferon.response.beta <- read_csv(here("annotations/interferon_response_beta.txt"))  %>% mutate(Interferon.Beta = TRUE)
interferon.response.gamma <-  read_csv(here("annotations/interferon_response_gamma.txt"))  %>% mutate(Interferon.Gamma = TRUE)
interferon.regulation.type1 <-  read_csv(here("annotations/regulation_of_type_I_interferon_mediated_signaling_pathway.txt"))  %>% mutate(Interferon.TypeI = TRUE)
interferon.regulation.type2.immune.response <-  read_csv(here("annotations/regulation_type_2_immune_response.txt"))  %>% mutate(Interferon.TypeII = TRUE)
antigen.processing.presentation <- read_csv(here("annotations/antigen_processing_and_presentation.txt"))  %>% mutate(Interferon.Antigen.Processing = TRUE)

proteins <- left_join(proteins, interferon.response.alpha, by=c("Gene names" = "GO_CELLULAR_RESPONSE_TO_INTERFERON_ALPHA"))
proteins <- left_join(proteins, interferon.response.beta, by=c("Gene names" = "GO_CELLULAR_RESPONSE_TO_INTERFERON_BETA"))
proteins <- left_join(proteins, interferon.response.gamma, by=c("Gene names" = "GO_RESPONSE_TO_INTERFERON_GAMMA"))
proteins <- left_join(proteins, interferon.regulation.type1, by=c("Gene names" = "GO_REGULATION_OF_TYPE_I_INTERFERON_MEDIATED_SIGNALING_PATHWAY"))
proteins <- left_join(proteins, interferon.regulation.type2.immune.response, by=c("Gene names" = "GO_REGULATION_OF_TYPE_2_IMMUNE_RESPONSE"))
proteins <- left_join(proteins, antigen.processing.presentation, by=c("Gene names" = "GO_ANTIGEN_PROCESSING_AND_PRESENTATION"))

proteins$Inteferon_Response <- proteins$Interferon.Alpha == TRUE | proteins$Interferon.Beta == TRUE | proteins$Interferon.Gamma == TRUE | proteins$Interferon.TypeI == TRUE | proteins$Interferon.TypeII == TRUE | proteins$Interferon.Antigen.Processing == TRUE

# is it mutated?
proteins <- left_join(proteins, H522.mutations, by=c("Gene names" = "GeneName"))
# is it a cell surface or plasma membrane protein?
proteins <- left_join(proteins, human_surface_and_plasma_membrane_protLevel, by=c("Leading canonical" = "UniProt"))
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

diff.proteins <- proteins %>% filter(abs(`logFC`) >= proteomics.log.fc.threshold & adj.P.Val <= proteomics.stats.threshold)

# calc enrichment of interactors
num.proteins <- sum(proteins$SARS_CoV_2 == F)
num.interactors <- sum(!is.na(proteins$`SARS2 Interactor`))
num.diff <- sum(diff.proteins$SARS_CoV_2 == F)
num.diff.interactors <- sum(diff.proteins$SARS_CoV_2 == F & !is.na(diff.proteins$`SARS2 Interactor`))
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
                                          maxK=8,
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

row.clusters <- clustering.results[[proteomics.num.clusters]][["consensusClass"]]

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
row.clusters[pos.2] <- 2
row.clusters[pos.4] <- 3
row.clusters[pos.1] <- 4
row.clusters[pos.7] <- 5
row.clusters[pos.3] <- 6
row.clusters[pos.5] <- 7
#row.clusters[pos.7] <- 8

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
                                          mutate(Mock = str_detect(Condition, "Mock")) %>%
                                          mutate(Time = str_extract(Condition, "\\d+")))

diff.proteins.averaged.condition <- (diff.proteins %>%
                                       pivot_longer(all_of(quant.colnames), names_to="Condition", values_to="FC") %>%
                                       mutate(Condition = str_extract(Condition, ".*(?=( Rep ))")) %>%
                                       group_by(`Gene names`, Condition) %>%
                                       summarise_at(c("FC","cluster","SARS_CoV_2"), mean, na.rm=T))

diff.proteins.averaged.condition.fc.mock <- (diff.proteins.averaged.condition %>%
                                               left_join(filter(., Condition == "Mock - 4 hr"), by="Gene names", suffix=c("",".mock")) %>%
                                               mutate(FC = FC - FC.mock) %>%
                                               mutate(Condition = factor(Condition, levels=c("Mock - 4 hr", "4 hr", "12 hr", "24 hr", "48 hr", "72 hr", "96 hr", "Mock - 96 hr"))) %>%
                                               mutate(Mock = str_detect(Condition, "Mock")) %>%
                                               mutate(Time = str_extract(Condition, "\\d+")))

proteins.fc.mock <- (proteins %>%
                       pivot_longer(all_of(quant.colnames), names_to="Condition", values_to="FC") %>%
                       mutate(Replicate = str_extract(Condition, "\\d+$")) %>%
                       mutate(Condition = str_extract(Condition, ".*(?=( Rep ))")) %>%
                       left_join(filter(., Condition == "Mock - 4 hr"), by=c("Gene names", "Replicate"), suffix=c("",".mock")) %>%
                       mutate(FC = FC - FC.mock) %>%
                       mutate(Condition = factor(Condition, levels=c("Mock - 4 hr", "4 hr", "12 hr", "24 hr", "48 hr", "72 hr", "96 hr", "Mock - 96 hr"))) %>%
                       mutate(Mock = str_detect(Condition, "Mock")) %>%
                       mutate(Time = str_extract(Condition, "\\d+")))
################################################################################



################################################################################
# Write processed data
################################################################################
write_tsv(normalized.data, here("data_processed/proteinsQuantNormalizedToBridge.txt"))
write_tsv(proteins, here("data_processed/proteins.txt"))
write_tsv(as_tibble(proteins.z.scored), here("data_processed/proteinsZscored.txt"))
write_tsv(diff.proteins, here("data_processed/diffProteins.txt"))
write_tsv(proteins.averaged.condition.fc.mock, here("data_processed/proteinsFCvsMockAveraged.txt"))
write_tsv(diff.proteins.averaged.condition.fc.mock, here("data_processed/diffProteinsFCvsMockAveraged.txt"))
write_tsv(proteins.fc.mock, here("data_processed/proteinsFCvsMock.txt"))
################################################################################



################################################################################
# Cytoscape output
################################################################################

createDir(here("tables/networks"))

# integrate CORUM
corum.matches <- corum %>% 
  select("ComplexName", "subunits(Entrez IDs)") %>%
  separate_rows("subunits(Entrez IDs)", sep=";") %>%
  filter(`subunits(Entrez IDs)` != "None") %>%
  mutate("subunits(Entrez IDs)" = as.numeric(`subunits(Entrez IDs)`)) %>%
  distinct() %>%
  group_by(ComplexName) %>%  mutate(num_subunits=n()) %>% 
  left_join(diff.proteins, by = c("subunits(Entrez IDs)" = "First_GeneID")) %>%
  group_by(ComplexName) %>%
  mutate(n=sum(!is.na(logFC))) %>% mutate(ratio = n/num_subunits) %>%
  filter((n > 1 & ratio >= 0.51 & num_subunits > 2)) %>%
  select("ComplexName", "n", "num_subunits", "ratio", "cluster", 
         "subunits(Entrez IDs)", "Gene names", "SARS2 Interactor", 
         "is_unique_mutation", "Num cell surface evidence") %>%
  mutate("Gene names" = getSYMBOL(as.character(`subunits(Entrez IDs)`), data='org.Hs.eg')) %>%
  distinct()

corum.not.matched <- corum.matches %>% filter(is.na(cluster)) %>%
  left_join(select(proteins, "First_GeneID", "is_unique_mutation"), by = c("subunits(Entrez IDs)" = "First_GeneID")) %>%
  mutate(is_unique_mutation = is_unique_mutation.y) %>% select(-c(is_unique_mutation.x,is_unique_mutation.y)) %>%
  distinct()

# Global node file
CoV2.proteins <- c("E","M","N","S",
                   str_c("Nsp", 1:16),
                   str_c("Orf", c("3a", "3b", "6", "7a", "7b", "8", "9b", "10")))
CoV2.proteins.not.diff <- CoV2.proteins[!CoV2.proteins %in% diff.proteins$`Gene names`]

nodes <- select(diff.proteins, "Gene names", "First_GeneID", "SARS_CoV_2", "cluster", "is_unique_mutation", "Num cell surface evidence") %>% 
  mutate(quantified=T)
nodes <- rbind(nodes, tibble("Gene names" = CoV2.proteins.not.diff, "First_GeneID"=NA, "SARS_CoV_2"=T, 
                             "cluster"=NA, "is_unique_mutation"=NA, "Num cell surface evidence"=NA, "quantified"=CoV2.proteins.not.diff %in% proteins$`Gene names`))
nodes <- rbind(nodes, tibble("Gene names" = corum.not.matched$"Gene names", "First_GeneID"=corum.not.matched$"subunits(Entrez IDs)", 
                             "SARS_CoV_2"=F, "cluster"=NA, "is_unique_mutation"=corum.not.matched$"is_unique_mutation",
                             "Num cell surface evidence"=NA, "quantified" = corum.not.matched$"subunits(Entrez IDs)" %in% proteins$First_GeneID) %>% distinct())
write_tsv(nodes, here("tables/networks/nodes.tsv"), na="")

bg_get_key("Dennis","Goldfarb","d.goldfarb@wustl.edu","interactome")



# edge file for each CORUM complex
for (complex in unique(corum.matches$ComplexName)) {
  complex.matches <- filter(corum.matches, ComplexName == complex)
  
  # make a clique for the CORUM complex
  complex.interactions <- tibble()
  for (interactor_a in unique(complex.matches$`Gene names`)) {
    for (interactor_b in unique(complex.matches$`Gene names`)) {
      if (interactor_a < interactor_b) {
        complex.interactions <- rbind(complex.interactions, 
                                      tibble(official_symbol_for_interactor_a=interactor_a,
                                             official_symbol_for_interactor_b=interactor_b,
                                             n=1))
      }
    }
  }
  
  
  
  interactions <- complex.interactions
  interactions <- interactions %>% mutate(source = official_symbol_for_interactor_a,
                                          target = official_symbol_for_interactor_b) %>%
    select(-official_symbol_for_interactor_a, -official_symbol_for_interactor_b)
  
  # Add SARS-CoV-2 interactions
  for (bait in CoV2.proteins) {
    num_evidence <- rowSums(tibble(B = str_detect(regex(complex.matches$`SARS2 Interactor`, ignore_case = T), bait)),
                            na.rm=T)
    int.pos <- which(num_evidence > 0)
    
    interactions <- rbind(interactions,
                          tibble(source=bait,
                                 target=complex.matches$`Gene names`[int.pos],
                                 n=0))#num_evidence[int.pos]))
  }
  # Edge file
  write_tsv(interactions, here(str_c("tables/networks/edges_", str_replace(complex, "/", "_"), ".tsv")), na="")
}


# get BioGRID interactions between diff proteins and save as edges tibble for tidygraph
graph_edges <- biogrid %>%
  filter(`Entrez Gene Interactor A` %in% diff.proteins$First_GeneID & `Entrez Gene Interactor B` %in% diff.proteins$First_GeneID) %>%
  left_join(diff.proteins %>% select(First_GeneID, `Gene names`), by = c(`Entrez Gene Interactor A` = "First_GeneID")) %>%
  left_join(diff.proteins %>% select(First_GeneID, `Gene names`), by = c(`Entrez Gene Interactor B` = "First_GeneID")) %>%
  select(from = `Gene names.x`, to = `Gene names.y`) %>%
  bind_rows(BioGRID.SARS2.interactions %>% 
              filter(humanGene %in% diff.proteins$First_GeneID) %>%
              left_join(diff.proteins %>% select(First_GeneID, `Gene names`), by = c("humanGene" = "First_GeneID")) %>%
              ungroup() %>%
              mutate(`from` = `Gene names`, `to` = `SARS2 Interactor`) %>%
              select(`from`, `to`)) %>%
  unique() %>%
  add_column(type = "BioGRID")
  









# get BioGRID interactions with the complex members
# geneList = c(as.vector(diff.proteins$`First_GeneID`[!is.na(diff.proteins$`First_GeneID`)]),
#              complex.matches$`subunits(Entrez IDs)`)
# interactions <- bg("interactions") %>%
#   bg_constrain(taxId = 9606, 
#                geneList = str_c(geneList, collapse="|"),
#                searchIds = T,
#                includeInteractors = F) %>%
#   bg_get_results()
# 
# # Check if any diff proteins interact with >50% of the complex
# interactions <- interactions %>%
#   filter(entrez_gene_id_for_interactor_a != "-" & entrez_gene_id_for_interactor_b != "-") %>%
#   filter(entrez_gene_id_for_interactor_a != entrez_gene_id_for_interactor_b) %>%
#   filter(entrez_gene_id_for_interactor_a %in% complex.matches$`subunits(Entrez IDs)` | entrez_gene_id_for_interactor_b %in% complex.matches$`subunits(Entrez IDs)`) %>%
#   mutate(official_symbol_for_interactor_a = getSYMBOL(as.character(entrez_gene_id_for_interactor_a), data='org.Hs.eg'),
#          official_symbol_for_interactor_b = getSYMBOL(as.character(entrez_gene_id_for_interactor_b), data='org.Hs.eg')) %>%
#   mutate(interactor_min = pmin(official_symbol_for_interactor_a, official_symbol_for_interactor_b),
#          interactor_max = pmax(official_symbol_for_interactor_a, official_symbol_for_interactor_b)) %>%
#   mutate(official_symbol_for_interactor_a = interactor_min, 
#          official_symbol_for_interactor_b = interactor_max) %>%
#   group_by(official_symbol_for_interactor_a, official_symbol_for_interactor_b) %>%
#   summarise(n = n()) %>% filter(n > 1)
# 
# interactions.A <- interactions %>% 
#   group_by(official_symbol_for_interactor_a) %>% 
#   filter(!(official_symbol_for_interactor_a %in% complex.matches$`Gene names`)) %>%
#   summarise(c.n = n())
# interactions.B <- interactions %>% 
#   group_by(official_symbol_for_interactor_b) %>% 
#   filter(!(official_symbol_for_interactor_b %in% complex.matches$`Gene names`)) %>%
#   summarise(c.n = n())
# interactions.tot <- full_join(interactions.A, interactions.B, by=c("official_symbol_for_interactor_a" = "official_symbol_for_interactor_b")) %>%
#   mutate(c.n = rowSums(select(., c.n.x, c.n.y), na.rm=T)) %>%
#   filter(c.n > 0.67*complex.matches$num_subunits[1])
# 
# if (nrow(interactions.tot) > 0) {
#   interactions.with.complex <- interactions %>%
#     filter(official_symbol_for_interactor_a %in% interactions.tot$official_symbol_for_interactor_a | official_symbol_for_interactor_b %in% interactions.tot$official_symbol_for_interactor_a)
# 
#   interactions <- rbind(complex.interactions, interactions.with.complex)
#   
#   if (nrow(interactions.tot) > 1) {
#     interactions.between <- bg("interactions") %>%
#       bg_constrain(taxId = 9606, 
#                    geneList = str_c(interactions.tot$official_symbol_for_interactor_a, collapse="|"),
#                    searchIds = F,
#                    includeInteractors = F) %>%
#       bg_get_results() 
#     
#     if (nrow(interactions.between) > 0) {
#       interactions.between <- interactions.between %>% 
#         filter(entrez_gene_id_for_interactor_a != entrez_gene_id_for_interactor_b) %>%
#         mutate(interactor_min = pmin(official_symbol_for_interactor_a, official_symbol_for_interactor_b),
#                interactor_max = pmax(official_symbol_for_interactor_a, official_symbol_for_interactor_b)) %>%
#         mutate(official_symbol_for_interactor_a = interactor_min, 
#                official_symbol_for_interactor_b = interactor_max) %>%
#         group_by(official_symbol_for_interactor_a, official_symbol_for_interactor_b) %>%
#         summarise(n = n()) %>% filter(n > 1)
#       
#       interactions <- rbind(interactions, interactions.between)
#     }
#   }
# } else {
#   interactions <- complex.interactions
# }



# interactions <- bg("interactions") %>%
#   bg_constrain(taxId = 9606, 
#                geneList = str_c(as.vector(complex.matches$`subunits(Entrez IDs)`[!is.na(complex.matches$`subunits(Entrez IDs)`)]), collapse="|"),
#                searchIds = T,
#                includeInteractors = F) %>%
#   bg_get_results()
# 
# if (nrow(interactions) > 0) {
#   interactions <- interactions %>%
#     filter(entrez_gene_id_for_interactor_a != "-" & entrez_gene_id_for_interactor_b != "-") %>%
#     filter(entrez_gene_id_for_interactor_a != entrez_gene_id_for_interactor_b) %>%
#     mutate(official_symbol_for_interactor_a = getSYMBOL(as.character(entrez_gene_id_for_interactor_a), data='org.Hs.eg'),
#            official_symbol_for_interactor_b = getSYMBOL(as.character(entrez_gene_id_for_interactor_b), data='org.Hs.eg')) %>%
#     mutate(interactor_min = pmin(official_symbol_for_interactor_a, official_symbol_for_interactor_b),
#            interactor_max = pmax(official_symbol_for_interactor_a, official_symbol_for_interactor_b)) %>%
#     mutate(official_symbol_for_interactor_a = interactor_min, 
#            official_symbol_for_interactor_b = interactor_max) %>%
#     group_by(official_symbol_for_interactor_a, official_symbol_for_interactor_b) %>%
#     summarise(n = n())
# }


#interactions <- interactions %>% filter(experimental_system_name != "Two-hybrid" & experimental_system_type == "physical") %>%
#  filter(entrez_gene_id_for_interactor_a != "-" & entrez_gene_id_for_interactor_b != "-") %>%
#  filter(entrez_gene_id_for_interactor_a != entrez_gene_id_for_interactor_b) %>%
#  mutate(official_symbol_for_interactor_a = getSYMBOL(as.character(entrez_gene_id_for_interactor_a), data='org.Hs.eg'),
#         official_symbol_for_interactor_b = getSYMBOL(as.character(entrez_gene_id_for_interactor_b), data='org.Hs.eg')) %>%
#  group_by(official_symbol_for_interactor_a, official_symbol_for_interactor_b) %>%
#  summarise(n = n()) %>% filter(n > 1)







# Edges SARS-CoV-2
#CoV2.interactions <- tibble()
#for (bait in CoV2.proteins) {
#  num_evidence <- rowSums(tibble(K=str_detect(diff.proteins$Bait.Krogan, bait),
#                                 M=str_detect(diff.proteins$Bait.Mann, bait),
#                                 L=str_detect(diff.proteins$Bait.Liang, bait)), 
#                          na.rm=T)
#  int.pos <- which(num_evidence > 0)
#  
#  CoV2.interactions <- rbind(CoV2.interactions, tibble(official_symbol_for_interactor_a=bait, 
#                                                       official_symbol_for_interactor_b=diff.proteins$`Gene names`[int.pos], 
#                                                       n=num_evidence[int.pos]))
#}
# Edge file
#edges <- rbind(interactions, CoV2.interactions)
#write_tsv(edges, here("tables/edges.tsv"), na="")




