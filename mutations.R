library(CellSurfaceTools)
library(tidyverse)
library(here)

# Read in cell line mutations
# H522, H1299, KYSE30, PC9, SCC25, OE21, A427, H596
data.path <- "~/Box/CellBio-MajorLab/UNC/P Drive/Users/Former students and postdocs/Emily/V Foundation/mutation and expression data/Vfoundation_TO_mut_calls 2/"
H522 <- read_tsv(str_c(data.path, "55828830_T_00001ADEFL_final_VCs_w_type.tsv"))
H1299 <- read_tsv(str_c(data.path, "55828831_T_00001ADF2L_final_VCs_w_type.tsv"))
KYSE30 <- read_tsv(str_c(data.path, "54972491_T_000019FD9L_final_VCs_w_type.tsv"))
PC9 <- read_tsv(str_c(data.path, "55828838_T_00001AE08L_final_VCs_w_type.tsv"))
SCC25 <- read_tsv(str_c(data.path, "54972494_T_000019FDCL_final_VCs_w_type.tsv"))
OE21 <- read_tsv(str_c(data.path, "54972492_T_000019FDAL_final_VCs_w_type.tsv"))
A427 <- read_tsv(str_c(data.path, "55828826_T_00001ADE4L_final_VCs_w_type.tsv"))
H596 <- read_tsv(str_c(data.path, "54972495_T_000019FDDL_final_VCs_w_type.tsv"))
# Read in expression
models <- read_tsv("~/Box/CellBio-GoldfarbLab/Projects/Current/COVID-19/RNAseq models/models_Jun1.txt")
# Read in proteomics
proteome <- read_tsv(here("data_processed/proteinsNormedToBridge.txt"))


# get all and unique H552 mutations
H522.all <- left_join(H522, H1299, by = c("Position", "Alt"), suffix=c("",".H1299")) %>%
  left_join(KYSE30, by = c("Position", "Alt"), suffix=c("",".KYSE")) %>%
  left_join(PC9, by = c("Position", "Alt"), suffix=c("",".PC9")) %>%
  left_join(SCC25, by = c("Position", "Alt"), suffix=c("",".SCC25")) %>%
  left_join(OE21, by = c("Position", "Alt"), suffix=c("",".OE21")) %>%
  left_join(A427, by = c("Position", "Alt"), suffix=c("",".A427")) %>%
  left_join(H596, by = c("Position", "Alt"), suffix=c("",".H596")) %>%
  filter(snpEff != "SYNONYMOUS_CODING") %>%
  mutate("is_unique_mutation" = is.na(GeneName.H1299) & is.na(GeneName.KYSE) & is.na(GeneName.PC9) & 
           is.na(GeneName.SCC25) & is.na(GeneName.OE21) & is.na(GeneName.A427) & is.na(GeneName.H596)) %>%
  select(1:19, "is_unique_mutation") %>%
  group_by(GeneName) %>%
  summarize(position = str_c(str_c(Ref, Position, Alt, sep="_"), collapse=";"), 
            prob=str_c(prob, collapse=";"),
            snpEff=str_c(snpEff, collapse=";"),
            snpEffImpact=str_c(snpEffImpact, collapse=";"),
            type=str_c(type, collapse=";"),
            is_unique_mutation=min(is_unique_mutation))

write_tsv(H522.all, here("annotations/H522_mutations.tsv"))




# annotate and filter for surface and plasma membrane proteins
H522.surf <- left_join(H522, human_surface_and_plasma_membrane_geneLevel, by = c("GeneName" = "Gene Name")) %>% filter(!is.na(GeneID))
H1299.surf <- left_join(H1299, human_surface_and_plasma_membrane_geneLevel, by = c("GeneName" = "Gene Name")) %>% filter(!is.na(GeneID))
KYSE30.surf <- left_join(KYSE30, human_surface_and_plasma_membrane_geneLevel, by = c("GeneName" = "Gene Name")) %>% filter(!is.na(GeneID))
PC9.surf <- left_join(PC9, human_surface_and_plasma_membrane_geneLevel, by = c("GeneName" = "Gene Name")) %>% filter(!is.na(GeneID))
SCC25.surf <- left_join(SCC25, human_surface_and_plasma_membrane_geneLevel, by = c("GeneName" = "Gene Name")) %>% filter(!is.na(GeneID))
OE21.surf <- left_join(OE21, human_surface_and_plasma_membrane_geneLevel, by = c("GeneName" = "Gene Name")) %>% filter(!is.na(GeneID))
A427.surf <- left_join(A427, human_surface_and_plasma_membrane_geneLevel, by = c("GeneName" = "Gene Name")) %>% filter(!is.na(GeneID))
H596.surf <- left_join(H596, human_surface_and_plasma_membrane_geneLevel, by = c("GeneName" = "Gene Name")) %>% filter(!is.na(GeneID))

# get unique H552 mutations
H522.surf.unique <- left_join(H522.surf, H1299.surf, by = c("Position", "Alt"), suffix=c("",".H1299")) %>%
  left_join(KYSE30.surf, by = c("Position", "Alt"), suffix=c("",".KYSE")) %>%
  left_join(PC9.surf, by = c("Position", "Alt"), suffix=c("",".PC9")) %>%
  left_join(SCC25.surf, by = c("Position", "Alt"), suffix=c("",".SCC25")) %>%
  left_join(OE21.surf, by = c("Position", "Alt"), suffix=c("",".OE21")) %>%
  left_join(A427.surf, by = c("Position", "Alt"), suffix=c("",".A427")) %>%
  left_join(H596.surf, by = c("Position", "Alt"), suffix=c("",".H596")) %>%
  filter(is.na(GeneName.H1299) & is.na(GeneName.KYSE) & is.na(GeneName.PC9) & 
           is.na(GeneName.SCC25) & is.na(GeneName.OE21) & is.na(GeneName.A427) & is.na(GeneName.H596)) %>%
  filter(snpEff != "SYNONYMOUS_CODING") %>%
  select(1:27) %>%
  left_join(select(models, "Gene", "H522", "r2", "slope", "intercept"), by=c("GeneName" = "Gene")) %>%
  filter(H522 > 10) %>%
  left_join(proteome, by=c("GeneName"="Gene names"))

H522.surf.only <- filter(H522.surf.unique, `Num cell surface evidence` > 0)

H522.surf.only.ms <- filter(H522.surf.only, !is.na(`Protein IDs`))
# unique(H522.surf.unique$GeneName[which(H522.surf.unique$`Num cell surface evidence` > 0)])



