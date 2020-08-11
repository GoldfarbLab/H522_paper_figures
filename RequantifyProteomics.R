library(tidyverse)

data.path <- "~/Box/CellBio-GoldfarbLab/Data/Mass Spec/Analysis/H522 SARS-CoV-2/txt/"


################################################################################
# Process evidence
################################################################################
evidence <- read_tsv(here('data_raw/MS/evidence.zip'), guess_max=10000)
intensity.col.names <- grep("Reporter intensity corrected", colnames(evidence), value=T)

evidence <- filter(evidence, `Reporter intensity corrected 9` > 0)
evidence <- filter(evidence, `Reporter intensity corrected 1` + `Reporter intensity corrected 2` + `Reporter intensity corrected 3`
                    + `Reporter intensity corrected 4` + `Reporter intensity corrected 5` + `Reporter intensity corrected 6`
                    + `Reporter intensity corrected 7` + `Reporter intensity corrected 8` + `Reporter intensity corrected 9` >= 4e3 )

# roll-up evidence quantifications to mod peptide level (peptide ID'ed in multiple fractions or charge states)
evidence.aggregate <- evidence %>% group_by(Experiment, `Mod. peptide ID`) %>% summarise_at(intensity.col.names, sum)




################################################################################
# Process modification specific peptides
################################################################################
mod.peptides <- read_tsv(here('data_raw/MS/modificationSpecificPeptides.zip'), guess_max=10000)
intensity.col.rep.names <- grep("Reporter intensity corrected.*R", colnames(mod.peptides), value=T)
mod.peptides <- filter(mod.peptides, `Acetyl (Protein N-term)` == 0 & `Deamidation (N)` == 0)

# clear old quantification
mod.peptides[, grep("Reporter intensity corrected ", colnames(mod.peptides))] <- 0

# get replacement values
tmp.1 <- left_join(mod.peptides, filter(evidence.aggregate, `Experiment` == "Rep 1"), by=c("id" = "Mod. peptide ID"), suffix=c("x",""))
tmp.2 <- left_join(mod.peptides, filter(evidence.aggregate, `Experiment` == "Rep 2"), by=c("id" = "Mod. peptide ID"), suffix=c("x",""))
tmp.3 <- left_join(mod.peptides, filter(evidence.aggregate, `Experiment` == "Rep 3"), by=c("id" = "Mod. peptide ID"), suffix=c("x",""))

# replace with new quant
mod.peptides[,  grep("Reporter intensity corrected.*Rep 1", colnames(mod.peptides), value=T)] <- tmp.1[, intensity.col.names]
mod.peptides[,  grep("Reporter intensity corrected.*Rep 2", colnames(mod.peptides), value=T)] <- tmp.2[, intensity.col.names]
mod.peptides[,  grep("Reporter intensity corrected.*Rep 3", colnames(mod.peptides), value=T)] <- tmp.3[, intensity.col.names]

# replace NAs
mod.peptide.quant <- mod.peptides[, intensity.col.rep.names]
mod.peptides[, intensity.col.rep.names] <- replace(mod.peptide.quant, is.na(mod.peptide.quant), 0)

# roll-up peptide quantifications to protein level
mod.peptides.aggregate <- mod.peptides %>% group_by(`Peptide ID`) %>% summarise_at(intensity.col.rep.names, sum)


################################################################################
# Process peptides
################################################################################
peptides <- read_tsv(here('data_raw/MS/peptides.zip'), guess_max=10000)
#peptides <- filter(peptides, `Unique (Groups)` == "yes")

# clear old quantification
peptides[, intensity.col.rep.names] <- 0

# get replacement values
tmp <- left_join(peptides[,"id"], mod.peptides.aggregate, by=c("id" = "Peptide ID"), suffix=c("","y"))

# replace NAs
tmp <- tmp[, intensity.col.rep.names]
tmp <- replace(tmp, is.na(tmp), 0)

# replace with new quant
peptides[, intensity.col.rep.names] <- tmp[, intensity.col.rep.names]

# roll-up peptide quantifications to protein level
peptides.aggregate <- peptides %>% group_by(`Leading razor protein`) %>% summarise_at(intensity.col.rep.names, sum)

# roll up peptide quantification counts to protein level
peptides.aggregate.count <- peptides %>% group_by(`Leading razor protein`) %>% summarise_at(intensity.col.rep.names, funs(sum(. > 0)))

# rename count columns
colnames(peptides.aggregate.count) <- str_replace(colnames(peptides.aggregate.count), "corrected", "count")

# join results
peptides.aggregate <- left_join(peptides.aggregate, peptides.aggregate.count)


################################################################################
# Process protein groups
################################################################################
proteins <- read_tsv(here('data_raw/MS/proteinGroups.zip'), guess_max=10000)
count.col.rep.names <- grep("Reporter intensity count.*R", colnames(proteins), value=T)
# extract leading razor protein
proteins <- separate(proteins, `Majority protein IDs`, c("Leading razor protein"), sep=";", remove=F, extra="drop")

# clear old quantification
proteins[, intensity.col.rep.names] <- 0
proteins[, count.col.rep.names] <- 0

# get replacement values
tmp <- left_join(proteins[,"Leading razor protein"], peptides.aggregate, "Leading razor protein", suffix=c("","y"))

# replace NAs
tmp <- select(tmp, -`Leading razor protein`)
tmp <- replace(tmp, is.na(tmp), 0)

# replace with new quant
proteins[, intensity.col.rep.names] <- tmp[, intensity.col.rep.names]
proteins[, count.col.rep.names] <- tmp[, count.col.rep.names]

# write new proteinGroups file
write_tsv(proteins, here("data_processed/proteins.txt"))
