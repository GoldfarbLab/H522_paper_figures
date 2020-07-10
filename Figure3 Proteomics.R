library(tidyverse)


#read file 
#path <- "/Users/riajasuja/Box/CellBio-GoldfarbLab/Users/Ria Jasuja/SARS CoV 2 Paper Data/proteinGroups1.txt"

#another function that reads Protein Groups file : readProteinGroups(file, meta, measure.cols = NULL, data.cols = proteinColumns)

data <- read_tsv(path)

#filter data for contaminants, cols ES (Only ID'd by site) ET (Reverse)  EU (Potential Contaminant) and K (at least 1 unique peptide)
# do i need to rename these cols to access them in with filter fx or will this work? also like should this be a fx of like organizing these files? 

filtered.data <- filter(data, is.na(data$"Only identified by site")) # data$"Reverse" != "+", data$"Potential contaminant" != "+", data$"Razor + unique peptides" > 1)

#condense data to just the useful columns, ie Gene Names and Reporter Intensities 

reporter.intensities <- grep("Reporter intensity corrected", colnames(data), ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE, value = TRUE)
filtered.data <- select(data, "Protein names", "Gene names", reporter.intensities)

#PCA plot: figure out what Reporter Intensity Cols correspond to which samples 

#rename and reorganize table to make sense 
#1- figure out how to rename reporter intensities to sample names (probably use Grep?)
#2- reorganize cols to combine time points 

# combine time points 

#run time course? 
