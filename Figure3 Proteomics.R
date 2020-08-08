library(timecourse)
library(tidyverse)

#READ FILE
path <- "/Users/riajasuja/Box/CellBio-GoldfarbLab/data/mass spec/analysis/H522 SARS-CoV-2/txt/proteinGroups.txt"  #ACTUAL PATH
  #practice data path: path <-  "/Users/riajasuja/Box/CellBio-GoldfarbLab/Users/Ria Jasuja/proteinGroups1.txt"
data <- read_tsv(path)
design <- read_csv("data/Experimental Design H522 Paper.csv") #ACTUAL DESIGN FILE
  #practice design file: design <- read_csv("data/Experimental Design Proteomics.csv")
  #TMT9: Reference Channel 
  #TMT10: remove

golden.ratio <- 1/1.618
grey <- "#333333"
light.grey <- "#AAAAAA"
theme.basic <- (theme_minimal() 
                + theme(axis.line = element_line(colour = grey, size = 1, linetype = "solid"), 
                        panel.grid = element_blank(),
                        axis.ticks = element_line(size = 0.5, colour = grey),
                        axis.ticks.length = unit(.25, "cm"),
                        aspect.ratio = golden.ratio,
                        plot.title = element_text(size=11, hjust=0.5))
)
#FILTERING DATA FOR NAs, CONTAMINANTS, NONZERO ROWS
data$`Gene names`[is.na(data$`Gene names`)] <- data$`Majority protein IDs`[is.na(data$`Gene names`)] #replace NA gene names with Protein ID column (Majority protein IDs)
filtered.data <- filter(data, is.na(data$"Only identified by site"), is.na(data$"Reverse"), is.na(data$"Potential contaminant"), data$"Razor + unique peptides" > 1) #filter data for contaminants, cols ES (Only ID'd by site) ET (Reverse)  EU (Potential Contaminant) and K (at least 1 unique peptide)

#condense data to just the useful columns, ie Gene Names and Reporter Intensities (AND REMOOVE TMT 10)
rep.10 <- grep("Reporter intensity corrected 10", colnames(filtered.data), ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE, value = TRUE)
filtered.data <- filtered.data[ , !names(filtered.data) %in% rep.10]
reporter.intensities <- grep("Reporter intensity corrected", colnames(filtered.data), ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE, value = TRUE)
filtered.data <- dplyr::select(filtered.data, "Gene names", reporter.intensities)


#At least 1 value for each Rep (aka not all 0s)
intensities.rep1 <- grep("Rep 1", colnames(filtered.data), ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE, value = TRUE) #filter data into Reps 
intensities.rep2 <- grep("Rep 2", colnames(filtered.data), ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE, value = TRUE)
intensities.rep3 <- grep("Rep 3", colnames(filtered.data), ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE, value = TRUE)

rep1.valid <- rowSums(dplyr::select(filtered.data, all_of(intensities.rep1)) > 0) > 0 #finding at least 1 nonzero value in each Rep 
rep2.valid <- rowSums(dplyr::select(filtered.data, all_of(intensities.rep2)) > 0) > 0
rep3.valid <- rowSums(dplyr::select(filtered.data, all_of(intensities.rep3)) > 0) > 0

valid.data <- filter(filtered.data, (rep1.valid + rep2.valid + rep3.valid) == 3)

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
normalization.data <- dplyr::select(valid.data, reporter.intensities)
normalization.data[normalization.data == 0] <- min(normalization.data[normalization.data > 0])

#REP1
normalization.data.rep1 <- dplyr::select(normalization.data, intensities.rep1)
normalization.factor.rep1 <- colSums(normalization.data.rep1)
normalization.factor.rep1 <- normalization.factor.rep1/max(normalization.factor.rep1)
normalization.data.rep1 <- normalization.data.rep1/normalization.factor.rep1
normalized.data.rep1 <- sweep(normalization.data.rep1, 1, normalization.data.rep1[, 9], "/") #dividing all rows by reference channel (TMT9)

#REP2
normalization.data.rep2 <- dplyr::select(normalization.data, intensities.rep2)
normalization.factor.rep2 <- colSums(normalization.data.rep2)
normalization.factor.rep2 <- normalization.factor.rep2/max(normalization.factor.rep2)
normalization.data.rep2 <- normalization.data.rep2/normalization.factor.rep2
normalized.data.rep2 <- sweep(normalization.data.rep2, 1, normalization.data.rep2[, 9], "/")

#REP3
normalization.data.rep3 <- dplyr::select(normalization.data, intensities.rep3)
normalization.factor.rep3 <- colSums(normalization.data.rep3)
normalization.factor.rep3 <- normalization.factor.rep3/max(normalization.factor.rep3)
normalization.data.rep3 <- normalization.data.rep3/normalization.factor.rep3
normalized.data.rep3 <- sweep(normalization.data.rep3, 1, normalization.data.rep3[, 9], "/")

normalized.data <- cbind(normalized.data.rep1, normalized.data.rep2, normalized.data.rep3)

#NaN and Inf Values? 
normalized.data <- log2(normalized.data)

#remove Ref Channel (TMT9) after normalization (after normalization: all 0s) 
ref.channel <- grep("Reporter intensity corrected 9", colnames(normalized.data), ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE, value = TRUE)
normalized.data <- normalized.data[ , !names(normalized.data) %in% ref.channel]
#normalized.data <- cbind(gene.names, normalized.data)

#PCA plot: figure out what Reporter Intensity Cols correspond to which samples 
plotPCA <- function(data, design)
{
  design <- design[-c(9,18,27),] #have to remove Ref Channel in design file too 
  data.t <- as.data.frame(t(data))
  pca <- prcomp(data.t)
  eigen <- pca$sdev^2
  variance <- (eigen/sum(eigen)) * 100
  variance <- format(round(variance, 2), nsmall=2) # show 2 digits after decimal
  
  plotting.data <- cbind(data.frame(pca$x), design)
  y.range <- max(plotting.data$PC2) - min(plotting.data$PC2)
  min.x <- min(plotting.data$PC1)
  max.x <- max(plotting.data$PC1)
  x.range <- max.x - min.x
  
  p <- (ggplot(plotting.data, aes(x=PC1, y=PC2, shape=as.factor(Infected), color=Time, label=Time)) 
        + geom_point(size=3)
        + geom_text(size=3, nudge_y = y.range/20)
        
        + scale_x_continuous(limits=c(min.x - x.range/10, max.x + x.range/10))
        
        + xlab(paste("PC1 (", variance[1], "% explained variance)", sep=""))
        + ylab(paste("PC2 (", variance[2], "% explained variance)", sep=""))
        + ggtitle("Principal Component Analysis (PCA) Proteomics")
        
        + theme.basic
  )
  print(p)
}
plotPCA(normalized.data, design)

#RUNNNIG TIMECOURSE 

#First run TC with just the INFECTED time points 
  #remove TMT 1,8 (Mock sample time points 4h 96h)
timecourse.infected.data <- dplyr::select(normalized.data, -c(1,8,9,16,17,24)) 
assay.1 <- rep(c("1", "2", "3"), each = 6)
time.grp.1 <- rep(1:6, 3)
reps <- rep(3, nrow(timecourse.infected.data))
timecourse.infected.results <- mb.long(timecourse.infected.data, method="1D", times=6, reps=reps, rep.grp = assay.1, time.grp = time.grp.1) 
print(timecourse.infected.results)
plotProfile(timecourse.infected.results, type="b", gnames=gnames, legloc=c(2,15), pch=c("1","2","3"), xlab="Time Point", gid = "MX1", col = c("pink", "black", "green")) 
#use gid= to look for a specific protein's profile 

#Then run TC with MOCK time points and INFECTED samples at those 2 time points 
  #remove TMT 3,4,5,6 (Infected middle time points)
timecourse.mock.data <- dplyr::select(normalized.data, -c(3,4,5,6,11,12,13,14,19,20,21,22)) 
trt <- rep(c("mock", "infected", "infected", "mock"), 3)
time.grp.2 <- rep(rep(1:2, each = 2),3)
assay.2 <- rep(c("1", "2", "3"), each = 4)
reps.2 <- matrix(3, nrow = nrow(timecourse.infected.data), ncol = 2)
timecourse.mock.results <-  mb.long(timecourse.mock.data, method="2D", times=2, reps=reps.2, rep.grp = assay.2, time.grp = time.grp.2, condition.grp = trt)
print(timecourse.mock.results)
plotProfile(timecourse.mock.results, type="b", gnames=gnames, legloc=c(2,15), pch=c("1","2","3"), xlab="Time Point", gid = "MX1", col = c("pink", "green")) 

#Finally make a matrix with gene names, then T2 stat for infected timecourse, then T2 stat for infected/mock timecourse
infected.data.results <- timecourse.infected.results[["HotellingT2"]]
mock.data.results <- timecourse.mock.results[["HotellingT2"]]
tc.results <- tibble(gene.names, infected.data.results, mock.data.results)





