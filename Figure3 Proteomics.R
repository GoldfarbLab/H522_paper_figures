library(tidyverse)

#read file 
#path <- "/Users/riajasuja/Box/CellBio-GoldfarbLab/Users/Ria Jasuja/proteinGroups.txt" 

#another function that reads Protein Groups file : readProteinGroups(file, meta, measure.cols = NULL, data.cols = proteinColumns)

data <- read_tsv(path)
design <- read_csv("data/Experimental Design Proteomics.csv")

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
#filter data for contaminants, cols ES (Only ID'd by site) ET (Reverse)  EU (Potential Contaminant) and K (at least 1 unique peptide)

filtered.data <- filter(data, is.na(data$"Only identified by site"), is.na(data$"Reverse"), is.na(data$"Potential contaminant"), data$"Razor + unique peptides" > 1)

#condense data to just the useful columns, ie Gene Names and Reporter Intensities 

reporter.intensities <- grep("Reporter intensity corrected", colnames(data), ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE, value = TRUE)
filtered.data <- select(filtered.data, "Protein names", "Gene names", reporter.intensities)

#At least 1 value for each Rep (aka not all 0s)
#filter data into Reps 
intensities.rep1 <- grep("Rep 1", colnames(filtered.data), ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE, value = TRUE)
intensities.rep2 <- grep("Rep 2", colnames(filtered.data), ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE, value = TRUE)
intensities.rep3 <- grep("Rep 3", colnames(filtered.data), ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE, value = TRUE)
#At least 1 nonzero value / Rep (idea: sum each row for each rep and if its more than 0 its good) 
rep1.valid <- rowSums(select(filtered.data, intensities.rep1) > 0) > 0
rep2.valid <- rowSums(select(filtered.data, intensities.rep2) > 0) > 0
rep3.valid <- rowSums(select(filtered.data, intensities.rep3) > 0) > 0

valid.data <- filter(filtered.data, (rep1.valid + rep2.valid + rep3.valid) == 3)

#PCA plot: figure out what Reporter Intensity Cols correspond to which samples 
plotPCA <- function(data, design)
{
  data <- select(data, -c(1,10,11,12,21,22,23,32,33)) #removes TMT1 (light standard) and TMT11 (heavy standard) (AND BRIDGE)
  design <- design[-c(1,10,11,12,21,22,23,32,33),]
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
  
  p <- (ggplot(plotting.data, aes(x=PC1, y=PC2, shape=as.factor(Replicate), color=Condition, label=Condition)) 
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
plotPCA(select(valid.data, reporter.intensities), design)

#rename and reorganize table to make sense 
#1- figure out how to rename reporter intensities to sample names (probably use Grep?)
#2- reorganize cols to combine time points 

# combine time points 

#run time course? 
