library(tidyverse)
library(scales)
library(RColorBrewer)
library(here)

source(here("common.R"))

################################################################################
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
################################################################################
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

################################################################################
# Figure 1A
################################################################################
plotViralLoadCellLines <- function(data, title, min.y) {
  dfwc_between <- summarySE(data=data, measurevar="Viral Load", groupvars=c("Cell.line", "Time"), na.rm=FALSE, conf.interval=.95)
  dfwc_between$logVL <- log10(dfwc_between$`Viral Load`)
  dfwc_between$logVL[which(is.infinite(dfwc_between$logVL))] <- 0
  
  p <- (ggplot(dfwc_between, aes(x=Time, 
                                 y=logVL, 
                                 group=Cell.line,
                                 color=Cell.line))
        
        + geom_errorbar(width=3, 
                        aes(ymin=log10(`Viral Load` - se),
                            ymax=log10(`Viral Load` + se)))
        + geom_line(size=0.5)
        + geom_point(size=1.5)
        
        
        + scale_x_continuous(name = "Hours post-infection", 
                             breaks = c(4, 72),
                             labels = c("4","72"),
                             expand = c(0.25,0))
        + scale_y_continuous(name = expression("Viral RNA (copies/"*mu*"L)"),
                             breaks = c(4, 5, 6, 7, 8, 9, 10, 11, 12),
                             labels = c("1e4", "1e5", "1e6", "1e7", "1e8", "1e9", "1e10", "1e11", "1e12"),
                             limits = c(min.y, 12))
        + scale_color_manual(name = "MOI", values=colors.Cell.line)
        
        + ggtitle(title)
        
        + theme.basic
        + theme(legend.position = "none")
  )
}



################################################################################
# Figure 1B
################################################################################
plotViralLoad <- function(data, title, y.axis.label) {
  dfwc_between <- summarySE(data=data, measurevar="Viral Load", groupvars=c("MOI","Time"), na.rm=FALSE, conf.interval=.95)
  dfwc_between$logVL <- log10(dfwc_between$`Viral Load`)
  dfwc_between$logVL[which(is.infinite(dfwc_between$logVL))] <- 0
  
  p <- (ggplot(dfwc_between, aes(x=Time, 
                                 y=logVL, 
                                 group=MOI,
                                 color=as.factor(MOI)))
        
        + geom_errorbar(width=3, 
                        aes(ymin=log10(`Viral Load` - se),
                            ymax=log10(`Viral Load` + se)))
        + geom_line(size=0.5)
        + geom_point(size=1.5)
        
        
        + scale_x_continuous(name = "Hours post-infection", 
                             breaks = c(0, 4, 24, 48, 72, 96),
                             labels = c("Mock", "4", "24", "48", "72", "96"))
        + scale_y_continuous(name = y.axis.label,
                             breaks = c(0, 2, 4, 6, 8, 10),
                             labels = c("0", "1e2", "1e4", "1e6", "1e8", "1e10"),
                             limits = c(0, 10.8))
        + scale_color_manual(name = "MOI", values=colors.MOI)
        
        + ggtitle(title)
        
        + theme.basic
        + theme(legend.position = "none")
  )
}


################################################################################
# Figure 1D
################################################################################
plotFACS <- function(data) {
  dfwc_between <- summarySE(data=data, measurevar="Infected Cells", groupvars=c("MOI","Time"), na.rm=FALSE, conf.interval=.95)
  
  p <- (ggplot(dfwc_between, aes(x=Time, 
                                 y=`Infected Cells`, 
                                 group=MOI,
                                 color=as.factor(MOI)))
        
        + geom_errorbar(width=3, 
                        aes(ymin=`Infected Cells`-se,
                            ymax=`Infected Cells`+se))
        + geom_line(size=0.5)
        + geom_point(size=1.5)
        
        + scale_x_continuous(name = "Hours post-infection", 
                             breaks = c(0, 4, 24, 48, 72, 96),
                             labels = c("Mock", "4", "24", "48", "72", "96"))
        + scale_y_continuous(name = "SARS-CoV-2 (+) cells",
                             breaks = c(0, 10, 20, 30, 40),
                             labels = c("0%", "10%", "20%", "30%", "40%"))
        + scale_color_manual(name = "MOI", values=colors.MOI)
        
        + ggtitle("H522 - FACS")
        
        + theme.basic
        + theme(legend.position = "none")
  )
}





################################################################################
# Read data
################################################################################
# Figure 1
FACS.H522.EXP08.09 <- read_csv(here("data/FACS H522_EXP08 and EXP09.csv"))
Viral.load.H522.EXP08.09 <- read_csv(here("data/Viral Load H522_EXP08 and EXP09.csv"))
Viral.load.insups.H522.EXP08.09 <- read_csv(here("data/Viral Load-insups H522 EXP08 and EXP09.csv"))
Viral.load.cell.lines <- read_csv(here("data/Viral Load-cell lines.csv"))
# Figure S1
Viral.load.cell.lines.low.MOI <- read_csv(here("data/Viral Load-cell lines MOI 015.csv"))

################################################################################
# Generate figures
################################################################################
# Figure 1
p1a <- plotViralLoadCellLines(Viral.load.cell.lines, "Cell-associated Viral RNA - All lines", 5)
p1b <- plotViralLoad(Viral.load.H522.EXP08.09, "Cell-associated Viral RNA - H522", expression("Viral RNA (copies/"*mu*"L)"))
p1c <- plotViralLoad(Viral.load.insups.H522.EXP08.09, "Viral RNA in media - H522", "Viral RNA (copies/mL)" )
p1d <- plotFACS(FACS.H522.EXP08.09)

F1 <- arrangeGrob(p1a, p1b, p1c,
             p1d,
             nrow = 2,
             ncol = 3)

#grid.draw(F1) # to view the plot
saveFig(F1, "Figure1", 7.5, 4)

# Figure S1
ps1a <- plotViralLoadCellLines(Viral.load.cell.lines.low.MOI, "Cell-associated Viral RNA (MOI 0.015) - All lines", 4)
FS1 <- arrangeGrob(ps1a,
             nrow = 2,
             ncol = 3)

saveFig(FS1, "FigureS1", 7.5, 4)

