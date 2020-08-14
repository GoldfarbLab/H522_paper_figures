library(here)
library(scales)
library(RColorBrewer)

source(here("common.R"))

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
                             expand = c(0.25, 0))
        + scale_y_continuous(name = "Viral RNA (copies/cell)",
                             breaks = c(2, 3, 4, 5, 6, 7, 8),
                             labels = c("1e2", "1e3", "1e4", "1e5", "1e6", "1e7", "1e8"),
                             limits = c(min.y, 8))
        + scale_color_manual(name = "MOI", values=colors.Cell.line)
        
        + ggtitle(title)
        
        + theme.basic
        + theme(legend.position = "none")
  )
}



################################################################################
# Figure 1B
################################################################################
plotViralLoad <- function(data) {
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
                             breaks = c(4, 24, 48, 72, 96),
                             labels = c("4", "24", "48", "72", "96"))
        + scale_y_continuous(name = "Viral RNA (copies/cell)",
                             breaks = c(0, 1, 2, 3, 4, 5, 6),
                             labels = c("0", "1e1", "1e2", "1e3", "1e4", "1e5", "1e6"),
                             limits = c(0, 6))
        + scale_color_manual(name = "MOI", values=colors.MOI)
        
        + ggtitle("Viral load in cells - H522")
        
        + theme.basic
        + theme(legend.position = "none")
  )
}

################################################################################
# Figure 1C
################################################################################
plotViralLoadSups <- function(data, y.axis.label) {
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
                             breaks = c(4, 24, 48, 72, 96),
                             labels = c("4", "24", "48", "72", "96"))
        + scale_y_continuous(name = "Viral RNA (copies/mL)",
                             breaks = c(6, 7, 8, 9, 10, 11),
                             labels = c("1e6", "1e7", "1e8", "1e9", "1e10", "1e11"),
                             limits = c(6, 11))
        + scale_color_manual(name = "MOI", values=colors.MOI)
        
        + ggtitle("Viral RNA in media - H522")
        
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
                             breaks = c(4, 24, 48, 72, 96),
                             labels = c("4", "24", "48", "72", "96"))
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
p1a <- plotViralLoadCellLines(Viral.load.cell.lines, "Viral load in cells - All lines", 2)
p1b <- plotViralLoad(Viral.load.H522.EXP08.09)
p1c <- plotViralLoadSups(Viral.load.insups.H522.EXP08.09)
p1d <- plotFACS(FACS.H522.EXP08.09)

F1 <- arrangeGrob(p1a, p1b, p1c,
             p1d,
             nrow = 2,
             ncol = 3)

#grid.draw(F1) # to view the plot
saveFig(F1, "Figure1", 4, 6.85)

# Figure S1
#ps1a <- plotViralLoadCellLines(Viral.load.cell.lines.low.MOI, "Cell-associated Viral RNA (MOI 0.015) - All lines", 4)
#FS1 <- arrangeGrob(ps1a,
#             nrow = 2,
#             ncol = 3)

#saveFig(FS1, "FigureS1", 4, 6.85)

