library(here)
library(scales)
library(RColorBrewer)

source(here("common.R"))

################################################################################
# Figure 1B ACE2
################################################################################
plotPCR.ACE2 <- function(data, title) {
  data.Vero <- data %>% filter(Cell.line == "Vero E6")
  data.woVero <- data %>% filter(Cell.line != "Vero E6")
  dfwc_between <- summarySE(data=data.woVero, measurevar="Norm.expression", groupvars=c("Cell.line"), na.rm=FALSE, conf.interval=.95)
  dfwc_between <- dfwc_between %>% add_row(Cell.line = data.Vero$Cell.line[1], Norm.expression = data.Vero$Norm.expression[1], se=data.Vero$sd[1])
  dfwc_between$Cell.line <- factor(dfwc_between$Cell.line, levels= order.Cell.line)
  
  p <- (ggplot(dfwc_between, aes(x=Cell.line,
                                 y=Norm.expression,
                                 fill=Cell.line,
                                 color=Cell.line))
        
        + geom_errorbar(width=0.33, size=0.3, 
                        aes(ymin=Norm.expression - se,
                            ymax=Norm.expression + se,
                            color=Cell.line))
        + geom_bar(stat="identity", width=0.5, position = position_dodge(width=1))
        
        + geom_vline(xintercept=1.5, linetype="dashed", size=0.25)
        
        + scale_y_continuous(name = "Normalized Expression")
        + scale_fill_manual(name = "Cell line", values=colors.Cell.line)
        + scale_color_manual(name = "Cell line", values=colors.Cell.line)
        
        + ggtitle(title)
        
        + theme.basic
        + theme(legend.position = "none", 
                axis.title.x = element_blank(), 
                axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  )
}

################################################################################
# Figure 1B TMPRSS2
################################################################################
plotPCR.TMPRSS2 <- function(data, title) {
  data.Vero <- data %>% filter(Cell.line == "Vero E6")
  data.woVero <- data %>% filter(Cell.line != "Vero E6")
  dfwc_between <- summarySE(data=data.woVero, measurevar="Norm.expression", groupvars=c("Cell.line"), na.rm=FALSE, conf.interval=.95)
  dfwc_between <- dfwc_between %>% add_row(Cell.line = data.Vero$Cell.line[1], Norm.expression = data.Vero$Norm.expression[1], se=data.Vero$sd[1])
  dfwc_between$Cell.line <- factor(dfwc_between$Cell.line, levels= order.Cell.line)
  
  trans <- function(x){pmin(x,8e-10) + 0.05*pmax(x-8e-10,0)}
  yticks <- c(0, 2e-10, 4e-10, 6e-10, 8e-10, 5e-9, 1e-8)
  
  dfwc_between$Norm.expression_t <- trans(dfwc_between$Norm.expression)
  dfwc_between$se_up_t <- trans(dfwc_between$Norm.expression + dfwc_between$se)
  dfwc_between$se_low_t <- trans(dfwc_between$Norm.expression - dfwc_between$se)
  
  
  p <- (ggplot(dfwc_between, aes(x=Cell.line,
                                 y=Norm.expression_t,
                                 fill=Cell.line,
                                 color=Cell.line))
        
        + geom_errorbar(width=0.33, size=0.3, 
                        aes(ymin=se_low_t,
                            ymax=se_up_t,
                            color=Cell.line),
                        position = "dodge")
        + geom_bar(stat="identity", width=0.5, position = "dodge")
        
        + geom_vline(xintercept=1.5, linetype="dashed", size=0.25)
        
        + geom_rect(aes(xmin=0, xmax=11, ymin=8.4e-10, ymax=9.2e-10), fill="white", color="white")
        
        + scale_y_continuous(name = "Normalized Expression", limits=c(0,NA), breaks=trans(yticks), labels=yticks)
        
        + scale_fill_manual(name = "Cell line", values=colors.Cell.line)
        + scale_color_manual(name = "Cell line", values=colors.Cell.line)
        
        + ggtitle(title)
        
        + theme.basic
        + theme(legend.position = "none", 
                axis.title.x = element_blank(), 
                axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  )
}

################################################################################
# Figure 1A
################################################################################
plotViralLoadCellLines <- function(data, title) {
  dfwc_between <- summarySE(data=data, measurevar="Viral Load", groupvars=c("Cell.line", "Time", "MOI"), na.rm=FALSE, conf.interval=.95)
  dfwc_between$logVL <- log10(dfwc_between$`Viral Load`)
  dfwc_between$logVL[which(is.infinite(dfwc_between$logVL))] <- 0
  
  p <- (ggplot(dfwc_between, aes(x=Time, 
                                 y=logVL, 
                                 group=Cell.line,
                                 color=Cell.line,
                                 facet=MOI))
        
        + geom_errorbar(width=3, size=0.3,
                        aes(ymin=log10(`Viral Load` - se),
                            ymax=log10(`Viral Load` + se)))
        + geom_line(size=0.5)
        + geom_point(size=1.5)
        
        
        + scale_x_continuous(name = "Hours post-infection", 
                             breaks = c(4, 72),
                             labels = c("4","72"),
                             expand = c(0.15, 0))
        + scale_y_continuous(name = "Viral RNA (copies/cell)",
                             breaks = seq(1, 8),
                             labels = paste("1e", seq(1, 8), sep=""),
                             limits = c(1, 8))
        + scale_color_manual(name = "Cell line", values=colors.Cell.line)
        
        + facet_wrap(~ MOI, nrow=1)
        
        + ggtitle(title)
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 1.3,
                strip.text.x = element_text(size = 6))
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
        
        + geom_errorbar(width=3, size=0.3, 
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
        
        + geom_errorbar(width=3, size=0.3,
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
  dfwc_between <- summarySE(data=data, measurevar="Infected Cells", groupvars=c("MOI","Time", "Cell.line"), na.rm=FALSE, conf.interval=.95)
  
  p <- (ggplot(dfwc_between, aes(x=Time, 
                                 y=`Infected Cells`, 
                                 group=interaction(MOI,Cell.line),
                                 color=as.factor(MOI),
                                 linetype=Cell.line))
        
        + geom_errorbar(width=3, size=0.3,
                        aes(ymin=`Infected Cells`-se,
                            ymax=`Infected Cells`+se))
        + geom_line(size=0.5)
        + geom_point(size=1.5)
        
        + scale_x_continuous(name = "Hours post-infection", 
                             breaks = c(4, 24, 48, 72, 96),
                             labels = c("4", "24", "48", "72", "96"))
        + scale_y_continuous(name = "SARS-CoV-2 (+) cells (%)",
                             breaks = c(0, 10, 20, 30, 40, 50, 60),
                             labels = c("0", "10", "20", "30", "40", "50", "60"),
                             limits = c(0, 60))
        + scale_color_manual(name = "MOI", values=colors.MOI)
        
        + ggtitle("Infected Cells")
        
        + theme.basic
        + theme(legend.position = "none")
  )
}





################################################################################
# Read data
################################################################################
# Figure 1
ACE2.PCR <- read_tsv(here("data/Figure1B_ACE2_PCR.txt"))
TMPRSS2.PCR <- read_tsv(here("data/Figure1B_TMPRSS2_PCR.txt"))
FACS.H522.EXP08.09 <- read_csv(here("data/FACS H522_EXP08 and EXP09.csv"))
Viral.load.H522.EXP08.09 <- read_csv(here("data/Viral Load H522_EXP08 and EXP09.csv"))
Viral.load.insups.H522.EXP08.09 <- read_csv(here("data/Viral Load-insups H522 EXP08 and EXP09.csv"))
Viral.load.cell.lines <- read_csv(here("data/Viral Load-cell lines - all HCC827.csv"))
# Figure S1
Viral.load.cell.lines.low.MOI <- read_csv(here("data/Viral Load-cell lines MOI 015.csv"))

################################################################################
# Generate figures
################################################################################
# Figure 1
p1b <- plotPCR.ACE2(ACE2.PCR, "ACE2 mRNA")
p1b2 <- plotPCR.TMPRSS2(TMPRSS2.PCR, "TMPRSS2 mRNA")


p1d <- plotViralLoadCellLines(Viral.load.cell.lines, "Viral load in cells - All lines")
#p1d2 <- plotViralLoadCellLines(Viral.load.cell.lines.low.MOI, "Viral load in cells - All lines", 1)
p1e <- plotViralLoad(Viral.load.H522.EXP08.09)
p1f <- plotViralLoadSups(Viral.load.insups.H522.EXP08.09)
p1g <- plotFACS(FACS.H522.EXP08.09)

F1.top <- arrangeGrob(
  p1b, p1b2,
  nrow = 1,
  ncol = 2)

F1.bottom <- arrangeGrob(
  p1d, p1e, p1f,
  p1g,
  nrow = 2,
  ncol = 3)

#grid.draw(F1) # to view the plot
saveFig(F1.top, "Figure1_top", 2, 4)
saveFig(F1.bottom, "Figure1_bottom", 4, 6.85)

# Figure S1
#ps1a <- plotViralLoadCellLines(Viral.load.cell.lines.low.MOI, "Cell-associated Viral RNA (MOI 0.015) - All lines", 4)
#FS1 <- arrangeGrob(ps1a,
#             nrow = 2,
#             ncol = 3)

#saveFig(FS1, "FigureS1", 4, 6.85)

