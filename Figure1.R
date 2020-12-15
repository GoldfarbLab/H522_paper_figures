library(here)
source(here("common.R"))

################################################################################
# Plotting functions
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
                                 color=Cell.line,
                                 ymin=Norm.expression - se,
                                 ymax=Norm.expression + se))
        
        + geom_errorbar(width=0.33, size=0.3)
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


################################################################################
plotPCR.TMPRSS2 <- function(data, title) {
  data.Vero <- data %>% filter(Cell.line == "Vero E6")
  data.woVero <- data %>% filter(Cell.line != "Vero E6")
  dfwc_between <- summarySE(data=data.woVero, measurevar="Norm.expression", groupvars=c("Cell.line"), na.rm=FALSE, conf.interval=.95)
  dfwc_between <- dfwc_between %>% add_row(Cell.line = data.Vero$Cell.line[1], Norm.expression = data.Vero$Norm.expression[1], se=data.Vero$sd[1])
  dfwc_between$Cell.line <- factor(dfwc_between$Cell.line, levels=order.Cell.line)
  
  trans <- function(x){pmin(x,8e-10) + 0.05*pmax(x-8e-10,0)}
  yticks <- c(0, 2e-10, 4e-10, 6e-10, 8e-10, 5e-9, 1e-8)
  
  dfwc_between$Norm.expression_t <- trans(dfwc_between$Norm.expression)
  dfwc_between$se_up_t <- trans(dfwc_between$Norm.expression + dfwc_between$se)
  dfwc_between$se_low_t <- trans(dfwc_between$Norm.expression - dfwc_between$se)
  
  p <- (ggplot(dfwc_between, aes(x=Cell.line,
                                 y=Norm.expression_t,
                                 fill=Cell.line,
                                 color=Cell.line,
                                 ymin=se_low_t,
                                 ymax=se_up_t))
        
        + geom_errorbar(width=0.33, size=0.3)
        + geom_bar(stat="identity", width=0.5)
        
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


################################################################################
plotViralLoadCellLines <- function(data, title) {
  dfwc_between <- summarySE(data=data, measurevar="Viral Load", groupvars=c("Cell.line", "Time", "MOI"), na.rm=FALSE, conf.interval=.95)
  dfwc_between$logVL <- log10(dfwc_between$`Viral Load`)
  dfwc_between$logVL[which(is.infinite(dfwc_between$logVL))] <- 0
  
  data.other <- filter(dfwc_between, !(Cell.line %in% c("H522", "Vero E6")))
  data.VeroE6 <- filter(dfwc_between, Cell.line == "Vero E6")
  data.H522 <- filter(dfwc_between, Cell.line == "H522")
  
  p <- (ggplot(dfwc_between, aes(x=Time, 
                                 y=logVL, 
                                 group=Cell.line,
                                 color=Cell.line,
                                 facet=MOI,
                                 ymin=log10(`Viral Load` - se),
                                 ymax=log10(`Viral Load` + se)))
        
        # separated to draw them in the desired z-order
        + geom_errorbar(data=data.other, width=4, size=0.3)
        + geom_line(data=data.other, size=0.5, alpha=0.5)
        + geom_point(data=data.other, size=1.5)
        
        + geom_errorbar(data=data.VeroE6, width=4, size=0.3)
        + geom_line(data=data.VeroE6, size=0.5)
        + geom_point(data=data.VeroE6, size=1.5)
        
        + geom_errorbar(data=data.H522, width=4, size=0.3)
        + geom_line(data=data.H522, size=0.5)
        + geom_point(data=data.H522, size=1.5)

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


################################################################################
plotViralLoad <- function(data) {
  dfwc_between <- summarySE(data=data, measurevar="Viral Load", groupvars=c("MOI","Time"), na.rm=FALSE, conf.interval=.95)
  dfwc_between$logVL <- log10(dfwc_between$`Viral Load`)
  dfwc_between$logVL[which(is.infinite(dfwc_between$logVL))] <- 0
  
  p <- (ggplot(dfwc_between, aes(x=Time, 
                                 y=logVL, 
                                 group=MOI,
                                 color=as.factor(MOI),
                                 ymin=log10(`Viral Load` - se),
                                 ymax=log10(`Viral Load` + se)))
        
        + geom_errorbar(width=3, size=0.3)
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


################################################################################
plotViralLoadSups <- function(data, y.axis.label) {
  dfwc_between <- summarySE(data=data, measurevar="Viral Load", groupvars=c("MOI","Time"), na.rm=FALSE, conf.interval=.95)
  dfwc_between$logVL <- log10(dfwc_between$`Viral Load`)
  dfwc_between$logVL[which(is.infinite(dfwc_between$logVL))] <- 0
  
  p <- (ggplot(dfwc_between, aes(x=Time, 
                                 y=logVL, 
                                 group=MOI,
                                 color=as.factor(MOI),
                                 ymin=log10(`Viral Load` - se),
                                 ymax=log10(`Viral Load` + se)))
        
        + geom_errorbar(width=3, size=0.3)
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



################################################################################
# Read data
################################################################################
ACE2.PCR <- read_tsv(here("data/Figure 1/ACE2_PCR_all_lines.txt"))
TMPRSS2.PCR <- read_tsv(here("data/Figure 1/TMPRSS2_PCR_all_lines.txt"))
FACS.H522 <- read_csv(here("data/Figure 1/FACS H522_EXP08 and EXP09.csv"))
Viral.load.H522 <- read_csv(here("data/Figure 1/Viral Load H522_EXP08 and EXP09.csv"))
Viral.load.insups.H522 <- read_csv(here("data/Figure 1/Viral Load-insups H522 EXP08 and EXP09.csv"))
Viral.load.cell.lines <- read_csv(here("data/Figure 1/Viral Load-cell lines - all HCC827.csv"))
################################################################################



################################################################################
# Generate figures
################################################################################
panel.PCR.ACE2 <- plotPCR.ACE2(ACE2.PCR, "ACE2 mRNA")
panel.PCR.TMPRSS2 <- plotPCR.TMPRSS2(TMPRSS2.PCR, "TMPRSS2 mRNA")
panel.ViralLoad.AllLines <- plotViralLoadCellLines(Viral.load.cell.lines, "Viral load in cells - All lines")
panel.ViralLoad.H522 <- plotViralLoad(Viral.load.H522)
panel.ViralLoad.H522.insups <- plotViralLoadSups(Viral.load.insups.H522)
panel.FACS.H522 <- plotFACS(FACS.H522)

arranged.PCR <- arrangeGrob(
  panel.PCR.ACE2, panel.PCR.TMPRSS2,
  nrow = 1,
  ncol = 2)

arranged.ViralLoad <- arrangeGrob(
  panel.ViralLoad.AllLines, panel.ViralLoad.H522, panel.ViralLoad.H522.insups,
  panel.FACS.H522,
  nrow = 2,
  ncol = 3)

#grid.draw(F1) # to view the plot
saveFig(arranged.PCR, "Figure1_B", 2, 4)
saveFig(arranged.ViralLoad, "Figure1_DEFG", 4, 6.85)
################################################################################
