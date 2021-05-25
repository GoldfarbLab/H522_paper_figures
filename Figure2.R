library(here)
source(here("common.R"))

################################################################################
# Plotting functions
################################################################################
plotAB <- function(data, title, x.axis, y.axis) {
  data$Cell.line <- factor(data$Cell.line, levels=c("Vero E6", "H522", "H522-ACE2"))
  dfwc_between <- summarySE(data=data, measurevar="Viral RNA", groupvars=c("Concentration", "Cell.line"), na.rm=FALSE, conf.interval=.95)
  
  p <- (ggplot(data, aes(x=as.factor(Concentration),
                         y=`Viral RNA`,
                         color=Cell.line,
                         fill=Cell.line,
                         facet=Cell.line))

        #+ stat_summary(fun = "mean", size= 0.15, geom = "bar", width=0.75)
        #+ geom_boxplot(width=0.65, outlier.shape=NA, size=0.4)
        + geom_errorbar(data=dfwc_between, aes(ymin = `Viral RNA` - se,
                                               ymax = `Viral RNA` + se), width=0.33, size=0.3)
        + stat_summary(fun = "mean", size= 0.15, geom = "bar", width=0.75, fill=lighter.grey, color=medium.grey)
        #+ stat_summary(fun = "mean", size= 0.15, geom = "crossbar", width=0.75)
        + geom_jitter(size=0.75, shape=21, stroke = 0.4, width=0.2, fill="white")
        

        + ggtitle(title)
        
        + scale_x_discrete(name = x.axis)
        + scale_y_continuous(name = y.axis)

        + scale_color_manual(name = "Cell.line", values=colors.Cell.line)
        + scale_fill_manual(name = "Cell.line", values=colors.Cell.line)
        
        + facet_wrap(~ Cell.line, scales="free")
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 0.7,
                strip.text.x = element_text(size = 6))#angle=45, vjust=1, hjust=1))
  )
}
################################################################################


################################################################################
plotIncucyte <- function(data, doBreak=F) {
  data$Cell.line <- factor(data$Cell.line, levels=c("Vero E6", "H522", "H522-ACE2", "Basal AEC", "Basal AEC-ACE2"))
  data$Cells <- pmin(data$Cells, 100)
  dfwc_between <- summarySE(data=data, measurevar="Cells", groupvars=c("MOI", "Time", "Cell.line"), na.rm=FALSE, conf.interval=.95)
  
  trans <- function(x){x}
  yticks <- c(0, 25, 50, 75, 100)
  
  if (doBreak) {
    trans <- function(x){pmin(x,0.75) + 0.01*pmax(x-0.75,0)}
    yticks <- c(0, 0.25, 0.5, 0.75, 50, 75, 100)
  }
  
  dfwc_between$se_up <- dfwc_between$Cells + dfwc_between$se
  dfwc_between$se_low <- dfwc_between$Cells - dfwc_between$se
  
  dfwc_between$Cells_t <- trans(dfwc_between$Cells)
  dfwc_between$se_up_t <- trans(dfwc_between$se_up)
  dfwc_between$se_low_t <- trans(dfwc_between$se_low)
  
  p <- (ggplot(dfwc_between, aes(x=Time,
                                 y=Cells_t,
                                 fill=Cell.line,
                                 color=as.factor(MOI),
                                 group=MOI,
                                 facet=Cell.line))
        + geom_ribbon(aes(ymin=se_up_t, ymax=se_low_t),
                      fill=lightest.grey, 
                      color=lightest.grey)
        
        + geom_line(size=0.5)
        
        #+ geom_rect(aes(xmin=0, xmax=48, ymin=1.1, ymax=49), fill="white", color="white")
        
        + scale_x_continuous(name = "Hours post-infection", 
                             breaks = c(0, 12, 24, 36, 48),
                             limits = c(0,48))
        + scale_y_continuous(name = "GFP+ cells (%)", limits=c(0,NA), breaks=trans(yticks), labels=yticks)

        + scale_color_manual(name = "MOI", values=colors.MOI.6)
        
        + facet_wrap(Cell.line ~ ., scales = "free_x", ncol=1)
        
        + theme.basic
        + theme(#legend.position = "none",
                strip.background = element_blank(),
                strip.text.x = element_blank(),
                aspect.ratio=1/2)
  )
}
################################################################################

################################################################################
plotMutantIncucyte <- function(data, doBreak=F) {
  data <- data %>% filter(MOI == 20, Time <= 48)
  data$Cell.line <- factor(data$Cell.line, levels=c("WT", "E484D", "E484K/R685S"))
  #data$Cells <- pmin(data$Cells, 100)
  dfwc_between <- summarySE(data=data, measurevar="Cells", groupvars=c("MOI", "Time", "Cell.line"), na.rm=FALSE, conf.interval=.95)
  
  trans <- function(x){x}
  yticks <- c(0, 25, 50, 75, 100)
  
  if (doBreak) {
    trans <- function(x){pmin(x,0.75) + 0.01*pmax(x-0.75,0)}
    yticks <- c(0, 0.25, 0.5, 0.75, 50, 75, 100)
  }
  
  dfwc_between$se_up <- dfwc_between$Cells + dfwc_between$se
  dfwc_between$se_low <- dfwc_between$Cells - dfwc_between$se
  
  dfwc_between$Cells_t <- trans(dfwc_between$Cells)
  dfwc_between$se_up_t <- trans(dfwc_between$se_up)
  dfwc_between$se_low_t <- trans(dfwc_between$se_low)
  
  p <- (ggplot(dfwc_between, aes(x=Time,
                                 y=Cells_t,
                                 fill=Cell.line,
                                 color=Cell.line,
                                 group=Cell.line))
                                 #facet=Cell.line))
        #+ geom_ribbon(aes(ymin=se_up_t, ymax=se_low_t),
        #              fill=lightest.grey, 
        #              color=lightest.grey)
        
        + geom_smooth(size=0.5, fill=light.grey)
        
        + scale_x_continuous(name = "Hours post-infection", 
                             breaks = c(0, 12, 24, 36, 48),
                             limits = c(0,48))
        
        + scale_y_continuous(name = "Normalized GFP intensity",
                             breaks = c(0, 2e5, 4e5, 6e5, 8e5, 1e6),
                             labels = c("0", "2e5", "4e5", "6e5", "8e5", "1e6"),
                             limits = c(-1e5, 1.1e6))
        
        #+ scale_y_continuous(name = "GFP intensity / cell confluence", limits=c(0,NA), breaks=trans(yticks), labels=yticks)
        
        + scale_color_manual(name = "Cell line", values=colors.mutant.lines)
        
        #+ facet_wrap(Cell.line ~ ., scales = "free_x", ncol=1)
        
        + theme.basic
        + theme(#legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          aspect.ratio=0.6)
  )
}
################################################################################


################################################################################
plotLentiMutants <- function(data) {
  data$Condition <- factor(data$Condition, levels=c("WT", "E484D", "R682W", "E484D/R682W"))
  dfwc_between <- summarySE(data=data, measurevar="GFP", groupvars=c("Condition"), na.rm=FALSE, conf.interval=.95)
  
  p <- (ggplot(data, aes(x=Condition,
                         y=`GFP`))
        
        + geom_errorbar(data=dfwc_between, aes(ymin = GFP - se,
                                               ymax = GFP + se), width=0.33, size=0.3, color=colors.Cell.line["H522"])
        + stat_summary(fun = "mean", size= 0.15, geom = "bar", width=0.65, fill=lighter.grey, color=medium.grey)
        
        #+ stat_summary(fun = "mean", size= 0.15, geom = "crossbar", width=0.75, color=colors.Cell.line["H522"])
        + geom_jitter(size=0.75, shape=21, stroke = 0.4, width=0.2, fill="white", color=colors.Cell.line["H522"])
        
        
        + ggtitle("Lenti Spike mutants in H522 cells")
        
        + scale_y_continuous(name = "% GFP positive")
        
        #+ scale_color_manual(name = "Condition", values=colors.Cell.line)
        #+ scale_fill_manual(name = "Condition", values=colors.Cell.line)
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 0.7,
                axis.title.x = element_blank(),
                axis.text.x = element_text(size = 6))#, angle=45, vjust=1, hjust=1))
  )
}
################################################################################


################################################################################
plot293TLentiMutants <- function(data) {
  data$Condition <- factor(data$Condition, levels=c("WT", "E484D", "R682W", "E484D/R682W"))
  dfwc_between <- summarySE(data=data, measurevar="GFP", groupvars=c("Condition", "Amount"), na.rm=FALSE, conf.interval=.95)
  
  p <- (ggplot(data, aes(x=as.factor(Amount),
                         y=`GFP`))
        
        + geom_errorbar(data=dfwc_between, aes(ymin = GFP - se,
                                               ymax = GFP + se), width=0.33, size=0.3, color=colors.Cell.line["293T"])
        + stat_summary(fun = "mean", size= 0.15, geom = "bar", width=0.6, fill=lighter.grey, color=medium.grey)
        
        #+ stat_summary(fun = "mean", size= 0.15, geom = "crossbar", width=0.55, color=colors.Cell.line["H522"])
        + geom_jitter(size=0.75, shape=21, stroke = 0.4, width=0.2, fill="white", color=colors.Cell.line["293T"])
        
        
        + ggtitle("Lenti Spike mutants in 293T cells")
        
        + xlab("Volume (ul)")
        + scale_y_continuous(name = "% GFP positive", limits=c(0,70))
        
        + facet_wrap(~ Condition, nrow=1, scales="free_x")
        
        #+ scale_color_manual(name = "Condition", values=colors.Cell.line)
        #+ scale_fill_manual(name = "Condition", values=colors.Cell.line)
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 1.2,
                strip.text.x = element_blank(), #element_text(size = 6),
                #axis.title.x = element_blank(),
                axis.text.x = element_text(size = 6))#, angle=45, vjust=1, hjust=1))
  )
}
################################################################################


################################################################################
plot293TVSVMutants <- function(data) {
  data$Condition <- factor(data$Condition, levels=c("WT", "E484D", "E484K/R685S"))
  dfwc_between <- summarySE(data=data, measurevar="Intensity", groupvars=c("Condition"), na.rm=FALSE, conf.interval=.95)
  
  p <- (ggplot(data, aes(x=Condition,
                         y=Intensity))
        
        + geom_errorbar(data=dfwc_between, aes(ymin = Intensity - se,
                                               ymax = Intensity + se), width=0.33, size=0.3, color=colors.Cell.line["293T"])
        + stat_summary(fun = "mean", size= 0.15, geom = "bar", width=0.5, fill=lighter.grey, color=medium.grey)
        
        #+ stat_summary(fun = "mean", size= 0.15, geom = "crossbar", width=0.55, color=colors.Cell.line["H522"])
        + geom_jitter(size=0.75, shape=21, stroke = 0.4, width=0.2, fill="white", color=colors.Cell.line["293T"])
        
        
        + ggtitle("VSV Spike mutants in 293T cells")
        
        + scale_y_continuous(name = "Relative infectivity", 
                             breaks = seq(0,2,0.5),
                             limits=c(0,2))
       
        
        #+ scale_color_manual(name = "Condition", values=colors.Cell.line)
        #+ scale_fill_manual(name = "Condition", values=colors.Cell.line)
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 0.7,
                strip.text.x = element_blank(), #element_text(size = 6),
                #axis.title.x = element_blank(),
                axis.text.x = element_text(size = 6))#, angle=45, vjust=1, hjust=1))
  )
}
################################################################################



################################################################################
plotH522TVSVMutants <- function(data) {
  data$Condition <- factor(data$Condition, levels=c("WT", "E484D", "E484K/R685S"))
  dfwc_between <- summarySE(data=data, measurevar="GFP", groupvars=c("Condition", "MOI"), na.rm=FALSE, conf.interval=.95)
  
  p <- (ggplot(data, aes(x=as.factor(MOI),
                         y=GFP))
        
        + geom_errorbar(data=dfwc_between, aes(ymin = GFP - se,
                                               ymax = GFP + se), width=0.33, size=0.3, color=colors.Cell.line["H522"])
        + stat_summary(fun = "mean", size= 0.15, geom = "bar", width=0.75, fill=lighter.grey, color=medium.grey)
        
        #+ stat_summary(fun = "mean", size= 0.15, geom = "crossbar", width=0.55, color=colors.Cell.line["H522"])
        + geom_jitter(size=0.75, shape=21, stroke = 0.4, width=0.2, fill="white", color=colors.Cell.line["H522"])
        
        
        + ggtitle("VSV-Spike mutants in H522 cells")
        
        + scale_y_continuous(name = "GFP+ cells (%)", 
                             breaks = seq(0,6,2),
                             limits=c(0,6))
        + xlab("MOI")
        
        #+ scale_color_manual(name = "Condition", values=colors.Cell.line)
        #+ scale_fill_manual(name = "Condition", values=colors.Cell.line)
        
        + facet_wrap(~ Condition, nrow=1)
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 1.0,
                strip.text.x = element_blank(), #element_text(size = 6),
                #axis.title.x = element_blank(),
                axis.text.x = element_text(size = 6))#, angle=45, vjust=1, hjust=1))
  )
}
################################################################################



################################################################################
plotH522VirusStocks  <- function(data)
{
  dfwc_between <- summarySE(data=data, measurevar="Copies", groupvars=c("Time", "Condition"), na.rm=FALSE, conf.interval=.95, removeOutliers=F, logt=T)
  dfwc_between$logVL <- log10(dfwc_between$Copies)
  dfwc_between$logVL[which(is.infinite(dfwc_between$logVL))] <- 0
  
  p <- (ggplot(dfwc_between, aes(x=Time, 
                                 y=logVL, 
                                 group=Condition,
                                 color=colors.Cell.line["H522"]))
        
        + geom_errorbar(width=7, size=0.3,
                        aes(ymin=log10(`Copies` - se),
                            ymax=log10(`Copies` + se)))
        + geom_line(aes(linetype = Condition), size=0.5)
        + geom_point(size=0.75)
        
        
        + scale_x_continuous(name = "Hours post-infection", 
                             breaks = c(4, 96),
                             labels = c("4","96"),
                             expand = c(0.25, 0))
        + scale_y_continuous(name = "Viral RNA (copies/cell)",
                             breaks = c(2, 3, 4, 5, 6),
                             labels = c("1e2", "1e3", "1e4", "1e5", "1e6"),
                             limits = c(2, 6))
        
        + scale_linetype_manual(name = "Condition", values=c("E484D" = "longdash", "WT" = "solid"))
        
        + ggtitle("Virus Stocks in H522 cells")
        
        + theme.basic
        + theme(legend.position = "none",
                strip.text.x = element_text(size = 6),
                aspect.ratio=golden.ratio)
  )
}
################################################################################

################################################################################
plotVeroNAB <- function(data) {
  data$Condition <- factor(data$Condition, levels=c("Mock", "2B04", "2H04", "SARS2-02", "SARS2-38", 
                                                    "SARS2-71", "SARS2-31", "SARS2-57", "SARS2-11"))
  
  
  dfwc_between <- summarySE(data=data, measurevar="Copies", groupvars=c("Condition"), na.rm=FALSE, conf.interval=.95)
  dfwc_between$logVL <- log10(dfwc_between$Copies)
  dfwc_between$logVL[which(is.infinite(dfwc_between$logVL))] <- 0
  
  p <- (ggplot(dfwc_between, aes(x=as.factor(Condition),
                                 y=logVL))
        
        + geom_errorbar(aes(ymin=log10(`Copies` - se),
                            ymax=log10(`Copies` + se)), width=0.33, size=0.3, color=colors.Cell.line["Vero E6"])
        + stat_summary(fun = "mean", size= 0.15, geom = "bar", width=0.65, fill=lighter.grey, color=medium.grey)
        
        + geom_jitter(data=data, aes(y=log10(Copies)), size=0.75, shape=21, stroke = 0.4, width=0.25, fill="white", color=colors.Cell.line["Vero E6"])
        
        + ggtitle("mAB neutralization in Vero E6 cells")
        
        + scale_y_continuous(name = "Viral RNA (copies/cell)",
                             breaks = c(0, 2, 4, 6, 8),
                             labels = c("1e0", "1e2", "1e4", "1e6", "1e8"),
                             limits = c(0, 8))
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 0.4,
                strip.text.x = element_text(size = 6),
                axis.title.x = element_blank(),
                plot.subtitle = element_text(size = 6, hjust=0.5))
  )
}
################################################################################

################################################################################
plotH522NAB <- function(data) {
  data$Condition <- factor(data$Condition, levels=c("Mock", "2B04", "2H04", "SARS2-02", "SARS2-38", 
                                                    "SARS2-71", "SARS2-31", "SARS2-57", "SARS2-11"))
  
  
  dfwc_between <- summarySE(data=data, measurevar="Copies", groupvars=c("Condition"), na.rm=FALSE, conf.interval=.95)
  dfwc_between$logVL <- log10(dfwc_between$Copies)
  dfwc_between$logVL[which(is.infinite(dfwc_between$logVL))] <- 0
  
  p <- (ggplot(dfwc_between, aes(x=as.factor(Condition),
                                 y=logVL))
        
        + geom_errorbar(aes(ymin=log10(`Copies` - se),
                            ymax=log10(`Copies` + se)), width=0.33, size=0.3, color=colors.Cell.line["H522"])
        + stat_summary(fun = "mean", size= 0.15, geom = "bar", width=0.65, fill=lighter.grey, color=medium.grey)
        
        + geom_jitter(data=data, aes(y=log10(Copies)), size=0.75, shape=21, stroke = 0.4, width=0.25, fill="white", color=colors.Cell.line["H522"])
        
        + ggtitle("mAB neutralization in Vero E6 cells")
        
        + scale_y_continuous(name = "Viral RNA (copies/cell)",
                             breaks = c(0, 1, 2, 3, 4, 5),
                             labels = c("1e0", "1e1", "1e2", "1e3", "1e4", "1e5"),
                             limits = c(0, 5.4))
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 0.4,
                strip.text.x = element_text(size = 6),
                axis.title.x = element_blank(),
                plot.subtitle = element_text(size = 6, hjust=0.5))
  )
}
################################################################################







################################################################################
# Read data
################################################################################
data.neutralizing <- read_csv(here("data/Figure 2/S_neutralizing.csv"))
data.fc <- read_csv(here("data/Figure 2/ACE2_Fc.csv"))
data.incucyte <- read_csv(here("data/Figure 2/incucyte.csv"))
data.mutant.incucyte <- read_csv(here("data/Figure 2/mutant_incucyte.csv"))
data.H522.lenti.mutants <- read_csv(here("data/Figure 2/H522_lenti_mutants.csv"))
data.293T.lenti.mutants <- read_csv(here("data/Figure 2/293T_lenti_mutants.csv"))
data.293T.ACE2.VSV <- read_csv(here("data/Figure 2/293-ACE2-VSV.csv"))
data.H522.VSV <- read_csv(here("data/Figure 2/H522_VSV.csv"))
data.H522.virus.stocks <- read_csv(here("data/Figure 2/H522_virus_stocks.csv"))
data.Vero.nAB <- read_csv(here("data/Figure 2/Vero_nAB.csv"))
data.H522.nAB<- read_csv(here("data/Figure 2/H522_nAB.csv"))
################################################################################



################################################################################
# Generate figures
################################################################################
panel.S.AB <- plotAB(data.neutralizing, "S neutralizing antibody", "2B04 (ug/mL)", "Viral RNA (%)")
panel.ACE2.Fc <- plotAB(data.fc, "ACE2 Fc", "hACE2 Fc (ug/mL)", "Viral RNA (%)")
# Execute this with doBreak=T to get y-axis break
panel.incucyte <- plotIncucyte(data.incucyte, doBreak=F) 
panel.mutant.incucyte <- plotMutantIncucyte(data.mutant.incucyte, doBreak=F) 
panel.H522.lenti.mutants <- plotLentiMutants(data.H522.lenti.mutants)
panel.293T.lenti.mutants <- plot293TLentiMutants(data.293T.lenti.mutants)
panel.293T.ACE2.VSV <- plot293TVSVMutants(data.293T.ACE2.VSV)
panel.H522.VSV <- plotH522TVSVMutants(data.H522.VSV)
panel.H522.virus.stocks <- plotH522VirusStocks(data.H522.virus.stocks)
panel.Vero.nAB <- plotVeroNAB(data.Vero.nAB)
panel.H522.nAB <- plotH522NAB(data.H522.nAB)

arranged.AB.Fc <- arrangeGrob(
  panel.S.AB, panel.ACE2.Fc,
  nrow = 2,
  ncol = 1)

arranged.nAB <- arrangeGrob(
  panel.Vero.nAB,
  panel.H522.nAB,
  nrow=1,
  ncol=2)

saveFig(arranged.AB.Fc, "Figure2_BC", 3, 3.7)
saveFig(panel.incucyte, "Figure2_D", 4.1, 6.85)
saveFig(panel.H522.lenti.mutants, "Figure2_E", 1.1, 6.85)
saveFig(panel.293T.lenti.mutants, "Figure2_F", 1.2, 6.85)
saveFig(panel.mutant.incucyte, "Figure2_G", 1.0, 6.85)
saveFig(panel.293T.ACE2.VSV, "Figure2_H", 1.2, 6.85)
saveFig(panel.H522.VSV, "Figure2_I", 1.2, 6.85)
saveFig(panel.H522.virus.stocks, "Figure2_J", 1.2, 6.85)
saveFig(arranged.nAB, "Figure2_KL", 1.15, 5.00)
################################################################################
