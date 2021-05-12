library(here)
source(here("common.R"))

################################################################################
# Plotting functions
################################################################################
plotBlocking <- function(data, title, x.axis, y.axis) {
  data$Condition <- factor(data$Condition, levels=c("Mock", "Anti-ACE2", "Anti-DC SIGN", "Anti-GFP"))
  data$Cell.line <- factor(data$Cell.line, levels=c("PgsA", "H522", "Calu-3"))
  dfwc_between <- summarySE(data=data, measurevar="Copies", groupvars=c("Condition", "Cell.line"), na.rm=FALSE, conf.interval=.95)
  
  anno <- tibble(x1=c(1,1,1, 1,1,1),
                 x2=c(2,3,4, 2,3,4),
                 y1=c(3.5,2.7,1.9, 1,4.7,5.5),
                 y2=c(3.5,2.7,1.9, 1,4.7,5.5),
                 xstar=c(1.5,2,2.5, 1.5,2,2.5),
                 ystar=c(3.2,2.4,1.6, 0.6,5.2,6),
                 lab=c(rep("ns",3), "***", rep("ns",2)),
                 Cell.line=factor(c(rep("H522",3),rep("Calu-3",3)), levels=c("PgsA", "H522", "Calu-3")))
  
  p <- (ggplot(data, aes(x=as.factor(Condition),
                         y=log10(Copies),
                         color=Cell.line,
                         fill=Cell.line,
                         facet=Cell.line))
        
        #+ stat_summary(fun = "mean", size= 0.15, geom = "bar", width=0.5)
        #+ geom_boxplot(width=0.65, outlier.shape=NA, size=0.4)
        + stat_summary(fun = "mean", size= 0.15, geom = "crossbar", width=0.75)
        + geom_jitter(size=0.75, shape=21, stroke = 0.4, width=0.25, fill="white")
        #+ geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab), size=2, color=medium.grey, family = "Arial")
        #+ geom_segment(data = anno, aes(x = x1, xend = x2, y = y2, yend = y2), colour=medium.grey)
        
        
        + ggtitle(title)
        
        + scale_y_continuous(name = y.axis,
                             breaks = c(0, 1, 2, 3, 4, 5, 6),
                             labels = c("1e0", "1e1", "1e2", "1e3", "1e4", "1e5", "1e6"),
                             limits = c(1,6))
        
        + scale_color_manual(name = "Cell.line", values=colors.Cell.line)
        + scale_fill_manual(name = "Cell.line", values=colors.Cell.line)
        
        + facet_wrap(~ Cell.line, scales="free")
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 0.7,
                strip.text.x = element_text(size = 6),
                axis.title.x = element_blank(),
                plot.subtitle = element_text(size = 6, hjust=0.5),
                axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  )
}
################################################################################


################################################################################
plotViralLoadCRISPR <- function(data, title) {
  #data <- data %>% filter(is.na(Outlier))
  dfwc_between <- summarySE(data=data, measurevar="Copies", groupvars=c("Cell.line", "Time", "KO"), na.rm=FALSE, conf.interval=.95, removeOutliers=F, logt=T)
  dfwc_between$logVL <- log10(dfwc_between$Copies)
  dfwc_between$logVL[which(is.infinite(dfwc_between$logVL))] <- 0
  
  p <- (ggplot(dfwc_between, aes(x=Time, 
                                 y=logVL, 
                                 group=interaction(Cell.line, KO),
                                 color=Cell.line))
        
        + geom_errorbar(width=4, size=0.3,
                        aes(ymin=log10(`Copies` - se),
                            ymax=log10(`Copies` + se)))
        + geom_line(aes(linetype = KO), size=0.5)
        + geom_point(size=0.75)
        
        
        + scale_x_continuous(name = "Hours post-infection", 
                             breaks = c(4, 72),
                             labels = c("4","72"),
                             expand = c(0.25, 0))
        + scale_y_continuous(name = "Viral RNA (copies/cell)",
                             breaks = c(3, 4, 5, 6),
                             labels = c("1e3", "1e4", "1e5", "1e6"),
                             limits = c(2.6, 6.3))
        
        + scale_color_manual(name = "Cell line", values=colors.Cell.line)
        + scale_linetype_manual(name = "KO", values=c("ACE2" = "longdash", "WT" = "solid"))
        
        + ggtitle(title)
        
        + theme.basic
        + theme(legend.position = "none")
  )
}
################################################################################


################################################################################
plotViralLoadCRISPRMono <- function(data, title) {
  dfwc_between <- summarySE(data=data, measurevar="Copies", groupvars=c("Cell.line", "Time", "KO", "Clone"), na.rm=FALSE, conf.interval=.95, removeOutliers=F, logt=T)
  dfwc_between$logVL <- log10(dfwc_between$Copies)
  dfwc_between$logVL[which(is.infinite(dfwc_between$logVL))] <- 0
  
  p <- (ggplot(dfwc_between, aes(x=Time, 
                                 y=logVL, 
                                 group=interaction(Clone, KO),
                                 color=Cell.line))
        
        + geom_errorbar(width=7, size=0.3,
                        aes(ymin=log10(`Copies` - se),
                            ymax=log10(`Copies` + se)))
        + geom_line(aes(linetype = KO), size=0.5)
        + geom_point(size=0.75)
        
        
        + scale_x_continuous(name = "Hours post-infection", 
                             breaks = c(4, 72),
                             labels = c("4","72"),
                             expand = c(0.25, 0))
        + scale_y_continuous(name = "Viral RNA (copies/cell)",
                             breaks = c(3, 4, 5, 6, 7),
                             labels = c("1e3", "1e4", "1e5", "1e6", "1e7"),
                             limits = c(2.5, 7.2))
        
        + scale_color_manual(name = "Cell line", values=colors.Cell.line)
        + scale_linetype_manual(name = "KO", values=c("ACE2" = "longdash", "WT" = "solid", "Mixed" = "dotdash"))
        
        + facet_grid(~ KO)
        
        + ggtitle(title)
        
        + theme.basic
        + theme(legend.position = "none",
                strip.text.x = element_text(size = 6),
                aspect.ratio=0.8)
  )
}
################################################################################


################################################################################
plotCRISPRBlocking <- function(data, title, x.axis, y.axis) {
  data$Condition <- factor(data$Condition, levels=c("Mock", "Anti-ACE2", "Anti-GFP"))
  data$Cell.line <- factor(data$Cell.line, levels=c("H522", "H522 ACE2 KO", "Calu-3", "Calu-3 ACE2 KO"))
  dfwc_between <- summarySE(data=data, measurevar="Copies", groupvars=c("Condition", "Cell.line"), na.rm=FALSE, conf.interval=.95)
  
  p <- (ggplot(data, aes(x=as.factor(Condition),
                         y=log10(Copies),
                         color=Cell.line,
                         fill=Cell.line,
                         facet=Cell.line))
        
        #+ stat_summary(fun = "mean", size= 0.15, geom = "bar", width=0.5)
        #+ geom_boxplot(width=0.65, outlier.shape=NA, size=0.4)
        + stat_summary(fun = "mean", size= 0.15, geom = "crossbar", width=0.75)
        + geom_jitter(size=0.75, shape=21, stroke = 0.4, width=0.25, fill="white")
        #+ geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab), size=2, color=medium.grey, family = "Arial")
        #+ geom_segment(data = anno, aes(x = x1, xend = x2, y = y2, yend = y2), colour=medium.grey)
        
        
        + ggtitle(title)
        
        + scale_y_continuous(name = y.axis,
                             breaks = c(2, 3, 4, 5, 6),
                             labels = c("1e2", "1e3", "1e4", "1e5", "1e6"),
                             limits = c(2,6))
        
        + scale_color_manual(name = "Cell.line", values=colors.Cell.line)
        + scale_fill_manual(name = "Cell.line", values=colors.Cell.line)
        
        + facet_wrap(~ Cell.line, scales="free_x", nrow=1)
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 0.8,
                strip.text.x = element_text(size = 6),
                axis.title.x = element_blank(),
                plot.subtitle = element_text(size = 6, hjust=0.5),
                axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  )
}
################################################################################


################################################################################
plotHeparanSulfate  <- function(data, title)
{
  dfwc_between <- summarySE(data=data, measurevar="Copies", groupvars=c("Condition", "Time"), na.rm=FALSE, conf.interval=.95, removeOutliers=F, logt=T)
  dfwc_between$logVL <- log10(dfwc_between$Copies)
  dfwc_between$logVL[which(is.infinite(dfwc_between$logVL))] <- 0
  
  p <- (ggplot(dfwc_between, aes(x=Time, 
                                 y=logVL, 
                                 group=interaction(Condition),
                                 color=Condition))
        
        + geom_errorbar(width=4, size=0.3,
                        aes(ymin=log10(`Copies` - se),
                            ymax=log10(`Copies` + se)))
        + geom_line(size=0.5)
        + geom_point(size=0.75)
        
        
        + scale_x_continuous(name = "Hours post-infection", 
                             breaks = c(4, 72),
                             labels = c("4","72"),
                             expand = c(0.25, 0))
        + scale_y_continuous(name = "Viral RNA (copies/cell)",
                             breaks = c(3, 4, 5, 6, 7),
                             labels = c("1e3", "1e4", "1e5", "1e6", "1e7"),
                             limits = c(2.8, 7.5))
        
        + scale_color_manual(name = "Condition", values=colors.heparan)
        
        + ggtitle(title)
        
        + theme.basic
        + theme(legend.position = "none")
  )
}
################################################################################

################################################################################
plotNRP1Monoclones  <- function(data, title)
  dfwc_between <- summarySE(data=data, measurevar="Copies", groupvars=c("Cell.line", "Time", "KO", "Clone"), na.rm=FALSE, conf.interval=.95, removeOutliers=F, logt=T)
  dfwc_between$logVL <- log10(dfwc_between$Copies)
  dfwc_between$logVL[which(is.infinite(dfwc_between$logVL))] <- 0
  
  p <- (ggplot(dfwc_between, aes(x=Time, 
                                 y=logVL, 
                                 group=interaction(Clone, KO),
                                 color=Cell.line))
        
        + geom_errorbar(width=7, size=0.3,
                        aes(ymin=log10(`Copies` - se),
                            ymax=log10(`Copies` + se)))
        + geom_line(aes(linetype = KO), size=0.5)
        + geom_point(size=0.75)
        
        
        + scale_x_continuous(name = "Hours post-infection", 
                             breaks = c(4, 72),
                             labels = c("4","72"),
                             expand = c(0.25, 0))
        + scale_y_continuous(name = "Viral RNA (copies/cell)",
                             breaks = c(3, 4, 5, 6, 7),
                             labels = c("1e3", "1e4", "1e5", "1e6", "1e7"),
                             limits = c(2.5, 7.2))
        
        + scale_color_manual(name = "Cell line", values=colors.Cell.line)
        + scale_linetype_manual(name = "KO", values=c("ACE2" = "longdash", "WT" = "solid", "Mixed" = "dotdash"))
        
        + facet_grid(~ KO)
        
        + ggtitle(title)
        
        + theme.basic
        + theme(legend.position = "none",
                strip.text.x = element_text(size = 6),
                aspect.ratio=0.8)
  )
}
################################################################################



################################################################################
# Read data
################################################################################
data.blockingAB <- read_csv(here("data/Figure 3/ACE2_blocking.csv"))
data.ace2.CRISPR <- read_csv(here("data/Figure 3/ACE2_CRISPR_BULK.csv"))
data.ace2.CRISPR.mono <- read_csv(here("data/Figure 3/ACE2_CRISPR_Monoclonals.csv"))
data.ace2.CRISPR.blockingAB <- read_csv(here("data/Figure 3/ACE2_CRISPR_blocking.csv"))
data.NRP1.AXL.siRNA <- read_csv(here("data/Figure 3/AXL_NRP1_siRNA.csv"))
data.NRP1.monoclones <- read_csv(here("data/Figure 3/NRP1_clones.csv"))
data.heparan.sulfate <- read_csv(here("data/Figure 3/Heparan_sulfate.csv"))
################################################################################



################################################################################
# Generate figures
################################################################################
panel.blockingAB <- plotBlocking(data.blockingAB, "Blocking antibody", "3 days post-infection", "Viral RNA (copies/cell)")
panel.ViralLoad.CRISPR.bulk <- plotViralLoadCRISPR(data.ace2.CRISPR, "ACE2 CRISPR KO cells")
panel.ViralLoad.CRISPR.monoclonal <- plotViralLoadCRISPRMono(data.ace2.CRISPR.mono, "Monoclonal ACE2 CRISPR KO cells")
panel.blockingAB.CRISPR <- plotCRISPRBlocking(data.ace2.CRISPR.blockingAB, "Blocking antibody in ACE2 CRISPR KO cells", "3 days post-infection", "Viral RNA (copies/cell)")
#panel.NRP1.AXL.siRNA <- plotSiRNA()
panel.NRP1.monoclones <- plotNRP1Monoclones(data.NRP1.monoclones, "Monoclonal NRP1 CRISPR KO cells")
panel.heparan.sulfate <- plotHeparanSulfate(data.heparan.sulfate, "Growth media on H552 cells")

arranged.blockingAB <- arrangeGrob(
  panel.blockingAB, panel.blockingAB.CRISPR,
  nrow = 2,
  ncol = 1)

arranged.alternative.receptors <- arrangeGrob(
  panel.NRP1.monoclones, panel.heparan.sulfate,
  nrow = 1,
  ncol = 2)

saveFig(panel.ViralLoad.CRISPR.bulk, "Figure3_B", 1.2, 3.7)
saveFig(panel.ViralLoad.CRISPR.monoclonal, "Figure3_D", 1.38, 3.7)
saveFig(arranged.blockingAB, "Figure3_AC", 3.5, 3.7)
saveFig(arranged.alternative.receptors, "Figure3_EFG", 1.2, 4.7)
################################################################################
