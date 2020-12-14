library(here)
library(scales)
library(RColorBrewer)

source(here("common.R"))


################################################################################
# Read data
################################################################################
# Figure 2
data.neutralizing <- read_csv(here("data/Figure2A_S_neutralizing.csv"))
data.fc <- read_csv(here("data/Figure2B_ACE2_Fc.csv"))
data.blockingAB <- read_csv(here("data/Figure2C_ACE2_blocking.csv"))
data.ace2.CRISPR <- read_csv(here("data/Figure2D_ACE2_CRISPR_BULK.csv"))
data.ace2.CRISPR.blockingAB <- read_csv(here("data/Figure2_ACE2_CRISPR_blocking.csv"))
data.incucyte <- read_csv(here("data/Figure2F_incucyte.csv"))

################################################################################
# Figure 2A / B
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
        + stat_summary(fun = "mean", size= 0.15, geom = "crossbar", width=0.75)
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
# Figure 2C
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
# Figure 2D
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
# CRISPR blocking antibody
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
# Figure 2F
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
# Generate figures
################################################################################
# Figure 2
p2a <- plotAB(data.neutralizing, "S neutralizing antibody", "2B04 (ug/mL)", "Viral RNA (%)")
p2c <- plotAB(data.fc, "ACE2 Fc", "hACE2 Fc (ug/mL)", "Viral RNA (%)")
p2d <- plotBlocking(data.blockingAB, "Blocking antibody", "3 days post-infection", "Viral RNA (copies/cell)")
p2e <- plotViralLoadCRISPR(data.ace2.CRISPR, "Viral load in ACE2 CRISPR KO cells")
p2f <- plotCRISPRBlocking(data.ace2.CRISPR.blockingAB, "Blocking antibody in ACE2 CRISPR KO cells", "3 days post-infection", "Viral RNA (copies/cell)")
p2g <- plotIncucyte(data.incucyte, doBreak=F)

F2.top.left <- arrangeGrob(
  p2a, p2c, p2d, p2f,
  nrow = 4,
  ncol = 1)

F2.e <- arrangeGrob(
  p2e,
  nrow = 1,
  ncol = 2)

F2.bottom <- arrangeGrob(
  p2g,
  nrow = 1,
  ncol = 1)

saveFig(F2.e, "Figure2_E", 1.2, 3.7)
saveFig(F2.top.left, "Figure2_top_left", 7, 3.7)
saveFig(F2.bottom, "Figure2_bottom", 4.1, 6.85)