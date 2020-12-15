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
# Read data
################################################################################
data.neutralizing <- read_csv(here("data/Figure 2/S_neutralizing.csv"))
data.fc <- read_csv(here("data/Figure 2/ACE2_Fc.csv"))
data.incucyte <- read_csv(here("data/Figure 2/incucyte.csv"))
################################################################################



################################################################################
# Generate figures
################################################################################
panel.S.AB <- plotAB(data.neutralizing, "S neutralizing antibody", "2B04 (ug/mL)", "Viral RNA (%)")
panel.ACE2.Fc <- plotAB(data.fc, "ACE2 Fc", "hACE2 Fc (ug/mL)", "Viral RNA (%)")
# Execute this with doBreak=T to get y-axis break
panel.incucyte <- plotIncucyte(data.incucyte, doBreak=F) 

arranged.AB.Fc <- arrangeGrob(
  panel.S.AB, panel.ACE2.Fc,
  nrow = 2,
  ncol = 1)

saveFig(arranged.AB.Fc, "Figure2_BC", 3, 3.7)
saveFig(panel.incucyte, "Figure2_D", 4.1, 6.85)
################################################################################
