library(here)
source(here("common.R"))

################################################################################
# Plotting functions
################################################################################
plotKD.efficiency <- function(data) {
  data <- data %>% filter(Gene != "TLR8")
  data$Gene <- factor(data$Gene, levels= c("TLR3", "TLR7","TLR9", "TRIF", "MyD88", "RIG-I", "MDA5", "LGP-2", "MAVS"))
  
  p <- (ggplot(data = data, aes(x=Time, y=KD, group=Gene))
        + geom_point(color = colors.Cell.line["H522"])
        + stat_summary(fun = "mean", size= 0.15, geom = "crossbar")
        
        + facet_wrap(~ Gene, scales="free", nrow=2)
        
        + scale_x_continuous(name = "Hours post-infection", 
                             breaks = c(24, 96, 120),
                             labels = c("24", "96", "120"),
                             limits = c(12, 132))
        + scale_y_continuous(name = "Relative copies (siRNA / siNT)",
                             breaks = seq(0,1,0.2),
                             labels = c("0", "0.2","0.4", "0.6", "0.8", "1.0"),
                             limits = c(0,1.05))
        
        + ggtitle("siRNA knockdown efficiency")
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 1.2,
                strip.text.x = element_text(size = 6))
        )
}

plotISG.induction <- function(data) {
  
  data$siRNA <- factor(data$siRNA, levels=c("siTLR3", "siTLR7", "siTLR8", "siTLR9",
                                            "siMyD88", "siTRIF (TICAM)", "siRIG-I (DDX58)",
                                            "siMDA5 (IFIH1)", "siLPG2 (DHX58)", "siMAVS"))
  data$gene <- factor(data$gene, levels=c("ISG15", "IFIT1", "IFIT2", "MX1"))
  
  data <- data %>% filter(is.na(outlier))
  
  p <- (ggplot(data = data, aes(x=gene, y=log2(abundance)))
        
        + geom_hline(yintercept=0, linetype="dotted", color=light.grey)
        
        + geom_point(color = colors.Cell.line["H522"], size=0.75, shape=21, stroke = 0.4, fill="white")
        + stat_summary(fun = "mean", size= 0.15, width=0.65, geom = "crossbar")
        
        + facet_wrap(~ siRNA, nrow=1, scales="free_x")
        
        + scale_y_continuous(name = "log2 (mRNA / NT infected)",
                             breaks = c(-8,-6,-4,-2,0,2),
                             #labels = c("-2", "0", "2", "4", "6", "8", "10"),
                             limits = c(-9, 2.5))

        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 1/golden.ratio,
                strip.text.x = element_text(size = 6),
                axis.title.x = element_blank(),
                axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  )
}

plotsiRNA.ViralRNA <- function(data) {
  
  data$siRNA <- factor(data$siRNA, levels=c("NT uninfected", "NT infected",
                                            "Mock uninfected", "Mock infected",
                                            "siTLR3", "siTLR7", "siTLR8", "siTLR9",
                                            "siMyD88", "siTRIF (TICAM)", "siRIG-I (DDX58)",
                                            "siMDA5 (IFIH1)", "siLPG2 (DHX58)", "siMAVS"))
  data$gene <- factor(data$gene, levels=c("ISG15", "IFIT1", "IFIT2", "MX1"))
  
  data <- data %>% filter(is.na(outlier))
  
  dfwc_between <- summarySE(data=data, measurevar="abundance", groupvars=c("siRNA"), na.rm=FALSE, conf.interval=.95)
  dfwc_between$logVL <- log10(dfwc_between$abundance)
  dfwc_between$logVL[which(is.infinite(dfwc_between$logVL))] <- 0
  
  #data$abundance <- log10(data$abundance)
  #data$abundance[which(is.infinite(data$abundance))] <- 0

  p <- (ggplot(data = dfwc_between, aes(x=siRNA, y=logVL))
        + geom_errorbar(aes(ymin=log10(`abundance` - se),
                            ymax=log10(`abundance` + se)), width=0.33, size=0.3, color=colors.Cell.line["H522"])
        + stat_summary(fun = "mean", size=0.15, width=0.65, geom = "bar", fill=lighter.grey, color=medium.grey)
        + geom_jitter(data=data, aes(y=log10(abundance)), size=0.75, shape=21, stroke = 0.4, width=0.25, fill="white", color=colors.Cell.line["H522"])
        
        + geom_hline(yintercept=0, linetype="dotted", color=light.grey)
        
        + scale_y_continuous(name = "Viral RNA (copies/cell)",
                             breaks = c(0,2,4,6,8),
                             labels = c("1e0", "1e2", "1e4", "1e6", "1e8"),
                             limits = c(0, 8))
        
        + ggtitle("Viral load after siRNA treatment")
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 0.3,
                strip.text.x = element_text(size = 6),
                axis.title.x = element_blank(),
                axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  )
}

plotsiRNA.genes.control <- function(data) {
  
  data$siRNA <- factor(data$siRNA, levels=c("NT uninfected", "NT infected",
                                            "Mock uninfected", "Mock infected"))
  data$gene <- factor(data$gene, levels=c("ISG15", "IFIT1", "IFIT2", "MX1"))
  
  data <- data %>% filter(is.na(outlier))
  
  p <- (ggplot(data = data, aes(x=gene, y=log2(abundance)))
        
        + geom_hline(yintercept=0, linetype="dotted", color=light.grey)
        
        + geom_point(color = colors.Cell.line["H522"], size=0.75, shape=21, stroke = 0.4, fill="white")
        + stat_summary(fun = "mean", size=0.16, width=0.65, geom = "crossbar")
        
        + facet_wrap(~ siRNA, nrow=1, scales="free_x")
        
        + scale_y_continuous(name = "log2 (mRNA / NT uninfected)",
                             breaks = c(-2,0,2,4,6,8,10),
                             labels = c("-2", "0", "2", "4", "6", "8", "10"),
                             limits = c(-2, 11.1))
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 1/golden.ratio,
                strip.text.x = element_text(size = 6),
                axis.title.x = element_blank(),
                axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  )
}
################################################################################




################################################################################
# Read data
################################################################################
KD.efficiency <- read_csv(here("data/Figure 7/kd-efficiency.csv"))
ISG.induction <- read_csv(here("data/Figure 7/ISG_induction.csv"))
siRNA.viralRNA <- read_csv(here("data/Figure 7/siRNA_viralRNA.csv"))
siRNA.genes.control <- read_csv(here("data/Figure 7/siRNA_genes_control.csv"))
################################################################################



################################################################################
# Generate figures
################################################################################
panel.KD.efficiency <- plotKD.efficiency(KD.efficiency)
panel.ISG.induction <- plotISG.induction(ISG.induction)
panel.siRNA.viralRNA <- plotsiRNA.ViralRNA(siRNA.viralRNA)
panel.siRNA.genes.control <- plotsiRNA.genes.control(siRNA.genes.control)


#print(panel.siRNA.genes.control)

saveFig(panel.siRNA.viralRNA, "Figure7Ci", 3, 3.6)
saveFig(panel.siRNA.genes.control, "Figure7Cii", 3, 2.9)
saveFig(panel.ISG.induction, "Figure7Ciii", 3, 6.85)
saveFig(panel.KD.efficiency, "FigureS4_A", 3, 5)
################################################################################
