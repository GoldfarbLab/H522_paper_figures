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
  data <- data %>% filter(is.na(outlier))
  data$siRNA <- factor(data$siRNA, levels=c("siScramble uninfected", "siScramble Infected",
                                            "Mock Uninfected", "Mock Infected",
                                            "siTLR3", "siTLR7", "siTLR8", "siTLR9",
                                            "siMyD88", "siTRIF - TICAM", "siRIG-I_DDX58",
                                            "siMDA5_IFIH1", "siLPG2_DHX58", "siMAVS"))
  
  p <- (ggplot(data = data, aes(x=siRNA, y=abundance))
        + geom_point(color = colors.Cell.line["H522"])
        
        + facet_wrap(~ gene, ncol=1, scales="free")

        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 0.2,
                strip.text.x = element_text(size = 6))
  )
}
################################################################################




################################################################################
# Read data
################################################################################
KD.efficiency <- read_csv(here("data/Figure 7/kd-efficiency.csv"))
ISG.induction <- read_csv(here("data/Figure 7/ISG_induction.csv"))
################################################################################



################################################################################
# Generate figures
################################################################################
panel.KD.efficiency <- plotKD.efficiency(KD.efficiency)
panel.ISG.induction <- plotISG.induction(ISG.induction)

print(panel.ISG.induction)
#saveFig(panel.KD.efficiency, "FigureS4_A", 3, 5)
################################################################################
