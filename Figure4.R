library(here)
source(here("common.R"))

################################################################################
# Plotting functions
################################################################################
plotInhibitors <- function(data) {
  data <- data %>% filter(Condition != "DMSO") %>%
    filter(!(Condition %in% c("EIPA", "Dynasore", "Latrinculin B")))
  #data <- data %>% filter(Infection < 500)
  data$Cell.line <- factor(data$Cell.line, levels=c("Vero E6", "H522", "H522-ACE2"))
  data$Condition <- factor(data$Condition, 
                           levels=c("DMSO", 
                                    "Bafilomycin A",
                                    "AAK1 inhibitor",
                                    "E64D", 
                                    "Apilimod", 
                                    "EIPA",
                                    "Dynasore",
                                    "Latrinculin B",
                                    "Camostat mesylate"),
                           labels=c("DMSO",
                                    "Bafilomycin A\n10uM",
                                    "AAK1 inhibitor\n5uM",
                                    "E64D\n5uM", 
                                    "Apilimod\n0.1uM", 
                                    "EIPA\n25uM",
                                    "Dynasore\n10uM",
                                    "Latrinculin B\n0.5uM",
                                    "Camostat mesylate\n50uM"))

  p <- (ggplot(data, aes(x=Cell.line,
                         y=Infection,
                         color=Cell.line,
                         group=interaction(Cell.line, Condition),
                         facet=Cell.line,
                         ))
        + geom_boxplot(width=0.65, outlier.shape=NA, size=0.3)
        + geom_point(size=0.5, shape=21, stroke = 0.5, fill="white", alpha=1.0, position = position_jitterdodge(jitter.width=0.25))
        
        + scale_y_continuous(name = "RNA copies relative to DMSO (%)",
                             breaks = c(0, 50, 100, 150),
                             limits = c(0,165))
        
        + scale_color_manual(name = "Cell.line", values=colors.Cell.line)
        
        + facet_wrap_custom(~ Condition, nrow=1, scales="free_y", scale_overrides = list(
          scale_override(5, scale_y_continuous(limits=c(0,315), breaks = c(0, 50, 100, 150, 200, 250, 300)))
        ))
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 0.8,
                strip.text.x = element_text(size=7),
                axis.title.x = element_blank(),
                axis.text.x = element_text(size = 6, angle=45, vjust=1, hjust=1))
  )
}
################################################################################


################################################################################
plotDoseResponse <- function(data) {
  data$Condition <- factor(data$Condition, levels=c("Bafilomycin A", "AAK1 inhibitor", "E64D", "Apilimod", "Dynasore", "Camostat mesylate"))
  data <- data %>%
    filter(!(Condition %in% c("EIPA", "Dynasore", "Latrinculin B")))
  dfwc_between <- summarySE(data=data, measurevar="Infection", groupvars=c("Concentration", "Condition"), na.rm=FALSE, conf.interval=.95)
  
  # adjust 0's so they can be plotted in log space
  dfwc_between$Concentration[which(dfwc_between$Concentration == 0 & dfwc_between$Condition == "Bafilomycin A")] <- 0.08
  dfwc_between$Concentration[which(dfwc_between$Concentration == 0 & dfwc_between$Condition == "AAK1 inhibitor")] <- 0.04
  dfwc_between$Concentration[which(dfwc_between$Concentration == 0 & dfwc_between$Condition == "E64D")] <- 0.04
  dfwc_between$Concentration[which(dfwc_between$Concentration == 0 & dfwc_between$Condition == "Apilimod")] <- 0.004
  dfwc_between$Concentration[which(dfwc_between$Concentration == 0 & dfwc_between$Condition == "Camostat mesylate")] <- 0.4
  
  p <- (ggplot(dfwc_between, aes(x=Concentration,
                         y=Infection,
                         facet=Condition))
  + geom_errorbar(color="#ff7f00",
                   size=0.3,
                  aes(ymin=Infection - se,
                      ymax=Infection + se))
  + geom_line(color="#ff7f00", size=0.5)
  + geom_point( color="#ff7f00", size=0.5)
  
  + scale_x_continuous(name = "Concentration (uM)", trans="log10")
  + scale_y_continuous(name = "Relative infection (%)", limits=c(0,180), breaks = c(0, 50, 100, 150))
  
  + facet_wrap_custom(~ Condition, scales = "free", nrow = 1, scale_overrides = list(
      scale_override(1, scale_x_continuous(trans="log10", limits=c(0.08, 20), breaks = c(0.08, 0.4, 2, 10, 20), labels=c("0", "0.4", "2", "10", "20"))), 
      scale_override(2, scale_x_continuous(trans="log10", limits=c(0.04, 10), breaks = c(0.04, 0.2, 1, 5, 10), labels=c("0", "0.2", "1", "5", "10"))),
      scale_override(3, scale_x_continuous(trans="log10", limits=c(0.04, 10), breaks = c(0.04, 0.2, 1, 5, 10), labels=c("0", "0.2", "1", "5", "10"))),
      scale_override(4, scale_x_continuous(trans="log10", limits=c(0.004, 1), breaks = c(0.004, 0.02, 0.1, 0.5, 1), labels=c("0", "0.02", "0.1", "0.5", "1"))),
      scale_override(5, scale_x_continuous(trans="log10", limits=c(0.4, 100), breaks = c(0.4, 2, 10, 50, 100), labels=c("0", "2", "10", "50", "100"))),

      scale_override(5, scale_y_continuous(limits=c(0,300), breaks = c(0, 50, 100, 150, 200, 250, 300)))
    ))
  
  + theme.basic
  + theme(legend.position = "none",
          strip.text.x = element_text(size=7))
  )
}
################################################################################


################################################################################
plotDonors <- function(data) {
  
  # adjust 0's so they can be plotted in log space
  #data$Concentration[which(data$Concentration == 0)] <- 0.2

  p <- (ggplot(data, aes(x=Concentration,
                         y=Copies,
                         color=Donor,
                         group=Donor
  ))

  + geom_line()
  + geom_point(size=0.5, stroke = 0.75)
  
  + ggtitle("AAK1 inhibitor")
  
  + scale_y_continuous(name = "Relative RNA copies (%)")
  + scale_x_continuous(name = "Concentration (uM)",
                       #trans="log10",
                       breaks = c(0, 1, 5, 10),
                       limits = c(0, 10))#,
                      # labels=c("0", "1", "5", "10"))
  

  
  + scale_color_manual(name = "Donors", values=colors.donors)
  
  + theme.basic
  + theme(legend.position = "none")
  )
}
################################################################################



################################################################################
# Read data
################################################################################
data.inhibitors <- read_csv(here("data/Figure 4/inhibitors_noOutliers.csv"))
data.dose.response <- read_csv(here("data/Figure 4/dose_response.csv"))
data.donors <- read_csv(here("data/Figure 4/donors.csv"))
################################################################################



################################################################################
# Generate figures
################################################################################
panel.inhibitors <- plotInhibitors(data.inhibitors)
panel.doseResponse <- plotDoseResponse(data.dose.response)
panel.donors <- plotDonors(data.donors)

arranged.inhibitors <- arrangeGrob(
  panel.inhibitors, 
  panel.doseResponse,
  nrow = 2,
  ncol = 1)

saveFig(arranged.inhibitors, "Figure4_AB", 3.5, 6.85)
saveFig(panel.donors, "Figure4_D", 1.35, 2)
################################################################################
