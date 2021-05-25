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
plotVeroConvalescent <- function(data) {
  data$Condition <- factor(data$Condition, levels=c("Mock infected", "095-02", "102-02", "113-02", "117-02", "119-02"))
  dfwc_between <- summarySE(data=data, measurevar="Copies", groupvars=c("Condition"), na.rm=FALSE, conf.interval=.95)
  dfwc_between$logVL <- log10(dfwc_between$Copies)
  dfwc_between$logVL[which(is.infinite(dfwc_between$logVL))] <- 0
  
  p <- (ggplot(dfwc_between, aes(x=as.factor(Condition),
                                 y=logVL))
        
        + geom_errorbar(aes(ymin=log10(`Copies` - se),
                            ymax=log10(`Copies` + se)), width=0.33, size=0.3, color=colors.Cell.line["Vero E6"])
        + stat_summary(fun = "mean", size= 0.15, geom = "bar", width=0.65, fill=lighter.grey, color=medium.grey)
        
        + geom_jitter(data=data, aes(y=log10(Copies)), size=0.75, shape=21, stroke = 0.4, width=0.25, fill="white", color=colors.Cell.line["Vero E6"])

        + ggtitle("Convalescent sera / Vero E6 cells")
        
        + scale_y_continuous(name = "Viral RNA (copies/cell)",
                             breaks = c(0, 2, 4, 6, 8),
                             labels = c("1e0", "1e2", "1e4", "1e6", "1e8"),
                             limits = c(0,8.2))
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 0.5,
                strip.text.x = element_text(size = 6),
                axis.title.x = element_blank(),
                plot.subtitle = element_text(size = 6, hjust=0.5))
  )
}
################################################################################

################################################################################
plotH522Convalescent <- function(data) {
  data$Condition <- factor(data$Condition, levels=c("Mock", "1", "2", "3", "4", "5"))
  dfwc_between <- summarySE(data=data, measurevar="Copies", groupvars=c("Condition"), na.rm=FALSE, conf.interval=.95)
  dfwc_between$logVL <- log10(dfwc_between$Copies)
  dfwc_between$logVL[which(is.infinite(dfwc_between$logVL))] <- 0
  
  p <- (ggplot(dfwc_between, aes(x=as.factor(Condition),
                                 y=logVL))
        
        + geom_errorbar(aes(ymin=log10(`Copies` - se),
                            ymax=log10(`Copies` + se)), width=0.33, size=0.3, color=colors.Cell.line["H522"])
        + stat_summary(fun = "mean", size= 0.15, geom = "bar", width=0.65, fill=lighter.grey, color=medium.grey)
        
        + geom_jitter(data=data, aes(y=log10(Copies)), size=0.75, shape=21, stroke = 0.4, width=0.25, fill="white", color=colors.Cell.line["H522"])
        
        + ggtitle("Convalescent sera / H522 cells")
        
        + scale_y_continuous(name = "Viral RNA (copies/cell)",
                             breaks = c(0, 2, 4, 6),
                             labels = c("1e0", "1e2", "1e4", "1e6"),
                             limits = c(0, 6.25))
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 0.5,
                strip.text.x = element_text(size = 6),
                axis.title.x = element_blank(),
                plot.subtitle = element_text(size = 6, hjust=0.5))
  )
}
################################################################################

################################################################################
plotVeroVaccinated <- function(data) {
  #data$Condition <- factor(data$Dilution, levels=c("1:320", "1:160", "1:80", "1:40"))
  
  p <- (ggplot(data, aes(x=Dilution,
                         y=log10(Copies),
                         group=interaction(Condition, Day),
                         color=as.factor(Day)))

        + geom_line()
        + geom_point(size=0.5)
        
        + ggtitle("Vaccinated sera / Vero E6 cells")
        
        + scale_y_continuous(name = "Viral RNA (copies/cell)",
                             breaks = c(2, 4, 6, 8),
                             labels = c("1e2", "1e4", "1e6", "1e8"),
                             limits = c(1.7, 8))
        
        + scale_x_continuous(name = "Dilution",
                            breaks = c(320, 160, 80, 40),
                            labels = c("1:320", "1:160", "1:80", "1:40"),
                            expand = c(0.1, 0),
                            trans = trans_reverser('log10'))
        
        + scale_color_manual(name = "Days", values=c("0" = "#dea4ca",
                                                     "35" = "#984ea3",
                                                     "-1" = "#666666"))
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 0.6,
                strip.text.x = element_text(size = 6),
                axis.title.x = element_blank(),
                plot.subtitle = element_text(size = 6, hjust=0.5))
  )
}
################################################################################


################################################################################
plotH522Vaccinated <- function(data) {
  dfwc_between <- summarySE(data=data, measurevar="Copies", groupvars=c("Condition", "Day", "Dilution"), na.rm=FALSE, conf.interval=.95)
  dfwc_between$logVL <- log10(dfwc_between$Copies)
  dfwc_between$logVL[which(is.infinite(dfwc_between$logVL))] <- 0
  
  p <- (ggplot(dfwc_between, aes(x=Dilution,
                         y=logVL,
                         group=interaction(Condition, Day),
                         color=as.factor(Day)))
        
        + geom_errorbar(aes(ymin=log10(`Copies` - se),
                            ymax=log10(`Copies` + se)), width=0.05, size=0.3)
        + geom_line()
        + geom_point(size=0.5)
        
        + ggtitle("Vaccinated sera / H522 cells")
        
        + scale_y_continuous(name = "Viral RNA (copies/cell)",
                             breaks = c(2, 4, 6),
                             labels = c("1e2", "1e4", "1e6"),
                             limits = c(1.7, 6.25))
        
        + scale_x_continuous(name = "Dilution",
                             breaks = c(320, 160, 80, 40),
                             labels = c("1:320", "1:160", "1:80", "1:40"),
                             expand = c(0.1, 0),
                             trans = trans_reverser('log10'))
        
        + scale_color_manual(name = "Days", values=c("0" = "#f2ceb6",
                                                     "35" = "#ff7f00",
                                                     "-1" = "#666666"))
        
        + theme.basic
        + theme(legend.position = "none",
                aspect.ratio = 0.6,
                strip.text.x = element_text(size = 6),
                axis.title.x = element_blank(),
                plot.subtitle = element_text(size = 6, hjust=0.5))
  )
}
################################################################################










################################################################################
# Read data
################################################################################
data.inhibitors <- read_csv(here("data/Figure 4/inhibitors_noOutliers.csv"))
data.dose.response <- read_csv(here("data/Figure 4/dose_response.csv"))
data.donors <- read_csv(here("data/Figure 4/donors.csv"))
data.vero.convalescent <- read_csv(here("data/Figure 4/Vero_convalescent.csv"))
data.H522.convalescent <- read_csv(here("data/Figure 4/H522_convalescent.csv"))
data.vero.vaccinated <- read_csv(here("data/Figure 4/Vero_vaccinated.csv"))
data.H522.vaccinated <- read_csv(here("data/Figure 4/H522_vaccinated.csv"))
################################################################################



################################################################################
# Generate figures
################################################################################
panel.inhibitors <- plotInhibitors(data.inhibitors)
panel.doseResponse <- plotDoseResponse(data.dose.response)
panel.donors <- plotDonors(data.donors)
panel.vero.convalescent <- plotVeroConvalescent(data.vero.convalescent)
panel.H522.convalescent <- plotH522Convalescent(data.H522.convalescent)
panel.vero.vaccinated <- plotVeroVaccinated(data.vero.vaccinated)
panel.H522.vaccinated <- plotH522Vaccinated(data.H522.vaccinated)

arranged.inhibitors <- arrangeGrob(
  panel.inhibitors, 
  panel.doseResponse,
  nrow = 2,
  ncol = 1)

arranged.lower <- arrangeGrob(
  panel.H522.convalescent, 
  panel.vero.vaccinated,
  panel.H522.vaccinated,
  nrow = 1,
  ncol = 3)

saveFig(arranged.inhibitors, "Figure4_AB", 3.5, 6.85)
saveFig(panel.donors, "Figure4_C", 1.35, 2)
saveFig(panel.vero.convalescent, "Figure4_D", 1.35, 2)
saveFig(arranged.lower, "Figure4_EFG", 1.25, 6.85)
################################################################################
