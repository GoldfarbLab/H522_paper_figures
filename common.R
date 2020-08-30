library(tidyverse)
library(RColorBrewer)
library(circlize)
library(egg)

golden.ratio <- 1/1.618
grey <- "#333333"
medium.grey <- "#666666"
light.grey <- "#AAAAAA"
COV2.color <- "#fb8072"

colors.MOI <- brewer.pal(5, 'YlGnBu')[2:5]

colors.Cell.line <- c("H522" = "#ff7f00",
                      "Vero E6" = "#984ea3", 
                      "Detroit562" = light.grey, 
                      "SCC25" = light.grey,
                      "OE21" = light.grey,
                      "KYSE30" = light.grey,
                      "H596" = light.grey,
                      "H1299" = light.grey,
                      "HCC827" = light.grey,
                      "PC-9" = light.grey,
                      "A427" = light.grey)

theme.basic <- (theme_minimal() 
                + theme(axis.line = element_line(colour = grey, size = 0.35, linetype = "solid"), 
                        panel.grid = element_blank(),
                        axis.ticks = element_line(size = 0.35, colour = grey),
                        axis.ticks.length = unit(.1, "cm"),
                        aspect.ratio = golden.ratio,
                        plot.title = element_text(size=7, hjust=0.5),
                        text = element_text(size=6),
                        axis.title = element_text(size=6)
                        )
)

saveFig <- function(p, filename, height, width)
{
  createDir(here("figures"))
  ggsave(paste(here("figures/"),filename, ".pdf", sep=""), width=width, height=height, compress=F, p)
}

# create directory if it doesn't exist
createDir <- function(path)
{
  if (!file.exists(path))
  {
    dir.create(path)
  }
}

################################################################################
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
################################################################################
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                       .fun = function(xx, col) {
                         c(N    = length2(xx[[col]], na.rm=na.rm),
                           mean = mean   (xx[[col]], na.rm=na.rm),
                           sd   = sd     (xx[[col]], na.rm=na.rm)
                         )
                       },
                       measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}