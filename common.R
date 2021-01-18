library(tidyverse)
library(RColorBrewer)
library(circlize)
library(egg)
library(extrafont)
library(outliers)
library(scales)
library(grid)
#library(ggsignif)

#font_import()
loadfonts()

################################################################################
# Proteomics parameters and thresholds
################################################################################
proteomics.stats.threshold = 0.05
proteomics.log.stats.threshold = -log10(proteomics.stats.threshold)
proteomics.log.fc.threshold = log2(1.3)
proteomics.num.clusters = 7
################################################################################


################################################################################
# Visualization
################################################################################
golden.ratio <- 1/1.618
grey <- "#333333"
medium.grey <- "#666666"
light.grey <- "#AAAAAA"
lightest.grey <- "#DDDDDD"
COV2.color <- "#fb8072"

colors.MOI <- brewer.pal(5, 'YlGnBu')[2:5]
colors.MOI.6 <- brewer.pal(7, 'YlGnBu')[2:7]

colors.donors <- brewer.pal(5, 'Dark2')

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
                      "A427" = light.grey,
                      "H522-ACE2" = "#377eb8", #"#fdbf6f",#4daf4a
                      "PgsA" = medium.grey,
                      "Calu-3" = "#4daf4a",
                      "Calu-3 ACE2 KO" = "#4daf4a",
                      "H522 ACE2 KO" = "#ff7f00")  #"#377eb8"

color.structural <-  "#1b9e77" #"#e41a1c"
color.orf <- "#d95f02" #"#377eb8"
color.nsp <- "#7570b3" #"#4daf4a"

colors.SARS.proteins <- c("M" = color.structural, 
                      "N" = color.structural, 
                      "S" = color.structural, 
                      "Orf7a" = color.orf, 
                      "Orf8" = color.orf,
                      "Orf3a" = color.orf,
                      "Orf3b" = color.orf,
                      "Orf9b" = color.orf,
                      "Nsp1" = color.nsp,
                      "Nsp3" = color.nsp,
                      "Nsp4" = color.nsp,
                      "Nsp2" = color.nsp,
                      "Nsp8" = color.nsp)

order.Cell.line <- c("Vero E6", "H522", "H596", "H1299", "A427", "HCC827", "PC-9", "Detroit562", "SCC25", "KYSE30", "OE21")

theme.basic <- (theme_minimal(base_family = "Arial")
                + theme(axis.line = element_line(colour = grey, size = 0.5, linetype = "solid"), 
                        panel.grid = element_blank(),
                        axis.ticks = element_line(size = 0.5, colour = grey),
                        axis.ticks.length = unit(.1, "cm"),
                        aspect.ratio = golden.ratio,
                        plot.title = element_text(size=7, hjust=0.5),
                        axis.text = element_text(size=6), #sans = arial
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
                      conf.interval=.95, .drop=TRUE, removeOutliers=F, logt=F) {
  
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                       .fun = function(xx, col) {
                         
                         if (removeOutliers) {
                           print('pre')
                           print(xx[[col]])
                           #xx[[col]] <- boxB(xx[[col]])
                           print('post')
                           #print(exp(rm.outlier(log(xx[[col]]), fill = FALSE, median = FALSE, opposite = FALSE)))
                           print(scores(log(xx[[col]]), type="t", prob=0.90) )
                           #print(xx[[col]][boxB(xx[[col]], method="resistant", logt=logt, k=2.5)[["outliers"]]])
                         }
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


# used to override scales on facet'ed plots
scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}

CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)

facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}