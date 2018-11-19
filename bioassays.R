# ########################################################################
# Generate a latex table of bioassay data for P. australia, P. australiana
# K. Perry, 23/5/2017 
#
# Notes:
# Data collected by SARDI, Kevin Powis and Greg Baker
# The table is designed for publication as a (probably) supplementary table in 
# ... a manuscript for submission to a journal
#
# Additional notes;:
# This is draft working code. It can be made more efficient, and object names are best not to include `.`,
# but they do for now.
# Whe revising code later, watch use of ED function. For log models (LL2), you must use interval = "fls" to get back-transformed estimates.
#
# #############################################

library(tidyverse)
library(readxl)
library(xtable)
library(gridGraphics)
library(gridExtra)
library(grid)
library(drc)
library(RColorBrewer)
library(ggplot2)


# ##################
# define functions #
# ##################

# a function to summarise the data across reps 
# ... for plotting averages of observed points only
# ... make sure you use the full raw dataset for running models 
aveReps <- function(df){
  df  %>%
    group_by(insecticide, strain, dose, dose0) %>%
    mutate(propResp = resp / n) %>%
    summarise(
      n  = mean(n),
      resp = mean(resp))
}


# define a function to run a drc::drm logit model
# ... This function runs drm model using the log(ED50) approach and outputs a list of:
# ... ... drm model 
# ... ... input data with zero doses upward adjusted, 
# ... ... new predicted data for plotting a fitted curve
runLogit <- function(inpdata) {
  m <- drm(resp/n ~ dose, weights = n, data = inpdata, 
           fct = LL2.2(), type = "binomial")
  # ... new dose levels as support for the line
  newdata <- expand.grid(conc = exp(seq(log(1e-04), log(1000), length=1000)))
  # ... predictions and confidence intervals
  pm <- predict(m, newdata = newdata, interval = "confidence", reference = "control")
  # ... new data with predictions
  newdata$p    <- pm[,1]
  newdata$pmin <- pm[,2]
  newdata$pmax <- pm[,3]
  # ... shift conc == 0 a bit up, otherwise there are problems with coord_trans
  inpdata$dose0 <- inpdata$dose
  inpdata$dose0[inpdata$dose0 == 0] <- min(inpdata$dose0[inpdata$dose0 != 0])/10
  
  # if true, average the observed resp. across reps. for plotting 
  list(inpdata, m, newdata) %>%
    setNames(c("inpdata", "m", "newdata"))
}

# function to extract drm results into a table
extractParams <- function(x = drmFits, # a name of drmFits
                          lcLev = c(50, 99)) {
  
  #  extract regression coefficient and SE
  # ... note: in drm models, b = regres. coef (slope) and e is the LD50 estimate
  # ... note, drm puts out a negative slope, which I need to check (possibly models survivorship)
  b   <- -1 * coef(summary(x))["b:(Intercept)", "Estimate"]
  bSE <- coef(summary(x))["b:(Intercept)", "Std. Error"]
  
  # extract estimated effective dose and lower, upper, 95% cis
  # ... use interval = "fls" to backtranform estimate/cis from a log(dose) model (LL2.2)
  # ... if using dose model (LL.2), use interval = "delta"
  ed <-  ED(x, lcLev, interval = "fls", reference = "control")
  
  # for each estimate, create new columns for estimates and cis
  lapply(row.names(ed), function(y){
    
    e <- ed[y, ]
    lev <- gsub("e\\:1\\:", "", y)
    nms <- c(paste0("LC", lev), paste0("lci", lev), paste0("uci", lev))
    data.frame(est = e[names(e) == "Estimate"], # if this fails, use colnames() instead
               lci = e[names(e) == "Lower"],
               uci = e[names(e) == "Upper"]) %>%
      setNames(nms)
    
  }) %>%
    bind_cols() %>%
    cbind(data.frame(b = b,
                     bSE = bSE))
}


# add the registered field rates for each chemical (in mL/ha a.i.)
fieldUseRate <- function(x){
  gsub("dom", 40, x) %>%
    gsub("pro", 13.2, .) %>%
    gsub("cor", 20, .) %>%
    gsub("suc", 48, .) %>% as.numeric()
}

# #################
# read in data
# #################

rawDf <- read_excel("C:/UserData/Kym/PhD/Data/P.australiana/bioassayRawDataPowisMay2017.xlsx",
                    sheet = "allDataKym", skip = 1) %>%
  rename(dose = `dose_mg/L_ai`) %>%
  mutate(strain = paste(species, strain, sep = "_") %>%
           gsub("P.x_field", "P.x (SA)", .) %>%
           gsub("P.x_WS", "P.x (WS)", .) %>%
           gsub("P.aus_field", "P.aus (SA)", .))

# subset the data into strains and pesticides
pxWS.dom <- filter(rawDf, strain == "P.x (WS)", insecticide == "dominex") # Px Waite Susceptible
pxWS.cor <- filter(rawDf, strain == "P.x (WS)", insecticide == "coragen") 
pxWS.pro <- filter(rawDf, strain == "P.x (WS)", insecticide == "proclaim")
pxWS.suc <- filter(rawDf, strain == "P.x (WS)", insecticide == "successneo")

pxSA.dom <- filter(rawDf, strain == "P.x (SA)", insecticide == "dominex") # Px field strain from SA
pxSA.cor <- filter(rawDf, strain == "P.x (SA)", insecticide == "coragen")
pxSA.pro <- filter(rawDf, strain == "P.x (SA)", insecticide == "proclaim")
pxSA.suc <- filter(rawDf, strain == "P.x (SA)", insecticide == "successneo")

paSA.dom <- filter(rawDf, strain == "P.aus (SA)", insecticide == "dominex") # P australiana field strain from SA
paSA.cor <- filter(rawDf, strain == "P.aus (SA)", insecticide == "coragen")
paSA.pro <- filter(rawDf, strain == "P.aus (SA)", insecticide == "proclaim")
paSA.suc <- filter(rawDf, strain == "P.aus (SA)", insecticide == "successneo")


# #########################################################################
# plot fitted and observed dose response curves from the raw bioassay data
# ... use r package drc
# ... use logit regression
# ... extract output for ggplot2
# ... need to generate separate fitted model for each strain x pesticide combo
# ###############################################

# *****************************************
# run model, plot fitted curves for DOMINEX
# *****************************************
pxWS.dom.m1 <- runLogit(inpdata = pxWS.dom)
pxSA.dom.m1 <- runLogit(inpdata = pxSA.dom)
paSA.dom.m1 <- runLogit(inpdata = paSA.dom)

# set plotting parameters
ciAlpha <- 0.1 # alpha (transparency) level for the 95% CL limit ribbon
edLevel <- 99 # effective dose level for plotting each pesticide (geom_vline)
edColour  <- "black" # colour for the edLeveL line on plots for each strain  
fldSize <- 0.7       # size for the field dose rate line on plots
#fldColour <- "blue"  # colour for the field dose rate line on plots
fldColour <- brewer.pal(9, "Reds")[7]

(domPl <- ggplot() +
    coord_trans(x = "log") +
    xlab("Log dose") + 
    ylab("Proportion response") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text   = element_text(size = 10),
          axis.title = element_blank()) +
    scale_x_continuous(limits = c(0.0001, 1000),
                       breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                       labels = c(0,      0.001, 0.01, 0.1, 1, 10, 100, 1000)) +
    # first put all the vertical segments underneath
    geom_segment(aes(x = fieldUseRate("dom"), xend = fieldUseRate("dom"), 
                     y = -Inf, yend = 1),
                 lty = 1, size = fldSize, colour = fldColour) +
    geom_segment(aes(x    = ED(pxWS.dom.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     xend = ED(pxWS.dom.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     y = -Inf, yend = edLevel/100),
                 lty = 1, colour = edColour) +
    geom_segment(aes(x    = ED(pxSA.dom.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     xend = ED(pxSA.dom.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     y = -Inf, yend = edLevel/100),
                 lty = 2, colour = edColour) +
    geom_segment(aes(x    = ED(paSA.dom.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     xend = ED(paSA.dom.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     y = -Inf, yend = edLevel/100),
                 lty = 3, colour = edColour) +
    # P.x (WS)
    geom_point(data = aveReps(pxWS.dom.m1$inpdata),
               aes(x = dose0, y = resp/n), pch = 1, size = 2) +
    geom_line(data = pxWS.dom.m1$newdata,
              aes(x = conc, y = p), lty = 1) +
    geom_ribbon(data  = pxWS.dom.m1$newdata,
                aes(x = conc, y = p, ymin = pmin, ymax = pmax),
                alpha = ciAlpha) + 
    # P.x (SA)
    geom_point(data = aveReps(pxSA.dom.m1$inpdata),
               aes(x = dose0, y = resp/n), pch = 2, size = 2) +
    geom_line(data = pxSA.dom.m1$newdata,
              aes(x = conc, y = p), lty = 2) +
    geom_ribbon(data = pxSA.dom.m1$newdata,
                aes(x = conc, y = p, ymin = pmin, ymax = pmax),
                alpha = ciAlpha) +
    # P.aus (SA)
    geom_point(data = aveReps(paSA.dom.m1$inpdata),
               aes(x = dose0, y = resp/n), pch = 3, size = 2) +
    geom_line(data = paSA.dom.m1$newdata,
              aes(x = conc, y = p), lty = 3) +
    geom_ribbon(data  = paSA.dom.m1$newdata, aes(x = conc, y = p, ymin = pmin, ymax = pmax),
                alpha = ciAlpha))

# *****************************************
# run model, plot fitted curves for CORAGEN
# *****************************************
pxWS.cor.m1 <- runLogit(inpdata = pxWS.cor)
pxSA.cor.m1 <- runLogit(inpdata = pxSA.cor)
paSA.cor.m1 <- runLogit(inpdata = paSA.cor)

(corPl <- ggplot() +
    coord_trans(x = "log") +
    xlab("Log dose") + 
    ylab("Proportion response") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text   = element_text(size = 10),
          axis.title = element_blank()) +
    scale_x_continuous(limits = c(0.0001, 1000),
                       breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                       labels = c(0,      0.001, 0.01, 0.1, 1, 10, 100, 1000)) +
    # first put all the vertical segments underneath 
    geom_segment(aes(x    = fieldUseRate("cor"), 
                     xend = fieldUseRate("cor"), 
                     y = -Inf, yend = 1),
                 lty = 1, size = fldSize, colour = fldColour) +
    geom_segment(aes(x    = ED(pxWS.cor.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     xend = ED(pxWS.cor.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     y = -Inf, yend = edLevel/100),
                 lty = 1, colour = edColour) +
    geom_segment(aes(x    = ED(pxSA.cor.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     xend = ED(pxSA.cor.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     y = -Inf, yend = edLevel/100),
                 lty = 2, colour = edColour) +
    geom_segment(aes(x    = ED(paSA.cor.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     xend = ED(paSA.cor.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     y = -Inf, yend = edLevel/100),
                 lty = 3, colour = edColour) +
    # P.x (WS)
    geom_point(data = aveReps(pxWS.cor.m1$inpdata), 
               aes(x = dose0, y = resp/n), pch = 1, size = 2) +
    geom_line(data = pxWS.cor.m1$newdata, 
              aes(x = conc, y = p), lty = 1) +
    geom_ribbon(data  = pxWS.cor.m1$newdata, aes(x = conc, y = p, ymin = pmin, ymax = pmax), 
                alpha = ciAlpha) +
    # P.x (SA)
    geom_point(data = aveReps(pxSA.cor.m1$inpdata), 
               aes(x = dose0, y = resp/n), pch = 2, size = 2) +
    geom_line(data = pxSA.cor.m1$newdata, 
              aes(x = conc, y = p), lty = 2) +
    geom_ribbon(data  = pxSA.cor.m1$newdata, aes(x = conc, y = p, ymin = pmin, ymax = pmax), 
                alpha = ciAlpha) +
    # P.aus (SA)
    geom_point(data = aveReps(paSA.cor.m1$inpdata), 
               aes(x = dose0, y = resp/n), pch = 3, size = 2) +
    geom_line(data = paSA.cor.m1$newdata, 
              aes(x = conc, y = p), lty = 3) +
    geom_ribbon(data  = paSA.cor.m1$newdata, aes(x = conc, y = p, ymin = pmin, ymax = pmax), 
                alpha = ciAlpha))

# *******************************************
# run model, plot fitted curves for PROCLAIM
# *******************************************
pxWS.pro.m1 <- runLogit(inpdata = pxWS.pro)
pxSA.pro.m1 <- runLogit(inpdata = pxSA.pro)
paSA.pro.m1 <- runLogit(inpdata = paSA.pro)

(proPl <- ggplot() +
    coord_trans(x = "log") +
    # xlab("Log dose") + 
    # ylab("Proportion response") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text   = element_text(size = 10),
          axis.title = element_blank()) + # remove axis titles
    scale_x_continuous(limits = c(0.0001, 1000),
                       breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                       labels = c(0,      0.001, 0.01, 0.1, 1, 10, 100, 1000)) +
    # first put all the vertical segments underneath
    geom_segment(aes(x    = fieldUseRate("pro"), 
                     xend = fieldUseRate("pro"), 
                     y = -Inf, yend = 1),
                 lty = 1, size = fldSize, colour = fldColour) +
    geom_segment(aes(x    = ED(pxWS.pro.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     xend = ED(pxWS.pro.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     y = -Inf, yend = edLevel/100),
                 lty = 1, colour = edColour) +
    geom_segment(aes(x    = ED(pxSA.pro.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     xend = ED(pxSA.pro.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     y = -Inf, yend = edLevel/100),
                 lty = 2, colour = edColour) +
    geom_segment(aes(x    = ED(paSA.pro.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     xend = ED(paSA.pro.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     y = -Inf, yend = edLevel/100),
                 lty = 3, colour = edColour) +
    # P.x (WS)
    geom_point(data = aveReps(pxWS.pro.m1$inpdata), 
               aes(x = dose0, y = resp/n), pch = 1, size = 2) +
    geom_line(data = pxWS.pro.m1$newdata, 
              aes(x = conc, y = p), lty = 1) +
    geom_ribbon(data  = pxWS.pro.m1$newdata, aes(x = conc, y = p, ymin = pmin, ymax = pmax), 
                alpha = ciAlpha) +
    # P.x (SA)
    geom_point(data = aveReps(pxSA.pro.m1$inpdata), 
               aes(x = dose0, y = resp/n), pch = 2, size = 2) +
    geom_line(data = pxSA.pro.m1$newdata, 
              aes(x = conc, y = p), lty = 2) +
    geom_ribbon(data  = pxSA.pro.m1$newdata, aes(x = conc, y = p, ymin = pmin, ymax = pmax), 
                alpha = ciAlpha) +
    # P.aus (SA)
    geom_point(data = aveReps(paSA.pro.m1$inpdata), 
               aes(x = dose0, y = resp/n), pch = 3, size = 2) +
    geom_line(data = paSA.pro.m1$newdata, 
              aes(x = conc, y = p), lty = 3) +
    geom_ribbon(data  = paSA.pro.m1$newdata, aes(x = conc, y = p, ymin = pmin, ymax = pmax), 
                alpha = ciAlpha))

# *********************************************
# run model, plot fitted curves for SUCCESS NEO
# *********************************************
pxWS.suc.m1 <- runLogit(inpdata = pxWS.suc)
pxSA.suc.m1 <- runLogit(inpdata = pxSA.suc)
paSA.suc.m1 <- runLogit(inpdata = paSA.suc)

(sucPl <- ggplot() +
    coord_trans(x = "log") +
    # xlab("Log dose") +
    # ylab("Proportion response") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text   = element_text(size = 10),
          axis.title =  element_blank()) + # remove axis titles
    scale_x_continuous(limits = c(0.0001, 1000),
                       breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                       labels = c(0,      0.001, 0.01, 0.1, 1, 10, 100, 1000)) +
    # first put all the vertical segments underneath
    geom_segment(aes(x    = fieldUseRate("suc"), 
                     xend = fieldUseRate("suc"), 
                     y = -Inf, yend = 1),
                 lty = 1, size = fldSize, colour = fldColour) +
    geom_segment(aes(x    = ED(pxWS.suc.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     xend = ED(pxWS.suc.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     y = -Inf, yend = edLevel/100),
                 lty = 1, colour = edColour) +
    geom_segment(aes(x    = ED(pxSA.suc.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     xend = ED(pxSA.suc.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     y = -Inf, yend = edLevel/100),
                 lty = 2, colour = edColour) +
    geom_segment(aes(x    = ED(paSA.suc.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     xend = ED(paSA.suc.m1$m, edLevel, interval = "fls")[,"Estimate"],
                     y = -Inf, yend = edLevel/100),
                 lty = 3, colour = edColour) +
    # P.x (WS)
    geom_point(data = aveReps(pxWS.suc.m1$inpdata), 
               aes(x = dose0, y = resp/n), pch = 1, size = 2) +
    geom_line(data = pxWS.suc.m1$newdata, 
              aes(x = conc, y = p), lty = 1) +
    geom_ribbon(data  = pxWS.suc.m1$newdata, aes(x = conc, y = p, ymin = pmin, ymax = pmax), 
                alpha = ciAlpha) +
    # P.x (SA)
    geom_point(data = aveReps(pxSA.suc.m1$inpdata), 
               aes(x = dose0, y = resp/n), pch = 2, size = 2) +
    geom_line(data = pxSA.suc.m1$newdata, 
              aes(x = conc, y = p), lty = 2) +
    geom_ribbon(data  = pxSA.suc.m1$newdata, aes(x = conc, y = p, ymin = pmin, ymax = pmax), 
                alpha = ciAlpha) +
    # P.aus (SA)
    geom_point(data = aveReps(paSA.suc.m1$inpdata), 
               aes(x = dose0, y = resp/n), pch = 3, size = 2) +
    geom_line(data = paSA.suc.m1$newdata, 
              aes(x = conc, y = p), lty = 3) +
    geom_ribbon(data  = paSA.suc.m1$newdata, aes(x = conc, y = p, ymin = pmin, ymax = pmax), 
                alpha = ciAlpha))

# #########################################################################################
# draw the 4 plots on a single grid
# https://stackoverflow.com/questions/11076567/plot-a-legend-and-well-spaced-universal-y-axis-and-main-titles-in-grid-arrange
# #########################################################################################


# first, make a dummy legend and save as a grob
dummyData <- expand.grid(Strain = c("P. xylostella (S)", "P. xylostella", "P. australiana"),
                         dose   = c(0, 1, 10)) %>%
  mutate(resp = seq(1:nrow(.)),
         # reorder variables to show P. australiana first in lower legend 
         Strain = factor(Strain, levels = c("P. australiana", "P. xylostella", "P. xylostella (S)")))
gPlot <- ggplot(dummyData, aes(x = dose, y = resp, shape = Strain, lty = Strain)) +
  #geom_line() +
  geom_point() +
  theme_bw() + 
  theme(legend.title = element_blank(),
        legend.text = element_text(face = "italic", size = 12)) +
  guides(shape = guide_legend(override.aes = list(size=2))) + # to manually adjust the point size 
  scale_linetype_manual(values  = c(1, 3, 2)) +
  theme(legend.position = "bottom") +
  scale_shape_manual(values = c(3, 2, 1)) 

saveLegend <- function(gplot){
  tmp <- ggplot_gtable(ggplot_build(gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
legend <- saveLegend(gPlot)

# make function to draw an x-axis break on plot using grid.segments()
drawAxisBreak <- function(xStart = 0.113, yStart = 0.1375, 
                          segLength = 0.008, sepWidth = 0.008, lineWidth = 1,
                          adjustAngle = 1, units = "npc") {
  joiningSeg <- grid.segments(
    x0 = unit(xStart + (segLength / 2), units), 
    y0 = unit(yStart + (segLength / 2), units),
    x1 = unit(xStart + sepWidth + (segLength / 2), units),
    y1 = unit(yStart + (segLength / 2), units),
    draw = TRUE, gp = gpar(col = "white", lwd = 1.5 * lineWidth)
  )
  seg1 <- grid.segments(
    x0 = unit(xStart, units), y0 = unit(yStart, units),
    x1 = unit(xStart + segLength, units),
    y1 = unit(yStart + segLength, units),
    gp = gpar(col = "black", lwd = lineWidth),
    draw = TRUE
  )
  seg2 <- grid.segments(
    x0 = unit(xStart + sepWidth, units), y0 = unit(yStart, units),
    x1 = unit(xStart + sepWidth + segLength, units),
    y1 = unit(yStart + segLength, units),
    gp = gpar(col = "black", lwd = lineWidth),
    draw = TRUE
  )
}




ylab <- textGrob("Proportion response", rot = 90, vjust = 1)
#xlab <- textGrob("Log dose")
xlab <- textGrob(expression('Dose (mg L' ^ -1 ~ a.i.*')'))
ywidth  <- unit(2, "lines")  # width  to allow for global y lab
xheight <- unit(2, "lines")  # height to allow for global x lab
legheight <- unit(2, "lines") # height to allow space to plot legend below
p <- arrangeGrob(domPl, corPl, sucPl, proPl, ncol = 2)


old <- file.path("C:", "UserData", "Kym", "PhD", "Data", "P.australiana")
outPath <- file.path("C:", "UserData", "Kym", "PhD", "Data", "P.australiana")
#outPath <- file.path("C:", "UserData", "Kym", "PhD", "thesis", "images")
setwd(outPath)
plotWidth <- 11
pdf("doseResponseCurvesFinalRevisedProof.pdf", width = plotWidth, height = 2 * (plotWidth / 3))
## the below grid.arrange code for legend position on the graph (right)). Use the grid.draw(legend) line below
# grid.arrange(ylab, arrangeGrob(p, 
#                                xlab, 
#                                heights = unit.c(unit(1, "npc") - xheight, xheight), 
#                                nrow = 2),
#              widths = unit.c(unit(3, "lines"), unit(1, "npc") - unit(3, "lines")),
#              nrow = 1)
grid.arrange(ylab, arrangeGrob(p, 
                               xlab,
                               legend,
                               heights = unit.c(unit(1, "npc") - xheight - legheight, xheight, legheight), 
                               nrow = 3),
             widths = unit.c(unit(3, "lines"), unit(1, "npc") - unit(3, "lines")),
             nrow = 1)
vp1 <- viewport(x = 0.93, y = 0.77, width = 0.25, height = 0.25)
pushViewport(vp1)
#grid.daw(legend)
upViewport()
# plot labels
vp2 <- viewport(x = 0.105, y = 0.9655, width = 0.2, height = 0.2)
pushViewport(vp2)
grid.text("Dominex", just = "left")
upViewport()#
vp3 <- viewport(x = 0.578, y = 0.9655, width = 0.2, height = 0.2)
pushViewport(vp3)
grid.text("Coragen", just = "left")
upViewport()
vp4 <- viewport(x = 0.105, y = 0.515, width = 0.2, height = 0.2)
pushViewport(vp4)
grid.text("Success Neo", just = "left")
upViewport()
vp5 <- viewport(x = 0.578, y = 0.515, width = 0.2, height = 0.2)
pushViewport(vp5)
grid.text("Proclaim", just = "left")
upViewport()
# draw x-axis breaks between 0 and 0.001 for each plot.
drawAxisBreak(xStart = 0.12, yStart = 0.1391, segLength = 0.008, sepWidth = 0.008)
drawAxisBreak(xStart = 0.12, yStart = 0.5842, segLength = 0.008, sepWidth = 0.008)
drawAxisBreak(xStart = 0.592, yStart = 0.1391, segLength = 0.008, sepWidth = 0.008)
drawAxisBreak(xStart = 0.592, yStart = 0.5842, segLength = 0.008, sepWidth = 0.008)

dev.off()
setwd(old)



# now plot on a 2 x 2 grid
# vp1 <- viewport(x = 0,   y = 0.5, width = 0.5, height = 0.5, just = c(0,0))
# vp2 <- viewport(x = 0.5, y = 0.5, width = 0.5, height = 0.5, just = c(0,0))
# vp3 <- viewport(x = 0,   y = 0, width = 0.5, height = 0.5, just = c(0,0))
# vp4 <- viewport(x = 0.5, y = 0, width = 0.5, height = 0.5, just = c(0,0))

# #pdf("testGrid.pdf")
# grid.newpage()
# print(domPl, vp = vp1)
# print(corPl, vp = vp2)
# print(proPl, vp = vp3)
# print(sucPl, vp = vp4)
# print(legend, vp = vp4)
# grid.text(label = "A", 0.1, 0.9, gp = gpar(fontsize = 12, fontface = "bold"))
# grid.text(label = "Coragen", 0.1, 0.9, gp = gpar(fontsize = 12))
# #dev.off()
# 
# # plot dose-response curves
# old <- getwd()
# setwd("C:/UserData/Kym/PhD/thesis/images")
# pdf("doseResponseCurvesDraft.pdf", height = 7, width = 10)
# grid.newpage()
# print(domPl, vp = vp1)
# print(corPl, vp = vp2)
# print(proPl, vp = vp3)
# print(sucPl, vp = vp4)
# grid.text(label = "Dominex",  0.125,  0.95, gp = gpar(fontsize = 12, fontface = "bold"))
# grid.text(label = "Coragen",  0.6125, 0.95, gp = gpar(fontsize = 12, fontface = "bold"))
# grid.text(label = "Proclaim", 0.125,  0.45, gp = gpar(fontsize = 12, fontface = "bold"))
# grid.text(label = "Success",  0.6125, 0.45, gp = gpar(fontsize = 12, fontface = "bold"))
# dev.off()
# setwd(old)


# ########################################################
# Extract logit analysis parameters for publication table
# ########################################################

# the task: 
# (i) from model outputs, extract:
#... strain, product
# ... n
# ... the regression coefficient and SE
# ... ed50 +- 95% ci
# ... ed99 +- 95% ci
# (ii) calculate resistance ratios
# (iii) possibly backtransform the coefficients
# ()


# create a list of the model outputs (4 pesticides x 3 strains) 
# ... note: for drm models 2 param log-logistic models, \
# ... \ b is the regression coefficient, e is the ED50 estimate 
drmFits <- list(pxWS.dom = pxWS.dom.m1$m,
                pxSA.dom = pxSA.dom.m1$m, 
                paSA.dom = paSA.dom.m1$m,
                pxWS.cor = pxWS.cor.m1$m,
                pxSA.cor = pxSA.cor.m1$m, 
                paSA.cor = paSA.cor.m1$m,
                pxWS.pro = pxWS.pro.m1$m,
                pxSA.pro = pxSA.pro.m1$m, 
                paSA.pro = paSA.pro.m1$m,
                pxWS.suc = pxWS.suc.m1$m,
                pxSA.suc = pxSA.suc.m1$m,
                paSA.suc = paSA.suc.m1$m) 


# extract results into a table
drmTable <- lapply(X = drmFits, FUN = extractParams) %>%
  bind_rows %>%
  mutate(fitNms = names(drmFits),
         chemical = gsub("^[a-zA-Z]{4}\\.", "", fitNms),
         strain   = gsub("\\.[a-z]{3}$", "", fitNms),
         fieldRate = fieldUseRate(chemical),
         fieldToLC99Ratio = fieldRate / LC99)

# get baseline susceptibility (lc estimates for the Px (WS) strain
# ... for each insecticide, for caculating resistance ratios
baselineDf <- drmTable %>%
  filter(strain == "pxWS") %>%
  dplyr::select(chemical, baseLC50 = LC50, baseLC99 = LC99)



# get sample sizes for each bioassay group (chemical x strain) from the raw data 
nDf <- rawDf %>%
  filter(dose != 0) %>% # don't count the zero dose controls
  group_by(strain, insecticide) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  # match the strain and chemical names to enable merging
  mutate(strain = gsub("P\\.x \\(WS\\)", "pxWS", strain) %>%
           gsub("P\\.x \\(SA\\)", "pxSA", .) %>%
           gsub("P\\.aus \\(SA\\)", "paSA", .),
         chemical = gsub("coragen", "cor", insecticide) %>%
           gsub("dominex", "dom", .) %>%
           gsub("proclaim", "pro", .) %>%
           gsub("successneo", "suc", .)) %>%
  dplyr::select(-insecticide)


# format for latex
# export Latex table code
rnd <- 3
drmLatex <- drmTable %>%
  # add group sample sizes
  merge(nDf, by = c("chemical", "strain")) %>%
  # calculate resistance ratios
  merge(baselineDf, by = "chemical") %>%
  mutate(b   = round(b, rnd),
         bSE = round(bSE, rnd), 
         slopeSE = paste(b, bSE, sep = " \\pm "),
         Chemical = gsub("dom", "Dominex", chemical) %>%
           gsub("cor", "Coragen", .) %>%
           gsub("pro", "Proclaim", .) %>%
           gsub("suc", "Success Neo", .),
         Strain  = gsub("pxWS", "\\\\textit{P. x} (S)", strain) %>%
           gsub("pxSA", "\\\\textit{P. x}", .) %>%
           gsub("paSA", "\\\\textit{P. aus}", .),
         lci50 = round(lci50, rnd),
         uci50 = round(uci50, rnd),
         LC50range = paste0("(", lci50, "--", uci50, ")"),
         LC50RR = LC50 / baseLC50,
         lci99 = round(lci99, rnd),
         uci99 = round(uci99, rnd),
         LC99range = paste0("(", lci99, "--", uci99, ")"),
         LC99RR = LC99 / baseLC99) %>%
  dplyr::select(Chemical, 
                Strain, n,
                `{Slope (SE)}` = slopeSE,
                `{$LC_{50}$}` = LC50, 
                LC50range, 
                LC50RR,
                `{$LC_{99}$}` = LC99, 
                LC99range,
                LC99RR)
                #fieldToLC99Ratio)

drmLatex$Chemical[duplicated(drmLatex$Chemical)] <- NA

# set up multicolumn headers
# set up add.to.rows for the table sub-headers
drmNms <- list()
drmNms$pos <- list(0, 0)
drmNms$command <- c(
  paste(c(
  "Product & Strain & $n$ & {Slope (SE)} & ",
  "\\multicolumn{2}{c}{$LC_{50}$ (95\\% CL)} & ",
  "$RR_{LC50}$ & ",
  "\\multicolumn{2}{c}{$LC_{99}$ (95\\% CL)} & ",
  "$RR_{LC99}$ \\\\"), collapse = ""),
  # 2nd header line
  paste("&", "&", "&", "&",
        "\\multicolumn{2}{c}{[\\si{\\milli\\gram\\per\\liter} a.i.]} & & ",
        "\\multicolumn{2}{c}{[\\si{\\milli\\gram\\per\\liter} a.i.]} & \\\\\n",
        collapse = ""))

drmAlgn <- c("l", "l", "l",
             "S[table-number-alignment=center,table-figures-decimal=1,table-figures-integer=4]",
             "S[separate-uncertainty,table-figures-uncertainty=1,table-number-alignment=center,table-figures-integer=1,table-figures-decimal=3]",
             "S[round-mode=places,round-precision=3]", 
             "@{ }l",
             "S[round-mode=places,round-precision=2]",
             "S[round-mode=places,round-precision=3]", 
             "@{ }l",
             "S[round-mode=places,round-precision=2]")

drmDig <- c(0, 0, 0, 0, 3, 3, 3, 2, 3, 3, 2)
drmCap <- "Log-logistic regression statistics for dose-response bioassays on \\textit{P. australiana} (\\textit{P. aus}) and \\textit{P. xylostella} (\\textit{P.x}) field strains and the \\textit{P. xylostella} (S) reference strain
exposed to four commercial insecticides. 
Statistics presented include the number of insects tested (\\textit{n}), $LC_{50}$ and $LC_{99}$ estimates with 95\\% confidence limits, 
and resistance ratios ($RR$) at each $LC$ level."

xtable(drmLatex,
       caption = drmCap,
       align  = drmAlgn,
       digits = drmDig,
       lab = "tab:doseResponseLogit") %>% 
  print.xtable(floating = TRUE,
               include.rownames = FALSE,
               include.colnames = FALSE,
               add.to.row = drmNms,
               table.placement = "p",
               caption.placement = "top",
               NA.string = "",
               scalebox = 0.8,
               booktabs = TRUE,
               sanitize.text.function = function(x){x})

# ###############################################################
# summarise the control mortality across insecticides and strains
# ###############################################################

controlMortDf <- rawDf %>% 
  filter(dose == 0) %>%
  group_by(strain, insecticide) %>%
  summarise(totResp = sum(resp),
            totN = sum(n)) %>%
  mutate(mnPcResp = 100* totResp/totN)

(controlMortStrain <- controlMortDf %>%
  group_by(strain) %>%
  summarise(cMort = mean(mnPcResp)))

  
# End script
# ######################################


# OLD CODE BELOW

# ########################################
# full dose response probit analysis table
# ########################################

# read in the tabulated data
# bPath <- "C:/UserData/Kym/PhD/Data/P.australiana/bioassayDataTablePowisMay2017.xlsx"
# bioassayDf<- read_excel(bPath, skip = 1, sheet = "fullDoseResponse")
# 
# # export Latex table code
# bioassayLatex <- bioassayDf %>%
#   mutate(slopeSE = paste(Slope, SE, sep = " \\pm ")) %>%
#   mutate(Population = gsub("P. xylostella", "\\\\textit{P. xylostella}", Population) %>%
#            gsub("P. australiana", "\\\\textit{P. australiana}", .),
#          LC50range = paste0("(", LC50lcl, "--", LC50ucl, ")"),
#          LC99range = paste0("(", LC99lcl, "--", LC99ucl, ")")) %>%
#   dplyr::select(Population, Chemical, n,
#                 `{Slope (SE)}` = slopeSE,
#                 `{$LC_{50}$}` = LC50, LC50range, LC50RR,
#                 `{$LC_{99}$}` = LC99, LC99range, LC99RR)
# 
# bioassayLatex$Population[duplicated(bioassayLatex$Population)] <- NA

# set up multicolumn headers
# set up add.to.rows for the table sub-headers
# colNms <- list()
# colNms$pos <- list(0, 0)
# colNms$command <- c(paste(c(
#   "Population & ",
#   "Product &",
#   "$n$ & ",
#   "{Slope (SE)} & ",
#   "\\multicolumn{2}{c}{$LC_{50}$ (95\\% FL)} & ",
#   "$RR_{LC50}$ & ",
#   "\\multicolumn{2}{c}{$LC_{99}$ (95\\% FL)} & ",
#   "$RR_{LC99}$ \\\\"), collapse = ""),
#   # 2nd header line
#   paste("&", "&", "&", "&",
#         "\\multicolumn{2}{c}{[\\si{\\milli\\gram\\per\\liter} a.i.]} &",
#         " & ",
#         "\\multicolumn{2}{c}{[\\si{\\milli\\gram\\per\\liter} a.i.]} &", 
#         "\\\\", collapse = ""))


# bioAlgn <- c("l", "l", "l",
#              "S[table-number-alignment=center,table-figures-decimal=1,table-figures-integer=4]",
#              "S[separate-uncertainty,table-figures-uncertainty=1,table-number-alignment=center,table-figures-integer=1,table-figures-decimal=3]",
#              "S[round-mode=places,round-precision=3]", 
#              "@{ }l",
#              "S[round-mode=places,round-precision=2]",
#              "S[round-mode=places,round-precision=3]", 
#              "@{ }l",
#              "S[round-mode=places,round-precision=2]")
# 
# bioDig <- c(0, 0, 0, 0, 3, 3, 3, 2, 3, 3, 2)
# bioCap <- "Dose-response bioassay probit analysis for field strains of \\textit{P. australiana} and \\textit{P. xylostella} and the Waite susceptible \\textit{P. xylostella} (WS) strain
# exposed to four commercial insecticides."
# 
# xtable(bioassayLatex,
#        caption = bioCap,
#        align  = bioAlgn,
#        digits = bioDig,
#        lab = "tab:bioassay") %>% 
#   print.xtable(floating = TRUE,
#                include.rownames = FALSE,
#                include.colnames = FALSE,
#                add.to.row = colNms,
#                table.placement = "p",
#                caption.placement = "top",
#                NA.string = "",
#                scalebox = 0.8,
#                booktabs = TRUE,
#                sanitize.text.function = function(x){x})

# #######
# note: to keep precision of to the LC50 and LC99 ranges consistent (at 3 digits) 
# ... add the trailing zeros range by hand in latex 
# ... coding this in R is complex
# ######

# #########################
# discriminating dose table
# #########################

# ddDf <- read_excel(bPath, sheet = "discDose")
# ddLatex <- ddDf %>%
#   mutate(aiUnits = gsub("[0-9 ]+", "", aiConc),
#          aiConc  = as.numeric(gsub("[/a-zA-Z]", "", aiConc)),
#          aiUnits = gsub("g/L", "\\\\gram\\\\per\\\\liter", aiUnits) %>%
#            gsub("g/Kg", "\\\\gram\\\\per\\\\kilogram", .)) %>%
#   dplyr::select(Product, 
#                 `Active ingredient` = ai,
#                 `{a.i.}` = aiConc,
#                 `{a.i. units}` = aiUnits,
#                 `{Field use conc.}` = fieldUseConc,
#                 `{DD}` = DD, 
#                 `{Ratio (Field:DD)}` = fieldUseDDRatio,
#                 `{\\textit{P. xyl}}` = mortPx,
#                 `{\\textit{P. aus}}` = mortPaus)

# set up multicolumn headers
# set up add.to.rows for the table sub-headers
# ddNms <- list()
# ddNms$pos <- list(0, 0)
# ddNms$command <- c(
#   paste(c("Product & ", "Active ingredient & ", "\\multicolumn{2}{c}{a.i. conc.} & ",
#   "{Field use conc.} & ", "{DD} & ", "{Ratio} & ", 
#   "\\multicolumn{2}{c}{Percent mortality} \\\\\n"), collapse = "") %>%
#     paste0("\\cmidrule(lr){8-9}\n"),
#   paste(c("& & & & ", "{(\\si{\\milli\\gram\\per\\liter} a.i.)} & ", "{(\\si{\\milli\\gram\\per\\liter} a.i.)} & ",
#   "{(Field:DD)} & ","{\\textit{P. xyl}} & ", "{\\textit{P. aus}} \\\\\n"), collapse = "")
# )
# 
# 
# ddAlgn <- c("l", 
#             "l", "l", 
#             "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=3]", 
#             "@{\\space}s[table-unit-alignment=left]", "S", "S", "S", 
#             "S", "S")
# ddDig <- c(0,
#            0, 0, 0, 2, 2, 2, 1, 1, 1)
# ddCap <- "Mean percentage mortality of $3^{rd}$ instar larvae of \\textit{P. australiana} (SA) and \\textit{P. xylostella} (SA)
# field strains exposed to a discriminating dose (DD) of five commercial insecticides. The DD is the $LC_{99.9}$ for 
# the susceptible \\textit{P. xylostella} (WS) strain. A comparison of the DD to the commercial field use concentration
# is also presented."
# 
# xtable(ddLatex,
#        caption = ddCap,
#        align  = ddAlgn,
#        digits = ddDig,
#        lab = "tab:discriminatingDose") %>% 
#   print.xtable(floating = TRUE,
#                include.rownames = FALSE,
#                include.colnames = FALSE,
#                add.to.row = ddNms,
#                table.placement = "p",
#                caption.placement = "top",
#                NA.string = "",
#                scalebox = 0.8,
#                booktabs = TRUE,
#                sanitize.text.function = function(x){x})


# #################################
# plot the LC50s, LC99s with 95% FL
# #################################

# # lc50
# ggplot(bioassayDf, aes(x = Population, y = LC50,
#                        ymin = LC50 - LC50lcl, 
#                        ymax = LC50 + LC50ucl)) +
#   theme_bw() +
#   geom_point(position    = position_dodge(width = 0.25)) +
#   facet_wrap(~Chemical, scales = "free_y", ncol = 2) +
#   geom_errorbar(position = position_dodge(width = 0.25), 
#                 width = 0.2)
# 
# # lc99
# ggplot(bioassayDf, aes(x = Population, y = LC99,
#                        ymin = LC99 - LC99lcl, 
#                        ymax = LC99 + LC99ucl)) +
#   theme_bw() +
#   geom_point(position    = position_dodge(width = 0.25)) +
#   facet_wrap(~Chemical, scales = "free_y", ncol = 2) +
#   geom_errorbar(position = position_dodge(width = 0.25), 
#                 width = 0.2)


# read in the formatted raw dose response data


# # calculate mean response across replicates and plot
# sumDf <- rawDf %>%
#   group_by(insecticide, strain, dose) %>%
#   mutate(propResp = resp / n) %>%
#   summarise(
#     # test whether using the total rather than means matters
#     nTot = sum(n),
#             n  = mean(n),
#             resp = mean(resp))



# end snippet
# #######


