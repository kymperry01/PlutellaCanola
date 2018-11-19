# ################################################################################
# Analysis of climatic conditions leading up to crop colonisation
# Kym D Perry, 4/10/2016
#
# This script does:
# (1) Extracts rainfall data for all field sites (corresponding grid cells)
# (2) Plots mean cumulative rainfall for all field sites
# (3) Plots mean CLIMEX weekly Growth Index (GI) for field sites

# This script uses the workspace in the project `climateData.R`
#
# Notes:
# The time period of interest is Feb to May
# The sites of interest are colonisation sites and autumn survey sites (78 cells)
# 
# ################################################################################

library(tidyverse)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

# standard error
se <- function(x, na.rm = TRUE){
  if(na.rm) x <- x[!is.na(x)]
  sd(x) / sqrt(length(x))
  }


# #################
# Graphics set up #
# #################

# For setting common x-axis as day of year, but labelled by month name  
# Use 2014 as a `template` sequence, 
dateSeq <- seq(as_date("2014-02-01"), as_date("2014-6-01"), by = "month")
ydayBreaks <- dateSeq %>% yday() 
ydayLabels <- dateSeq %>% format("%b")

# # Function to offset 2016 by 1 day (leap year) so that `date` v `day of year` match across years 
# # Can probably ignore *weekly measures only)
# ydayAdjust <- function(x, adjustYear = 2016, adjustment = 1){
#   y  <- year(x)
#   yd <- yday(x) 
#   yd[which(y == adjustYear)] <- yd[which(y == adjustYear)] + adjustment
#   yd
# }

# Colour palettes
oranges <- brewer.pal(9, "Oranges")
greys   <- brewer.pal(9, "Greys")
greens  <- brewer.pal(9, "Greens")


# ###############################################
# Read in climate data extracted for grid cells #
# ###############################################

# Read in metadata for all field sites (colonisation, autumn survey, light traps)
# This file was made in R script `dataDrill.R`
climatePath <- file.path("C:", "UserData", "Kym", "PhD", "Data", "climateData")
fieldSites <- file.path(climatePath, "fieldSitesRastCells.csv") %>% 
  read_csv()

# Read in rainfall data for the correponding grid cells for all field sites
gridCells <- fieldSites %>%
  filter(studyType %in% c("colonisation", "survey", "lightTrap")) %>%
  distinct(cellNr) %>%
  unlist() %>% 
  as.numeric() %>%
  sort() # to ensure names in correct order below

fileNameToGridCell <- function(x){
  x %>% 
    str_extract("_\\d{1,8}") %>% str_extract("\\d{1,8}") %>%
    as.numeric()
}

addGridCell <- function(x, y){
  x %>% mutate(cellNr = y)
  }

climateDailyDf <- file.path(climatePath, "dailyDataDrill2016") %>%
  list.files(pattern = "dd_", full.names = TRUE) %>%
  magrittr::extract(which(fileNameToGridCell(.) %in% gridCells)) %>%
  map(read.table, skip = 8, header = TRUE, stringsAsFactors = FALSE) %>%
  # annoyingly, these files don't have the grid cell number internally, so add
  map2(gridCells, addGridCell) %>% 
  bind_rows() %>%
  mutate(date  = lubridate::ymd(date),
         year  = lubridate::year(date),
         month = lubridate::month(date)) 

# For each year (2014-2016), calculate a mean of cumulative rainfall across grid cells.
# We're interested in Feb-May (pre-colonisation) from 2014-2016 (study years)
cumRainDf <- climateDailyDf %>%
  filter(year %in% 2014:2016, month %in% 2:5) %>%
  split(list(.$year, .$cellNr)) %>%
  map(mutate, cumRain = cumsum(dRai)) %>%
  bind_rows() %>%
  group_by(date) %>%  # now calculate mean rainfall across cells
  summarise(nCells = length(date),
            meanCumRain = mean(cumRain),
            sdCumRain   = sd(cumRain),
            seCumRain   = se(cumRain)) %>%
  ungroup() %>%
  mutate(year = factor(lubridate::year(date)),
         yearDay = yday(date))
           
# =============================================================
# Calculate total monthly rainfall each for the manuscript text
# =============================================================
climateDailyDf %>%
  filter(year %in% 2014:2016, month %in% 2:5) %>% 
  # first calculate total monthly rainfall for each grid cell
  group_by(cellNr, year, month) %>%
  summarise(nDays = length(cellNr),
            rainMonthlyTotal = sum(dRai, na.rm = TRUE)) %>% 
  # now average across grid cells
  group_by(year, month) %>%
  summarise(nCell = length(cellNr),
            meanRainMonthlyTotal = mean(rainMonthlyTotal),
            sdRainMonthlyTotal   = sd(rainMonthlyTotal))


# ########################################
# Read in CLIMEX results for field sites #
# ########################################

climexPath <- file.path("C:", "UserData", "Kym", "PhD", "Data", "climexOutput")
climexWeeklyGI <- file.path(climexPath, "climexCompareYearsPointLocations") %>%
  list.files(pattern = "csv", full.names = TRUE) %>%
  magrittr::extract(which(fileNameToGridCell(.) %in% gridCells)) %>% # these are in ascending order
  map(read_csv, skip = 3, col_names = TRUE) %>%
  map2(gridCells, addGridCell) %>%
  bind_rows() %>%
  mutate(date  = lubridate::dmy(SimulaDate),
         month = lubridate::month(date),
         year  = lubridate::year(date))
         
# For each year and date, calculate a mean of the Weekly GI and stress indices across sites
climexIndices <- c("GI(w)", "DS(w)", "WS(w)", "CS(w)", "HS(w)", "HWS(w)", "DS(w/acc)", 
                   "WS(w/acc)", "CS(w/acc)", "HS(w/acc)", "HWS(w/acc)")
climexSiteMeans <- climexWeeklyGI %>%
  split(.$year) %>%
  map(group_by, date, month, year) %>% # now calculate means of all growth & stress indices  
  #map(summarise_at, "GI(w)", mean)
  map(summarise_at, climexIndices, funs(mean, se, length)) %>% # length, just to check number of rows summarised
  map(dplyr::select, c("date", "month", "year", nCells = "GI(w)_length", contains("_mean"), contains("_se"))) %>%
  bind_rows() %>%
  ungroup() %>%
  # keep Feb to May, 2014-2016
  filter(month %in% 2:5, year %in% 2014:2016)

# convert to long, but with mean and se data in separate columns for geom_ribbon()
climexSiteMeansMelt <- climexSiteMeans %>%
  dplyr::select(-nCells, -month) %>%
  reshape2::melt(
    id.vars = c("date", "year"),
    variable.name = "index",
    value.name = "value"
    ) %>% 
  mutate(metric = str_extract(index, "_[a-z]{1,4}"),
         metric = str_replace(metric, "_", ""),
         index  = str_replace(index, "_[a-z]{1,4}", "")) %>%
  reshape2::dcast(date + year + index ~ metric) %>% 
  dplyr::rename(meanValue = mean, seValue = se) %>%
  arrange(index, year, date) %>%
  mutate(yearDay = lubridate::yday(date),
         year = factor(year))


# plot climex indices
climexSiteMeansMelt %>%
  ggplot(aes(x = date, y = meanValue, colour = index)) +
  geom_line() +
  facet_wrap(~year, scales = "free_x", ncol = 1)

#> Of the weekly stress indices, only accumulate heat stress (HS(w/acc)) is limiting.
#> This suggests P. xylostella may start to build up in autumn when HS declines,
#> and moisture provides host plants.

# ===========================================================
# Therefore: Plot GIw and HS(w/acc) to compare between years
# ===========================================================

# ##########################  
# plot together using grid #
# ##########################

# plot rain
(plotRain <- cumRainDf %>%
   ggplot(aes(x = yearDay, y = meanCumRain, colour = year,
              ymin = meanCumRain - seCumRain, ymax = meanCumRain + seCumRain)) +
   geom_ribbon(aes(group = year), fill = alpha("grey", 0.7), colour = "transparent") + 
   geom_line() +
   theme_bw(base_size = 11) +
   theme(legend.direction = "vertical",
         legend.position = c(0.1, 0.85),
         legend.justification = "centre",
         panel.grid.major.y = element_blank(),
         panel.grid.minor   = element_blank(),
         plot.margin = unit(c(2, 2, 2, 5), "mm")) + # t r b l
   scale_x_continuous(breaks = ydayBreaks, labels = ydayLabels) +
   scale_colour_manual(values = c("violetred1", "black", "steelblue")) +
   labs(x = NULL, y = "Cumulative rainfall (mm)", colour = NULL))
#ggsave("cumRainfall.pdf", width = 7, height = 7)

# plot climex
# Keep as separate panels because it's nice and easy to read
plotClimex <- climexSiteMeansMelt %>%
  filter(index %in% c("GI(w)", "HS(w/acc)")) %>%
  mutate(index = str_replace_all(
    index, c("GI\\(w\\)" = "Growth index", "HS\\(w/acc\\)" = "Heat stress index"))
  ) %>%
  ggplot(aes(x = yearDay, y = meanValue, colour = index,
             ymin = meanValue - seValue, ymax = meanValue + seValue)) +
  geom_ribbon(aes(group = index), fill = alpha("grey", 0.5), colour = "transparent") +
  geom_line() + 
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        #strip.text.x = element_blank(),
        strip.background = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.5, 0.95),
        legend.justification = "centre",
        plot.margin = unit(c(2, 2, 2, 5), "mm")) +
  labs(x = NULL, y = "CLIMEX weekly index", colour = NULL) +
  scale_x_continuous(breaks = ydayBreaks, labels = ydayLabels) +
  scale_y_continuous(breaks = seq(0, 100, by = 20) / 100) +
  scale_colour_manual(values = c(oranges[6], greys[6])) +
  facet_wrap(~year, ncol = 1)
#ggsave("climexFieldSites.pdf", width = 8, height = 6)

# plot together
g <- gridExtra::arrangeGrob(plotRain, plotClimex, ncol = 2)
ggsave("RainfallAndClimex181cellsDRAFT.pdf", g, width = 10, height = 5)






        
  
