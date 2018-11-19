# ######################################################################
# Analysis of datasets on timing of canola colonisation by P. xylostella
# K. Perry, 17/7/2018 

# Background:
# ... 2014-2016

# This script does:
# (1) plots histograms of P. xylostella head capsule width data
# (2) Back-predicts the timing of egg lay in each sentinel canola crop
# (3) Produces static plots of colonisation timing, for publication 
# (4) Plots a time series movie of geographic maps showing spatiotemporal
# ... colonisation dynamics
# 
# Notes:
# This script uses the workspace in `colonisation.Rproj`.

# TO DO:
# Change naming of `sitesData`. Code is run for this object 5 times through the script, meaning all need to be re-run every time something changes. Just save 5 objects.
# 
# NEXT: Check the site names/locations for 2015 data (iron out any bugs)

# #######################################################################

library(readxl)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(RColorBrewer)
library(geosphere)
library(raster)
library(dbmdev) # my package for modelling P. xylostella development
#> https://community.rstudio.com/t/latest-rtools-seems-to-be-incompatible-with-the-latest-r/7657

# ###########
# Functions #
# ###########

# a function to average trap catches over a time period with missing values and an accumulated value at the end
imputeMeanValues <- function(x) {
  
  # group vec into chunks containing consecutive NAs and the next number
  group <- cumsum(!is.na(x)) + 1 # direct piping doesn't work!
  group <- lag(group, 1, default = 1)
  
  # split-apply-combine
  split(x, group) %>%
    lapply(
      function(x){
        # impute mean values in a vector containing NAs and a single integer
        if (length(which(!is.na(x))) > 0) {
          ave <- x[!is.na(x)]/length(x)
          rep(ave, length(x))
        } else {
          # in a vector without an integer (trailing NAs), keep NAs to maintain length(x)
          rep(NA, length(x))
        }
      }) %>%
    unlist() %>%
    as.numeric()
}

yearFromSiteCode <- function(x){
  x %>% str_extract("\\d{4}")
}

# standard error
se <- function(x, na.rm = TRUE){
  if(na.rm) x <- x[!is.na(x)]
  sd(x) / sqrt(length(x))
}

# ###########################
# Read in and organise data #
# ###########################

# ==============================
# head capsule data for 3 years
# ==============================

colonisationPath <- file.path(
  "C:", "UserData", "Kym", "PhD", "Data", "colonisation"
  )
headCapsule2014 <- file.path(
  colonisationPath, "data2014", "colonisationData2014.xlsx"
  ) %>%
  read_excel(sheet = "headCapsuleWidthData2014")
headCapsule2015 <- file.path(
  colonisationPath, "data2015", "colonisationData2015.xlsx"
) %>%
  read_excel(sheet = "headCapsuleWidthData2015")
headCapsule2016 <- file.path(
  colonisationPath, "data2016", "colonisationData2016.xlsx"
) %>%
  read_excel(sheet = "headCapsuleWidthData2016")

headCapsuleData <-headCapsule2014 %>%
  bind_rows(headCapsule2015) %>%
  bind_rows(headCapsule2016)

# ------------------------------------------------------------------------------
# determine the instar/lifestage of each individual based on head capsule widths
# ------------------------------------------------------------------------------


# First plot a frequency distribution of head capsule widths
# ... as a kernel-smoothed density curve overlaid on a histogram (uses "Gaussian" as default)
binWidth <- 1
plotHeadCap <- headCapsuleData %>%
  mutate(year = yearFromSiteCode(siteCode),
         year = factor(year)) %>%
  ggplot(aes(x = headWidthEPU)) +
  #ggplot(aes(x = headWidthEPU, y = ..density..)) + # for use with geom_density
  geom_histogram(binwidth = binWidth, alpha = 1, 
                 fill = "lightgrey", colour = "darkgrey") +
  #geom_line(stat = "density", adjust = 0.25, colour = "black", size = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 0.6) +
  labs(x = "Head capsule width (Eyepiece Units)", y = "Density")

# visually determine (iteratively) the size range break points for larval instars
instarBreaks <- c(14, 21, 31, 38, 50) # range of data is 14-50. Include the lower and upper.
instarLabels <- 1:4
blue <- brewer.pal(9, "Blues")[5]
plotHeadCap +
  geom_vline(aes(xintercept = instarBreaks[2] - 0.5), linetype = 2, size = 0.7, colour = blue) +
  geom_vline(aes(xintercept = instarBreaks[3] - 0.5), linetype = 2, size = 0.7, colour = blue) +
  geom_vline(aes(xintercept = instarBreaks[4] - 0.5), linetype = 2, size = 0.7, colour = blue) +
  geom_text(aes(x = 17.5, y = 35, label = "I"),   size = 8) +
  geom_text(aes(x = 26,   y = 35, label = "II"),  size = 8) +
  geom_text(aes(x = 34.5, y = 35, label = "III"), size = 8) +
  geom_text(aes(x = 47,   y = 35, label = "IV"),  size = 8) +
  theme_bw() +
  theme(panel.grid = element_blank())
  
ggsave("plotHeadCapsuleWidthEPUDraft.pdf", height = 4.5, width = 7.5)


# determine the development stage of each instar from the raw hc width data
headCapsuleDataLifeStages <- headCapsuleData %>%
  mutate(instar = cut(headWidthEPU, include.lowest = TRUE,
                       breaks = instarBreaks, 
                       labels = instarLabels),
         instar = as.numeric(instar),
         # assign "pupa" to lifestage 5
         instar = ifelse(lifeStage == "pupa", 5, instar))

# ---------------------------------------------
# determine the modelling biofix for each site
# ---------------------------------------------

# At sites where juveniles were detected, collaborators were asked to 
# collect juvenile samples on an extra date to increase the sample (and confidence).

# Determining the biofix would involve, for each site:
# (i)   Taking the oldest instar for each date where larvae were detected
# (ii)  Back-predicting development until egg lay
# (iii) Determining the earliest median egg lay date, and
# (iv)  ... reporting those starting params (date, lifestage) as the biofix.

# This has been done semi-manually for all sites, and the biofixes used for the 
# analysis included in the `sitesData` metadata in the code below.
# Coding this step in (to automate the analysis) can be achieved but 
# would require some additional work. Move on for now.

# --- draft code ---
# headCapsuleDataLifeStages %>%
#   # filter(firstDetectionFlag == 1) %>% 
#   group_by(siteCode, dateCollection) %>%
#   summarise(nIndividuals = length(instar),
#             oldestInstar = max(instar))
# --- end draft code


# ======================================
# read in the sites metadata for 3 years
# ======================================

# these files contain the metadata and the biofixes for modelling

sites2014 <- file.path(
  colonisationPath, "data2014", "colonisationData2014.xlsx"
) %>%
  read_excel(sheet = "sitesMetaData2014")
sites2015 <- file.path(
  colonisationPath, "data2015", "colonisationData2015.xlsx"
) %>%
  read_excel(sheet = "sitesMetaData2015")
sites2016 <- file.path(
  colonisationPath, "data2016", "colonisationData2016.xlsx"
) %>%
  read_excel(sheet = "sitesMetaData2016")

sitesData <- sites2014 %>%
  bind_rows(sites2015) %>%
  bind_rows(sites2016)

# write a local file with sites metadata for use in other scripts (e.g. `mapFieldSites.R`)
sitesData %>%
  write_csv("colonisationSitesMetaData.csv")


# =========================================
# read in the colonisation data for 3 years
# =========================================

colonisation2014 <- file.path(
  colonisationPath, "data2014", "colonisationData2014.xlsx"
) %>%
  read_excel(sheet = "colonisationData2014")
colonisation2015 <- file.path(
  colonisationPath, "data2015", "colonisationData2015.xlsx"
) %>%
  read_excel(sheet = "colonisationData2015")
colonisation2016 <- file.path(
  colonisationPath, "data2016", "colonisationData2016.xlsx"
) %>%
  read_excel(sheet = "colonisationData2016")

# combine years and take a mean pheromone trap moth counts and larval counts per date per site
colonisationData <- colonisation2014 %>%
  bind_rows(colonisation2015) %>%
  bind_rows(colonisation2016) %>%
  group_by(name, siteLocation, siteCode, date) %>%
  summarise(nTraps = length(trapNumber),
            growthStage      = growthStage[1],
            mnMothCount      = mean(mothCount, na.rm = TRUE),
            mnLarvaeCount    = mean(larvaeCount, na.rm = TRUE),
            samplingEffort   = samplingEffort[1],
            samplingDistance = samplingDistance[1]) %>%
  mutate(year = lubridate::year(date),
         date = as.Date(date)) %>%
  ungroup()


# ========================================================
# Aggregate colonisation data into daily and weekly counts
# ========================================================

# First, add our predicted dates of initial egg lay in each crop

# ########################################################################
# Back-predict the timing of initial colonisation (egg lay) in each crop #
# ########################################################################

# For this you need for each site:
# (1) Lat/long
# (2) Biofix for modelling (age of the 1st juvenile cohort)
# (2) Daily max/min temperature data time series

# As output you want for each site: 
# (1) Date range for initial egg lay (median, earliest, latest)

# steps:
# (1) Extract a daily max/min temperature time series from point locations for each site 
# (2) Read time series data into a list (named using the SILO grid cell number)
# (3) Interpolate daily to hourly temperatures
# (4) Run models for all sites simultaneously

# ==========================================================================
# Extract daily temperature time series for point locations from raster data
# ==========================================================================

# We are using daily gridded temperature data from SILO for 2014-2016

# --------------------------------------------------
# Identify the raster grid cells for all field sites
# --------------------------------------------------

# read in a reference raster
rasterPath <- file.path("C:", "UserData", "Kym", "PhD", "Data", "climateData", 
                        "siloDailyGrid2016", "dailyMax_2013-2016")
rasterTemplate <- raster(file.path(rasterPath, "20130101_max.txt"))

# add the silo grid cell for each site
sitesData <- sitesData %>%
  mutate(siloGridCell = cellFromXY(
    rasterTemplate, matrix(c(longitude, latitude), ncol = 2)
  ))

# --------------------------------------------
# read climate time series data into a list
# --------------------------------------------

# Time series data were previously extracted and written to text files
# using a custom R script `dataDrill.R` by K Perry.

# read daily temperature files for grid cells of interest
dataDrillPath  <- file.path(
  "C:", "UserData", "Kym", "PhD", "Data", "climateData", "dailyDataDrill2016"
)

parseGridCellNumber <- function(x){
  x %>% str_extract("\\d{6}") 
}

dataDrillFilePaths <- dataDrillPath %>%
  list.files(pattern = "dd", full.names = TRUE) %>%
  .[parseGridCellNumber(.) %in% unique(sitesData$siloGridCell)] # keep only grid cells we want
 
dailyObsList <- dataDrillFilePaths %>%
  map(read.table, skip = 8, header = TRUE) %>%
  set_names(parseGridCellNumber(dataDrillFilePaths))

# ------------------------------------------------
# interpolate daily max/min to hourly temperatures
# ------------------------------------------------

# Note for below:
# The data drill files did not have a column containing the grid cell number.
# I needed to add this column to each df
# It was not possible to extract the grid value directly from the list names in map(),
# so I created a function that takes two arguments, and interated over
# the list and the names in parallel using map2().

# Function with 2 arguments.
addGridCellColumn <- function(df, value){
  df %>% mutate(siloGridCell = value)
}

# TURNED OFF SO IT DOESN'T ACCIDENTALLY RUN
# (Takes about 12 minutes to run)
hourlyObsList <- dailyObsList %>%
  map(dplyr::select, date, min = dMin, max = dMax) %>% # presently must rename vars. Fix this!
  map2(addGridCellColumn, .x =  ., .y = names(dailyObsList) %>% as.numeric()) %>%
  map(mutate, # need lats/long to calculate sunrise/set
      latitude  = yFromCell(rasterTemplate, siloGridCell),
      longitude = xFromCell(rasterTemplate, siloGridCell)) %>%
  map(msunriset) %>% # need sunrise/set times to interpolate hourly temperatures
  map(dplyr::rename, station = siloGridCell) %>% # presently must rename to station. Fix!!
  map(hourlyTemperatures)

# write an external file just in case you don't want to wait again
# hourlyObsList %>%
#   bind_rows() %>%
#   write.csv("hourlyTemperatureObsDataframe.csv")

# Note: 
# Above, you cannot obtain lats/longs for each grid cell by joining with sites metadata.
# This is because multiple lat/longs can map to grid cells, creating duplicate dates!
# hourlyTemperatures() is designed to throw an error if there are duplicate dates in the input df.  

# Interesting, there are a few overlapping datetimes in the data.
# because these were daylight savings!
# dupDates7 <- hourlyObsList[[7]] %>% 
#   filter(duplicated(datetime)) %>%  
#   dplyr::select(datetime) %>% unlist() as.numeric()
# hourlyObsList[[7]] %>% filter(datetime %in% dupDates7) %>% View()
# (this may explain the tryCatch warning! Possibly need to set timezone as CST, 
# (as daylight savings time don't match SA DS times. they prob match UTC zone.)


# ===============================================================
# Run P. xylostella development models for each site using biofix
# ===============================================================

# Here, we use Briere's model formula II (Briere, 1999)

# Set model parameters for P. xylostella
# We parameterised the model using a separate R script in package `nlsreg`
# note, the code presently uses the rownames for stage, but doesn't require the colnames for each param 
# (params are internally matched by position within the workhorse function fwdBriere)
devParams <- data.frame(
  egg        = c(0.0003592,   1.754,    35.08,     4.347),
  instar_1_2 = c(0.00036,    -2.9122,   32.4565,  60.5092),
  instar_3   = c(0.0009688,   0.7893,   32.04,    14.53),
  instar_4   = c(0.000674,   -1.22561,  32.04873, 15.27334),
  prePupa    = c(0.00181533,  3.963006, 33.04467,  4.317786),
  pupaMF     = c(0.000396921, 2.417172, 32.44556, 11.99131)
) %>%
  magrittr::set_rownames(c("a", "Tmin", "Tmax", "m")) %>%
  t()

# -----------------------------------------------------------
# Notes: The below code presently runs models three times for each site
# to predict the date range of egg lay based on three biofixes:
# startDev = 1 (stage was fully developed), startDev = 0 (stage dev was just commencing)
# startDev = 0.5 (mid stage)
# note that instars 1 & 2 are combined, so beginning mid-2nd instar needs startDev = 0.75 (an approximation)

# The modelling package needs some dev upgrades so that:
# A list of biofixes can be input and models run with one call to mlocDev.
# The input data.frame is kept, with a single column added for each input biofix 
# (possibly add a function argument to specify a single column, or all cols for each biofix input)
# Some notes on upgrades have been added to the jobsList in the package dir.
# ---------------------------------------------------------------------------

# run model three times to predict the date range of egg lay at each site 
# run models only for sites where larvae were detected
eggLayPredictionsMedian <- sitesData %>%
  filter(!is.na(startDate)) %>% # remove sites with no larvae
  dplyr::rename(site = siteCode,     # presently, need to rename variable to these names. Fix!! should be able to tell it what the names are.
                station = siloGridCell) %>%
  mutate(gens     = 1,        # at present, need a column for gens. Ok if want to specify different gens for diff sites.
         startDev = 0.5) %>%  # at present, need a column for startDev. Would better if `missing` to be able set a global value for all sites, without so many resrictions on the input data
  mlocDev(tempObsDfList = hourlyObsList,
          devParamsDf   = devParams,
          timedir = "rev",
          output  = "generations") %>%
  bind_rows() %>%
  # rename variables again for joining (clunky! - that's why model code should not prescribe col names)
  dplyr::select(siteCode = site, eggLayMedian = dev0, totalDaysMedian = totDays)

eggLayPredictionsEarliest <- sitesData %>%
  filter(!is.na(startDate)) %>% # remove sites with no larvae
  dplyr::rename(site = siteCode,     # presently, need to rename variable to these names. Fix!! should be able to tell it what the names are.
                station = siloGridCell) %>%
  mutate(gens     = 1,        # at present, need a column for gens. Ok if want to specify different gens for diff sites.
         startDev = 1) %>%  # at present, need a column for startDev. Would better if `missing` to be able set a global value for all sites, without so many resrictions on the input data
  mlocDev(tempObsDfList = hourlyObsList,
          devParamsDf   = devParams,
          timedir = "rev",
          output  = "generations") %>%
  bind_rows() %>%
  # rename variables again for joining 
  dplyr::select(siteCode = site, eggLayEarliest = dev0, totalDaysEarliest = totDays)

# (note, at the moment, you can't input startDev = 0 when direction = reverse! Fix!)
# there has to be a single dev increment within the start stage for it to work.
# so, we have given it startDev = 0.025.
eggLayPredictionsLatest <- sitesData %>%
  filter(!is.na(startDate)) %>% # remove sites with no larvae
  dplyr::rename(site = siteCode,     # presently, need to rename variable to these names. Fix!! should be able to tell it what the names are.
                station = siloGridCell) %>%
  mutate(gens     = 1,        # at present, need a column for gens. Ok if want to specify different gens for diff sites.
         startDev = 0.025) %>%  # at present, need a column for startDev. Would better if `missing` to be able set a global value for all sites, without so many resrictions on the input data
  mlocDev(tempObsDfList = hourlyObsList,
          devParamsDf   = devParams,
          timedir = "rev",
          output  = "generations") %>%
  bind_rows() %>%
  # rename variables again for joining 
  dplyr::select(siteCode = site, eggLayLatest = dev0, totalDaysLatest = totDays)

# join predictions with original data to restore the sites without larvae
sitesDataEggLay <- eggLayPredictionsEarliest %>%
  left_join(eggLayPredictionsMedian) %>% 
  left_join(eggLayPredictionsLatest) %>% 
  left_join(sitesData, .)


# ######################################################################################
# Aggregate colonisation data into weekly counts of moth capture and juvenile presence #
# ###################################################################################### 

# Fit a full daily time series to each site, then aggregate to weekly counts
# Find the full date range across all years
colonisationData %>%
  split(.$year) %>%
  lapply(function(x){range(x$date)}) %>%
  bind_rows() %>%
  mutate(measure = c("min", "max")) %>%
  reshape2::melt(variable.name = "year", value.name = "date") %>%
  mutate(week = cut(date, breaks = "week", start.on.monday = FALSE)) # use cut so we know when sequence starts on sunday
#> full time sequence should be 19th April to 14th October

# create the full time sequence for each year
dateSeq2014 <- seq(as.Date("2014-04-19"), as.Date("2014-10-14"), by = "day")
dateSeq2015 <- seq(as.Date("2015-04-19"), as.Date("2015-10-14"), by = "day")
dateSeq2016 <- seq(as.Date("2016-04-19"), as.Date("2016-10-14"), by = "day")

# ---------------------------
# disaggregate to daily data 
# ---------------------------

# Add a variable describing the presence or absence of 
# juveniles in the crop from our modelling results
# note: this is the back-predicted first egg lay date, NOT the first detection of larvae! 
juvenilesPresent <- function(date, eggLayDate, larvaeCount){
  
  startedSampling <- min(which(!is.na(larvaeCount)))
  out <- rep(0, length(date))
  out[seq(startedSampling) - 1] <- NA  # return NA values until the first larval sampling date
  
  # if no larvae were detected, return zero counts starting from first sampling date
  if(all(is.na(eggLayDate))) {return(out)}
  out[which(date >= eggLayDate)]  <- 1
  out
}

dupRows <- function(df, n){
  nObs <- nrow(df)
  ind <- rep(1:nObs, each = n)
  df[ind, ]
  }

colonisationDaily <- colonisationData %>%
  split(.$siteCode) %>%
  map(function(x){
    siteData <- x %>%
      dplyr::select("name", "siteLocation", "siteCode", "year") %>%
      distinct()
    year <- unique(x$year)
    if(year == 2014){
      dateSeq <- dateSeq2014
      }
    else if(year == 2015){
      dateSeq <- dateSeq2015
    }
    else if(year == 2016){
      dateSeq <- dateSeq2016
    }
    siteData %>% 
      dupRows(n = length(dateSeq)) %>%
      mutate(date = dateSeq) %>%
      left_join(x)
  }) %>%
  map(left_join, sitesDataEggLay %>% 
        dplyr::select(siteCode, eggLayMedian, dateTrappingStart)) %>%
  map(arrange, year, siteCode, date) %>%
  map(mutate, juvPresence = juvenilesPresent(date, eggLayMedian, mnLarvaeCount)) %>%
  bind_rows()


# ---------------------------------------------------------------------
# Aggregate to weekly data and average moth counts across missing weeks
# ---------------------------------------------------------------------

# Function to sum non-NA values but return NA if there are no integers
sumValues <- function(x){
  if(all(is.na(x))){return(NA)}
  sum(x, na.rm = TRUE)
}

maxValues <- function(x){
  if(all(is.na(x))){return(NA)}
  max(x, na.rm = TRUE)
}

# Function to average pheromone trap counts across missing values,
# taking into account the date trapping commenced
averageMothCounts <- function(date, dateTrappingStart, mothCount) {
  # for sites with no trapping start date shown, assume 1 week earlier
  if(all(is.na(dateTrappingStart))){
    firstSamplingDate <- as.Date(date[min(which(!is.na(mothCount)))])
    dateTrappingStart <- firstSamplingDate - 7
  }
  mothCount[which(date <= dateTrappingStart)] <- 0 # a place `stopper` for imputeMeanValues
  meanCount <- imputeMeanValues(mothCount)
  meanCount[which(date <= dateTrappingStart)] <- NA # remove stopper, restore NAs.
  meanCount
}

colonisationWeekly <- colonisationDaily %>%
  # cut dates at the same 7-day intervals across years 
  split(.$year) %>%
  map(arrange, year, siteCode, date) %>%
  map(mutate, 
      week  = cut(date, breaks = "7 days"),
      week  = as.Date(week),
      day   = lubridate::day(week),
      month = lubridate::month(week),
      # create a `day-month` variable for mapping the same weeks across years
      weekOfYear = format(week, "%d-%b")
      ) %>%
  map(group_by, 
      year, name, siteCode, siteLocation, week, weekOfYear,
      eggLayMedian, dateTrappingStart) %>% # we've dropped growth stage, samplingEffort, samplingDistance vars
  map(summarise,
      nSamples = length(siteCode),
      nMoths  = sumValues(mnMothCount),
      nLarvae = sumValues(mnLarvaeCount),
      juvPresence = maxValues(juvPresence)) %>%
  map(ungroup) %>%
  bind_rows() %>%
  split(.$siteCode) %>%
  map(mutate, aveMoths = averageMothCounts(week, dateTrappingStart, nMoths)) %>%
  map(arrange, year, siteCode, week) %>%
  bind_rows() %>%
  left_join(sitesDataEggLay) # join the sites data for plotting geographic maps 

# write a file of weekly colonisation data (for use in other scripts - eg plotting on top of climex rasters)
colonisationWeekly %>%
  write_csv("colonisationWeekly.csv")

# TO DO:
# possibility to select fewer columns in sitesDataEggLay when join to colonisationWeekly
# for sites without larvae detected, currently filled zeros for juvenile presence
# starting fro the first sampling date. could make this fill NAs after last sampling date
# (edit the juvenilesPresent() function)


# #################################################################
# Create static plots to compare colonisation timing in each year #
# #################################################################

# ========================================================
# Plot the average weekly moth abundance across all sites
# ========================================================

# colonisationWeekly %>%
#   filter(!is.na(aveMoths),
#          lubridate::month(week) %in% c(5, 6, 7)) %>%
#   group_by(year, weekOfYear) %>%
#   summarise(week = unique(week), # keep the week variable to allow correct date sorting
#             nSitesWithData = length(weekOfYear),
#             meanMothsPerSite = mean(aveMoths),
#             seMothsPerSite   = se(aveMoths)) %>% 
#   ungroup() %>%
#   arrange(year, week) %>% 
#   mutate(year = factor(year),
#          weekOfYear = reorder(weekOfYear, week)) %>%
#   ggplot(aes(x = weekOfYear, y = meanMothsPerSite, linetype = year, group = year,
#              ymin = meanMothsPerSite - seMothsPerSite, 
#              ymax = meanMothsPerSite + seMothsPerSite)) +
#   geom_line() +
#   geom_errorbar(width = 0.1, size = 0.5) +  
#   geom_point() +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1),
#         panel.grid  = element_blank(),
#         aspect.ratio = 0.6) +
#   labs(x = "Week of year", y = "Mean +- SEM weekly moth count per site", linetype = NULL)
# ggsave("plotMeanMothsPerSiteDRAFT.pdf", width = 11, height = 7)

# =============
# Crossbar plot
# =============

# Here we need to use a numeric axis based on `day of year` to compare dates on a common axis across years.
# Then, manually annotate the date label on the y-axs. Do this by 
# (i)   checking the dates vs day of year for 2014, 2015, 2016
# (ii)  there was a leap year in 2016, so manually offset the number for 2016 (ie yday = yday + 1)
# (iii) set breaks and labels manually on the y-axis using scale_y_xcontinuous() (make labels into dates)

# set breaks and labels
# determine the offset for `day of year` for the same date each year
data_frame(date = as_date(c("2014-02-01", "2015-02-01", "2016-02-01")),
           yday = yday(date))
# # A tibble: 3 x 2
# date        yday
# <date>     <dbl>
#   1 2014-05-01   121
# 2 2015-05-01   121
# 3 2016-05-01   122 # leap year in 2016
#> Shows that 2016 is one day later.
# Therefore, to generate a common `date` x-axis labels across years from day of year, offset 2016 by 1 day 

# now set breaks as 1st day of the month, then set labels
# ok to use 2014 as a sequence, as yday has been adjusted for 2016 so that dates match across years. 
dateSeq <- seq(as_date("2014-03-01"), as_date("2014-11-01"), by = "month")
ydayBreaks <- dateSeq %>% yday() 
ydayLabels <- dateSeq %>% format("%d-%b")

# function to adjust (add 1) day of year for 2016 dates only
ydayAdjust <- function(x, adjustYear = 2016, adjustment = 1){
  y  <- year(x)
  yd <- yday(x) 
  yd[which(y == adjustYear)] <- yd[which(y == adjustYear)] + adjustment
  yd
}

# to add rows `NA` row for each year as space for a label
# siteSpacer <- colonisationWeekly %>%
#   group_by(year) %>%
#   summarise_all(function(x){x[1]}) %>%  #take the first 2 rows for each year
#   #dplyr::select(contains("eggLay"))
#   mutate_at(vars(-contains("year")), .funs = function(x){NA}) %>%
#   # bind_rows(., .) %>% # duplicate rows to give two sites spacing per year
#   # convert egg lay cols to NA
#   mutate(siteCode = paste("dummy", 2014:2016))
#   # two rows mutate(siteCode = paste("dummy", rep(2014:2016, 2), c(rep("a", 3), rep("b", 3)), sep = "_"))
# #   

stageLabels <- c("instar_1_2" = "1st & 2nd instar",
                 "instar_3" = "3rd instar",
                 "instar_4" = "4th instar",
                 "pupaMF"   = "Pupa")

colonisationWeekly %>%
  filter(!is.na(eggLayMedian)) %>% 
  # bind_rows(siteSpacer) %>%
  arrange(desc(year), eggLayMedian) %>%
  mutate(plotOrder = 1:nrow(.),
         siteCode  = reorder(siteCode, plotOrder),
         year = factor(year),
         startStage = str_replace_all(startStage, stageLabels)) %>% 
  ggplot(aes(x = siteCode, y = ydayAdjust(eggLayMedian),
             ymin = ydayAdjust(eggLayEarliest), ymax = ydayAdjust(eggLayLatest))) +
  geom_boxplot(#aes(group = year),
               outlier.shape = NA,
               fill = alpha("lightgrey", 0.5)) +
  # add the back-predicted median and date range (earliest and latest) egg lay 
  geom_crossbar(fill = "white") + # FOR SOME REASON IT NEED THIS EARLIER CALL!
  # split the years using a vertical line
  # for categorical data (site), map the line to factor index to plot at the correct point)
  geom_vline(aes(xintercept = 43.5), colour = "black") +
  geom_vline(aes(xintercept = 16.5), colour = "black") +
  geom_text(aes(x = 46, y = 280, label = "2014"), size = 6) +
  geom_text(aes(x = 18, y = 280, label = "2015"), size = 6) +
  geom_text(aes(x =  2, y = 280, label = "2016"), size = 6) +
  # add a line linking the dot and the crossbar (for points over-plotted below)
  geom_segment(aes(x = siteCode, xend = siteCode,
                   y = ydayAdjust(eggLayLatest), yend = ydayAdjust(startDate)),
               linetype = 1, size = 0.2) +
  # add the first detection date and lifestage
  geom_point(aes(x = siteCode, y = ydayAdjust(startDate), fill = startStage),
             pch = 21, colour = "black") + 
  # # add another line linking the sowing date with the crossbar (for points overplotted below)
  # geom_segment(aes(x = siteCode, xend = siteCode,
  #                  y = ydayAdjust(dateSowing), yend = ydayAdjust(eggLayEarliest)),
  #              linetype = 1, size = 0.2) +
  # add a symbol showing the crop sowing date
  geom_point(aes(x = siteCode, y = ydayAdjust(dateSowing)), 
             size = 0.7, colour = "black") +
  theme_bw(base_size = 16) +
  theme(axis.text.y = element_blank(), # turn off site labels
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = c(0.9, 0.8)) +
  labs(x = "Site", y = NULL, fill = "Life stage biofix") +
  #scale_colour_discrete(limits = c("instar_1_2", "instar_3", "instar_4", "pupaMF")) +
  scale_y_continuous(breaks = ydayBreaks, labels = ydayLabels) +
  coord_flip()
ggsave("plotColonisationCrossbarDraft.pdf", width = 12, height = 8) 
#thesisPathChap5 <- file.path("C:", "UserData", "Kym", "PhD", "thesis", "Chapter5", "Figs")
#file.path(thesisPathChap5, "plotColonisationCrossbarDraft.pdf") %>%
 # ggsave(width = 9, height = 6) 

# Make a box plot
colonisationWeekly %>%
  filter(!is.na(eggLayMedian)) %>%
  group_by(siteCode, year) %>%
  summarise(eggLayMedian = unique(eggLayMedian)) %>%
  ungroup() %>%
  mutate(doy = yday(eggLayMedian)) %>%
  ggplot(aes(x = siteCode, y = doy, colour = factor(year))) +
  geom_boxplot() +
  coord_flip()
#> Later...


## TO FIX!!!
# Check your instar_1_2 models - these need a different biofix! median egg_lay for 1_2 should start with a biofix at 0.75!!!



# ##########################################################################
# Create static panel plots for each site describing dynamics at each site #
# ##########################################################################

# For each site, we want to produce a plot that shows:
# (a) Climex Weekly Growth Index for several generations preceding the first colonisation (egg lay) date
# (b) Mean weekly moth counts in pheromone traps
# (c) Crop sowing date (or emergence date)
# (d) First colonisation date

# ============================================================================
# Back-predict the dates for 5 generations before colonisation (first egg lay)
# ============================================================================

# In 2014, we may not have a long enough time series of data to back-predict > 3 generations

# set the development parameters for P. xylostella for the entire lifecycle
# We parameterised Briere's model using R package `nlsreg`.
devParamsTotal <- data.frame(
  lifecycle = c(0.0000938, 1.042, 32.73, 9.943)) %>%
  magrittr::set_rownames(c("a", "Tmin", "Tmax", "m")) %>%
  t()

# Back-predict development rates for 3 generations at sites where larvae were detected
# (note: the sites without larvae do not have a biofix at this stage - we could add one to each site as a guide.)
backPredictions5Generations <- sitesDataEggLay %>%
  filter(!is.na(startDate)) %>%
  dplyr::rename(site    = siteCode,
                station = siloGridCell) %>%
  mutate(gens = 5,
         startStage = "lifecycle",  # set biofix to start at the end of lifecycle (must match rownames(devParamsTotal))
         startDate  = eggLayMedian, # set biofix to back-predict from the date of first egg lay
         startDev   = 1) %>%        # set biofix to begin at the completion of the lifecycle (adult emergence, approximating egg lay)
  mlocDev(tempObsDfList = hourlyObsList,
          devParamsDf   = devParamsTotal,
          timedir = "rev",
          output  = "generations") %>%
  bind_rows()
#> A good idea to visually check that back predictions are sensible (total days v aveTemp is a decent guide)

# cast each generation to a separate variable and rejoin with the original data
sitesDataEggLay5GensBack <- backPredictions5Generations %>%
  reshape2::dcast(site ~ gen, value.var = "dev0") %>% # dev0 is the commencement date for each generation 
  # convert numeric data to a datetime again
  map_dfr(function(x){
    if(!is.numeric(x)){return(x)}
    as_datetime(x, origin = lubridate::origin, tz = "UTC")
  }) %>%
  set_names(c("site", "back1Gen", "back2Gens", "back3Gens", "back4Gens", "back5Gens")) %>%
  dplyr::rename(siteCode = site) %>%
  left_join(sitesDataEggLay, ., by = "siteCode")


# ===========================================================================================
# Read in the output of Climex Compare Years `Weekly Growth Index` predictions for each site
# ===========================================================================================

climexDir  <- file.path(
  "C:", "UserData", "Kym", "PhD", "Data", "climexOutput", "climexCompareYearsPointLocations"
  )
climexFilePaths <- climexDir %>%
  list.files(pattern = ".csv", full.names = TRUE) %>%
  .[parseGridCellNumber(.) %in% unique(sitesData$siloGridCell)] # keep only grid cells we want

gridCellNumbers <- parseGridCellNumber(climexFilePaths)

climexResults <- climexFilePaths %>%
  map(read_csv, skip = 3) %>% 
  set_names(gridCellNumbers) %>%
  # add a column with the grid cell number
  map2(addGridCellColumn, .x =  ., .y = gridCellNumbers %>% as.numeric()) %>%
  map(mutate, 
      date = lubridate::dmy(SimulaDate),
      date = lubridate::as_datetime(date)) %>% # check that we still need this extra step to correctly offset the seconds
  map(dplyr::select, siloGridCell, date, weeklyGI = `GI(w)`) %>%
  bind_rows() %>%
  # add a site code to the data: this will increase number of rows as some grid cells contain more than site (this is fine!)
  left_join(sitesData %>% dplyr::select(siloGridCell, siteCode), .) %>%
  arrange(siteCode, date)


# ===============================================================================
# Find the Climex weekly GI value on the back-predicted dates for each generation
# ===============================================================================

# This will allow us to visualise the generation times more clearly against climex values

# First, diaggregate climex weekly index to daily values
climexResultsDaily <- climexResults %>%
  split(.$siteCode) %>%
  # for now, just keep the time series for the same year as the site-code (this can be altered if multi-year analysis becomes necessary)
  map(function(x){
    siteYear <- yearFromSiteCode(x$siteCode) %>% unique()
    x %>% filter(lubridate::year(date) == siteYear)
  }) %>%
  map(function(x){
    dateSeq <- seq(min(x$date), max(x$date), by = "day")
    x %>%
      dplyr::select(siloGridCell, siteCode) %>%
      distinct() %>%
      dupRows(n = length(dateSeq)) %>%
      mutate(date = dateSeq) %>%
      full_join(x, by = c("siloGridCell", "siteCode", "date")) %>%
      arrange(date) %>%
      mutate(weeklyGIInterpolated = zoo::na.approx(weeklyGI))
  }) %>%
  map(mutate, date = as_date(date)) %>% # for convenient joining below
  bind_rows()


# make an additional data frame with (for each site) the the egg lay predictions for 5 gens
# and the (daily interpolated) climex value for those dates.
sitesDataEggLay5GensBackWithClimex <- sitesDataEggLay5GensBack %>%
  reshape2::melt(
    id.vars = c("siteCode"),
    measure.vars = c("eggLayMedian", # include the egg lay date as well 
                     "back1Gen", "back2Gens", "back3Gens", "back4Gens", "back5Gens"),
    variable.name = "gensBack",
    value.name = "date"
  ) %>%
  mutate(date = as_date(date)) %>%
  split(.$siteCode) %>%
  map(function(x){
    tempDfClimex <- climexResultsDaily %>% 
      filter(siteCode %in% x$siteCode, 
             date     %in% lubridate::date(x$date)) %>%
      dplyr::select(siteCode, date, weeklyGIInterpolated)
    left_join(x, tempDfClimex) # my function does *not* fail for sites without predicted dates (no larvae) detected, which is great.
  }) %>%
  map(mutate, date = as_datetime(date, tz = "UTC")) %>% # convert back to posixct date times for correct plotting along x-axis (to understand this step, check as.numeric() with and without this conversion.)
  bind_rows()

# =========================================
# Plot each site individually as a `story`
# =========================================

# Normalise moth counts to a proportion of maximum moth counts. This will:
# (i)  allow plotting on the same axis as climex Weekly GI, and 
# (ii) remove the need for a `free_y` scale, which would require extra space for y-labels on each panel
# Note that normalising moth counts means that you lose the abundance data, therefore
# the spatio-tenporal time series movie and the line graph of average abundance across sites above are important figures.


normaliseMothCounts <- function(x){
  maxCount <- max(x, na.rm = TRUE)
  x[!is.na(x)] <- x[!is.na(x)] / maxCount
  x
}

orange <- rgb(1, 0.5, 0, 1)
red  <- "red" #brewer.pal(9, "Reds")[9]
blue <- brewer.pal(9, "Blues")[5]

# ----------
# 2014 sites
# ----------
colonisationWeekly %>% 
  split(.$siteCode) %>%
  map(mutate, 
      mothCountNormalised = normaliseMothCounts(aveMoths),
      week = as_datetime(week)) %>% # weeks have been stored incorrectly! (they encode time near 1970!)
  bind_rows() %>%# get the moth counts from weekly data 
  filter(year == 2014) %>%
  ggplot(aes(x = week, y = mothCountNormalised)) +
  # add moth counts
  geom_line()  +
  # add climex weekly growth index
  geom_line(data = climexResults %>% # keep the sites and dates relevant to 2014 (otherwise everything will be plotted)
              filter(yearFromSiteCode(siteCode) == 2014, # keep only sites monitored in 2014
                     lubridate::year(date) == 2014,      # keep only time series for 2014 (each yearly site has multi-year climex predictions)
                     lubridate::month(date) >= 1, lubridate::month(date) < 10), # Keep Jan to September
            aes(x = date, y = weeklyGI), linetype = 3, colour = alpha("black", 0.8)) +
  # add the crop sowing date
  geom_vline(aes(xintercept = dateSowing),   linetype = 1, colour = blue) +
  # add the median date of initial egg lay, and plot in red dots
  geom_point(data = sitesDataEggLay5GensBackWithClimex %>% 
               filter(yearFromSiteCode(siteCode) == 2014,
                      gensBack == "eggLayMedian"), # plot just the eggLay data in red dots dots
             aes(x = date, y = weeklyGIInterpolated), colour = red) +
  #geom_vline(aes(xintercept = eggLayMedian), linetype = 1, colour = alpha("red", 1)) +
  #plot gen 1 back, placing point on the predicted climex weekly GI for that date
  geom_point(data = sitesDataEggLay5GensBackWithClimex %>%
               filter(yearFromSiteCode(siteCode) == 2014,
                      gensBack != "eggLayMedian"), # just plot the previous generations in black dots
             aes(x = date, y = weeklyGIInterpolated), colour = "black", size = 0.7) +
  # now add the vertical segments back to the x-axis for each generation
  geom_segment(data = sitesDataEggLay5GensBackWithClimex %>% 
                 filter(yearFromSiteCode(siteCode) == 2014,
                        gensBack != "eggLayMedian"), # just plot the previous generations in black dots
               aes(x = date, xend = date,
                   y = -Inf, yend = weeklyGIInterpolated), colour = "black", linetype = 2) +
  # now add the vertical segment for the egg lay date
  geom_segment(data = sitesDataEggLay5GensBackWithClimex %>% 
                 filter(yearFromSiteCode(siteCode) == 2014,
                        gensBack == "eggLayMedian"), # plot just the eggLay data in red dots dots
             aes(x = date, xend = date,
                 y = -Inf, yend = weeklyGIInterpolated), colour = red, linetype = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.7) +
  labs(x = NULL, y = "Moths trapped per week (solid line) scaled to the CLIMEX Weekly Growth Index (dotted line)") +
  facet_wrap(~siteCode, ncol = 6) +
  scale_x_datetime(breaks = scales::date_breaks(width = "1 month"),
                   labels = scales::date_format(format = "%b-%Y"))
ggsave("plots2014Draft.pdf", width = 10.25, height = 7.75)

# ----------
# 2015 sites
# ----------
colonisationWeekly %>% 
  split(.$siteCode) %>%
  map(mutate, 
      mothCountNormalised = normaliseMothCounts(aveMoths),
      week = as_datetime(week)) %>% # weeks have been stored incorrectly! (they encode time near 1970!)
  bind_rows() %>%# get the moth counts from weekly data 
  filter(year == 2015) %>%
  ggplot(aes(x = week, y = mothCountNormalised)) +
  # add moth counts
  geom_line()  +
  # add climex weekly growth index
  geom_line(data = climexResults %>% # keep the sites and dates relevant to 2014 (otherwise everything will be plotted)
              filter(yearFromSiteCode(siteCode) == 2015, # keep only sites monitored in 2014
                     lubridate::year(date) == 2015,      # keep only time series for 2014 (each yearly site has multi-year climex predictions)
                     lubridate::month(date) >= 1, lubridate::month(date) < 10), # Keep Jan to September
            aes(x = date, y = weeklyGI), linetype = 3, colour = alpha("black", 0.8)) +
  # add the crop sowing date
  geom_vline(aes(xintercept = dateSowing),   linetype = 1, colour = blue) +
  # add the median date of initial egg lay, and plot in red dots
  geom_point(data = sitesDataEggLay5GensBackWithClimex %>% 
               filter(yearFromSiteCode(siteCode) == 2015,
                      gensBack == "eggLayMedian"), # plot just the eggLay data in red dots dots
             aes(x = date, y = weeklyGIInterpolated), colour = red) +
  #geom_vline(aes(xintercept = eggLayMedian), linetype = 1, colour = alpha("red", 1)) +
  #plot gen 1 back, placing point on the predicted climex weekly GI for that date
  geom_point(data = sitesDataEggLay5GensBackWithClimex %>%
               filter(yearFromSiteCode(siteCode) == 2015,
                      gensBack != "eggLayMedian"), # just plot the previous generations in black dots
             aes(x = date, y = weeklyGIInterpolated), colour = "black", size = 0.7) +
  # now add the vertical segments back to the x-axis for each generation
  geom_segment(data = sitesDataEggLay5GensBackWithClimex %>% 
                 filter(yearFromSiteCode(siteCode) == 2015,
                        gensBack != "eggLayMedian"), # just plot the previous generations in black dots
               aes(x = date, xend = date,
                   y = -Inf, yend = weeklyGIInterpolated), colour = "black", linetype = 2) +
  # now add the vertical segment for the egg lay date
  geom_segment(data = sitesDataEggLay5GensBackWithClimex %>% 
                 filter(yearFromSiteCode(siteCode) == 2015,
                        gensBack == "eggLayMedian"), # plot just the eggLay data in red dots dots
               aes(x = date, xend = date,
                   y = -Inf, yend = weeklyGIInterpolated), colour = red, linetype = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.7) +
  labs(x = NULL, y = "Moths trapped per week (solid line) scaled to the CLIMEX Weekly Growth Index (dotted line)") +
  facet_wrap(~siteCode, ncol = 6) +
  scale_x_datetime(breaks = scales::date_breaks(width = "1 month"),
                   labels = scales::date_format(format = "%b-%Y"))
ggsave("plots2015Draft.pdf", width = 10.25, height = 7.75)

# ----------
# 2016 sites
# ----------
colonisationWeekly %>% 
  split(.$siteCode) %>%
  map(mutate, 
      mothCountNormalised = normaliseMothCounts(aveMoths),
      week = as_datetime(week)) %>% # weeks have been stored incorrectly! (they encode time near 1970!)
  bind_rows() %>%# get the moth counts from weekly data 
  filter(year == 2016) %>%
  ggplot(aes(x = week, y = mothCountNormalised)) +
  # add moth counts
  geom_line()  +
  # add climex weekly growth index
  geom_line(data = climexResults %>% # keep the sites and dates relevant to 2014 (otherwise everything will be plotted)
              filter(yearFromSiteCode(siteCode) == 2016, # keep only sites monitored in 2014
                     lubridate::year(date) == 2016,      # keep only time series for 2014 (each yearly site has multi-year climex predictions)
                     lubridate::month(date) >= 1, lubridate::month(date) < 10), # Keep Jan to September
            aes(x = date, y = weeklyGI), linetype = 3, colour = alpha("black", 0.8)) +
  # add the crop sowing date
  geom_vline(aes(xintercept = dateSowing),   linetype = 1, colour = blue) +
  # add the median date of initial egg lay, and plot in red dots
  geom_point(data = sitesDataEggLay5GensBackWithClimex %>% 
               filter(yearFromSiteCode(siteCode) == 2016,
                      gensBack == "eggLayMedian"), # plot just the eggLay data in red dots dots
             aes(x = date, y = weeklyGIInterpolated), colour = red) +
  #geom_vline(aes(xintercept = eggLayMedian), linetype = 1, colour = alpha("red", 1)) +
  #plot gen 1 back, placing point on the predicted climex weekly GI for that date
  geom_point(data = sitesDataEggLay5GensBackWithClimex %>%
               filter(yearFromSiteCode(siteCode) == 2016,
                      gensBack != "eggLayMedian"), # just plot the previous generations in black dots
             aes(x = date, y = weeklyGIInterpolated), colour = "black", size = 0.7) +
  # now add the vertical segments back to the x-axis for each generation
  geom_segment(data = sitesDataEggLay5GensBackWithClimex %>% 
                 filter(yearFromSiteCode(siteCode) == 2016,
                        gensBack != "eggLayMedian"), # just plot the previous generations in black dots
               aes(x = date, xend = date,
                   y = -Inf, yend = weeklyGIInterpolated), colour = "black", linetype = 2) +
  # now add the vertical segment for the egg lay date
  geom_segment(data = sitesDataEggLay5GensBackWithClimex %>% 
                 filter(yearFromSiteCode(siteCode) == 2016,
                        gensBack == "eggLayMedian"), # plot just the eggLay data in red dots dots
               aes(x = date, xend = date,
                   y = -Inf, yend = weeklyGIInterpolated), colour = red, linetype = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.7) +
  labs(x = NULL, y = "Moths trapped per week (solid line) scaled to the CLIMEX Weekly Growth Index (dotted line)") +
  facet_wrap(~siteCode, ncol = 6) +
  scale_x_datetime(breaks = scales::date_breaks(width = "1 month"),
                   labels = scales::date_format(format = "%b-%Y"))
ggsave("plots2016Draft.pdf", width = 10.25, height = 7.75)



# TO do:
# Checking sowing dates in th raw data (guess where absent)
# Estimate harvest dates for each site (perhaps by region - South East, typically xxx?)
# Cast forward x generations from first egg lay until harvest
# the package needs to be able to calculate the number of generations for a set time window!

# For the manuscript, do:
# Calculate grand median of egg lay dates +- median for each year? (try raw numbers, otherwise present a boxplot)
# Calculate number of days betweeen crop sowing and first egg lay 
# Graph of proportion of crops colonised over time?
# 

# Cumulative rainfall over (a) a set time period. 
# Perhaps overlay vertical lines for DBM generations for each year as well, to highlight temperature conditions in lead up to sowing?

# =======
# boxplot
# =======

sitesDataEggLay %>%
  mutate(year = yearFromSiteCode(siteCode)) %>%
  split(.$year) %>%
  map_dfr(function(df){
    data_frame(year = lubridate::year(median(df$eggLayMedian, na.rm = TRUE)),
               medianEggLay = median(df$eggLayMedian, na.rm = TRUE),
               sdEggLay     = as_datetime(sd(df$eggLayMedian, na.rm = TRUE), origin = lubridate::origin),
               meanEggLay   = mean(df$eggLayMedian, na.rm = TRUE))
  })

sitesDataEggLay %>%
  mutate(year = yearFromSiteCode(siteCode),
         year = factor(year),
         dayOfYear  = lubridate::yday(eggLayMedian)) %>%
  ggplot(aes(x = year, y = dayOfYear)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.2) +
  coord_flip()


  



# #####################################################
# Plot time series geographic maps of the moth counts #
# #####################################################


# ===================================================================
# Read in a large Australia polygon with states and local admin areas
# ===================================================================

geoSpatialDataPath <- file.path("C:", "UserData", "Kym", "PhD", "Data", "geoSpatialData")
ausPolygon <- file.path(geoSpatialDataPath, "ausAdminAreas_shp") %>%
  rgdal::readOGR(dsn = ., layer = "AUS_adm2")

# convert to a standard data.frame for ggplot2
ausPolygonFortified <- fortify(ausPolygon)

# ====================================================
# Read in a much smaller Australia polygon with states (faster processing) 
# ====================================================

ausPolygonSmall <- file.path(geoSpatialDataPath, "nsaasr9nnd_02211a04ec_alb132") %>%
  rgdal::readOGR(dsn = ., layer = "aust_cd66states")
ausPolygonSmallFortified <- fortify(ausPolygonSmall) # convert to df for ggplot


# ========================================================
# write plot frames in a loop (for every week of the year)
# ========================================================

myPal <- brewer.pal(9, "Greens")[8] # use green circle
yBreaks <- seq(-38, -32, by = 2)
yLabels <- parse(text = paste(abs(yBreaks), "*degree ~ S", sep = ""))
xBreaks <- seq(134, 142, by = 2)
xLabels <- parse(text = paste(xBreaks, "*degree ~ E", sep = ""))
sizeBreaks <- c(0, 1, 10, 50, 100, 250)

counter <- c(0)
for(w in unique(colonisationWeekly$weekOfYear)){ # 5:20 is a good subset for testing
  
  #tmpDf <- colonisationWeekly[colonisationWeekly$weekOfYear == w, ]
  tmpDf <- colonisationWeekly %>%
    filter(weekOfYear == w) %>%
    mutate(group = 1,
           juvPresence = factor(juvPresence))
  
  
  # # edit the panels so they don't jump
  # # see: https://github.com/tidyverse/ggplot2/issues/2288
  # temp <- editGrob(text_grob, label = "gjpqyQ")
  # descent <- descentDetails(temp)
  
  ausPolygonSmallFortified %>%
    ggplot(aes(x = long, y = lat, group = group)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          strip.background = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle(w) +
    labs(x = "Longitude", y = "Latitude") +
    geom_polygon(fill = "lightgrey", colour = "white") +
    geom_point(data = tmpDf %>% filter(nMoths > 0),
               aes(x = longitude, y = latitude, size = nMoths, alpha = nMoths),
               pch = 21, fill = myPal, colour = myPal) +
    # map the zero values as open circles
    geom_point(data = tmpDf %>% mutate(group = 1) %>% filter(nMoths == 0),
               aes(x = longitude, y = latitude),
               pch = 21, fill = "lightgrey", colour = myPal) +
    geom_point(data = tmpDf %>% 
                 filter(!is.na(as.numeric(juvPresence))), 
               aes(x = longitude, y = latitude, fill = juvPresence),
               pch = 21, colour = "black") +
    coord_map(ylim = c(-38, -32), xlim = c(133.5, 141.5)) +
    scale_x_continuous(breaks = xBreaks, labels = xLabels) +
    scale_y_continuous(breaks = yBreaks, labels = yLabels) +
    # set limits and breaks so they are consistent across plots
    scale_size(limits = c(0, 350), breaks = sizeBreaks, labels = sizeBreaks) +
    scale_fill_manual(limits = c(0, 1),
                      values = c("white", "black")) +
    scale_alpha_continuous(range = c(0.5, 0.8), limits = c(0, 350),
                           breaks = sizeBreaks, labels = sizeBreaks) +
    facet_grid(weekOfYear~year)
  
  # ausPolygonSmallFortified %>%
  #   ggplot(aes(x = long, y = lat, group = group)) +
  #   theme_bw() +
  #   theme(panel.grid.major = element_blank(),
  #         strip.background = element_blank(),
  #         plot.title = element_text(hjust = 0.5, face = "bold")) +
  #   ggtitle(w) +
  #   labs(x = "Longitude", y = "Latitude") +
  #   geom_polygon(fill = "lightgrey", colour = "white") +
  #   geom_point(data = tmpDf %>% mutate(group = 1) %>% filter(nMoths > 0), # needed to match the aes() call in the polygons
  #              aes(x = longitude, y = latitude, size = nMoths, alpha = nMoths),
  #              pch = 21, fill = myPal, colour = myPal) +
  #   # map the zero values as open circles
  #   geom_point(data = tmpDf %>% mutate(group = 1) %>% filter(nMoths == 0),
  #              aes(x = longitude, y = latitude),
  #              pch = 21, fill = "lightgrey", colour = myPal) +
  #   coord_map(ylim = c(-38, -32), xlim = c(133.5, 141.5)) +
  #   scale_x_continuous(breaks = xBreaks, labels = xLabels) +
  #   scale_y_continuous(breaks = yBreaks, labels = yLabels) +
  #   # set limits and breaks so they are consistent across plots
  #   scale_size(limits = c(0, 350), breaks = sizeBreaks, labels = sizeBreaks) +
  #   scale_alpha_continuous(range = c(0.5, 0.8), limits = c(0, 350),
  #                          breaks = sizeBreaks, labels = sizeBreaks) +
  #   facet_grid(weekOfYear~year)
  
  
  counter <- counter + 1
  plotNum  <- str_pad(counter, side = "left", width = 2, pad = 0)
  plotName <- paste0("plotMovie", plotNum, ".png")
  ggsave(plotName, width = 10, height = 6, dpi = 300) 
}

# now make a time series movie from the plot frames
movieMEDIUM <- 'ffmpeg -y -i plotMovie%02d.png -c:v libx264 -crf 10 -b:v 1M -bufsize 1M -vf "setpts=(15/1)*PTS" plotMovieMEDIUM.mp4'
system(movieMEDIUM)


# TO DO:
# map the NAs to an `X`
# map the zero moths as open circles
# check the mapping. Weekly data are not correct yet.

# NEXT:
# We might have
# Impute data across missing weeks.
# Add the trapping start data to the imputations
# For more accuracy, we might have to do the imputing/averagig on daily counts first (we'll see)


# End script
# ############

# OLD CODE BELOW
# # ####################################################################
# # Calculate a cumulative rainfall graph for colonisation field sites #
# # ####################################################################
# 
# # as a quick test, use rainfall data for the grid cells used across all 3 years
# # to make a cumulative rainfall graph.
# 
# x <- dailyObsList %>% 
#   bind_rows() %>% # simply ignore field site (grid cell) and average the rainfall data for each date
#   dplyr::select(date, dRai) %>%
#   mutate(date = as_date(date),
#          year = lubridate::year(date),
#          year = factor(year),
#          dayOfYear = lubridate::yday(date)) %>%
#   group_by(year, date, dayOfYear) %>%
#   summarise(meanDailyRain = mean(dRai, na.rm = TRUE)) %>%
#   ungroup() %>%
#   split(.$year) %>%
#   map(mutate, cumRain = cumsum(meanDailyRain)) %>%
#   bind_rows()
# 
# 
# x %>%
#   filter(#dayOfYear > 31,
#          dayOfYear < 120,
#          year %in% c("2014", "2015", "2016")) %>% # just keep approx the frst four months
#   ggplot(aes(x = dayOfYear, y = cumRain, colour = year)) +
#   geom_line(size = 1) +
#   theme_bw()

# ######################################
# Calculate metrics for the manuscript #
# ######################################

# ===================================================
# calculate some data metrics for the manuscript text
# ===================================================

# tabulate the median dates of sowing and colonisation
# We use weekly data, but fort each site, the actual colonisation date is stored in column `eggLayMedian`
colonisationWeekly %>%
  group_by(siteCode, year) %>%
  summarise(nWeeks = length(siteCode),
            dateSowing = unique(dateSowing),
            dateColon  = unique(eggLayMedian)) %>% 
  ungroup() %>% 
  mutate(colonDaysAfterSowing = yday(dateColon) - yday(dateSowing)) %>%
  group_by(year) %>%
  summarise(nSites = length(siteCode),
            mnDateColon  = mean(dateColon, na.rm = TRUE),
            sdDateColon  = sd(yday(dateColon), na.rm = TRUE),
            seDateColon  = se(yday(dateColon), na.rm = TRUE),
            mnDateSowing = mean(dateSowing, na.rm = TRUE),
            seDateSowing = se(yday(dateSowing), na.rm = TRUE),
            mnColonDaysAfterSowing = mean(colonDaysAfterSowing, na.rm = TRUE),
            seColonDaysAfterSowing = se(colonDaysAfterSowing, na.rm = TRUE),
            sdColonDaysAfterSowing = sd(colonDaysAfterSowing, na.rm = TRUE)) %>% View()

# Calculate mean moth counts for May and June across sites for each year
colonisationWeekly %>%
  mutate(month = lubridate::month(week)) %>%
  filter(month %in% 5:6, # May to July
         !is.na(aveMoths)) %>%
  group_by(year) %>%
  summarise(nObs = length(year),
            meanMothsPerSite = mean(aveMoths, na.rm = TRUE),
            sdMothsPerSite   = sd(aveMoths, na.rm = TRUE)) %>% View()

# Calculate proportion of crops colonised in May and June
colonisationWeekly %>% 
  group_by(siteCode) %>%
  summarise(date = unique(eggLayMedian)) %>%
  ungroup() %>%
  mutate(month = month(date),
         year  = year(date)) %>%
  group_by(month) %>%
  summarise(nSitesColonised = length(month)) %>%
  ungroup() %>%
  mutate(pSitesColonised = 100 * nSitesColonised / sum(nSitesColonised))


# plot the proportion of crops colonised over time by week as a LINE graph
colonisationWeekly %>%
  mutate(dayOfYear   = yday(week),
         month = month(week), 
         year = factor(year)) %>%
  group_by(year, dayOfYear, weekOfYear, month) %>% # keep `dayOfYear` just for sorting in date order
  summarise(nSites = length(unique(siteCode)),
            nSitesColonised = sumValues(juvPresence), 
            pSitesColonised = nSitesColonised / nSites) %>%
  ungroup() %>%
  arrange(dayOfYear, year) %>% 
  ggplot(aes(x = dayOfYear, y = pSitesColonised, colour = year)) +
  geom_point() +
  geom_line()


# Summarise the proportion of crops colonised by the end of each month
colonisationProportion %>%
  group_by(year, month) %>%
  summarise(nSites = unique(nSites),
            nSitesColonised = maxValues(nSitesColonised),
            pSitesColonised = maxValues(pSitesColonised)) %>% 
  ungroup() %>%
  filter(month %in% 5:9) # just keep May to September, as no additional info (no more crops colonised) beyond Sept. 


# porportion crops colonised overall
colonisationProportion %>% 
  group_by(year) %>%
  summarise(nSites = unique(nSites),
            nSitesColonised = maxValues(nSitesColonised),
            pSitesColonised = nSitesColonised / nSites)


# -----------------------------------------------------------------
# Compare colonisation timing in the South East with other regions 
# -----------------------------------------------------------------

# First, for all years combined, split sites into two groups: 
# (i) above (Eyre Peninsula / Mid North) and (i) below (South East, Kanagoor Island) -35o latitude

colonisationWeekly %>%
  # create a factor for above and below -35o N
  mutate(regionFactor = ifelse(latitude < -35, "South East", "Eyre and North")) %>%
  # first combine the weeks, just keep colonisation date
  group_by(siteCode, regionFactor, year) %>%
  summarise(nWeeks = length(siteCode),
            dateSowing = unique(dateSowing),
            dateColon  = unique(eggLayMedian)) %>% 
  ungroup() %>% 
  mutate(colonDaysAfterSowing = yday(dateColon) - yday(dateSowing)) %>%
  group_by(regionFactor, year) %>%
  summarise(nSites = length(siteCode),
            mnDateColon  = mean(dateColon, na.rm = TRUE),
            sdDateColon  = sd(yday(dateColon), na.rm = TRUE),
            seDateColon  = se(yday(dateColon), na.rm = TRUE),
            mnDateSowing = mean(dateSowing, na.rm = TRUE),
            seDateSowing = se(yday(dateSowing), na.rm = TRUE),
            mnColonDaysAfterSowing = mean(colonDaysAfterSowing, na.rm = TRUE),
            seColonDaysAfterSowing = se(colonDaysAfterSowing, na.rm = TRUE))
#> Colonisation timing wasn't different between South East and all other sites (Eyre and North), other than by 2 weeks in 2014

# Were the following metrics `later` in the South East: sowing dates, colonisation day of the year, colonisation days after sowing check whether days after sowing different
colonisationWeekly %>%
  filter(year == 2015) %>% # to check the pattern just for 2015, when colonisation was later (comment out to look at all years combined)
  # create a factor for above and below -35o N
  mutate(regionFactor = ifelse(latitude < -35, "South East", "Eyre and North")) %>%
  # first combine the weeks, just keep colonisation date
  group_by(siteCode, regionFactor, year) %>%
  summarise(nWeeks = length(siteCode),
            dateSowing = unique(dateSowing),
            dateColon  = unique(eggLayMedian)) %>% 
  ungroup() %>% 
  mutate(dayOfYearSowing = yday(dateSowing),
         dayOfYearColonisation = yday(dateColon),
         colonDaysAfterSowing = dayOfYearColonisation - dayOfYearSowing) %>%
  group_by(regionFactor) %>%
  summarise(nSites = length(siteCode),
            mnDayOfYearColon  = mean(dayOfYearColonisation, na.rm = TRUE),
            seDayOfYearColon  =   se(dayOfYearColonisation, na.rm = TRUE),
            mnDayOfYearSowing   = mean(dayOfYearSowing, na.rm = TRUE),
            seDayOfYearSowing   =   se(dayOfYearSowing, na.rm = TRUE),
            mnColonDaysAfterSowing = mean(colonDaysAfterSowing, na.rm = TRUE),
            seColonDaysAfterSowing =   se(colonDaysAfterSowing, na.rm = TRUE))
# This shows that sowing is 7 days later in the south east, but colonisation is only 5 days (no diff between regions)


# ==========================================
# Calculate some metrics for the manuscript
# ==========================================

# The geographic bounding box of the colonisation sites
sitesData %>%
  mutate(year = yearFromSiteCode(siteCode)) %>%
  group_by(year) %>%
  summarise(minLon = min(longitude), maxLon = max(longitude),
            minLat = min(latitude),  maxLat = max(latitude))

# distance between the further points
# 2014
geoDist2014 <- geosphere::distm(
  sitesData %>%
    mutate(year = year(dateTrappingStart)) %>%
    filter(year == 2014) %>% 
    dplyr::select(longitude, latitude) %>% as.matrix(ncol = 2))
geoDist2014 <- geoDist2014 / 1000
# exclude the zeros for your calculations
max(geoDist2014[geoDist2014 != 0])  #> 778.5928
mean(geoDist2014[geoDist2014 != 0]) #> 310.0351
sd(geoDist2014[geoDist2014 != 0])   #> 185.788
min(geoDist2014[geoDist2014 != 0])  #> 6.295426

# 2015
geoDist2015 <- geosphere::distm(
  sitesData %>%
    mutate(year = year(dateTrappingStart)) %>%
    filter(year == 2015) %>% 
    dplyr::select(longitude, latitude) %>% as.matrix(ncol = 2))
geoDist2015 <- geoDist2015 / 1000
# exclude the zeros for your calculations
max(geoDist2015[geoDist2015 != 0])  #> 815.2895
mean(geoDist2015[geoDist2015 != 0]) #> 281.7316
sd(geoDist2015[geoDist2015 != 0])   #> 164.5533
min(geoDist2015[geoDist2015 != 0])  #> 0.610207

# 2016
geoDist2016 <- geosphere::distm(
  sitesData %>%
    mutate(year = year(dateTrappingStart)) %>%
    filter(year == 2016) %>% 
    dplyr::select(longitude, latitude) %>% as.matrix(ncol = 2))
geoDist2016 <- geoDist2016 / 1000
# exclude the zeros for your calculations
max(geoDist2016[geoDist2016 != 0])  #> 796.6897
mean(geoDist2016[geoDist2016 != 0]) #> 302.9074
sd(geoDist2016[geoDist2016 != 0])   #> 165.4281
min(geoDist2016[geoDist2016 != 0])  #> 8.216268


#  sitesData metrics
sitesData %>%
  mutate(year = yearFromSiteCode(siteCode)) %>%
  group_by(year) %>%
  summarise(nSites = length(siteCode),
            nSitesJuvsDetected = length(which(larvaeDetected == "y")),
            pSitesJuvsDetected = round(nSitesJuvsDetected / nSites, 2),
            nSitesJuvSamplesProvided = length(which(larvalSamples == "y"))) # could also look at nrow() of headCapsuleData for each year

# headCapsuleData metrics
(headCapsuleMetrics <- headCapsuleDataLifeStages %>%
    filter(firstDetectionFlag == 1) %>% # turn off filter to jkeep all samples (filter keeps samples from first cohort only)
    group_by(instar) %>%
    count() %>%
    ungroup() %>%
    mutate(p = round(n / sum(n), 2)))
sum(headCapsuleMetrics$n) # total n samples


# ##################################################################################
# Determine generations possible per canola crop cycle, based on colonisation date #
# ##################################################################################


# Earlier colonisation could allow more generations to develop
# Calculate using mean hourly temperatures across the 66 grid cells
# and local temps and median colonisation date in 2014, 2015, 2016
# Report these metrics in the manuscript

# Assume germination date of 1-May, harvest date 30th October (6 months)
# Note Ridland and Endersby 2008 used 1st April to 30th Sept for their calculations based on degree days 

# First average hourly data across all grid cells (66)
hourlyObsSiteMeans <- hourlyObsList %>%
  bind_rows %>%
  group_by(datetime) %>%
  summarise(nCells = length(station),
            obs    = mean(obs, na.rm = TRUE)) %>%
  ungroup()

# Now run some predictions for each year using the average colonisation dates across sites

# 2014
fwdDev(tempObsDf   = hourlyObsSiteMeans, 
       devParamsDf = devParamsTotal, 
       startDate   = "2014-05-27", # the average predicted colonisation date across all sites
       startStage = "lifecycle", startDev = 0, 
       gens = 4, output = "increments") %>%
  bind_rows() %>% 
  mutate(cumDev = cumsum(dev)) %>% # For `total` lifecycle, calculating cumulative deveopment should come standard in the output 
  filter(date(datetime) == as_date("2014-10-31"), hour(datetime) == 12)
#> cumdev = 3.27

# 2015
fwdDev(tempObsDf   = hourlyObsSiteMeans, 
       devParamsDf = devParamsTotal, 
       startDate   = "2015-07-10", # the average predicted colonisation date across all sites
       startStage = "lifecycle", startDev = 0, 
       gens = 4, output = "increments") %>%
  bind_rows() %>% 
  mutate(cumDev = cumsum(dev)) %>% # For `total` lifecycle, calculating cumulative deveopment should come standard in the output 
  filter(date(datetime) == as_date("2015-10-31"), hour(datetime) == 12)
#> cumdev = 2.45

# 2016
fwdDev(tempObsDf   = hourlyObsSiteMeans, 
       devParamsDf = devParamsTotal, 
       startDate   = "2016-06-11", # the average predicted colonisation date across all sites
       startStage = "lifecycle", startDev = 0, 
       gens = 4, output = "increments") %>%
  bind_rows() %>% 
  mutate(cumDev = cumsum(dev)) %>% # For `total` lifecycle, calculating cumulative deveopment should come standard in the output 
  filter(date(datetime) == as_date("2016-10-31"), hour(datetime) == 12)
#> cumdev = 2.46



# End script
# ###################################################################################




