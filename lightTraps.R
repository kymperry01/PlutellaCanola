# ##########################################################################
# Analysis for light trap data from four locations
# K. Perry, 30/4/2017
# 
# This script does:
# 1. read in daily light trap Plutella moth time series count data for four locations
# 2. impute mean values to generate average nightly moth counts
# 3. Plot counts at daily, weekly, monthly intervals
# 4. Compare trap counts against CLIMEX weekly Growth Index on chronological time axis
# 5. Compare trap counts against Cx wGI on a phenological axis

# Notes:
# The analysis uses mean nightly catch, not raw
#
# To do:
# ... for each trap, calculate the date coverage of data ...
# ... change labels of facet plots
# ... plot a horizontal line on each trap to show date coverage (minnipa, two wells missing data)
# ... try over-plotting the sex ratio (transformed to half max count) for male positive, female negative plots
# ... add to job notes list for DBMdevmod2::mlocDev(): extra argument for each site, `gens`, to be able to model a different number of generations for each site

# ####################################################

library(tidyverse)
library(readxl)
library(reshape2)
library(lubridate)
library(ggplot2)
library(RColorBrewer)
library(DBMdevmod2)
library(raster)

# ##################
# Define functions #
# ##################

# ... a function to impute means over accumulated time series values
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
          # in a vector without an integer (trailing NAs) 
          # ... keep NAs to maintain length(x)
          rep(NA, length(x))
        }
      }) %>%
    unlist() %>%
    as.numeric()
} 

# ##########################
# Organise moth count data #
# ##########################

trapDataPath <- file.path("C:", "UserData", "Kym", "PhD", "Data", "lightTraps", "lightTrapsMaster.xlsx")
trapSheetNms <- c("Minnipa", "TwoWells", "Sherwood", "DrummondPoint")

# dates when the Minnipa trap did not operate
noTrapDates <- seq(as.Date("2014-11-10"), as.Date("2015-03-11"), by = "day")

mothData <- trapSheetNms %>% 
  map(read_excel, path = trapDataPath) %>% 
  bind_rows() %>%
  mutate(date = as_date(date),
         mnMale   = imputeMeanValues(male),
         mnFemale = imputeMeanValues(female),
         totalMoths = sum(mnMale, mnFemale),
         # impute NA values for the dates the Minnipa trap did not operate 
         mnMale   = ifelse(siteCode == "mac" & date %in% noTrapDates, NA, mnMale),
         mnFemale = ifelse(siteCode == "mac" & date %in% noTrapDates, NA, mnFemale),
         mnMale   = mnMale * -1) %>% # transform to negative for plotting below axis
  dplyr::select(siteCode, date, male = mnMale, female = mnFemale) %>%
  # interpolate the climex weekly to a daily interpolation 
  melt(id.vars       = c("siteCode", "date"),
       variable.name = "sex",
       value.name    = "count")


# #############################################################################
# Organise CLIMEX Weekly growth index data predicted for the 4 trap locations # 
# #############################################################################

# We want to plot moth trap data against predicted climate suitabilty for P. xylostella,
# to infer whether moths are likely to have been locally derived.
# Use time series climate data extracted from interpolated SILO gridded datasets.
# These were extracted using a custom script `dataDrill.R` by K. Perry

# Use a reference raster to determine the grid cell number for each light trap location.
refRaster <- file.path(
  "C:", "UserData", "Kym", "PhD", "Data", "climateData", "siloDailyGrid2016", "dailyMax_2013-2016", "20131201_max.txt"
) %>% raster()

# light trap lats and longs
trapMetaData <- tribble(
  ~siteCode, ~lat,      ~lon,
  "sh",  -36.01998, 140.5766,
  "tw",  -34.5859,  138.5127,
  "dp",  -34.1369,  135.3312,
  "mac", -32.8332,  135.1512
) %>%
  mutate(siloGridCell = cellFromXY(
    refRaster, xy = as.matrix(dplyr::select(., lon, lat), ncol = 2)
  ))

# read in the output of CLIMEX `compare years` model for the four grid cells.
gridCellFromFileName <- function(x){
  x %>% str_extract("\\d{1,10}\\.") %>% as.numeric()
}

climexFilePaths <- file.path(
  "C:", "UserData", "Kym", "PhD", "Data", "climexOutput", "climexCompareYearsPointLocations"
) %>%
  list.files(full.names = TRUE) %>%
  magrittr::extract(which(gridCellFromFileName(.) %in% trapMetaData$siloGridCell))

climexPredictions <- climexFilePaths %>%
  map(function(x){
    x %>% 
      read_csv(skip = 3, col_names = TRUE) %>%
      mutate(siloGridCell = gridCellFromFileName(x)) # Needed this extra step in the read process. to grab the cell number from filename.
  }) %>%
  bind_rows() %>%
  mutate(date = as_date(SimulaDate, format = "%d/%m/%Y", tz = "Australia/Adelaide")) %>%
  left_join(trapMetaData %>% dplyr::select(siteCode, siloGridCell), by = "siloGridCell") %>% # add the site code
  dplyr::select(siteCode, siloGridCell, date, weeklyGI = `GI(w)`)


# ############################################
# Plot moth data on a phenological timescale #
# ############################################

# Plot moth counts against P. xylostella generations predicted using a temperature-based model.

# =================================================
# Extract daily max and min temperature time series
# =================================================

# These time series were extracted from SILO daily gridded data for Australia using
# a custom R script `dataDrill.R` by K Perry.

dailyObsPaths <- file.path(
  "C:", "UserData", "Kym", "PhD", "Data", "climateData", "dailyDataDrill2016"
  ) %>%
  list.files(pattern = "^dd_", full.names = TRUE) %>%
  magrittr::extract(which(gridCellFromFileName(.) %in% trapMetaData$siloGridCell))

# dailyObsList <- dailyObsPaths %>%
#   map(function(x){
#     x %>%
#       read.table(skip = 8, header = TRUE) %>%
#       mutate(siloGridCell = gridCellFromFileName(x)) # needed this extra step in this first loop to add the grid cell
#     }) %>%
#   set_names(gridCellFromFileName(dailyObsPaths)) # modelling packages needs named list to match climate data to site.
#   
# # interpolate to hourly temperatures using custom package DBMdevmod2 by K. Perry
# hourlyObsList <- dailyObsList %>%
#   map(dplyr::select, siloGridCell, date, min = dMin, max = dMax) %>% # presently must rename vars to run models using DBMdevmod2. Fix this in the package!
#   # add the lats and long to calculate sunrise and sunset times
#   map(left_join, trapMetaData %>% dplyr::select(siloGridCell, latitude = lat, longitude = lon)) %>%
#   map(msunriset) %>% # need sunrise/set times to interpolate hourly temperatures
#   map(dplyr::rename, station = siloGridCell) %>% # presently must rename to station. Fix!! %>% 
#   map(hourlyTemperatures)


# ==============================================================================================================
# Run models to predict P. xylostella generation times for the four light trap locations across the study period
# ==============================================================================================================

# Set P. xylostella parameters for the entire lifecycle
# These params were determined by Perry and Keller, 2018 using R package `nlsreg`
dbmParams <- data.frame(
  lifecycle = c(9.38E-05, 1.042, 32.73, 9.943)
  ) %>%
  magrittr::set_rownames(c("a", "Tmin", "Tmax", "m")) %>%
  t()

# run models for all sites across the date range of interest
# first, set the number of generations to run for each site. We set these 
# by trial and error to cover the full tie range of interest. Currently, the package 
# does not calculate cumulative development for a given time period - a development task. 
 
gensPerSite <- tribble(
  ~siteCode, ~ gens,
  "sh", 32,
  "tw", 36,
  "dp", 33,
  "mac", 37
)

# run models to predict the generation time using the same biofix for each site
dbmDevPredictionsIncrements <- trapMetaData %>%
  mutate(startDate  = "2014-03-01", # start from autumn, 3 months before trapping
         startStage = "lifecycle",  # this value must match one of rownames(dbmParams)
         startDev   = 0) %>% 
  left_join(gensPerSite, by = "siteCode") %>%
  dplyr::rename( # Currently, mlocdev() takes set column names in the input DF. FIX
    site = siteCode, latitude = lat, longitude = lon, station = siloGridCell
    ) %>% 
  DBMdevmod2::mlocDev(
    tempObsDfList = hourlyObsList,
    devParamsDf  = dbmParams,
    timedir = "fwd",
    output  = "increments"
    )
# > in package, fix to suppress warning messages about binding unequal factors. Messy. 

# run models and output the `generations` predictions
dbmDevPredictionsGens <- trapMetaData %>%
  mutate(startDate  = "2014-03-01", # start from autumn, 3 months before trapping
         startStage = "lifecycle", # this value must match one of rownames(dbmParams)
         startDev   = 0) %>% 
  left_join(gensPerSite, by = "siteCode") %>%
  dplyr::rename( # Currently, mlocdev() takes set column names in the input DF. FIX
    site = siteCode, latitude = lat, longitude = lon, station = siloGridCell
    ) %>% 
  DBMdevmod2::mlocDev(
    tempObsDfList = hourlyObsList,
    devParamsDf   = dbmParams,
    timedir       = "fwd",
    output        = "generations"
    )

# ===============================================================================================
# Prepare data frames for plotting moth counts and climex predictions on a phenological timescale
# ===============================================================================================

# To do this, we need to plot the nightly catches against the cumulative 
# number of generations instead of the date. Hence, calculate the
# predicted cumulative development (in generations as units) for each date 
# of the entire date range. We'll use the cumulative development at 12pm midday
# as the `generations` value for each date, then join these values to the 
# main data to use as the time axis.

# extract the cumulative development at 12pm daily
dbmDevPredictionsDaily12pm <- dbmDevPredictionsIncrements %>%
  map(mutate, cumulativeGens = cumsum(dev)) %>% # must calculate cumul. total before filtering rows!
  map(filter, lubridate::hour(datetime) == 12) %>%
  map(mutate, date = lubridate::date(datetime)) %>% # variable for joining
  bind_rows() %>%
  dplyr::select(siteCode = site, siloGridCell = station, date, cumulativeGens)

# set the breaks and labels for adding chronological dates to the x-axis separately
breaksDates <- seq(as_date("2014-03-01"), as_date("2016-12-01"), by = "1 month")
breaksDBMGens <- dbmDevPredictionsDaily12pm %>%
  filter(date %in% breaksDates) %>%
  dplyr::select(siteCode, date, cumulativeGens)

# rescale the climex weekly growth index predictions [0, 1] against max moth counts for each site
# for plotting on the same y-axis, and add the DBM development for plotting on the timescale axis
rescaleToMaxCount <- function(x, y){
  x * max(y, na.rm = TRUE)
}

maxCounts <- mothData %>%
  group_by(siteCode, date) %>% 
  summarise(count = sum(abs(count))) %>%
  group_by(siteCode) %>%
  summarise(maxCount = max(count, na.rm = TRUE))

climexPredictionsDBMDevPredictions <- climexPredictions %>%
  split(.$siteCode) %>%
  map(left_join, maxCounts, by = "siteCode") %>%
  map(mutate, weeklyGIRescaled = rescaleToMaxCount(weeklyGI, y = maxCount)) %>%
  bind_rows() %>%
  # add the phenological timescale for the x-axis
  left_join(dbmDevPredictionsDaily12pm, by = c("siteCode", "siloGridCell", "date"))

# create a dataframe with the canola cropping windows for each site 
# for plotting against the timescale x-axis using main ggplot aes()
canolaWindows <- tribble(
  ~siteCode, ~datePlant,   ~dateHarvest,
  "dp"     , "2014-05-01", "2014-11-01",
  "dp"     , "2015-05-01", "2015-11-01",
  "dp"     , "2016-05-01", "2016-11-01",
  "mac"    , "2014-05-01", "2014-11-01",
  "mac"    , "2015-05-01", "2015-11-01",
  "mac"    , "2016-05-01", "2016-11-01",
  "tw"     , "2014-05-01", "2014-11-01",
  "tw"     , "2015-05-01", "2015-11-01",
  "tw"     , "2016-05-01", "2016-11-01",
  "sh"     , "2014-05-01", "2014-11-01",
  "sh"     , "2015-05-01", "2015-11-01",
  "sh"     , "2016-05-01", "2016-11-01" 
  ) %>%
  reshape2::melt(id.vars       = "siteCode",
       variable.name = "cropDates",
       value.name    = "date") %>%
  mutate(date = as_date(date),
         year = year(date)) %>%
  left_join(dbmDevPredictionsDaily12pm, by = c("siteCode", "date")) %>% # now add the DBM development predictions
  reshape2::dcast(siteCode + year ~ cropDates, value.var = "cumulativeGens") %>%
  rowwise() %>%
  mutate(midPoint = mean(c(datePlant, dateHarvest))) %>%
  data.frame() # needed this to override rowwise

# ==================================================================================================
# Now analyse peaks in flight activity against peaks in the prior generation for individual sites 
# ==================================================================================================

# Conduct individual analysis of each site to identify significant flights, and determine
# whether they coincide with prior peaks.

myArrow <- grid::arrow(angle = 20, length = unit(2, "mm"), ends = "both", type = "closed")
siteNames <- c("sh"  = "Sherwood, South East",
               "dp"  = "Drummond Point, Lower Eyre Peninsula",
               "mac" = "Minnipa, Upper Eyre Peninsula",
               "tw"  = "Two Wells, Northern Adelaide Plains")

# *********************
# Drummond point (dp) #
# *********************

# determine significant flight dates for drummond point
mothFlightsDrummondPoint <- mothData %>% 
  # first sum the male and female counts
  group_by(siteCode, date) %>% 
  summarise(count = sum(abs(count))) %>%
  filter(siteCode == "dp", count > 50) %>%
  # add the cumulative DBM generations
  left_join(dbmDevPredictionsDaily12pm, by = c("siteCode", "date")) %>%
  mutate(retrace1Gen = cumulativeGens - 1)

# manually identify some dates of `smaller` flights of interest (after visual examination)
# and create a df for plotting on the main ggplot() with the same aes()
smallerFlightDatesDrummondPoint <- c(
  "2015-12-20", "2016-01-10", "2016-02-03", "2016-02-24", "2016-03-14"
  ) %>% as_date()

mothFlightsDrummondPointSmaller <- mothData %>% 
  # first sum the male and female counts
  group_by(siteCode, date) %>%
  summarise(count = sum(abs(count))) %>%
  filter(siteCode == "dp", date %in% smallerFlightDatesDrummondPoint) %>%
  # add the cumulative DBM generations
  left_join(dbmDevPredictionsDaily12pm, by = c("siteCode", "date")) %>%
  mutate(retrace1Gen = cumulativeGens - 1)

# Now select the flight dates of interest for re-tracing and plotting in the main figure
# These were selected after visual examination and retracing of many peaks
selectedFlightDatesDrummondPoint <- c(
  "2015-05-08", "2015-05-26", "2015-09-12", "2015-10-06", "2015-10-25", 
  "2015-12-20", "2016-01-10", "2016-02-03", "2016-02-24", "2016-03-14", 
  "2016-04-11", "2016-09-24", "2016-10-02", "2016-10-15"
  ) %>% as_date()
selectedFlightsDrummondPoint <- mothData %>% # add to the main data for plotting with correct aes() 
  group_by(siteCode, date) %>%  # first sum the male and female counts
  summarise(count = sum(abs(count))) %>%  
  filter(siteCode == "dp", date %in% selectedFlightDatesDrummondPoint) %>%
  left_join(dbmDevPredictionsDaily12pm, by = c("siteCode", "date")) %>% # add the cumulative DBM generations
  mutate(retrace1Gen = cumulativeGens - 1)

# plot the significant flights and retrace 1 generation (Drummond Point)
(plotDrummond <- mothData %>%
  filter(siteCode == "dp") %>%
  group_by(siteCode, date) %>%
  summarise(count = sum(abs(count))) %>%
  left_join(dbmDevPredictionsDaily12pm, by = c("siteCode", "date")) %>%
  mutate(siteName = str_replace_all(siteCode, siteNames)) %>%
  ggplot(aes(x = cumulativeGens, y = count)) + # this just adds the siteName nicely above the plot
  facet_wrap(~siteName, ncol = 1, scales = "free") +
  theme_light() +
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black"), # theme_light has white text which is invisible.
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey"),
        panel.grid.minor.x = element_blank()) +
  labs(x = "Cumulative development (generations)", y = "Mean no. moths trapped per night") +
  # Add chronological dates as labels
  geom_text(data = breaksDBMGens %>% filter(siteCode == "dp"),
            aes(x = cumulativeGens , y = -90, label = format(date, "-%b-%Y")), # extra dash as tick mark
            angle = 90, size = 3, hjust = "left",
            vjust = "center") + # when angle = 90, vjust = "center" aligns in the correct place on the timescale x-axis
  scale_x_continuous(breaks = 1:35, labels = 1:35, expand = c(0, 0)) +
  scale_y_continuous(limits = c(-90, 425), expand = c(0, 0)) + # to align date labels on x-axis, set min y limit same y in aes()
  coord_cartesian(xlim = c(3.75, 29.5)) +
    # ------------------------------
  # Add CLIMEX Weekly Growth Index
  # ------------------------------
  geom_line(data = climexPredictionsDBMDevPredictions %>%
              filter(siteCode == "dp"),
            aes(x = cumulativeGens, y = weeklyGIRescaled), 
            linetype = 3, colour = "black") +
  # --------------------------------------------------------------
  # Add segments and text labels to highlight the cropping windows
  # --------------------------------------------------------------
  geom_segment(data = canolaWindows %>% filter(siteCode == "dp"), 
               aes(x = datePlant, xend = dateHarvest, y = -10, yend = -10),
               colour = brewer.pal(9, "Greens")[5], size = 3) +
  geom_label(data = canolaWindows %>% filter(siteCode == "dp"), 
             aes(x = midPoint, y = -10, label = "Canola")) +
  # # -------------------------------------------------------------------
  # # Draw segments retracing significant flights by one generation (red)
  # # -------------------------------------------------------------------
  # # (note the flight in May 2015! no prior pop.)
  # geom_segment(data = mothFlightsDrummondPoint,
  #              (aes(x = retrace1Gen, xend = retrace1Gen,
  #                   y = 0, yend = count)), linetype = 4, colour = "red") +
  # geom_segment(data = mothFlightsDrummondPoint,
  #              (aes(x = retrace1Gen, xend = cumulativeGens,
  #                   y = count, yend = count)), linetype = 1, colour = "red",
  #              arrow = myArrow) +
  # # ----------------------------------------------------------
  # # draw segments retracing smaller flights of interest (blue)
  # # ----------------------------------------------------------
  # # note how there appear to be a local cohort going through discrete generations during summer 2016!
  # geom_segment(data = mothFlightsDrummondPointSmaller,
  #              (aes(x = retrace1Gen, xend = retrace1Gen, y = 0, yend = count)), 
  #              linetype = 4, colour = "blue") +
  # geom_segment(data = mothFlightsDrummondPointSmaller,
  #              (aes(x = retrace1Gen, xend = cumulativeGens, y = count, yend = count)), 
  #              linetype = 1, colour = "blue", arrow = myArrow) +
  # ------------------------------------------------------------------
  # draw segments retracing the final SELECTED flights of interest (red): replaces the exploratory retracing above.
  # --------------------------------------------------------------------
  geom_segment(data = selectedFlightsDrummondPoint,
               (aes(x = retrace1Gen, xend = retrace1Gen, y = 0, yend = count)),
               linetype = 4, colour = brewer.pal(9, "Blues")[6]) +
  geom_segment(data = selectedFlightsDrummondPoint,
               (aes(x = retrace1Gen, xend = cumulativeGens, y = count, yend = count)), 
               linetype = 1, colour = brewer.pal(9, "Blues")[6], arrow = myArrow) +
  geom_line(aes(group = siteCode))) # overplot the actual moth trap data last
ggsave("lightTrapsFlightsDrummondPointDraftV2.pdf", width = 20, height = 8)

# make a zoomed plot of Sherwood for publication
plotDrummond +
  #theme_bw(base_size) +
  coord_cartesian(x = 8:28) +
  theme_classic(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "white"))
ggsave("zoomPlotDrummondDraft.pdf", width = 11, height = 7)

# ***********************************
# Minnipa Agricultural Centre (mac) #
# ***********************************

# determine significant flight dates for Minnipa
# mothFlightsMinnipa <- mothData %>% 
#   # first sum the male and female counts
#   group_by(siteCode, date) %>% 
#   summarise(count = sum(abs(count))) %>%
#   filter(siteCode == "mac", count > 50) %>%
#   # add the cumulative DBM generations
#   left_join(dbmDevPredictionsDaily12pm, by = c("siteCode", "date")) %>%
#   mutate(retrace1Gen = cumulativeGens - 1)
# 
# # manually identify some dates of `smaller` flights of interest (after visual 
# # examination of plots) and create a df with same aes() for adding to the main ggplot
# smallerFlightDatesMinnipa <- c(
#   "2015-12-20", "2016-01-10", "2016-02-03", "2016-02-24", "2016-03-14"
# ) %>% as_date()
# 
# mothFlightsMinnipaSmaller <- mothData %>% 
#   # first sum the male and female counts
#   group_by(siteCode, date) %>%
#   summarise(count = sum(abs(count))) %>%
#   filter(siteCode == "mac", date %in% smallerFlightDatesMinnipa) %>%
#   # add the cumulative DBM generations
#   left_join(dbmDevPredictionsDaily12pm, by = c("siteCode", "date")) %>%
#   mutate(retrace1Gen = cumulativeGens - 1)

# Now select the flight dates of interest for re-tracing and plotting in the main figure
# These were selected after visual examination and retracing of many peaks
selectedFlightDatesMinnipa <- c(
  "2014-09-20", "2014-09-23", "2015-08-19", "2015-08-31", "2015-09-29", "2016-08-16", "2016-09-02"
  ) %>% as_date()

selectedFlightsMinnipa <- mothData %>% # add to the main data for plotting with correct aes() 
  group_by(siteCode, date) %>%  # first sum the male and female counts
  summarise(count = sum(abs(count))) %>%  
  filter(siteCode == "mac", date %in% selectedFlightDatesMinnipa) %>%
  left_join(dbmDevPredictionsDaily12pm, by = c("siteCode", "date")) %>% # add the cumulative DBM generations
  mutate(retrace1Gen = cumulativeGens - 1)


# plot the significant flights and retrace 1 generation (Minnipa)
mothData %>%
  filter(siteCode == "mac") %>%
  group_by(siteCode, date) %>%
  summarise(count = sum(abs(count))) %>%
  left_join(dbmDevPredictionsDaily12pm, by = c("siteCode", "date")) %>%
  mutate(siteName = str_replace_all(siteCode, siteNames)) %>%
  ggplot(aes(x = cumulativeGens, y = count)) + # this just adds the siteName nicely above the plot
  facet_wrap(~siteName, ncol = 1, scales = "free") +
  theme_light() +
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black"), # theme_light has white text which is invisible.
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey"),
        panel.grid.minor.x = element_blank()) +
  labs(x = "Cumulative development (generations)", y = "Mean no. moths trapped per night") +
  # Add chronological dates as labels
  geom_text(data = breaksDBMGens %>% filter(siteCode == "mac"),
            aes(x = cumulativeGens , y = -40, label = format(date, "-%d-%b-%Y")), # extra dash as tick mark
            angle = 90, size = 3, hjust = "left",
            vjust = "center") + # when angle = 90, vjust = "center" aligns in the correct place on the timescale x-axis
  scale_x_continuous(breaks = 1:35, labels = 1:35, expand = c(0, 0)) +
  scale_y_continuous(limits = c(-40, 150), expand = c(0, 0)) + # to align date labels on x-axis, set min y limit same y in aes()
  coord_cartesian(xlim = c(4.25, 32.5)) +
  # ------------------------------
# Add CLIMEX Weekly Growth Index
# ------------------------------
  geom_line(data = climexPredictionsDBMDevPredictions %>%
            filter(siteCode == "mac"),
            aes(x = cumulativeGens, y = weeklyGIRescaled),
            linetype = 3, colour = "black") +
  # -----------------------------------------------------------------
  # Add a segment and label showing the non-trapping period at Minnipa
  # -----------------------------------------------------------------
  geom_segment(data = data_frame(siteCode = "mac", date = range(noTrapDates)) %>%
               left_join(dbmDevPredictionsDaily12pm),# grab the DBM cumulativeGens for plotting against x.
               aes(x = min(cumulativeGens), xend = max(cumulativeGens), y = 2.5, yend = 2.5),
               linetype = 2, arrow = myArrow) +
  geom_label(data = data_frame(siteCode = "mac", date = range(noTrapDates)) %>%
               left_join(dbmDevPredictionsDaily12pm),# grab the DBM cumulativeGens for plotting against x.
             aes(x = mean(cumulativeGens), y = 2.5, label = "No trapping")) +
  # --------------------------------------------------------------
  # Add segments and text labels to highlight the cropping windows
  # --------------------------------------------------------------
  geom_segment(data = canolaWindows %>% filter(siteCode == "mac"),
               aes(x = datePlant, xend = dateHarvest, y = -7.5, yend = -7.5),
               colour = "green", size = 3) +
  geom_label(data = canolaWindows %>% filter(siteCode == "mac"),
             aes(x = midPoint, y = -7.5, label = "Canola")) +
  # # -------------------------------------------------------------------
# # Draw segments retracing significant flights by one generation (red)
# # -------------------------------------------------------------------
# # (note the flight in May 2015! no prior pop.)
# geom_segment(data = mothFlightsMinnipa,
#              (aes(x = retrace1Gen, xend = retrace1Gen,
#                   y = 0, yend = count)), linetype = 4, colour = "red") +
# geom_segment(data = mothFlightsMinnipa,
#              (aes(x = retrace1Gen, xend = cumulativeGens,
#                   y = count, yend = count)), linetype = 1, colour = "red",
#              arrow = myArrow) +
# # ----------------------------------------------------------
# # draw segments retracing smaller flights of interest (blue)
# # ----------------------------------------------------------
# # note how there appear to be a local cohort going through discrete generations during summer 2016!
# geom_segment(data = mothFlightsMinnipaSmaller,
#              (aes(x = retrace1Gen, xend = retrace1Gen, y = 0, yend = count)), 
#              linetype = 4, colour = "blue") +
# geom_segment(data = mothFlightsMinnipaSmaller,
#              (aes(x = retrace1Gen, xend = cumulativeGens, y = count, yend = count)), 
#              linetype = 1, colour = "blue", arrow = myArrow) +
# ------------------------------------------------------------------
# draw segments retracing the final SELECTED flights of interest (blue): replaces the exploratory retracing above.
# --------------------------------------------------------------------
  geom_segment(data = selectedFlightsMinnipa,
               (aes(x = retrace1Gen, xend = retrace1Gen, y = 0, yend = count)),
               linetype = 4, colour = "blue") +
  geom_segment(data = selectedFlightsMinnipa,
               (aes(x = retrace1Gen, xend = cumulativeGens, y = count, yend = count)), 
               linetype = 1, colour = "blue", arrow = myArrow) +
  geom_line(aes(group = siteCode)) # overplot the actual moth trap data last
ggsave("lightTrapsFlightsMinnipaDraftV2.pdf", width = 20, height = 8)




# ***************
# Sherwood (sh) #
# ***************

# # determine significant flight dates for drummond point
# mothFlightsSherwood <- mothData %>% 
#   # first sum the male and female counts
#   group_by(siteCode, date) %>% 
#   summarise(count = sum(abs(count))) %>%
#   filter(siteCode == "sh", count > 20) %>%
#   # add the cumulative DBM generations
#   left_join(dbmDevPredictionsDaily12pm, by = c("siteCode", "date")) %>%
#   mutate(retrace1Gen = cumulativeGens - 1)

# Now select the flight dates of interest for re-tracing and plotting in the main figure
# These were selected after visual examination and retracing of many peaks
selectedFlightDatesSherwood <- c(
  "2014-09-24", "2014-10-06", "2014-10-25", "2015-09-13", "2015-11-08", "2016-10-07",
  "2016-10-21", "2016-10-30"
  ) %>% as_date()

selectedFlightsSherwood <- mothData %>% # add to the main data for plotting with correct aes() 
  group_by(siteCode, date) %>%  # first sum the male and female counts
  summarise(count = sum(abs(count))) %>%  
  filter(siteCode == "sh", date %in% selectedFlightDatesSherwood) %>%
  left_join(dbmDevPredictionsDaily12pm, by = c("siteCode", "date")) %>% # add the cumulative DBM generations
  mutate(retrace1Gen = cumulativeGens - 1)


# plot the significant flights and retrace 1 generation (Sherwood)
(plotSherwood <- mothData %>%
  filter(siteCode == "sh") %>%
  group_by(siteCode, date) %>%
  summarise(count = sum(abs(count))) %>%
  left_join(dbmDevPredictionsDaily12pm, by = c("siteCode", "date")) %>%
  mutate(siteName = str_replace_all(siteCode, siteNames)) %>%
  ggplot(aes(x = cumulativeGens, y = count)) + # this just adds the siteName nicely above the plot
  facet_wrap(~siteName, ncol = 1, scales = "free") +
  theme_light() +
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black"), # theme_light has white text which is invisible.
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey"),
        panel.grid.minor.x = element_blank()) +
  labs(x = "Cumulative development (generations)", y = "Mean no. moths trapped per night") +
  # Add chronological dates as labels
  geom_text(data = breaksDBMGens %>% filter(siteCode == "sh"),
            aes(x = cumulativeGens , y = -85, label = format(date, "-%d-%b-%Y")), # extra dash as tick mark
            angle = 90, size = 3, hjust = "left",
            vjust = "center") + # when angle = 90, vjust = "center" aligns in the correct place on the timescale x-axis
  scale_x_continuous(breaks = 1:35, labels = 1:35, expand = c(0, 0)) +
  scale_y_continuous(limits = c(-85, 350), expand = c(0, 0)) + # to align date labels on x-axis, set min y limit same y in aes()
  coord_cartesian(xlim = c(3.1, 28)) +
  # Add CLIMEX Weekly Growth Index # ------------------------------
  geom_line(data = climexPredictionsDBMDevPredictions %>%
              filter(siteCode == "sh"),
            aes(x = cumulativeGens, y = weeklyGIRescaled), 
            linetype = 3, colour = "black") +
  # --------------------------------------------------------------
  # Add segments and text labels to highlight the cropping window
  # --------------------------------------------------------------
  geom_segment(data = canolaWindows %>% filter(siteCode == "sh"),
               aes(x = datePlant, xend = dateHarvest, y = -15, yend = -15), 
               colour = "green", size = 3) +
  geom_label(data = canolaWindows %>% filter(siteCode == "sh"),
             aes(x = midPoint, y = -15, label = "Canola")) +
  # # -------------------------------------------------------------------
# # Draw segments retracing significant flights by one generation (red)
# # -------------------------------------------------------------------
# # (note the flight in May 2015! no prior pop.)
# geom_segment(data = mothFlightsSherwood,
#              (aes(x = retrace1Gen, xend = retrace1Gen,
#                   y = 0, yend = count)), linetype = 4, colour = "red") +
# geom_segment(data = mothFlightsSherwood,
#              (aes(x = retrace1Gen, xend = cumulativeGens,
#                   y = count, yend = count)), linetype = 1, colour = "red",
#              arrow = myArrow) +
# # ----------------------------------------------------------
# # draw segments retracing smaller flights of interest (blue)
# # ----------------------------------------------------------
# # note how there appear to be a local cohort going through discrete generations during summer 2016!
# geom_segment(data = mothFlightsSherwoodSmaller,
#              (aes(x = retrace1Gen, xend = retrace1Gen, y = 0, yend = count)), 
#              linetype = 4, colour = "blue") +
# geom_segment(data = mothFlightsSherwoodSmaller,
#              (aes(x = retrace1Gen, xend = cumulativeGens, y = count, yend = count)), 
#              linetype = 1, colour = "blue", arrow = myArrow) +
# ------------------------------------------------------------------
# draw segments retracing the final SELECTED flights of interest (red): replaces the exploratory retracing above.
# --------------------------------------------------------------------
  geom_segment(data = selectedFlightsSherwood,
               (aes(x = retrace1Gen, xend = retrace1Gen, y = 0, yend = count)),
               linetype = 4, colour = "blue") +
  geom_segment(data = selectedFlightsSherwood,
               (aes(x = retrace1Gen, xend = cumulativeGens, y = count, yend = count)), 
               linetype = 1, colour = "blue", arrow = myArrow) +
  geom_line(aes(group = siteCode))) # overplot the actual moth trap data last
ggsave("lightTrapsFlightsSherwoodDraftV2.pdf", width = 20, height = 8)


# ****************
# Two Wells (tw) #
# ****************

# # determine significant flight dates for drummond point
# mothFlightsTwoWells <- mothData %>% 
#   # first sum the male and female counts
#   group_by(siteCode, date) %>% 
#   summarise(count = sum(abs(count))) %>%
#   filter(siteCode == "tw", count > 50) %>%
#   # add the cumulative DBM generations
#   left_join(dbmDevPredictionsDaily12pm, by = c("siteCode", "date")) %>%
#   mutate(retrace1Gen = cumulativeGens - 1)

# Now select the flight dates of interest for re-tracing and plotting in the main figure
# These were selected after visual examination and retracing of many peaks
selectedFlightDatesTwoWells <- c(
  "2014-09-17", "2014-09-24", "2014-10-22", "2014-11-21", "2015-09-13", "2015-10-05", "2015-10-23"
  ) %>% as_date()

selectedFlightsTwoWells <- mothData %>% # add to the main data for plotting with correct aes() 
  group_by(siteCode, date) %>%  # first sum the male and female counts
  summarise(count = sum(abs(count))) %>%  
  filter(siteCode == "tw", date %in% selectedFlightDatesTwoWells) %>%
  left_join(dbmDevPredictionsDaily12pm, by = c("siteCode", "date")) %>% # add the cumulative DBM generations
  mutate(retrace1Gen = cumulativeGens - 1)

# plot the significant flights and retrace 1 generation (Two Wells)
mothData %>%
  filter(siteCode == "tw") %>%
  group_by(siteCode, date) %>%
  summarise(count = sum(abs(count))) %>%
  left_join(dbmDevPredictionsDaily12pm, by = c("siteCode", "date")) %>%
  mutate(siteName = str_replace_all(siteCode, siteNames)) %>%
  ggplot(aes(x = cumulativeGens, y = count)) + # this just adds the siteName nicely above the plot
  facet_wrap(~siteName, ncol = 1, scales = "free") +
  theme_light() +
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black"), # theme_light has white text which is invisible.
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey"),
        panel.grid.minor.x = element_blank()) +
  labs(x = "Cumulative development (generations)", y = "Mean no. moths trapped per night") +
  # Add chronological dates as labels
  geom_text(data = breaksDBMGens %>% filter(siteCode == "tw"),
            aes(x = cumulativeGens , y = -75, label = format(date, "-%d-%b-%Y")), # extra dash as tick mark
            angle = 90, size = 3, hjust = "left",
            vjust = "center") + # when angle = 90, vjust = "center" aligns in the correct place on the timescale x-axis
  scale_x_continuous(breaks = 1:35, labels = 1:35, expand = c(0, 0)) +
  scale_y_continuous(limits = c(-75, 300), expand = c(0, 0)) + # to align date labels on x-axis, set min y limit same y in aes()
  coord_cartesian(xlim = c(4, 32)) +
  # ------------------------------
  # Add CLIMEX Weekly Growth Index 
  geom_line(data = climexPredictionsDBMDevPredictions %>%
            filter(siteCode == "tw"),
            aes(x = cumulativeGens, y = weeklyGIRescaled),
            linetype = 3, colour = "black") +
  # --------------------------------------------------------------
  # Add segments and text labels to highlight the cropping window
  # --------------------------------------------------------------
  geom_segment(data = canolaWindows %>% filter(siteCode == "tw"),
               aes(x = datePlant, xend = dateHarvest, y = -15, yend = -15),
               colour = "green", size = 3) +
  geom_label(data = canolaWindows %>% filter(siteCode == "tw"),
             aes(x = midPoint, y = -15, label = "Canola")) +
  # # -------------------------------------------------------------------
# # Draw segments retracing significant flights by one generation (red)
# # -------------------------------------------------------------------
# # (note the flight in May 2015! no prior pop.)
# geom_segment(data = mothFlightsTwoWells,
#              (aes(x = retrace1Gen, xend = retrace1Gen,
#                   y = 0, yend = count)), linetype = 4, colour = "red") +
# geom_segment(data = mothFlightsTwoWells,
#              (aes(x = retrace1Gen, xend = cumulativeGens,
#                   y = count, yend = count)), linetype = 1, colour = "red",
#              arrow = myArrow) +
# # ----------------------------------------------------------
# # draw segments retracing smaller flights of interest (blue)
# # ----------------------------------------------------------
# # note how there appear to be a local cohort going through discrete generations during summer 2016!
# geom_segment(data = mothFlightsTwoWellsSmaller,
#              (aes(x = retrace1Gen, xend = retrace1Gen, y = 0, yend = count)), 
#              linetype = 4, colour = "blue") +
# geom_segment(data = mothFlightsTwoWellsSmaller,
#              (aes(x = retrace1Gen, xend = cumulativeGens, y = count, yend = count)), 
#              linetype = 1, colour = "blue", arrow = myArrow) +
# ------------------------------------------------------------------
# draw segments retracing the final SELECTED flights of interest (red): replaces the exploratory retracing above.
# --------------------------------------------------------------------
  geom_segment(data = selectedFlightsTwoWells,
               (aes(x = retrace1Gen, xend = retrace1Gen, y = 0, yend = count)),
               linetype = 4, colour = "blue") +
  geom_segment(data = selectedFlightsTwoWells,
               (aes(x = retrace1Gen, xend = cumulativeGens, y = count, yend = count)), 
               linetype = 1, colour = "blue", arrow = myArrow) +
  geom_line(aes(group = siteCode)) # overplot the actual moth trap data last
ggsave("lightTrapsFlightsTwoWellsDraftV2.pdf", width = 20, height = 8)



# ======================================================
# Plot raw moth data on a date axis showing sex ratios #
# ======================================================

# Supplementary figure:

# Set the labels and breaks for each site so that 0 line sits in the middle of each plots
# Different scales on different facets is not implemented in ggplot2, 
# So use scales = "free_y") with a hack: add a single white bar underneath the main data on the plot,
# With transparent bar magnitude of max male or female counts for that site.
# Under plot the white bar on the same date as the actual count occurred, so that main data will hide it. 
maxCountEitherSex <- mothData %>%
  group_by(siteCode) %>%
  summarise(count = max(abs(count), na.rm = TRUE)) %>%
  # now grab the dates the max counts occurred (so we can underplot the white bars on the same dates) 
  left_join(mothData %>% 
              mutate(count = abs(count)) %>% 
              dplyr::select(-sex), by = c("siteCode", "count")) %>%
  # now duplicate these max counts in positive (female) and negative (male) space
  # to ensure xaxis at zero count sits at centre when these data are plotted
  bind_rows(., .) %>% # duplicate the data_frame
  mutate(count = ifelse(duplicated(count), count * -1, count),
         sex   = ifelse(count < 1, "Male", "Female"),
         siteName = str_replace_all(siteCode, siteNames)) %>% # need siteName var here to mesh with main ggplot aes()
  arrange(siteCode)
  
# multipanel plot  
mothData %>%
  mutate(sex = str_replace_all(sex, c("^female" = "Female", "^male" = "Male")),
         sex = factor(sex, levels = c("Female", "Male")),
         siteName = str_replace_all(siteCode, siteNames)) %>%
  ggplot(aes(x = date, y = count, colour = sex)) +
  # first underplot a transparent bar to set the male/female y scale for each site
  # this horizontally centres the zero x-axis line for each plot
  geom_bar(data = maxCountEitherSex, stat = "identity",
           aes(x = date, y = count, colour = sex), colour = "transparent", fill = "transparent") +
  geom_line() +
  labs(x = NULL, y = "Mean no. moths trapped per night", colour = NULL) +
  scale_colour_manual(values = c("violetred1", "steelblue")) +
  scale_x_date(breaks = scales::date_breaks(width = "1 month"),
               labels = scales::date_format("%b-%Y")) +
  scale_y_continuous(labels = function(x){abs(x)}) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.position = c(0.94, 0.92),
        legend.title = element_blank()) +
  facet_wrap(~siteName, scales = "free_y")
ggsave("lightTrapsMultipanelNightlyDraft.pdf", width = 10, height = 6.5)

# =====================================
# summarise and plot monthly trap catch
# =====================================

# Possible supplementary figure, but
# mainly for the text, so that we can describe whether which months
# were detected in.
mothData %>% 
  mutate(month = cut(date, breaks = "1 month"),
         month = as_date(month),
         siteName = str_replace_all(siteCode, siteNames),
         sex = str_replace_all(sex, c("^female" = "Female", "^male" = "Male"))) %>%
  filter(abs(count) > 0) %>% # just to avoid a pink/blue line appearing for that month
  ggplot(aes(x = month, y = count, fill = sex)) +
  geom_bar(stat = "identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.minor = element_blank()) +
  labs(x = NULL, y = "Mean no. moths trapped per month", fill = NULL) +
  scale_x_date(breaks = scales::date_breaks(width = "1 month"),
               labels = scales::date_format("%b-%Y")) +
  scale_fill_manual(values = c("violetred1", "steelblue")) +
  facet_wrap(~siteName)
ggsave("lightTrapsMultipanelMonthlyDraft.pdf", width = 9, height = 5)

# summarise monthly trap catch by site for describing in the manuscript text
mothData %>% 
  mutate(month = cut(date, breaks = "1 month"),
         month = as_date(month)) %>%
  group_by(siteCode, month) %>%
  summarise(totalMoths = sum(abs(count), na.rm = TRUE)) %>%
  # just cast sites to separate columns to read the summary
  reshape2::dcast(month ~ siteCode, value.var = "totalMoths") %>%
  mutate_if(is.numeric, round, 2)

# =======================================
# summarise sex ratios in trap catch data
# =======================================

# overall sex ratio
mothData %>% 
  group_by(sex) %>%
  summarise(count = sum(abs(count), na.rm = TRUE))

# change in sex ratio by month
mothData %>% 
  mutate(month = lubridate::month(date)) %>%
  group_by(month, sex) %>%
  summarise(count = sum(abs(count), na.rm = TRUE)) %>%
  reshape2::dcast(month ~ sex, value.var = "count") %>%
  filter(male > 1, female > 1) %>%
  mutate(month = factor(month), 
         sexRatio = male / female) %>%
  ggplot(aes(x = month, y = sexRatio, group = 1)) +
  geom_line()

# ===============================================================
# Examine the nearest detection to dispersal flight in April 2015
# ===============================================================

# The Drummond Point light trap detected possible dispersal events in May 2015.
# Moths did not correspond with peaks a generation earlier in April.
# From autumn survey data from 2015, the nearest April detections were:
# -- Moths
# -- Larvae
lightTrapDP    <- c(135.33121, -34.13687)
mothsNearest  <- c(135.2684,  -33.92549) # Swamp, south of Sheringa (look at autumn survey map)
larvaeNearest <- c(134.8598,  -33.55526) # Walkers Beach, sea rocket
geosphere::distGeo(lightTrapDP, mothsNearest)  / 1000  # 24.2 km
geosphere::distGeo(lightTrapDP, larvaeNearest) / 1000  # 77.9 km

#> No conclusions can be drawn about the seasonal patterns of sex ratios
#> because low numbers of moths trapped in summer/autumn bias the results.

# End script
# ##############################################

# # --------
# # steve's example for automating many segments on plots 
# grawSgments <- function(plot, i){
#   for (val in 1){
#     plot <- plot + geom_segment()
#   }
#   plot
# }
# # end example 
# # ------------- 
#   














