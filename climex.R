# ######################################################################################################
# Plot time series rasters of Climex Weekly Growth Index (WGI), Annual Growth Index (AGI) 
# ... and Ecoclimatic Index (EI) for the diamondback moth, Plutella xylostella
# K. Perry, 1/7/2017

# This script uses the workspace in `climexOutput.Rproj`

# Notes:
# Climex is a bioclimatic modelling program that predicts the climate niche for a \
# \ species based on climate data and the known responses of the species to climate variables. 
# ... Climex outputs several indices:
# ... (i)   Weekly Growth Index (WGI): Suitability of a location for population growth 
# ... (ii)  Annual Growth Index: Yearly snapshot of suitability for population growth
# ... (iii) Ecoclimatic Index: Suitability of a location for year-round persistence
# ... (iv)  Several stress indices, taken into account in calculatiob of the EI.

# This script analyses and plots the results of `CLIMEX Compare Locations Years` \
# \ for the diamondback moth in Australia for three full seasons, summer 2013 to spring 2016.
# \ Climex was run using `grim` formatted files, created using a script `asciiToGrim.R`

# This script does:
# ... extracts climex run results from a (8Gb) netcdf (WGI) 
# ... creates movies of time series raster plots of WGI for Australia, Eyre Peninsula
# ... analyses the proportions of raster pixels in WGI classes across time for Eyre.
# ... creates static plots of annual AGI and EI indices
# Memory use is an issue. 

# Addditional note: 
# ... in the code, the stacked area time series graphs need to be plotted before the maps.
# ... this is because zero values are set to NA for maps to allow underplotting a grey polygon.

# ######################################################################################################

library(tidyverse)
library(ncdf4)
library(ggplot2)
library(RColorBrewer)
library(raster)
library(rasterVis)
library(rgdal)
library(lubridate)
library(parallel) # for mclapply when running on Clarence or other linux pc

# NOTES:
# CLARENCE: You need to create an output directory for SouthAus


# #######################
# Set global file paths #
# #######################

# ---------------
# To run locally:
# ---------------
basePath <- file.path("C:", "UserData", "Kym", "PhD", "Data")
old  <- file.path(basePath, "climexOutput") # home dir
ncWd <- file.path(basePath, "climexOutput", "compareLocationsYears") # path to netcdf files
ausPolyWd  <- file.path(basePath, "geoSpatialData", "ausAdminAreas_rds") # path to aus polygons
ibraPolyWd <- file.path(basePath, "geoSpatialData", "ibraBioRegions")    # path to ibra polygon
outEyre    <- file.path(basePath, "climexOutput", "compareLocationsYears", "eyrePlotsWGI") # outPath for Eyre Peninsula plots/movies
outAus     <- file.path(basePath, "climexOutput", "compareLocationsYears", "ausPlotsWGI")  # outPath for Australia plots/movies
outSouthAus   <- file.path(basePath, "climexOutput", "compareLocationsYears", "southAusPlotsWGI")

# -----------------------------
# To run on Clarence remote PC:
# -----------------------------
# This is not platform independent as I had to add the preceding `/` for the root paths to work

# old  <- file.path("/home", "perry")
# ncWd <- old # needed.
# ausPolyWd   <- file.path("/home", "perry", "ausPoly")
# ibraPolyWd  <- file.path("/home", "perry", "ibraPoly")
# outEyre     <- file.path("/home", "perry", "eyrePlotsWGI")
# outAus      <- file.path("/home", "perry", "ausPlotsWGI")
# outSouthAus <- file.path("/home", "perry", "southAusPlotsWGI")

# ###################
# Set global themes #
# ###################

# set consistent palette for rasters and stacked area time series plots
bl      <- brewer.pal(9, 'Blues')
yor     <- brewer.pal(9, 'YlOrRd')
greyPal <- brewer.pal(9, "Greys")
wyor  <- c("white", brewer.pal(9, 'YlOrRd')[2:9]) 

myPal <- yor # set palette here
myTheme <- rasterTheme(region = myPal) # for raster time series plot frames
myThemeNoBorder <- rasterTheme(region = myPal, axis.line = list(col = "transparent"))#,
                               #layout.heights = list(xlab.key.padding = 3)) # add space if colour key is at the bottom
myRamp  <- colorRampPalette(myPal)     # for stacked area time series plot
myBreaks <- seq(0, 100) / 100 # set legend breaks for raster::levelplot (make a smooth colour gradient. max length = 100) 



# ##################
# Define functions #
# ##################

# sum non-NA values but return NA if all are values are NA
sumValues <- function(x){
  if(all(is.na(x))){return(NA)}
  sum(x, na.rm = TRUE)
}

# ########################################################################
# Download polygons of administrative boundaries for regions of interest #
# ########################################################################


# latlong projection string
projLL <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Australian coastline and states SpatialPointsDataFrame
setwd(ausPolyWd)
ausPoly <- readRDS("AUS_adm1.rds") # coast and states
setwd(old)
proj4string(ausPoly) # check projection
# [1] "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# IBRA bioregions polygon shapefile for the Eyre Peninsula region, South Australia
ibraSub <- readOGR(dsn = ibraPolyWd, layer = "ibra7_subregions")
projection(ibraSub) # check, transform projection
# [1] "+proj=longlat +ellps=GRS80 +no_defs"
ibraSub <- spTransform(ibraSub, projLL) # transform to latlong
projection(ibraSub)
# [1] "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# extract ibra subregion polygons for Eyre Peninsula (Eyre Yorke block, 03, 04, 05)
eyrePoly <- ibraSub[ibraSub$SUB_CODE_7 %in% c("EYB03", "EYB04", "EYB05"), ]


# ##################################################
# Create Raster* data from CLIMEX model run output #
# ##################################################

# extract data from a netcdf file containing CLIMEX run output results  
# set time period to extract
origin <- as.Date("2013-01-01") # start of the time series (climex run start)
timePoint  <- 49           # begin at time step 49 (04-Dec-2013)
nTimePoints <- 208 - timePoint - 3 # num time steps. drop the last 3 layers to finish at 205 (04-Dec-2016). 

ncPath <- paste0(ncWd, "/climexPxylostellaCompareLocationsYearsWeeklyStep2013-2016.nc")
stopifnot(file.exists(ncPath))
nc <- nc_open(ncPath)
print(nc) # to get information about the netcdf, variables
long <- ncvar_get(nc, varid = "Longitude") # the unique longs
lat  <- ncvar_get(nc, varid = "Latitude")  # the unique lats
step <- ncvar_get(nc, varid = "Step", start = timePoint, count = nTimePoints) # dimension variable (single dimension)
days <- ncvar_get(nc, varid = "Days since Start", start = c(1, 1, timePoint), count = c(1, 1, nTimePoints))
simDate <- origin + days
WGI  <- ncvar_get(nc, varid = "Weekly Growth Index", start = c(1, 1, timePoint), count = c(-1, -1, nTimePoints))
nc_close(nc) # always close .nc connection to avoid data loss 

coords <- expand.grid(long, lat) %>%
  setNames(c("longX", "latY")) %>%
  mutate(longX = as.numeric(longX),
         latY  = as.numeric(latY))

# create a template raster
e <- extent(c(xmin = min(long), xmax = max(long),
              ymin = min(lat),  ymax = max(lat)))
r <- raster(e, ncol = length(long), nrow = length(lat))

# create a raster brick
# ... use `rasterize`, as coordinate points are not on a regular grid
# ... very slow on a single core ( ~ 45 mins for 208 layers). mclapply parallelises where available. 
# ausWGIBrick <- WGI %>%
#   plyr::alply(.margins = 3) %>%
#   mclapply(function(x) {
#     rasterize(coords, r, field = x, fun = mean)
#   }) %>%
#   brick() 

# set raster projection 
proj4string(ausWGIBrick)
# [1] NA
# . SILO gridded data (from which climex run was derived) are formatted in a geographic 
# . coordinate system, therefore a longlat projection is appropriate
projection(ausWGIBrick) <- projLL


# trim to the Australian coast
# ausWGIBrick <- ausWGIBrick %>%
#   mask(ausPoly) %>%
#   trim() %>%
#   setNames(as.character(format(simDate, "%d-%b-%Y"))) # set layer names


# trim to Eyre Peninsula poygons, IBRA block Eyre Yorke Block 03, 04, 05
eyreWGIBrick <- ausWGIBrick %>%
  crop(y = eyrePoly) %>% # crop needed for y = *SpatialPointsDataFrame
  mask(eyrePoly) %>%
  trim()

# #######################################################################################
# Create a static figure of MONTHLY mean rasters of CLIMEX weekly GI for the manuscript #
# #######################################################################################

# This section plots colonisation data on top of raster maps of CLIMEX Weekly Growth Index

# Aggregate weekly into monthly mean raster stacks 
timeZone <- "Australia/Adelaide" 
layerWeeks  <- names(ausWGIBrick) %>% str_replace("^X", "") %>% as_date(format = "%d.%b.%Y", tz = timeZone) # parse the date
layerMonths <- cut(layerWeeks, breaks = "1 month") %>% as.character() %>% as_date(tz = tzone)

# Define the geographic extent of interest for the colonisation study
colonisationExtent <- c(xmin = 133, xmax = 143, ymin = -32, ymax = -38) %>%
  matrix(byrow = TRUE, nrow = 2) %>%
  extent()

# Make a monthly raster brick, crop to South Australia, subset the months
southAusWGIBrickMonths <- ausWGIBrick %>%
  crop(colonisationExtent) %>%
  stackApply(layerMonths, mean) %>% # layerMonths is used as the unique index 
  setNames(unique(layerMonths)) %>%
  subset(c(3:9, 15:21, 27:33)) # keep Feb to Aug each year for the monthly plot (not much info in Sep)

# Grab the full sequence of months
monthSeq <- names(southAusWGIBrickMonths) %>% 
  str_replace("^X", "") %>% 
  as_date(tz = timeZone) %>%
  unique()

# Read in colonisation data (generated in R script `colonisation.R`), 
# then aggregate to monthly values for juvenile presence.
# To plot the points on the corresponding panels, we need to create a list 
# with same length as number of raster layers
colonisationWeekly <- file.path(
  "C:", "UserData", "Kym", "PhD", "Data", "colonisation", "colonisationWeekly.csv"
) %>% read_csv()

colonisationMonthly <- colonisationWeekly %>% 
  mutate(month = cut(week, breaks = "1 month"),
         month = as.character(month) %>% as_date(tz = "Australia/Adelaide")) %>%
  filter(month %in% monthSeq) %>%
  group_by(siteCode, latitude, longitude, month) %>% # carry lat/lon along for the ride
  summarise(juvPresence = max(juvPresence, na.rm = TRUE)) %>% # take the juvenile presence up to end of month
  mutate(juvPresence = ifelse(juvPresence == -Inf, NA, juvPresence)) %>%
  ungroup()

# create a dataframe with all the points for each year for all months in monthSeq so that
# (1) we can plot where the sites are in all months as needed,
# (2) we need lat/lons for all dfs to easily convert the df list to SP class
fitPointsEachMonth <- function(pointsDf, month){
  month %>% # for each month, fit a df with all the lat/lons
    map(function(x){
     pointsDf %>% mutate(month = x) 
    })
  }

fitMonthsDf <- colonisationMonthly %>%
  mutate(year = lubridate::year(month)) %>%
  split(.$year) %>%
  map(distinct, siteCode, .keep_all = TRUE) %>%
  map(mutate, juvPresence = NA) %>%
  # for each year, take the monthSeq *for that year* and fit a points df for each month
  map2(monthSeq %>% split(year(.)),  fitPointsEachMonth) %>%
  map(bind_rows) %>%
  bind_rows() %>%
  dplyr::select(-year)

# Split into monthly list elements and convert to SP class for plotting on the panels
colonisationMonthlySP <- colonisationMonthly %>%
  bind_rows(fitMonthsDf) %>% # to fit the dfs, stack and sum values for unique month x location (a join method would create dup rows)
  group_by(siteCode, latitude, longitude, month) %>%
  summarise(juvPresence = sumValues(juvPresence)) %>%
  split(.$month) %>% # these do split in correct date order
  map(function(df){
    df %>% SpatialPointsDataFrame(
      coords = df[, c("longitude", "latitude")],
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  })


out <- colonisationMonthlySP # should prob call these `points``
southAusWGIBrickMonths[southAusWGIBrickMonths == 0] <- NA # to underplot a grey polygon, make 0 values transparent

trellis.device(png, file = file.path(outSouthAus, "southAusMonthsClimexColonisationDRAFT.png"),
               res = 300, width = 3000, height = 1500)
rasterVis::levelplot(southAusWGIBrickMonths,
                     maxpixels = ncell(southAusWGIBrickMonths),
                     layout    = c(7, 3),
                     par.strip.text = list(cex = 0.9, lines = 2),
                     at = myBreaks,
                     scales = list(draw = TRUE, col = 'black',
                                   x  = list(cex = 0.7),
                                   y  = list(cex = 0.7)),
                     xlab = list(cex = 0.9, label = "Longitude"),
                     ylab = list(cex = 0.9, label = "Latitude"),
                     colorkey = list(height = 0.2, width = 1,
                                     labels = list(cex = 0.8, at = c(0, 1))),
                     par.settings = myThemeNoBorder,
                     names.attr = names(southAusWGIBrickMonths) %>% 
                       str_replace("^X", "") %>% 
                       parse_date_time("Ymd", tz = timeZone) %>% 
                       format("%b %Y")) +
  layer(sp.polygons(ausPoly, col = NA, fill = 'lightgrey'), under = TRUE) +
  layer(sp.lines(ausPoly, lwd = 0.4, col = 'white')) +
  layer(sp.points(out[[panel.number()]], # panel.number draws the point `packets` (list indices) on corresponding panels (index)
                  cex = 0.5, 
                  pch  = ifelse(is.na(out[[panel.number()]]$juvPresence), 3, 21),
                  # print non-NA values with a circle border
                  col  = ifelse(is.na(out[[panel.number()]]$juvPresence), "transparent", alpha(greyPal[8], 1)),
                  fill = ifelse(out[[panel.number()]]$juvPresence == 0, "transparent", 
                                ifelse(out[[panel.number()]]$juvPresence == 1, alpha(greyPal[8], 1), "transparent"))))
dev.off()

# ####################################################################################
# Create weekly raster time series movie of CLIMEX Weekly GI and colonisation timing #
# ####################################################################################

# This will be a Supplementary Figure.

# TO DO: Add a subcaption for the movie.

layerWeeks  <- names(ausWGIBrick) %>% 
  str_replace("^X", "") %>% as_date(format = "%d.%b.%Y", tz = timeZone) # parse the date

# Make a weekly raster brick, crop to South Australia, subset the months 
southAusWGIBrickWeeks <- ausWGIBrick %>%
  crop(colonisationExtent) %>%
  stackApply(layerWeeks, mean) %>% # layerMonths is used as the unique index 
  setNames(unique(layerWeeks)) %>%
  subset(c(10:43, 62:95, 114:147)) # keep weeks in Feb to Sept each year 

# Grab the full sequence of weeks from the raster stack
weekSeq <- names(southAusWGIBrickWeeks) %>% 
  str_replace("^X", "") %>% 
  as_date(tz = timeZone) %>%
  unique() %>%
  `+`(2) # Add 2 days so that weeks (from CLIMEX analysis) align with weeks from the colonisation analaysis.
# Explain this in the figure text.

colonisationWeeks <- colonisationWeekly %>% # colonisationWeekly has the colonisation data (rw)
  dplyr::select(siteCode, latitude, longitude, juvPresence, week) %>%
  filter(week %in% weekSeq) %>%
  arrange(week)

# create a dataframe with all the points for each year for all weeks in weekSeq, so that
# (1) we can plot where the sites are in all weeks if needed, and 
# (2) we need lat/lons for all dfs to be able to easily convert all elements of df list to SP class
fitPointsEachWeek <- function(pointsDf, week){
  week %>% # for each week, fit a df with all the lat/lons
    map(function(x){
      pointsDf %>% mutate(week = x) 
    })
  }

fitWeeksDf <- colonisationWeeks %>%
  mutate(year = lubridate::year(week)) %>%
  split(.$year) %>%
  map(distinct, siteCode, .keep_all = TRUE) %>%
  map(mutate, juvPresence = NA) %>%
  # for each year, take the weekSeq *for that year* and fit a points df for each week
  map2(weekSeq %>% split(year(.)),  fitPointsEachWeek) %>%
  map(bind_rows) %>%
  bind_rows() %>%
  dplyr::select(-year)

# Split into monthly list elements and convert to SP class for plotting on the panels
colonisationWeeksSP <- colonisationWeeks %>%
  bind_rows(fitWeeksDf) %>% # to fit the dfs, stack and sum values for unique week x location (a join method would create dup rows)
  group_by(siteCode, latitude, longitude, week) %>%
  summarise(juvPresence = sumValues(juvPresence)) %>%
  split(.$week) %>% # these do split in correct date order
  map(function(df){
    df %>% SpatialPointsDataFrame(
      coords = df[, c("longitude", "latitude")],
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  })

# ====================================================================================================
# define the panel order for rasters and points to plot corresponding weeks side-by side for each year  
# ====================================================================================================

# Define the panel order to plot corresponding weeks for each year side by side
panelOrderDf <- data_frame(
  simDate = as_date(simDate, tz = timeZone) + 2, # adjust 2 days to align weeks with the colonisation analysis
  month   = month(simDate),
  day     = day(simDate)
  ) %>%
  filter(simDate %in% weekSeq) %>%
  arrange(simDate) %>%
  mutate(panelIndex = 1:nrow(.),
         # order with the preceding December first (start of summer)
         monthDecFirst = str_replace(month, "12", "0") %>% as.numeric()) %>%
  # now arrange all the weeks in day-of-year order, plot the preceding December first
  arrange(monthDecFirst, day)

# Now define the points order to correspond with panel order
pointsWeeks <- colonisationWeeksSP[panelOrderDf$panelIndex] # reorder using the panel index

# To underplot a grey polygon, make 0 values transparent
southAusWGIBrickWeeks[southAusWGIBrickWeeks == 0] <- NA 

subCaptionColonisation  <- "CLIMEX Weekly Growth Index (map colours) for the diamondback moth in South Australia, and\n colonisation of winter canola crops by the diamondback moth, shown by the presence (filled circles) or absence (empty circles) of immature lifestages. \n Kym D. Perry and Michael A. Keller, 2018."

# trim white space and add padding at bottom to vertically centre plot accounting for sub-caption 
# (doing this via par.settings arg in the levelplot call didn't adjust the margins.)
latticeOptsDefaults <- lattice.options()
lattice.options(layout.heights = list(bottom.padding = list(x = 3), top.padding   = list(x = 0)))
trellis.device(png, file = file.path(outSouthAus, "weekColonisation%03d.png"), 
               res = 300, width = 3000, height = 1300)
rasterVis::levelplot(southAusWGIBrickWeeks,
                     maxpixels = ncell(southAusWGIBrickWeeks),
                     layout    = c(3, 1),
                     par.strip.text = list(cex = 0.9, lines = 2),
                     at = myBreaks,
                     index.cond = list(panelOrderDf$panelIndex), 
                     #index.cond = list(1:3), #this line for testing only
                     scales = list(draw = TRUE, col = 'black',
                                   x  = list(cex = 0.7),
                                   y  = list(cex = 0.7)),
                     xlab = list(cex = 0.9, label = "Longitude"),
                     ylab = list(cex = 0.9, label = "Latitude"),
                     sub  = list(cex = 0.7, label = subCaptionColonisation, hjust = 0.5, vjust = 1, lineheight = 1),
                     colorkey = list(height = 0.3, width = 1.25,
                                     labels = list(cex = 0.8, at = c(0, 1))),
                     par.settings = myThemeNoBorder,
                     names.attr = names(southAusWGIBrickWeeks) %>% 
                       str_replace("^X", "") %>% 
                       parse_date_time("Ymd", tz = timeZone) %>% 
                       format("%d %b %Y")) +
  layer(sp.polygons(ausPoly, col = NA, fill = "lightgrey"), under = TRUE) +
  layer(sp.lines(ausPoly, lwd = 0.4, col = "white")) +
  layer(sp.points(pointsWeeks[[panel.number()]], # panel.number draws the point `packets` (list indices) on corresponding panels (index)
                  cex = 0.6,
                  pch  = ifelse(is.na(pointsWeeks[[panel.number()]]$juvPresence), 3, 21),
                  # print non-NA values with a circle border
                  col  = ifelse(is.na(pointsWeeks[[panel.number()]]$juvPresence), "transparent", alpha(greyPal[8], 1)),
                  fill = ifelse(pointsWeeks[[panel.number()]]$juvPresence == 0, "transparent",
                                ifelse(pointsWeeks[[panel.number()]]$juvPresence == 1, alpha(greyPal[8], 1), "transparent"))))
dev.off()
lattice.options(latticeOptsDefaults) # reset to defaults

# Make movies
setwd(outSouthAus)
southAusWeeksMovie <- 'ffmpeg -y -i weekColonisation%03d.png -c:v libx264 -crf 10 -b:v 1M -bufsize 1M -vf "setpts=(30/1)*PTS" DiamondbackMoth_Climex_CanolaColonisation_SouthAustralia_KDPerry_Oct2018.mp4'
system(southAusWeeksMovie)
setwd(old)




# ######################################
# Plot time series movie for Australia #
# ######################################

# Plot 3 panels, corresponding to years, side by side for comparison

# Reorder panel plotting order so that corresponding weeks are plottedside by side across years
simDateDf <- data.frame(
  simDate = as.Date(simDate),
  month   = lubridate::month(simDate),
  day     = lubridate::day(simDate)
) %>%
  arrange(simDate) %>%
  mutate(panelIndex = 1:nrow(.),
         # order with the preceding December first (start of summer)
         monthDecFirst = gsub("12", "0", month) %>%
           as.numeric()) %>%
  # now arrange all the weeks in day-of-year order, plotting the preceding December first
  arrange(monthDecFirst, day)

# set zero values to NA to underplot a grey-filled polygon representing zero values
ausWGIBrick[ausWGIBrick == 0] <- NA

subCaptionAus  <- "CLIMEX Weekly Growth Index for the diamondback moth in Australia between summer 2013 and spring 2016. Kym D. Perry and Michael A. Keller, 2018."

# trim white space and add padding at bottom to vertically centre plot accounting for sub-caption 
# (doing this via par.settings arg in the levelplot call didn't adjust the margins.)
latticeOptsDefaults <- lattice.options()
lattice.options(layout.heights = list(bottom.padding = list(x = 3), top.padding   = list(x = 0)))

trellis.device(png, file = file.path(outAus, "Rplot%03d.png"), res = 300, width = 3000, height = 1300)
levelplot(ausWGIBrick, 
          maxpixels = ncell(ausWGIBrick),
          layout = c(3, 1), 
          index.cond = list(simDateDf$panelIndex), # set panel plotting order
          #index.cond = list(1:3), # this line just when testing (subset to three panels first)
          par.strip.text = list(cex = 0.9, lines = 2),
          at = myBreaks,
          scales = list(draw = TRUE, col = 'black',
                        x  = list(cex = 0.7),
                        y  = list(cex = 0.7)),  
          xlab = list(cex = 0.9, label = 'Longitude'),
          ylab = list(cex = 0.9, label = 'Latitude'),
          sub  = list(cex = 0.7, label = subCaptionAus, hjust = 0.5, vjust = 3),
          colorkey = list(height = 0.2, width = 1.25,
                          labels = list(cex = 0.8, at = c(0, 1))), 
          par.settings = myThemeNoBorder,
          names.attr = names(ausWGIBrick) %>% str_replace_all(c("^X" = "", "\\." = " "))) +
  layer(sp.polygons(ausPoly, col = NA, fill = 'lightgrey'), under = TRUE) +
  layer(sp.lines(ausPoly, lwd = 0.4, col = 'white'))
dev.off()
lattice.options(latticeOptsDefaults) # reset to defaults

# create 3 movies with different frame speeds
setwd(outAus)
ausMovie <- 'ffmpeg -y -i Rplot%03d.png -c:v libx264 -crf 10 -b:v 1M -bufsize 1M -vf "setpts=(30/1)*PTS" DiamondbackMoth_Climex_Australia_KDPerry_Oct2018.mp4'
system(ausMovie)
setwd(old)

# to change the frame speed of movie:
# ... change -vf "setpts=(5/1)*PTS" option (5/1 = fast, 20/1 = slow). 
# ... don't change framerate args; will drop or duplicate frames
# -y  = overwrite silently

# ####################################################################
# Create a time series movie for movie for the Eyre Peninsula region #
# ####################################################################

# Plot 3 panels side-by side plots corresponding to years

# set NA to zeros to allow underplotting thea grey polygon to represent the zero values
eyreWGIBrick[eyreWGIBrick == 0] <- NA

subCaptionEP  <- "CLIMEX Weekly Growth Index for the diamondback moth on Eyre Peninsula, South Australia, between summer 2013 and spring 2016.\n Kym D. Perry and Michael A. Keller, 2018."

# trim white space and add padding at bottom to vertically centre plot accounting for sub-caption 
# (doing this via par.settings arg in the levelplot call didn't adjust the margins.)
latticeOptsDefaults <- lattice.options()
lattice.options(layout.heights = list(bottom.padding = list(x = 3), top.padding   = list(x = 0)))

trellis.device(png, file = file.path(outEyre, "Rplot%03d.png"), res = 300, width = 3000, height = 1300)
levelplot(eyreWGIBrick, # for testing, use eyreWGIBrick, 10:12) here, and in names.attr arg.
          maxpixels = ncell(eyreWGIBrick),
          layout = c(3, 1),
          index.cond = list(simDateDf$panelIndex), # set panel order
          #index.cond = list(1:3), # this line for testing only, subset to thee panels first.
          par.strip.text = list(cex = 0.9, lines = 2),
          at = myBreaks,
          scales = list(draw = TRUE, col = 'black',
                        x  = list(cex = 0.7),
                        y  = list(cex = 0.7)),
          xlab = list(cex = 0.9, label = 'Longitude'),
          ylab = list(cex = 0.9, label = 'Latitude'),
          sub  = list(cex = 0.7, label = subCaptionEP, hjust = 0.5, vjust = 2),
          colorkey = list(height = 0.23, width = 1.25,
                          labels = list(cex = 0.8, at = c(0, 1))),
          par.settings = myThemeNoBorder,
          names.attr = names(eyreWGIBrick) %>% str_replace_all(c("^X" = "", "\\." = " "))) +
  layer(sp.polygons(eyrePoly, col = NA, fill = 'lightgrey'), under = TRUE) +
  layer(sp.lines(eyrePoly, lwd = 0.4, col = 'lightgrey'))
dev.off()
lattice.options(latticeOptsDefaults) # reset to defaults


# create movies
setwd(outEyre)
eyreMovie <- 'ffmpeg -y -i Rplot%03d.png -c:v libx264 -crf 10 -b:v 1M -bufsize 1M -vf "setpts=(30/1)*PTS" DiamondbackMoth_Climex_EyrePeninsula_KDPerry_Oct2018.mp4'
system(eyreMovie)
setwd(old)


# #######################################################################
# Plot annual snapshots of ecoclimatic index (EI) and growth index (GI) #
# #######################################################################

# extract data
# set time period to extract
annTp  <- 2 # Time points. Begin at 2014
annNtp <- 3 # Num time steps to extract. Finish at 2016
#
# # extract data for the time period
ncPathAnn <- file.path(ncWd, "climexPxylostellaCompareLocationsYearsYearlySnapshot2013-2016.nc")
stopifnot(file.exists(ncPathAnn))
ncAnn <- nc_open(ncPathAnn)
#print(ncAnn) # to get information about the netcdf
year <- ncvar_get(ncAnn, varid = "Year", start = annTp, count = annNtp)
# Growth Index and Ecoclimatic Index
GI   <- ncvar_get(ncAnn, varid = "GI", start = c(1, 1, annTp), count = c(-1, -1, annNtp))
EI   <- ncvar_get(ncAnn, varid = "EI", start = c(1, 1, annTp), count = c(-1, -1, annNtp))
# Stress indices
CS   <- ncvar_get(ncAnn, varid = "CS", start = c(1, 1, annTp), count = c(-1, -1, annNtp))
HS   <- ncvar_get(ncAnn, varid = "HS", start = c(1, 1, annTp), count = c(-1, -1, annNtp))
WS   <- ncvar_get(ncAnn, varid = "WS", start = c(1, 1, annTp), count = c(-1, -1, annNtp))
DS   <- ncvar_get(ncAnn, varid = "DS", start = c(1, 1, annTp), count = c(-1, -1, annNtp))
nc_close(ncAnn) # always close .nc connection to avoid data loss

# convert coordinate points to raster layers
# ... use `rasterize`, as coordinate points are not on a regular grid
bGI <- plyr::alply(
  GI, .margins = 3,
  .fun = function(x) {
    rasterize(coords, r, field = x, fun = mean)
  }) %>%
  brick()

bEI <- plyr::alply(
  EI, .margins = 3,
  .fun = function(x) {
    rasterize(coords, r, field = x, fun = mean)
  }) %>%
  brick()

# set raster projection
projection(bGI) <- projLL
projection(bEI) <- projLL
#
# Trim rasters to the Australian coastline
bGI <- trim(mask(bGI, ausPoly))
bEI <- trim(mask(bEI, ausPoly))
#
# # set layer index names
bGI[] <- getValues(bGI)
bEI[] <- getValues(bEI)
names(bGI) <- year
names(bEI) <- year

# Combine GI and EI into single stack
sEGI <- raster::stack(bGI, bEI)
sEGI[sEGI == 0] <- NA
#subCaptionAnnual  <- "CLIMEX Annual Growth Index (GI, upper panels) and Ecoclimatic Index (EI, lower panels) for the diamondback moth in Australia, 2014-2016.\nKym D. Perry and Michael A. Keller, 2018"

# trim white space (doing this via par.settings arg in the levelplot call didn't adjust the margins.)
latticeOptsDefaults <- lattice.options()
lattice.options(layout.heights = list(bottom.padding = list(x = 0), top.padding   = list(x = 0)),
                layout.widths  = list(left.padding   = list(x = 0), right.padding = list(x = 0)))

trellis.device(png, file = file.path(ncWd, "annualEIandGIAustralia2014-2016.png"),
               res = 300, width = 1500, height = 1000)
levelplot(sEGI,
          maxpixels = ncell(sEGI),
          layout = c(3, 2),
          par.strip.text = list(cex = 0.6, lines = 2),
          at = myBreaks * 100,
          scales = list(draw = TRUE, col = 'black',
                        x  = list(cex = 0.5),
                        y  = list(cex = 0.5)),
          xlab = list(cex = 0.6, label = 'Longitude'),
          ylab = list(cex = 0.6, label = 'Latitude'),
          #sub  = list(cex = 0.4, label = subCaptionAnnual, hjust = 0.5, vjust = 2),
          colorkey = list(height = 0.175, width = 1,
                          labels = list(cex = 0.5, at = c(0, 100))),
          par.settings = myThemeNoBorder,
          names.attr = names(sEGI) %>% str_replace_all(
            c("^X" = "", "2014.1" = "GI 2014", "2014.2" = "EI 2014", "2015.1" = "GI 2015", 
            "2015.2" = "EI 2015", "2016.1" = "GI 2016", "2016.2" = "EI 2016"))) +
  layer(sp.polygons(ausPoly, col = NA, fill = 'lightgrey'), under = TRUE) +
  layer(sp.lines(ausPoly, lwd = 0.4, col = 'white'))
dev.off()
lattice.options(latticeOptsDefaults) # reset to defaults



# #####################
# Plot stress indices #
# #####################

# Convert coordinate points to raster layers (`rasterize` as coordinate points are not on a regular grid)
bCS <- plyr::alply(
  CS, .margins = 3,
  .fun = function(x) {
    rasterize(coords, r, field = x, fun = mean)
  }) %>%
  brick()

bHS <- plyr::alply(
  HS, .margins = 3,
  .fun = function(x) {
    rasterize(coords, r, field = x, fun = mean)
  }) %>%
  brick()

bWS <- plyr::alply(
  WS, .margins = 3,
  .fun = function(x) {
    rasterize(coords, r, field = x, fun = mean)
  }) %>%
  brick()

bDS <- plyr::alply(
  DS, .margins = 3,
  .fun = function(x) {
    rasterize(coords, r, field = x, fun = mean)
  }) %>%
  brick()

# Set a Long lat raster projection
projection(bCS) <- projLL
projection(bHS) <- projLL
projection(bWS) <- projLL
projection(bDS) <- projLL
#
# Trim rasters to the Australian coastline
bCS <- trim(mask(bCS, ausPoly))
bHS <- trim(mask(bHS, ausPoly))
bWS <- trim(mask(bWS, ausPoly))
bDS <- trim(mask(bDS, ausPoly))

# Set layer index names
bCS[] <- getValues(bCS)
bHS[] <- getValues(bHS)
bWS[] <- getValues(bWS)
bDS[] <- getValues(bDS)
names(bCS) <- paste("CS", year)
names(bHS) <- paste("HS", year)
names(bWS) <- paste("WS", year)
names(bDS) <- paste("DS", year)

# Combine into a single stack
stressIndicesStack <- stack(bCS, bHS, bWS, bDS)

# Set zero values as NA to underplot the polygon
stressIndicesStack[stressIndicesStack == 0] <- NA

# Plot
# # trim white space (doing this via par.settings arg in the levelplot call didn't adjust the margins.)
# latticeOptsDefaults <- lattice.options()
# lattice.options(layout.heights = list(bottom.padding = list(x = 0), top.padding   = list(x = 0)),
#                 layout.widths  = list(left.padding   = list(x = 0), right.padding = list(x = 0)))

trellis.device(png, file = file.path(ncWd, "annualStressIndicesAustralia2014-2016NEW.png"),
               res = 300, width = 1500, height = 1200)
levelplot(stressIndicesStack,
          maxpixels = ncell(stressIndicesStack),
          layout = c(4, 3),
          par.strip.text = list(cex = 0.6, lines = 2),
          #at = myBreaks * 100,
          scales = list(draw = TRUE, col = 'black',
                        x  = list(cex = 0.5),
                        y  = list(cex = 0.5)),
          xlab = list(cex = 0.6, label = 'Longitude'),
          ylab = list(cex = 0.6, label = 'Latitude'),
          #sub  = list(cex = 0.4, label = subCaptionAnnual, hjust = 0.5, vjust = 2),
          colorkey = list(height = 0.175, width = 1,
                          labels = list(cex = 0.5, at = c(0, 100))),
          par.settings = myThemeNoBorder,
          names.attr = names(stressIndicesStack)) +
  layer(sp.polygons(ausPoly, col = NA, fill = 'lightgrey'), under = TRUE) +
  layer(sp.lines(ausPoly, lwd = 0.4, col = 'white'))
dev.off()
# lattice.options(latticeOptsDefaults) # reset to defaults




# End script
# ###############################################################################################
# ###############################################################################################

# ================================================================================
# # Code to plot stacked area time series figures for Australia and Eyre Peninsula
# ================================================================================
# 
# # Define a function to create a dataframe from a raster brick
# # Cuts raster values into classes for plotting a proportional stacked area time series plot 
# rasterToPixelClassDataframe <- function(rasterBrick,  
#                                         dates,         # a vector of layer names in date format
#                                         classBreaks) { # a numeric vector of class breaks including the lowest number bound (include a -ve value to encode zeros)
#   rasterCut <- rasterBrick %>%
#     cut(breaks = classBreaks) %>%
#     setNames(as.character(format(dates, "%d-%b-%Y"))) 
#   
#   rasterCut[] <- (rasterCut[] - 1) / 10 # offset classes to match the values
#   
#   # count the non-ocean (non-NA) pixels for pixel proportion calculations
#   refLayer <- subset(rasterBrick, 1)
#   numLandPixels <- length(refLayer[!is.na(refLayer)])
#   
#   # format time series for plotting
#   layersDf <- mclapply(seq(raster::nlayers(rasterCut)), function(x) {
#     lyr <- subset(rasterCut, x)
#     data.frame(
#       date = lubridate::parse_date_time(gsub("^X", "", names(lyr)), "dbY"),
#       classVal = lyr[]
#     )
#   }) %>%
#     bind_rows() %>%
#     filter(!is.na(classVal)) %>%
#     group_by(date, classVal) %>%
#     summarise(numPix = length(classVal)) %>%
#     arrange(date, classVal)
#   
#   # for each date, each class must have a value for the number of pixels
#   allDates <- lubridate::parse_date_time(gsub("^X", "", names(rasterCut)), "dbY")
#   fitDf <- data.frame(date = sort(rep(allDates, length(classes[-1]))), # [-1] corrects classes offset classes
#                       classVal = rep(classes[-1], nlayers(rasterCut)), 
#                       numPix   = 0)
#   fitLayers <- bind_rows(layersDf, fitDf) %>%
#     group_by(date, classVal) %>%
#     summarise(numPix = sum(numPix)) %>%
#     ungroup() %>%
#     mutate(date = as.Date(date),
#            `Weekly GI` = factor(classVal),
#            proPix = numPix / numLandPixels) 
# }
# 
# 
# 
# # set classes
# classes <- seq(-1, 10, by = 1) / 10
# 
# # --------------------------------------------
# # stacked area time series for Eyre Peninsula
# # --------------------------------------------
# 
# pixelClassDfEyre <- rasterToPixelClassDataframe(eyreWGIBrick, simDate, classes)
# nClasses <- pixelClassDfEyre$classVal %>% unique() %>% length()
# 
# pixelClassDfEyre <- rasterToPixelClassDataframe(epBrickWGI, simDate, classes) 
# # create an object with the dates and labels I need
# # We want to find the breaks for the first day of each month
# dateSeq <- seq(as_date("2014-01-01"), as_date("2014-12-01"), by = "month") # use 2014 as a model sequence
# dateBreaks <- dateSeq %>% yday()# plus one day so that the Jan label appears on the plot.
# dateLabels <- dateSeq %>% format("%b") # drop year - we were only using 2014 as a model
# 
# # We also want to add vertical lines to separate the seasons
# seasonSeq <- seq(as_date("2013-12-01"), as_date("2014-12-01"), by = "3 months")
# seasonBreaks <- seasonSeq %>% yday()
# 
# # add the cropping and non-cropping breaks
# canolaWindow <- c("2014-05-01", "2014-11-01")
# canolaBreaks <- yday(canolaWindow)
# 
# # function to adjust (add 1) day of year for 2016 dates only
# ydayAdjust <- function(x, adjustYear = 2016, adjustment = 1){
#   y  <- year(x)
#   yd <- yday(x) 
#   yd[which(y == adjustYear)] <- yd[which(y == adjustYear)] + adjustment
#   yd
# }
# 
# pixelClassDfEyre %>% 
#   mutate(month = month(date),
#          year  = year(date),
#          yearLagged = ifelse(month == 12, year + 1, year),
#          dayOfYear  = ydayAdjust(date)) %>% # I need to facet by season (starting summer, Dec), not year
#   ggplot(aes(x = dayOfYear, y = proPix, fill = `Weekly GI`, group = `Weekly GI`)) +
#   geom_area() +
#   coord_polar() +
#   theme_light() +
#   #theme_dark() +
#   theme(panel.grid.minor.x = element_blank()) +
#   labs(x = NULL, y = "Proportion of pixels") +
#   geom_vline(xintercept = canolaBreaks, 
#              linetype = 3, colour = alpha("black", 0.7), size = 0.5) +
#   scale_fill_manual(values = c(alpha("grey", 0.5), myRamp(nClasses - 1))) +
#   scale_x_continuous(breaks = dateBreaks, labels = dateLabels,
#                      limits = c(0, 366)
#                      ) +
#   facet_wrap(~yearLagged) # using yearLagged mean Dec is from the year before
# ggsave("climexStackedPolarEyreLightDraft.pdf")
# # TO DO:
# # YOu need to add interpolate the data for the last week of dec so the background plots correctly
# # YOu might need to run this for Dec of the same year, not the year before (confusing)
# # ... which means you need to read in a month extra data for end of 2016.
# # Because 2016 is our by one day, subtract the day from 2016 data.
# 
# 
# # NEXT: (27/8/2018)
# # (1) Resolve the gap in January 2016
# # (2) REsolve whether to include data from current year in circles, or previous year (seasons)
# # (3) Try continuous ramp
# # 
# # # plot a stacked area time series for Australia
# # pixelClassDfAus <- rasterToPixelClassDataframe(ausWGIBrick, simDate, classes)
# # nClasses <- length(unique(pixelClassDfAus$classVal))
# # setwd(outAus)
# # pdf("climexWGIStackedAreaAustralia.pdf", height = 7, width = 10)
# # ggplot(pixelClassDfAus, 
# #        aes(x = date, y = proPix, fill = `Weekly GI`, group = `Weekly GI`)) +
# #   geom_area() +
# #   theme_bw() +
# #   theme(axis.text.x  = element_text(angle = 60, hjust = 1),
# #         axis.title.x = element_blank(),
# #         axis.title.y = element_text(face = "bold"),
# #         panel.grid   = element_blank(),
# #         axis.text    = element_text(colour = 'black'),
# #         legend.title = element_text(face = "bold"),
# #         plot.title = element_text(hjust = 0.5)) +
# #   ggtitle("CLIMEX Weekly Growth Index for P. xylostella across Australia over time") +
# #   guides(fill = guide_legend(reverse = TRUE)) +
# #   ylab("Proportion of raster pixels in each Weekly GI class") +
# #   scale_x_date(breaks = date_breaks(width = "1 month"),
# #                labels = date_format("%b-%y")) +
# #   scale_fill_manual(values = myRamp(nClasses))
# # dev.off()
# # setwd(old)
# 
# # # plot a stacked area time series for the Eyre Peninsula region
# # pixelClassDfEyre <- rasterToPixelClassDataframe(eyreWGIBrick, simDate, classes)
# # nClasses <- length(unique(pixelClassDfEyre$classVal))
# # setwd(outEyre)
# # pdf("climexWGIStackedAreaEyre.pdf", height = 7, width = 10)
# # ggplot(pixelClassDfEyre, 
# #        aes(x = date, y = proPix, fill = `Weekly GI`, group = `Weekly GI`)) +
# #   geom_area() +
# #   theme_bw() +
# #   theme(axis.text.x  = element_text(angle = 60, hjust = 1),
# #         axis.title.x = element_blank(),
# #         axis.title.y = element_text(face = "bold"),
# #         panel.grid   = element_blank(),
# #         axis.text    = element_text(colour = 'black'),
# #         legend.title = element_text(face = "bold"),
# #         plot.title = element_text(hjust = 0.5)) +
# #   ggtitle("CLIMEX Weekly Growth Index for P. xylostella across Eyre Peninsula over time") +
# #   guides(fill = guide_legend(reverse = TRUE)) +
# #   ylab("Proportion of raster pixels in each Weekly GI class") +
# #   scale_x_date(breaks = date_breaks(width = "1 month"),
# #                labels = date_format("%b-%y")) +
# #   scale_fill_manual(values = myRamp(nClasses))
# # dev.off()
# # setwd(old)


