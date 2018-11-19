# #############################################################################
# Make a collections table and geographic map for the RAD2 P. xylostella study 
# KD Perry, 12/9/2017
#
# Notes: 
# This script uses the workspace in `genotypes.Rproj`
#
# This script does:
# (i)   makes population summary of collection details for P. xylostella
# ... ... and exports LaTex table code
# (ii)  makes a geographic map of the Px collection details
#
# #############################################################################

library(tidyverse)
library(readxl)
library(xtable)
library(raster)
library(ggmap)
library(ggplot2)


# ##################################################
# Define functions
# ##################################################

# function to rename the 13 mislabelled Paus samples (species changed after COI sequencing)
renamePxToPaus <- function(x) {
  gsub("gilN14swx-01f.26", "gilN14swa-01f.26", x) %>%
    gsub("gilN14swx-06m.26", "gilN14swa-06m.26", .) %>%
    gsub("gilN14swx-07m.26", "gilN14swa-07m.26", .) %>%
    gsub("gilN14swx-08m.33", "gilN14swa-08m.33", .) %>%
    gsub("gilN14swx-02f.33", "gilN14swa-02f.33", .) %>%
    gsub("gilN14swx-03f.33", "gilN14swa-03f.33", .) %>%
    gsub("gilN14swx-04f.40", "gilN14swa-04f.40", .) %>%
    gsub("gilN14swx-10m.40", "gilN14swa-10m.40", .) %>%
    gsub("boyW14scx-05m.03", "boyW14sca-05m.03", .) %>%
    gsub("espW14scx-01f.04", "espW14sca-01f.04", .) %>%
    gsub("espW14scx-06m.11", "espW14sca-06m.11", .) %>%
    gsub("espW14scx-07m.11", "espW14sca-07m.11", .) %>%
    gsub("espW14scx-04m.04", "espW14sca-04m.04", .)
}

pasteAngle <- function(x, dig = 2, abs = FALSE, # for converting latitude from negative N to positive S
                         direction = c("N, S, E, W")) {
  if(abs){
    x <- abs(x)
  }
  format(round(x, dig), nsmall = dig) %>%
    paste0("\\ang{", ., "}", direction)
}

wrapParenth <- function(x){
  paste0("{(}", x, "{)}")
}
getSex <- function(x){
  sexes <- c("^m$" = "Male", "^f$" = "Female", "^u$" = "Unknown")
  x %>%
    str_sub(13, 13) %>%
    str_replace_all(sexes)
}
naToZero <- function(x){
  x[is.na(x)] <- 0
  x
}

# ##################################################
# Collections table
# ##################################################

# list individuals names used in the P. xylostella study (written in R script `writeSampleNames.R`)
px833ind <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "px833Ind.txt") %>%
  readLines()

popsMasterRAD2 <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "popsMasterRAD2.csv") %>%
  read_csv() 

# create population summary table
pxCollections <- data_frame(
  ind = renamePxToPaus(px833ind),
  popString = str_sub(ind, 1, 9),
  sex = getSex(ind)
  ) %>%
  # calculate the number of samples of each sex in each pop
  group_by(popString, sex) %>%
  summarise(n = length(popString)) %>%
  reshape2::dcast(popString ~ sex, value.var = "n") %>%
  map_dfr(naToZero) %>%
  rowwise() %>%
  mutate(totalInd = sum(Male, Female, Unknown)) %>%
  # join the population metadata
  left_join(popsMasterRAD2, by = "popString") %>%
  arrange(State, Location, collectionDate) %>%
  mutate(locationState  = paste(Location, State),
         collectionDate = format(collectionDate, "%b-%Y"))

# ----------------------------------------------------------------
# Make a quick summary of number of pops by host for Px manuscript
pxCollections %>% 
  group_by(hostType) %>%
  count()
# # Groups:   hostType [4]
# hostType     n
# <chr> <int>
#   1     Canola    29
# 2     Forage     3
# 3 Vegetables    15
# 4       Wild    12
# now group by season
pxCollections %>% 
  group_by(Season) %>%
  count()
# ----------------------------------------------------------------

# ---- check the sex ratio ----
pxCollections %>%
  summarise(males = sum(Male),
            females = sum(Female),
            unknown = sum(Unknown))
# A tibble: 1 x 3
# males females unknown
# <dbl>   <dbl>   <dbl>
#   1   364     317     152


# format for Latex
pxCollectionsLatex <- pxCollections %>%
  mutate(Latitude  = pasteAngle(Latitude,  abs = TRUE, direction = "S"),
         Longitude = pasteAngle(Longitude, abs = TRUE, direction = "E"),
         latLong   = paste(Latitude, Longitude, sep = "  ")) %>%
  dplyr::select(`Location` = locationState,
                `Collection date` = collectionDate,
                latLong, 
                Host = host,
                `{Total}` = totalInd,
                `{\\male{ }}` = Male,
                `{\\female{ }}` = Female)


pxAlign <- c("l",
             "l", # location
             "l", # collection date
             "l", # lat/lon,
             "l\n", # host
             "S[table-number-alignment=center,table-text-alignment=center,table-figures-integer=2,table-figures-decimal=0]\n", # Females
             "S[table-number-alignment=center,table-text-alignment=center,table-figures-integer=2,table-figures-decimal=0]\n", # Males
             "S[table-number-alignment=center,table-text-alignment=center,table-figures-integer=2,table-figures-decimal=0]\n") # Total

pxDigits <- c(0, 0, 0, 0, 0, 0, 0, 0)

# Add header rows
pxAddRows <- list()
pxAddRows$pos <- list(0, 0)
pxAddRows$command <- c(
  "Location & Collection & Coordinates & Host plant & \\multicolumn{3}{c}{No. sequenced} \\\\\n \\cmidrule(lr){5-7}\n",
  "& date & & & Total & {\\male} & {\\female} \\\\\n"
  )
pxCaption <- "Summary of \\textit{P. xylostella} collections from Australia."
pxCollectionsLatex %>%
  xtable(align = pxAlign,
         digits = pxDigits,
         caption = pxCaption,
         lab = "tab:collectionsPx") %>% 
  print.xtable(include.rownames = FALSE,
               include.colnames = FALSE,
               table.placement = "p",
               caption.placement = "top",
               add.to.row = pxAddRows,
               NA.string = "",
               booktabs = TRUE,
               sanitize.text.function = function(x){x})

  
# ===============================================================
# create a geographic map of collection locations or publication
# ===============================================================
# 
# use ggplot2 to allow aesthetics for year/host

library(broom)
library(ggsn) # for adding scalebar to ggplot

# read in an Australia polygon
old <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "genotypes")
ausPolyWd  <- file.path("C:", "UserData", "Kym", "PhD", "Data", "geoSpatialData", "ausAdminAreas_rds")
setwd(ausPolyWd)
ausPoly <- readRDS("AUS_adm1.rds") # coast and states
setwd(old)

# convert SpatialPointsDataFrame object to dataframe for plotting in ggplot2
# # extent coords for the Australian coast, taken from rasters
# 112.945, 153.595, -43.59501, -10.10499
ausPolygon <- broom::tidy(ausPoly) %>% # convert to standard data.frame for ggplot (replaces fortify())
  filter(lat > -43.59501, lat < -10.10499,
         long > 112.945, long < 153.595)

pxCollections <- pxCollections %>%
  mutate(Year = factor(Year),
         legendOrder = gsub("Canola", 1, Host) %>%
           gsub("Vegetables", 2, .) %>%
           gsub("Wild", 3, .) %>%
           gsub("Forage", 4, . ) %>%
           as.factor(),
         Host = gsub("Forage", "Forage brassicas", Host) %>%
           gsub("Vegetables", "Vegetable brassicas", .) %>%
           gsub("Wild", "Wild brassicas", .),
         Host  = factor(Host)) %>%
  # now set legend plot order
  arrange(legendOrder) %>%
  mutate(plotOrder = 1:nrow(.))

# reorder for plotting
pxCollections$Host <- reorder(pxCollections$Host, pxCollections$plotOrder)
  
# plot. to do: 
# ... check whether projection is best (coord_map(projection = "..."))
# ... add states text to centroid for each polygon (maybe inkscape)
# ... add place names/codes? or just esperance / southend?
# ... add capital city place names?
# ... add shaded regions for Vegetable regions?
# ... add shapefile shading for wheat/canola production regions

pxCollections %>%
  ggplot() +
  geom_polygon(data = ausPolygon, aes(x = long, y = lat, group = group),
              fill = "lightgrey", colour = "white",
              lwd = 0.3) +
  geom_point(data = filter(pxCollections, Year == 2014),
             aes(x = Longitude, y = Latitude, shape = Host),
             fill = 'black', size = 4) +
  geom_point(data = filter(pxCollections, Year == 2015),
             aes(x = Longitude, y = Latitude, shape = Host),
             fill = 'white', size = 2.75, show.legend = FALSE) +
  scale_shape_manual(values = c(21, 23, 24, 22),
                     guide  = guide_legend(title = "Host type")) +
  ggsn::blank() +
  theme(legend.title = element_text(face = "bold")) +
  ggsn::scalebar(x.min = 112.945, x.max = 153.595, 
                 y.min = -43.59501, y.max = -10.10499,
                 dist = 500, dd2km = TRUE, model = 'WGS84', 
                 height = 0.01, st.size = 3,
                 location = 'bottomleft',
                 anchor = c(x = 115, y = -40)) +
  coord_map()

setwd("C:/UserData/Kym/PhD/thesis/images")
ggsave(filename = "collectionsMapPx.pdf", height = 8, width = 10)
setwd(old)

stop please

# ==============================================
# Make a raster map of P. xylostella collections
# ==============================================
# 6/5/2018

library(ggmap) # you need at least version 2.7
library(magrittr)
library(RgoogleMaps)
#citation("ggmap")

# =========================
# grab a stamen map terrain
# =========================

#ausBoundingBox <- c(112.945, -43.59501, 153.595, -10.10499)
ausBoundingBox <- c(110.945, -43.59501, 155.595, -10.10499) # add some width around plot
ausMap <- get_map(ausBoundingBox, 
                    maptype = "terrain-background", # uses stamen
                    zoom = 7) # maps 

# draft code ------
library(rgeos)
trueCentroids  <- rgeos::gCentroid(ausPoly, byid = TRUE) # https://gis.stackexchange.com/questions/43543/how-to-calculate-polygon-centroids-in-r-for-non-contiguous-shapes?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
stateCentroids <- coordinates(trueCentroids) %>% 
  as_data_frame()
# get coordinates of Australian captical cities
ausCities <- c(
  "Adelaide, South Australia", "Canberra, Australian Capital Territory",
  "Melbourne, Victoria", "Perth, Western Australia", "Sydney, New South Wales",
  "Brisbane, Queensland", "Darwin, Northern Territory", "Hobart, Tasmania"
  )
# ausCityCoords <- geocode(
#   location = ausCities, output = "latlon"
# )
ausCityCoords <- ausCityCoords %>%
  mutate(labels = c("Adelaide", "Canberra", "Melbourne", "Perth",
                    "Sydney", "Brisbane", "Darwin", "Hobart"))
# draft code -----


ausStates <- geocode(
  location = c(
    "South Australia", "Australian Capital Territory", "Victoria, Australia",
    "Western Australia", "New South Wales", "Queensland, Australia", 
    "Northern Territory, Australia", "Tasmania, Australia"
    ))
ausStates<- ausStates %>%
  mutate(states = c(
    "SOUTH AUSTRALIA", "AUSTRALIAN CAPITAL TERRITORY", "VICTORIA",
    "WESTERN AUSTRALIA", "NEW SOUTH WALES", "QUEENSLAND",
    "NORTHERN TERRITORY", "TASMANIA"
    ))

yBreaks <- seq(-40, -20, by = 10)
yLabels <- parse(text = paste(yBreaks, "*degree ~ N", sep = ""))
xBreaks <- seq(120, 150, by = 10)
xLabels <- parse(text = paste(xBreaks, "*degree ~ W", sep = ""))

# # =================================
# # read land use `cropping` layer in
# # ==================================
# r <- file.path("C:", "UserData", "Kym", "PhD", "Data", "geoSpatialData",
#           "landUseDataAustralia", "luav5g9abll20160704a02egigeo___", "lu10v5ug", "vat.adf") %>% 
#   raster()
# 
# # get the attribute for cropping
# `3.0.0 Production from dryland agriculture and plantations`
# rLevsDf <- levels(r)[[1]]
# rLevsDf$LU_CODEV7N %>% unique() # find all the classification levels
# #> "3.3.1 Cereals" # I want this one
# # find the integer value corresponding the level I want
# z <- rLevsDf %>%
#   filter(LU_CODEV7N == 331)
# # integer values are "ID" I think. There are multiple IDs for 331.
# keepValues <- rLevsDf %>%
#   filter(LU_CODEV7N == 331) %>%
#   dplyr::select(ID) %>%
#   unlist()
# 
# # now make a new a raster layer with these cells showing a value, and others NA
# r[!r %in% keepValues] <- NA
# r[!is.na(r)] <- 1000
# 
# # we need to convert out raster to mercator projection to overlay with ggmap
# projection(r)
# library(rgdal)
# # get mercator projection
# EPSG <- make_EPSG() %>% data.frame()
# mercPrj <- EPSG %>%
#   filter(code == 3857) %>%
#   dplyr::select(prj4) %>% unlist %>% as.character
# projection(r) # "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
# projection(ausMap)
# plot(r, add = TRUE, colour = "black")
# # =================================================

ausMap %>%
   ggmap() +
# ggplot() +
  # geom_polygon(data = ausPolygon, aes(x = long, y = lat, group = group),
  #              fill = "transparent", colour = "lightgrey",
  #              lwd = 0.3) +
  # need to nudge these labels
  geom_text(data = ausStates, 
            aes(x = lon, y = lat, label = states), size = 3) +
  geom_text(data = ausCityCoords,
            aes(x = lon, y = lat, label = labels), size = 3) +
  geom_point(data = ausCityCoords,
             aes(x = lon, y = lat), 
             pch = 22, colour = "black", fill = "white", size = 4) +
  # geom_point(data = stateCentroids, 
  #            aes(x = x, y = y),
  #            size = 5, colour = "red") +
  geom_point(data = filter(pxCollections, Year == 2014),
             aes(x = Longitude, y = Latitude, fill = Host),
             pch = 21, colour = "transparent",
             #fill = 'black', 
             size = 5) +
  geom_point(data = filter(pxCollections, Year == 2015),
             aes(x = Longitude, y = Latitude, fill = Host),
             pch = 24, colour = "white", show.legend = FALSE,
             #fill = 'white', 
             size = 3) +
  scale_fill_manual(values   = c("black", "blue", "red", "yellow")) +
  # scale_y_continuous(breaks = yBreaks, labels = yLabels, expand = c(0, 0)) +
  # scale_x_continuous(breaks = xBreaks, labels = xLabels, expand = c(0, 0)) +
  #ggsn::blank() +
  theme_bw(base_size = 14) +
  theme(legend.position = c(0.15, y = 0.17),
        legend.title    = element_text(face = "bold")) +
  labs(x = "Longitude", y = "Latitude", fill = "Host") +
  ggsn::scalebar(x.min = 112.945, x.max = 153.595, y.min = -43.59501, y.max = -10.10499,
                 dist = 500, dd2km = TRUE, model = 'WGS84',
                 height = 0.01, st.size = 3,
                 location = 'bottomleft', anchor = c(x = 115, y = -42))

  # geom_polygon(data  = rtp, aes(x = long, y = lat), fill = "black") 
  # 
  # inset_raster(raster = as.raster(r), 
  #              xmin = attributes(ausMap)$bb$ll.lon,
  #              xmax = attributes(ausMap)$bb$ur.lon,
  #              ymin = attributes(ausMap)$bb$ll.lat,
  #              ymax = attributes(ausMap)$bb$ur.lat)
ggsave("PxCollectionsMapDraft.pdf", width = 12, height = 10)


# TO DO:
# Label states usng centroids
# circles around the vege production region
# shapefile or polygon arond the winter cropping zone
# point aesthetic scheme




# OLD NOTES BELOW
# https://gis.stackexchange.com/questions/3083/seeking-examples-of-beautiful-maps/45518#45518
# stamenAuckland <- get_map(location='Auckland', source="stamen", maptype="terrain", zoom=13)
# stamenAuckland %>% ggmap()
# stamenAdelaide <- get_map(location='Adelaide', source="stamen", maptype="terrain", zoom=10)
# stamenAdelaide %>% ggmap()

# get australia tiles (takes a long time)
# FULL ausStamen <- get_stamenmap(c(112.945, -43.59501, 153.595, -10.10499),
 #                          source = "stamen", maptype = "terrain")
#> crashes due to RAM.
#> https://stackoverflow.com/questions/22433245/make-large-ggmap-map-without-running-out-of-ram

# ==========================
# try some google maps tiles
# ==========================

# loc <- c(mean(pxCollections$Longitude), mean(pxCollections$Latitude))#
# plotmap(lon = pxCollections$Longitude,
#         lat = pxCollections$Latitude, API = "OSM", zoom = 3)

# bBox <- matrix(c(112.945, -43.59501, 153.595, -10.10499), nrow = 2) %>%
#   set_colnames(c("min", "max")) %>%
#   set_rownames(c("long", "lat"))
# 
# ggMap <- get_map(c(112.945, -43.59501, 153.595, -10.10499), source = "stamen") # zoom = 3
# 
# get_map("Australia", maptype = "roadmap", zoom = 3) %>% 
#   ggmap() + 
#   coord_map(xlim = c(112.945, 153.595), ylim = c(-43.59501, -10.10499))
#   
# 
# # example from here:
# # http://www.molecularecologist.com/2012/09/making-maps-with-r/
# library(maps)
# library(mapdata)
# map("worldHires", "Canada", xlim=c(-141,-53), ylim=c(40,85), col="gray90", fill=TRUE)
# map("worldHires", "Australia", xlim = c(112.945, 153.595), ylim = c(-43.59501, -10.10499), col="gray90", fill=TRUE)
# 
# # stamen maps
# stamenWater <- get_map(location='Auckland', source="stamen", maptype="watercolor", zoom=13)
# ggmap(stamenWater)
# 
# # stamen terain looks good, but can I get it ...
# 
# 
# ggmap(ggMap, extent = "normal", maprange = FALSE) +
#   geom_point(data = pxCollections, 
#              aes(Longitude, Latitude, colour = Year, shape = Year)) =
# 
# 
#   coord_map(xlim = c(112.945, 153.595), ylim = c(-43.59501, -10.10499))
#   
#   geom_text(x = loc[1], y = -31.29, 
#             aes(label = "Flinders Ranges NP"),
#             colour = "grey10", size = 4) +
#   geom_point(data= pcaForPlot, aes(Longitude, Latitude, colour = Population),
#              size = 2.5) +
#   scale_colour_manual(values = c(rgb(0, 0, 0, 0.7),
#                                  rgb(0.8, 0.2, 0.2, 0.9),
#                                  rgb(0.2, 0.8, 0.2, 0.9),
#                                  rgb(0.5, 0.5, 1, 0.7))) +
# 
#   guides(colour = FALSE) +
#   labs(x = "Longitude",
#        y = "Latitude") +
#   theme_bw() +
#   theme(text = element_text(size = 11))```
# 
# ```{r fullAusMap, results='hide'}
# ausPolygon <- file.path("..", "60804_shp", "framework") %>%
#   readOGR(layer = "aus25fgd_r") %>%
#   subset(AREA > 6)
# ausPts <- SpatialPoints(cbind(x = loc[1], y = loc[2]))
# proj4string(ausPts) <- proj4string(ausPolygon)
# fullPlot <- ggplot()+
#   geom_polygon(data=ausPolygon, 
#                aes(long, lat, group = group), fill = "white", colour = "black") + 
#   geom_point(data = as.data.frame(ausPts), aes(x, y), size = 3) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())```
# 
# 
# ```{r plotWithInset, echo = FALSE, fig.cap="Collection points for all 2012 samples with colours showing sub-populations initially defined by PCA analysis and *k*-means clustering.", fig.height=7, fig.width=8}
# grid.newpage()
# vp1 <- viewport(0.5,0.5, 1, 1)
# vp2 <- viewport(0.25, 0.28, 0.3, 0.3)
# print(zoomPlot, vp = vp1)
# print(fullPlot, vp = vp2)```




