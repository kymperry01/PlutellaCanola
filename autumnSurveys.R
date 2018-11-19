# ####################################################################################
# Plot autumn survey data and make tables and figures for publication
# K. Perry, 6/3/2017
#
# This script does:
# 1. Maps autumn survey count data for moths and larvae on geographic maps using proportional symbols
# 2. Produces data summaries for each year, and a table for publication

# Datasets include 3 years, 2014, 2015, 2016, with two time points in March and April 
# #####################################################################################

library(readxl) 
library(tidyverse)
library(reshape2)
library(rgdal)
library(lattice)
library(latticeExtra)
library(ggplot2)
library(RColorBrewer)
library(xtable)


# ###########
# functions #
# ###########

# standard error
se <- function(x, na.rm = TRUE){
  if(na.rm) x <- x[!is.na(x)]
  sd(x) / sqrt(length(x))
}


# ##############################################
# Tidy autumn survey data for 2014, 2015, 2016 #
# ##############################################

# read in the autumn survey data for all three years
wdir  <- file.path("C:", "UserData", "Kym", "PhD", "Data", "autumnSurveys")
pathData2014 <- file.path(wdir, "autumnSurveys2014", "autumnSurveyData2014.xlsx")
pathData2015 <- file.path(wdir, "autumnSurveys2015", "autumnSurveyData2015.xlsx")
pathData2016 <- file.path(wdir, "autumnSurveys2016", "autumnSurveyData2016.xlsx")

# Read in the data and collapse to single row per site per year
# In 2014 we ran multiple traps at some sites. Take a mean of moth and larval counts.
meanNumeric <- function(x){
  if(is.numeric(x)){
    x <- mean(x, na.rm = TRUE)
  }
  unique(x) # output a value of length 1.
}

surveyDataRaw <- c(pathData2014, pathData2015, pathData2016) %>%
  map_dfr(read_excel) %>%
  group_by(year, region, consensusSiteID, siteNum) %>%
  # take a per-site average of trap counts for 2014, by just averaging all numeric variables (some sites had two traps)
  summarise_all(., .fun = meanNumeric) %>% 
  # create a unique site identifier within each year
  mutate(siteID = str_pad(siteNum, side = "left", width = 2, pad = 0),
         siteID = paste(siteID, year, sep = "-")) %>%
  dplyr::select(siteID, everything()) # move it to the first column

# ===========================================
# convert the time series data to long format
# ===========================================

# Raw data are in an inconvenient format for melting: not all the variables are associated with each other.   
# There are:
# -- multiple date columns associated with time0 (March) time1 (April) sampling
# -- separate variables associated with each time point.
# Therefore, we will (1) split out the time points and associated variables, (2) melt, (3) recombine

# Create a character variable of ID variables common across all time points
# We'll drop the vegetation assessment variables (can be joined later as required)
idVars <- c("sampled2014", "sampled2015", "sampled2016",
            "siteID", "consensusSiteID", "region", "location", "lat", "long", "host")


# Split data for time point 0 (March) and time point 1 (April)
time0 <- surveyDataRaw %>%
  melt(
    id.vars = c(idVars, "date0"), # grab the dates for time point 0
    measure.vars = c("larvae0"),
    variable.name = "lifestage",
    value.name = "count"
  ) %>%
  dplyr::rename(date = date0)

time1 <- surveyDataRaw %>%
  melt(
    id.vars = c(idVars, "date1"), # grab the dates for time point 1 
    measure.vars = c("larvae1", "moths1"),
    variable.name = "lifestage",
    value.name = "count"
    ) %>%
  dplyr::rename(date = date1)

# Re-combine data in long format
surveyDataTidy <- time0 %>%
  bind_rows(time1) %>%
  mutate(timePoint = str_replace_all(lifestage, "[a-z]", ""),
         month = format(date, "%B"),
         year  = format(date, "%Y"),
         lifestage = str_replace_all(lifestage, "\\d$", "")) %>%
  arrange(siteID)

# ################################################################################
# Publication table: P. xylostella incidence and abundance by brassicaceous host #
# ################################################################################

# ==========================================================================
# Summarise the indicence of Plutella at sampled sites among years and hosts
# ==========================================================================

# For each year (taking both months together), calculate: 
# (a) how many sites were sampled for both moths and larvae 
# note: `sea rocket` sites were not sampled for moths, and `nil` host sites were not sampled for larvae.
# (b) how many sites were positive for moths, larvae or either stage 

# First calculate sites how many sites positive for any life stage
# ie: Plutella moths/larvae present in March/April)
sitesPositiveAnyStage <- surveyDataTidy %>%
  group_by(year, siteID, host, lifestage) %>%
  summarise(nTimesSampled = length(which(!is.na(count))),
            sumCount = sum(count, na.rm = TRUE)) %>% 
  filter(nTimesSampled > 0) %>% # remove rows for lifestages not sampled at a site
  group_by(year, siteID, host) %>% 
  summarise(sumCountAnyStage = sum(sumCount, na.rm = TRUE)) %>% 
  group_by(year, host) %>%
  summarise(nSitesSampled = length(unique(siteID)),
            nSitesPositiveAnystage = length(which(sumCountAnyStage > 0)),
            pSitesPositiveAnyStage = round(nSitesPositiveAnystage / nSitesSampled, 2))


# Now calculate how many sites positive for moths or larvae, and join `anystage` data
# Group in 2 steps because these stats required different grouping for the calculations
sitesSampled <- surveyDataTidy %>%
  group_by(year, siteID, host, lifestage) %>%
  summarise(nTimesSampled = length(which(!is.na(count))),
            sumCount = sum(count, na.rm = TRUE)) %>% 
  filter(nTimesSampled > 0) %>%  # remove rows for lifestages not sampled at a site
  group_by(year, host) %>%
  summarise(nSitesSampled = length(unique(siteID)),
            nSitesSampledMoths = length(which(lifestage  == "moths")),
            nSitesSampledLarvae = length(which(lifestage == "larvae")),
            nSitesPositiveMoths = length(
              intersect(which(lifestage == "moths"), which(sumCount > 0))
              ),
            pSitesPositiveMoths = round(nSitesPositiveMoths / nSitesSampledMoths, 2),
            nSitesPositiveLarvae = length(
              intersect(which(lifestage == "larvae"), which(sumCount > 0))
              ),
            pSitesPositiveLarvae = round(nSitesPositiveLarvae / nSitesSampledLarvae, 2)) %>%
  full_join(sitesPositiveAnyStage)

# ==========================================================================
# Summarise the mean abundance of Plutella at sampled sites by year and host
# ==========================================================================

# CHANGE HERE (abundance just at positive sites)

# Note:
# Calculate abundance across all sites positive or negative for Plutella
# at each sampling time point
sitesPlutellaAbundance <- surveyDataTidy %>%
  # filter(year == 2014, host == "sea rocket", lifestage == "larvae", month == "April") %>%
  # filter(count > 0) %>%
  group_by(year, month, host, lifestage) %>%
  summarise(nSites = length(count),
            nSitesPositive = length(count[which(count > 0)]),
            # mnCount = mean(count, na.rm = TRUE), # across all sites. Not useful.
            # seCount =   se(count, na.rm = TRUE),
            mnCountPositive = mean(count[count > 0], na.rm = TRUE), # mean across positive sites only (DBM present)
            seCountPositive =   se(count[count > 0], na.rm = TRUE)  # se as above
            )

# Convert to a single row per site. 
# Impossible to dcast() everything in one step, so split moths/larvae, then join
mothAbundance <- sitesPlutellaAbundance %>%
  filter(lifestage == "moths") %>%
  dplyr::rename(mnCountMoths = mnCountPositive,
                seCountMoths = seCountPositive) %>%
  ungroup() %>%
  dplyr::select(-lifestage, -month, -nSites, -nSitesPositive) %>%
  arrange(host, year)
  

# Split larvae summary by month
larvaeAbundanceMarch <- sitesPlutellaAbundance %>%
  filter(lifestage == "larvae", month == "March") %>%
  dplyr::rename(mnCountLarvaeMarch = mnCountPositive,
                seCountLarvaeMarch = seCountPositive) %>%
  ungroup() %>%
  dplyr::select(-lifestage, -month, -nSites, -nSitesPositive)

larvaeAbundanceApril <- sitesPlutellaAbundance %>%
  filter(lifestage == "larvae", month == "April") %>%
  dplyr::rename(mnCountLarvaeApril = mnCountPositive,
                seCountLarvaeApril = seCountPositive) %>%
  ungroup() %>%
  dplyr::select(-lifestage, -month, -nSites, -nSitesPositive)


# Now join the full dataset
surveySummary <- sitesSampled %>%
  full_join(mothAbundance) %>%
  full_join(larvaeAbundanceMarch) %>%
  full_join(larvaeAbundanceApril)


# ====================================
# output a Latex table for publication
# ====================================

wrapParenth <- function(x){
  #x[!is.na(x)] <- paste0("{(}", sprintf("%.2f", x[!is.na(x)]), "{)}")
  x[!is.na(x)] <- paste0("{(}", format(x[!is.na(x)], nsmall = 2), "{)}")
  x
}

roundNumeric <- function(x, digits = 2){
    if(is.numeric(x)){
      x <- round(x, digits)
    }
  x
    #format(x, nsmall = 2) # output a summary value of length 1. Make this character format, to force it to retain two decimal places
  }

# a function to make all the minor hacks needed to the final table
cleanupLatex <- function(x){
  x[x == "NaN \\pm NA"] <- NA
  x[x == "7.00 \\pm NA"] <- "7.00"
  x[x == "1.00 \\pm NA"] <- "1.00"
  x[x == "2.00 \\pm NA"] <- "2.00"
  x[x == "NaN"] <- NA
  x
}

# for table display 
rowOrder <- c("buchan weed"  = "5", "sea rocket" = "2",
              "lincoln weed" = "1", "forage brassica" = "4",
              "nil" = "6", "vol\\. canola" = "3")

hostNames <- c("buchan weed"  = "Buchan weed, \\\\textit{Hirschfeldia incana} L., ", 
               "sea rocket" = "Sea rocket, \\\\textit{Cakile maritima} Scop.,",
               "lincoln weed" = "Lincoln weed, \\\\textit{Diplotaxis tenuifolia} L.DC", 
               "forage brassica" = "\\\\textit{Brassica} forage crops",
               "nil" = "Nil host", "vol\\. canola" = "Weedy canola, \\\\textit{Brassica napus}, L.")


dupToBlank <- function(x){
  x[duplicated(x)] <- ""
  x
}

surveySummaryLatex <- surveySummary %>%
  map_dfr(roundNumeric, 2) %>%
  mutate(
    rowOrder = str_replace_all(host, rowOrder), 
    host     = str_replace_all(host, hostNames),
    pSitesPositiveAnyStage = wrapParenth(pSitesPositiveAnyStage),
    nSitesWithLarvae = paste0(nSitesPositiveLarvae, "/", nSitesSampledLarvae),
    pSitesWithLarvae = wrapParenth(pSitesPositiveLarvae),
    nSitesWithMoths = paste0(nSitesPositiveMoths, "/", nSitesSampledMoths),
    pSitesWithMoths = wrapParenth(pSitesPositiveMoths),
    mnMothsSE       = paste(sprintf("%.2f", mnCountMoths), sprintf("%.2f", seCountMoths), sep = " \\pm "), # sprintf hack esure 2 decimals
    mnLarvaeMarchSE = paste(sprintf("%.2f", mnCountLarvaeMarch), sprintf("%.2f", seCountLarvaeMarch), sep = " \\pm "),
    mnLarvaeAprilSE = paste(sprintf("%.2f", mnCountLarvaeApril), sprintf("%.2f", seCountLarvaeApril), sep = " \\pm ")
  ) %>%
  arrange(rowOrder, host, year) %>%
  dplyr::select(-rowOrder) %>% # exclude cols and select cols must be in separate calls to dplyr::select
  dplyr::select(host, 
                year,
                nSitesSampled, nSitesPositiveAnystage, pSitesPositiveAnyStage,
                nSitesWithMoths, pSitesWithMoths, 
                nSitesWithLarvae, pSitesWithLarvae,
                mnLarvaeMarchSE, mnLarvaeAprilSE, mnMothsSE
                ) %>%
  mutate(host = dupToBlank(host)) %>%
  map_dfr(cleanupLatex) 


# output Latex code
colAlign <- c(
  "l",
  "l", # host 
  "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=4]\n", # year,
  "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2,round-mode=places,round-precision=0]\n", # sites sampled
  "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]\n", # n sites positive any stage
  "@{}S[table-align-text-pre=false,table-figures-integer=1,round-mode=places,round-precision=2]\n", # p sites positive any stage,
  "c", # n sites positive larvae
  "@{}S[table-align-text-pre=false,table-figures-integer=1,round-mode=places,round-precision=2]\n", # p sites positive larvae
  "c", # n sites positive moths
  "@{}S[table-align-text-pre=false,table-figures-integer=1,round-mode=places,round-precision=2]\n", # p sites positive moths
  "S[separate-uncertainty,table-figures-uncertainty=1,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=2,round-mode=places,round-precision=2]\n", # mean larvae march
  "S[separate-uncertainty,table-figures-uncertainty=1,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=2]\n", # mean larvae april
  "S[separate-uncertainty,table-figures-uncertainty=1,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=2]\n"  # mean moths
  )

colDigits <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

addHeader <- list()
addHeader$pos  <- list(0, 0)
addHeader$command <- c(
  paste(c(
  # topline
  " & & & \\multicolumn{6}{c}{\\textit{Plutella} incidence} & ",
  "\\multicolumn{3}{c}{\\textit{Plutella} abundance}\\\\\n",
  "\\cmidrule(lr){4-9}\n", "\\cmidrule(lr){10-12}\n"
  ), collapse = ""),
  # second line
  paste(c(
  "Host & ",
  "{Year} & ", 
  "{$N$} & ",
  "\\multicolumn{2}{c}{$N_{\\textsc{positive}}$} & ",
  "\\multicolumn{2}{c}{$N_{\\textsc{moths}}$} & ",
  "\\multicolumn{2}{c}{$N_{\\textsc{larvae}}$} & ",
  "{\\makecell{Mean $\\pm$ SEM\\\\larvae per site\\\\(March)}} & ",
  "{\\makecell{Mean $\\pm$ SEM\\\\larvae per site\\\\(April)}} & ",
  "{\\makecell{Mean $\\pm$ SEM\\\\moths per site\\\\(March-April)}}\\\\\n"
  ), collapse = "")
  )

captionSurveys <- c(
"The incidence and abundance of \\textit{Plutella} at sites sampled during field surveys 
conducted in South Australia during March and April in each year from 2014-2016. 
Presented are the numbers and proportion in parentheses of: 
Number of sites sampled ($N$), 
sites ($N_{\\textsc{positive}}$) where any lifestage was detected in either month, 
sites where moths ($N_{\\textsc{moths}}$) or larvae ($N_{\\textsc{larvae}}$) were detected in either month 
and proportion in parentheses based on sites sampled for the respective lifestages,
and the mean $\\pm$ standard error of the mean numbers of larvae collected from positive plants in each month, and moths trapped in 
pheromone traps during the trapping period at positive sites.
A dash indicates no sampling."
)


surveySummaryLatex %>% 
  xtable(align  = colAlign,
         digits = colDigits,
         caption = captionSurveys,
         lab = "tab:autumnSurveys") %>%
  print.xtable(include.rownames = FALSE,
               include.colnames = FALSE,
               table.placement = "p",
               caption.placement = "top",
               add.to.row = addHeader,
               NA.string = "\\textendash",
               #scalebox = 0.5,
               booktabs = TRUE,
               sanitize.text.function = function(x){x})



# ###################################################
# Create geographic bubble plots of the survey data #
# ###################################################

# ========================================
# download polygons for geographic context
# ========================================

geoSpatialDataPath <- file.path("C:", "UserData", "Kym", "PhD", "Data", "geoSpatialData")

# Large Australia polygon with states and local admin areas
ausPolygon <- file.path(geoSpatialDataPath, "ausAdminAreas_shp") %>%
  rgdal::readOGR(dsn = ., layer = "AUS_adm2")
ausPolygonFortified <- fortify(ausPolygon) # convert to a standard data.frame for ggplot2


# Much smaller Australia polygon with states (faster processing) 
ausPolygonSmall <- file.path(geoSpatialDataPath, "nsaasr9nnd_02211a04ec_alb132") %>%
  rgdal::readOGR(dsn = ., layer = "aust_cd66states")
ausPolygonSmallFortified <- fortify(ausPolygonSmall) # convert to df for ggplot


# IBRA regions and sub-regions polygons
ibraRegionsPolygon <- file.path(geoSpatialDataPath, "ibraBioRegions") %>%
  rgdal::readOGR(dsn = ., layer = "ibra7_regions")
ibraRegionsPolygonFortified <- fortify(ibraRegionsPolygon)

ibraSubRegionsPolygon <- file.path(geoSpatialDataPath, "ibraBioRegions") %>%
  rgdal::readOGR(dsn = ., layer = "ibra7_subregions")
ibraSubRegionsPolygonFortified <- fortify(ibraSubRegionsPolygon)

# Extract the Eyre Yorke block polygons (see pdf of the ibra sub regions)
eyre <- ibraSubRegionsPolygon[ibraSubRegionsPolygon$SUB_CODE_7 %in% c("EYB03", "EYB04", "EYB05"), ]
eyreFortified <- fortify(eyre)


# =============================
# plot all survey data together
# =============================

# Modify the survey data for ggplot
surveyDataTidyGGplot <- surveyDataTidy %>%
  mutate(group = 1, # for geom_polygon() so that aesthetics match the ggplot() call
         stageMonth = paste(lifestage, month, sep = "_"))

# Plot all data (slow. takes ~5 mins with the polygon included. turn that off for testing)
ausPolygonSmallFortified %>%
  ggplot(aes(x = long, y = lat, group = group)) +
  theme_bw() +
  geom_polygon(fill = "lightgrey", colour = "white") +
  geom_point(data = surveyDataTidyGGplot,
             aes(x = long, y = lat, colour = host, size = count,
                 group = group)) +
  # subset using coord_map rather than filter: keeps the whole polygon which you need for fill.
  coord_map(xlim = c(134, 141), ylim = c(-38, -32)) +
  facet_grid(year~stageMonth)
ggsave("autumnSurveyMapFullDRAFT.png")

# To DO:
# plot moth and larvae data separately (the counts are different metrics)
# decide on bubbledot scheme:
# ... set the colour scheme.
# degrees on plot axes
# make the labels and host names nice

# Add a scalebar, north arrow, and an australia inset map.
# (for international audience)

# ===================================================================
# plot the Eyre Peninsula data only (possibility to use as an inset)
# ===================================================================

# we will subset the point data to include just Eyre Peninsula data,
# so that fill/colour aesthetics are mapped to just the hosts we want
eyreFortified %>%
  ggplot(aes(x = long, y = lat, group = group)) +
  theme_bw() +
  geom_polygon(fill = "lightgrey", colour = "white") +
  geom_point(data = surveyDataTidyGGplot %>%
               filter(region == "Eyre Peninsula"),
             aes(x = long, y = lat, colour = host, fill = host,
                 size = count, group = group)) + # adding `pch = 21` outside the aes() will geive you legend with open circles`
  # plot the zero counts as open circles with white fill to overwrite the filled circles
  geom_point(data = surveyDataTidyGGplot %>%
               filter(region == "Eyre Peninsula",
                      count  == 0),
             aes(x = long, y = lat, colour = host),
             pch = 21, fill = "white",
             show.legend = FALSE) +
  # subset using coord_map rather than filter: keeps the whole polygon which you need for fill.
  #coord_map(xlim = c(134, 141), ylim = c(-38, -32)) +
  coord_map() +
  facet_grid(year~stageMonth)
ggsave("autumnSurveyMapEyrePeninsulaDRAFT.pdf", height = 10, width = 15)


# ============================================================================
# Publication figure: Plot moth counts and larval averaged for March and April
#=============================================================================

# Average the larval counts across March and April each year
# This has the advantage of showing presence in either month
hostNames <- c("^lincoln weed$" = "Lincoln weed", "nil" = "Nil host",
               "^sea rocket$" = "Sea rocket", "vol\\. canola" = "Weedy canola",
               "^forage brassica$" = "Brassica forage crops",
               "^vol\\. canola$" = "Weedy canola",
               "^buchan weed" = "Buchan weed")

surveyDataTidyGrouped <- surveyDataTidy %>% 
  group_by(siteID, lat, long, region, host, lifestage, year) %>%
  summarise(nRows = length(siteID), # always good to know how many rows were grouped
            count = mean(count, na.rm = TRUE)) %>%  
  ungroup() %>%
  mutate(group = 1, # for geom_polygon() so that aesthetics match the ggplot() call
         #presence  = factor(ifelse(count > 0, 1, count)), #  variable to map open versus filled circles 
         lifestage = str_replace_all(lifestage, c("larvae" = "Larvae", "moths" = "Moths")),
         host      = str_replace_all(host, hostNames))# reorder for plotting
         

# I added blank labels to the breaks to avoid overlap on the facet plots
yBreaks <- seq(-35, -32, by = 1)
yLabels <- c("", parse(text = paste(abs(yBreaks[2:4]), "*degree ~ S", sep = ""))) # no label for extremes, to avoid overlap
xBreaks <- seq(134, 138, by = 1)
xLabels <- c("", parse(text = paste(xBreaks[2:4], "*degree ~ E", sep = "")), "") # now labels for extremes

hostColours <- c(rgb(0, 0, 0, 0.7),
                 rgb(0.8, 0.2, 0.2, 0.9),
                 rgb(0.2, 0.8, 0.2, 0.9),
                 rgb(0.5, 0.5, 1, 0.7))
pdf("autumnSurveyMapPublicationDRAFT.pdf", width = 9, height = 5)
surveyMaps <- ausPolygonSmallFortified %>%
  ggplot(aes(x = long, y = lat, group = group)) +
  theme_light() +
  theme(panel.grid = element_blank()
        #strip.background = element_blank()
        ) +
  labs(x = "Longitude", y = "Latitude", colour = "Host", size = "Count") +
  geom_polygon(fill = "lightgrey", colour = "white") +
  geom_point(data = surveyDataTidyGrouped %>%
               filter(region == "Eyre Peninsula", count != 0),
             aes(x = long, y = lat, colour = host, size = count, group = group),
             pch = 19, alpha = 0.7) +
  geom_point(data = surveyDataTidyGrouped %>%
               filter(region == "Eyre Peninsula", count == 0),
             aes(x = long, y = lat, colour = host),
             pch = 21, fill = "white", alpha = 0.7, show.legend = FALSE) +
  facet_grid(lifestage ~ year) +
  # subset using coord_map rather than filter: keeps the whole polygon which you need for fill.
  #coord_map(xlim = c(134, 141), ylim = c(-38, -32)) + # for all South Australia
  #coord_map(xlim = c(134.1, 137.9), ylim = c(-34.9, -31.9)) + # a hack to remove outer limits to avoid overlap of axis labels
  coord_map(xlim = c(134, 138), ylim = c(-35, -32)) +
  scale_size_continuous(# limits = c(0, 125), Use this if you want to add an empty circle for zero in the legend (have to be done manually using grid)
                        breaks = c(0, 1, 10, 25, 50, 125), # you need to manually make the zeros an open circle!
                        labels = c(0, 1, 10, 25, 50, 125)) +  
  scale_colour_manual(values = hostColours,
                      limits = c("Lincoln weed", "Sea rocket", "Weedy canola", "Nil host")) + # to reorder for plotting
  scale_x_continuous(breaks = xBreaks, labels = xLabels) +
  scale_y_continuous(breaks = yBreaks, labels = yLabels) +
  ggsn::scalebar(data = surveyDataTidyGrouped %>% 
                   filter(region == "Eyre Peninsula", count != 0),
                 x.min = 134, x.max = 138, y.min = -35, y.max = -32,
                 dist = 75, dd2km = TRUE, model = 'WGS84',
                 height = 0.05, st.size = 2.5, st.dist = 0.1, 
                 location = "topright", anchor = c(x = 137.6, y = -32.3),
                 facet.var = c("lifestage", "year"), # you need these facet.var, otherwise throws error about aes() 
                 facet.lev = c("Larvae", "2016"))
ggsn::north2(surveyMaps, x = 0.7825, y = 0.46, scale = 0.07, symbol = 1)
dev.off()
#ggsave("autumnSurveyMapPublicationDRAFT.pdf", width = 8, height = 6)


# ============================================================================
# Make a supplementary figure of the 2015 data including the South East region
# ============================================================================

# I added blank labels to the breaks to avoid overlap on the facet plots
yBreaks <- seq(-38, -32, by = 1)
yLabels <- parse(text = paste(abs(yBreaks), "*degree ~ S", sep = ""))
xBreaks <- seq(135, 141, by = 2)
xLabels <- parse(text = paste(xBreaks, "*degree ~ E", sep = ""))

pdf("autumnSurveyPublication2015SouthEastDRAFT.pdf", width = 9, height = 4)
surveyMaps2015 <- ausPolygonSmallFortified %>%
  ggplot(aes(x = long, y = lat, group = group)) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", colour = "Host", size = "Count") +
  geom_polygon(fill = "lightgrey", colour = "white") +
  geom_point(data = surveyDataTidyGrouped %>% filter(year == 2015, count != 0),
             aes(x = long, y = lat, colour = host, size = count, group = group),
             pch = 19, alpha = 0.7) +
  geom_point(data = surveyDataTidyGrouped %>% filter(year == 2015, count == 0),
             aes(x = long, y = lat, colour = host),
             pch = 21, fill = "white", alpha = 0.7) +
  facet_wrap(~lifestage) +
  # subset using coord_map rather than filter: keeps the whole polygon which you need for fill.
  coord_map(xlim = c(134, 141.5), ylim = c(-38, -32)) +
  #coord_map(xlim = c(134, 141), ylim = c(-38, -32)) + # for all South Australia
  #coord_map(xlim = c(134.1, 137.9), ylim = c(-34.9, -31.9)) + # a hack to remove outer limits to avoid overlap of axis labels
  #coord_map(xlim = c(134, 138), ylim = c(-35, -32)) +
  scale_size_continuous(# limits = c(0, 125), Use this if you want to add an empty circle for zero in the legend (have to be done manually using grid)
    breaks = c(0, 1, 10, 25, 50), # you need to manually make the zeros an open circle!
    labels = c(0, 1, 10, 25, 50)) +  
  scale_colour_manual(values = c(rgb(0, 0, 0, 0.7),
                                 rgb(0.8, 0.2, 0.2, 0.9),
                                 rgb(0.2, 0.8, 0.2, 0.9),
                                 "magenta",
                                 "goldenrod",
                                 rgb(0.5, 0.5, 1,   0.7)),
                      limits = c("Lincoln weed", "Sea rocket", "Weedy canola", "Buchan weed",
                                 "Brassica forage crops", "Nil host")) +
  scale_x_continuous(breaks = xBreaks, labels = xLabels) +
  scale_y_continuous(breaks = yBreaks, labels = yLabels) +
  ggsn::scalebar(data = surveyDataTidyGrouped %>% 
                   filter(year == 2015, count != 0),
                 x.min = 134, x.max = 138, y.min = -35, y.max = -32,
                 dist = 100, dd2km = TRUE, model = 'WGS84',
                 height = 0.02, st.size = 2.5, st.dist = 0.05, 
                 location = "bottomleft", 
                 anchor = c(x = 134.5, y = -37.5),
                 facet.var = c("lifestage"), # you need these facet.var, otherwise throws error about aes() 
                 facet.lev = c("Larvae"))
ggsn::north2(surveyMaps2015, x = 0.74, y = 0.83, scale = 0.08, symbol = 4)
dev.off()

# To Do:
# Possibly, the legend for counts zero size to an open circle (use grid)
# Plot the South Australian major roads polygon on top
# Add a separate plot with all South Australia as a supplementary figure so that South East is included
# change red / green colours.

# BEWARE BELOW. OLD, BAD, BAD, CODE.
# 
# 
# # ############################################################
# # plot larval count data on geographic map using ggplot method
# # ############################################################
# 
# lvData <- autumnData %>%
#   filter(stage == "larvae")
# 
# # make a spatial points dataframe (sp class) to generate a spatial extent (bbox)
# lvDataSP <- SpatialPointsDataFrame(
#   data = lvData,
#   coords = lvData[ , c("long", "lat")],
#   proj4string = CRS("+proj=longlat +ellps=WGS84")
#   )
# 
# # set symbology for larval counts based on manually set fixed classes
# lvBrks <- c(0, 0.01, 5, 20, 50, 100, 250) # breaks for larval counts
# nClasses <- 6
# intervals <- classIntervals(lvData$count, n = nClasses, style = 'fixed', fixedBreaks = lvBrks)
# tab <- print(intervals)
# nClasses <- length(intervals$brks) - 1
# # op <- options(digits = 4)
# # options(op)
# 
# # --------------------------------------------
# # To set circle size using Borden Dent scheme: 
# # ... (ref book Displaying TIme series data in R, O. Perpinan, pg 95)
# # complete Dent set of circle radii (mm)
# # dent <- c(0.64, 1.14, 1.65, 2.79, 4.32, 6.22, 9.65, 12.95, 15.11)
# # lvCex <- dent[seq(nClasses)]
# # -----------------------------
# 
# # manually set circle class sizes
# # link size and class: findCols returns the class number of each point; 
# # cex is the vector of sizes for each data point
# lvCex <- c(0.64, 1.14, 1.5, 2, 2.5, 3)
# idx <- findCols(intervals)
# 
# # manually set class colour scheme (heat colours)
# heatPal <- rev(heat.colors(nClasses))
# heatPal[1] <- "black" # plot zero counts in black 
# 
# # set classes
# lvDf <- data.frame(lvDataSP) # fortify to a standard df (ggplot does not handle SP classes)
# lvDf$class <- factor(names(tab)[idx])
# # format and reorder class for improved plotting display;
# lvDf$class2 <- as.character(lvDf$class) %>% # a hack to ensure correct class sorting and better legend display
#   gsub("\\[0,0.01\\)", "0", .) %>%
#   gsub("\\[0.01,5\\)", "001-5", .) %>%
#   gsub("\\[5,20\\)", "005-20", .) %>%
#   gsub("\\[20,50\\)", "020-50", .) %>%
#   gsub("\\[50,100\\)", "050-100", .) %>%
#   gsub("\\[100,250\\]", "100-250", .)
# 
# # to change the text of facet labels (for using facet_wrap by timePoint)
# levels(lvDf$monYear)[levels(lvDf$monYear) == "03-2014"] <- "March, 2014"
# levels(lvDf$monYear)[levels(lvDf$monYear) == "03-2015"] <- "March, 2015"
# levels(lvDf$monYear)[levels(lvDf$monYear) == "03-2016"] <- "March, 2016"
# levels(lvDf$monYear)[levels(lvDf$monYear) == "04-2014"] <- "April, 2014"
# levels(lvDf$monYear)[levels(lvDf$monYear) == "04-2015"] <- "April, 2015"
# levels(lvDf$monYear)[levels(lvDf$monYear) == "04-2016"] <- "April, 2016"
# 
# # to change the text of facet labels (for using facet_grid by month and year) 
# lvDf$month[lvDf$month == 3] <- "March"
# lvDf$month[lvDf$month == 4] <- "April"
# 
# # reorder March and April for correct plotting order
# lvDf$month <- factor(lvDf$month)
# lvDf <- lvDf %>%
#   arrange(desc(month)) %>%
#   mutate(monthOrder = 1:nrow(.)) 
# lvDf$month <- reorder(lvDf$month, lvDf$monthOrder)
# 
# # get base maps for spatial context using the ggmap::get_map() method
# # use different bounding boxes as needed
# # Note: you can manually define a bounding box for base tiles by trial and error and using crop = TRUE
# 
# bbEyre <- bbox(lvDataSP[lvDataSP$year == 2016, ]) # use the bbox for larval sampling 2016
# bbSouthEast <- bbox(lvDataSP[lvDataSP$region == "South East", ])
# eyreGG <- get_map(c(bbEyre), 
#                   maptype = 'watercolor', source = 'stamen', crop = FALSE)
# seGG   <- get_map(c(bbSouthEast), 
#                   maptype = 'watercolor', source = 'stamen', crop = FALSE)
# 
# # plot larval count data for the Eyre region (all hosts) 
# setwd("C:/UserData/Kym/PhD/Data/autumnSurveys/")
# 
# 
# # It would be perfect if colours were set by classes, and 
# # size/area proportional to value!
# 
# #pdf("larvalCountsAllHosts3yrs.pdf", height = 7, width = 7)
# ggmap(eyreGG) +
#   geom_point(data = lvDf,
#              aes(long, lat, size = class2, fill = class2),
#              pch = 21, col = 'black', alpha = 0.8) +
#   # geom_point(data = filter(lvDf, host == "sea rocket"),
#   #            aes(long, lat),
#   #            pch = "|", col = 'black', size = 4) +
#   theme_bw() +
#   labs(size = "No. larvae \n per 20 sweeps",
#        fill = "No. larvae \n per 20 sweeps", 
#        x = "Longitude", y = "Latitude") +
#   theme(strip.text = element_text(face = "bold", colour = "black", size = rel(1)),
#         #strip.background = element_rect(fill = "white", size = 1),
#         axis.text = element_text(colour = "black", size = rel(1)),
#         axis.title = element_text(face = "bold")) +
#   scale_fill_manual(values = heatPal) +
#   scale_size_manual(values = lvCex*2) +
#   facet_grid(year ~ month) # this order make more sense visually with displaying time sequence heading left-to-right into cropping period 
# #dev.off()
# 
# # plot larval count data for the south East region in 2015
# 
# #pdf("larvalCountsSouthEast2015.pdf", height = 7, width = 7)
# ggmap(seGG) +
#   geom_point(data = filter(lvDf, region == "South East", 
#                            year == 2015,
#                            host == "sea rocket"),
#              aes(long, lat, size = class2, fill = class2),
#              pch = 21, col = 'black', alpha = 0.8) +
#   theme_bw() +
#   labs(size = "No. larvae \n per 20 sweeps",
#        fill = "No. larvae \n per 20 sweeps", 
#        x = "Longitude", y = "Latitude") +
#   theme(strip.text = element_text(face = "bold", colour = "black", size = rel(1)),
#         #strip.background = element_rect(fill = "white", size = 1),
#         axis.text = element_text(colour = "black", size = rel(1)),
#         axis.title = element_text(face = "bold")) +
#   scale_fill_manual(values = heatPal) +
#   scale_size_manual(values = lvCex*2) +
#   facet_grid(year ~ month) # this order make more sense visually with displaying time sequence heading left-to-right into cropping period 
# #dev.off()
# 
# # group buchans, vol canola, 
# 
# 
# # pdf("mothCounts3yrs.pdf", height = 10, width = 10)
# # ggmap(eyreGG) +
# #   geom_point(data=mdf, aes(long, lat, size=class2, fill=class2),
# #              pch=21, col='black', alpha = 0.8) +
# #   theme_bw() +
# #   labs(size = "No. moths",
# #        fill = "No. moths", # change legend title
# #        x = "Longitude", y = "Latitude") +
# #   theme(strip.text = element_text(face = "bold", colour = "black", size = rel(1)),
# #         #strip.background = element_rect(fill = "white", size = 1),
# #         axis.text = element_text(colour = "black", size = rel(1)),
# #         axis.title = element_text(face = "bold")) +
# #   scale_fill_manual(values=airPal) +
# #   scale_size_manual(values=dentAQ*2) + 
# #   facet_wrap(~ monYear)
# # dev.off()
# 
# # ##################################################3
# 
# mBreaks <- c(0, 0.01, 5, 20, 50, 107)
# mothData <- autumnData %>%
#   filter(stage == "moths") %>%
#   # add a sampling timepoint factor
#   mutate(monYear = paste(str_pad(month, width = 2, side = "left", pad = 0),
#                          year, sep = "-") %>% 
#            factor())
# 
# # make a spatial points dataframe (sp class)
# msp <- SpatialPointsDataFrame(data = mothData, 
#                               coords = mothData[ , c("long", "lat")], 
#                               proj4string = CRS("+proj=longlat +ellps=WGS84"))
# 
# # set symbology based on fixed classes
# library(classInt)
# nClasses <- 5
# intervals <- classIntervals(msp$count, n = nClasses, style = 'fixed', fixedBreaks = mBreaks)
# nClasses <- length(intervals$brks) - 1
# op <- options(digits=4)
# tab <- print(intervals)
# options(op)
# 
# ## Complete Dent set of circle radii (mm)
# dent <- c(0.64, 1.14, 1.65, 2.79, 4.32, 6.22, 9.65, 12.95, 15.11)
# ## Subset for our dataset
# dentAQ <- dent[seq_len(nClasses)]
# ## Link Size and Class: findCols returns the class number of each
# ## point; cex is the vector of sizes for each data point
# idx <- findCols(intervals)
# #cexNO2 <- dentAQ[idx]
# msp$class <- factor(names(tab)[idx])
# 
# #NO2df <- data.frame(NO2sp)
# library(RColorBrewer)
# brewer.pal.info
# heatPal <- rev(heat.colors(nClasses))
# greenPal <- brewer.pal(nClasses, "Greens")
# #ylOrRd <- brewer.pal(nClasses, "YlOrRd")
# #airPal <- colorRampPalette(c('springgreen1', 'sienna3', 'gray5'))(5)
# airPal <- heatPal
# airPal[1] <- "black"
# # get a base map
# #bbPoints <- bbox(l16Msp) # use bb for larvae march 2016
# #eyreGG <- get_map(c(bbPoints), 
# #                  maptype = 'watercolor', source = 'stamen', crop = FALSE)
# 
# # try rescaling points based on proportional area
# # scalePoints <- function(x) {sqrt(x/pi)} # calculate radius to scale points by area
# # magnitude<- c(1, 2.5, 12.5, 35, 78.5) # average value of class midpoints
# # dentAQ <- scalePoints(magnitude)
# dentAQ <- c(0.64, 1.14, 1.5, 2, 2.5)
# 
# # ggplot version
# mdf <- data.frame(msp)
# mdf$class2 <- as.character(mdf$class) %>% # a hack to ensure correct class sorting, and better legend display
#   gsub("\\[0,0.01\\)", "0", .) %>%
#   gsub("\\[0.01,5\\)", "01-05", .) %>%
#   gsub("\\[5,20\\)", "05-20", .) %>%
#   gsub("\\[20,50\\)", "20-50", .) %>%
#   gsub("\\[50,107\\]", "50-107", .)
# 
# 
# # to change the text of facet labels
# mdf$month[mdf$month == 4] <- "April"
# levels(mdf$monYear)[levels(mdf$monYear) == "04-2014"] <- "April, 2014"
# levels(mdf$monYear)[levels(mdf$monYear) == "04-2015"] <- "April, 2015"
# levels(mdf$monYear)[levels(mdf$monYear) == "04-2016"] <- "April, 2016"
# 
# # plot moth count data facted by timepoint (April 2014, 2015, 2016)
# setwd("C:/UserData/Kym/PhD/Data/autumnSurveys/")
# pdf("mothCounts3yrs_vert.pdf", height = 7, width = 7)
# ggmap(eyreGG) +
#   geom_point(data=mdf, aes(long, lat, size=class2, fill=class2),
#              pch=21, col='black', alpha = 0.8) +
#   theme_bw() +
#   labs(size = "No. moths",
#        fill = "No. moths", # change legend title
#        x = "Longitude", y = "Latitude") +
#   theme(strip.text = element_text(face = "bold", colour = "black", size = rel(1)),
#         #strip.background = element_rect(fill = "white", size = 1),
#         axis.text = element_text(colour = "black", size = rel(1)),
#         axis.title = element_text(face = "bold")) +
#   scale_fill_manual(values=airPal) +
#   scale_size_manual(values=dentAQ*2) + 
#   facet_grid(year ~ month)
# dev.off()
# 
# # #########################################################################
# # Summarise incidence and abundance data from autumn surveys from 2014-2016.
# # The aim here is to compare across years, not between hosts
# # 21/3/2017
# # ###########################################################
# 
# # first define some helper functions
# 
# # ... add the trapping dates for each month and year
# addSurveyDates <- function(x) {
#   gsub("3-2014", "10-14 Mar 2014", x) %>%
#     gsub("4-2014", "7-10 Apr 2014", .) %>%
#     gsub("3-2015", "6-18 Mar 2015", .) %>%
#     gsub("4-2015", "6-16 Apr 2015", .) %>%
#     gsub("3-2016", "7-10 Mar 2016", .) %>%
#     gsub("4-2016", "4-7 Apr 2016", .)
# } # for moth data, these days are converted to trapping periods below ...
# 
# # ... hack fixer function to restore NAs below, and any other manual fixes in the table ...
# fixNA <- function(x){ 
#   x[x == "NA (NA)"] <- NA
#   x[x == "NA NA"] <- NA
#   x[x == "2 NaN"] <- 2
#   x[x == "2 NA"]  <- 2
#   x[x == "0 NaN"] <- 0
#   x[x == "NaN"] <- NA
#   x
# }
# 
# # ... fill all rows with identical survey dates col to enable \multirow
# fillDates <- function(x){ 
#   val  <- unique(x[!is.na(x)])
#   rep(val, length(x))
# }
# 
# # ... a workaround to enable latex \multirow with xtable
# # ... adapted from: http://blogs.reed.edu/ed-tech/2015/10/creating-nice-tables-using-r-markdown/
# # ... note that \multirow verically aligns text to the middle. See here: http://tex.stackexchange.com/questions/74108/text-alignment-on-top-with-multirow-command
# make_multirow <- function(x, parbox = NULL){ # parbox is the latex \parbox command, to wrap column in a \parbox 
#   rle_len <- rle(x)$lengths
#   first   <- !duplicated(x)
#   x[!first] <- ""
# 
#     if(is.null(parbox)){
#     x[first]  <- paste0("\\multirow{", rle_len, "}{*}{", x[first],  "}")
#   } else {
#     x[first]  <- paste0("\\multirow{", rle_len, "}{*}{", parbox,  "}")
#     }
#   x
# }
# 
# # ... function to abbreviate host names for convenience
# abbrHost <- function(x) {
#   gsub("lincoln weed", "LW", x) %>%
#     gsub("sea rocket", "SR", .) %>%
#     gsub("buchan weed", "BW", .) %>%
#     gsub("forage brassica", "FB", .) %>%
#     gsub("nil", "NL", .) %>%
#     gsub("vol\\. canola", "VC", .) %>%
#     gsub("wards weed", "WW", .)
# }
# 
# # ... format zeros as hyphens for display in latex tables (
# # ... this is useful when table has missing data. If no missing data, you can format zeros as hyphens by changing zeros to Na, then using latex's na.string = "-")
# zeroToHyphen <- function(x) {
#   x <- str_trim(x)
#   x[x == 0 | x == "0" | x == "0.00"] <- "-"
#   x
# }
# 
# # ... paste together only the elements of two vectors that are not blacklisted ("hyphens, nas etc)
# pasteNonHyphens <- function(x, y, # two vectors or strings to paste together  
#                             blacklist = c("-"), # blacklisted symbols
#                             na.value = "-"){ # value to fill elements not pasted (those with blacklisted symbols)
#   hyphIdx <- union(which(x %in% blacklist), which(y %in% blacklist))
#   vec <- rep(na.value, length(x))
#   if (length(vec[-hyphIdx]) > 0){
#     vec[-hyphIdx] <- paste(x[-hyphIdx], y[-hyphIdx], sep = " ")
#     vec
#   } else {
#     paste(x, y)
#   }
# }
# 
# # summarise counts by season, year, stage (moth/larvae)
# byYearHostDf <- autumnData %>%
#   group_by(year, month, host, stage) %>%
#   summarise(n = length(host),
#             nPres = sum(count > 0), # number of sites with DBM stage present
#             pPres = nPres / n, # proportion of sites where DBM were present 
#             # calc abundance data across positive sites only or all sites as desired ...
#             mnDens  = mean(count),#mean(count[which(count > 0)]), # mean abundance of positive sites
#             sdDens = sd(count)) %>% #sd(count[which(count > 0)])) %>%
#   mutate(surveyDates = addSurveyDates(paste(month, year, sep = "-"))) %>% 
#   arrange(stage, host, year, month) %>%
#   ungroup()
# 
# 
# # #############################################
# # make latex table of LARVAL autumn survey data  
# # #############################################
# 
# lvWide <- byYearHostDf %>%
#   filter(stage == "larvae",
#          host != "nil") %>% # remove pheromone trap data from sea rocket. Traps were placed way back in dune and data not really informative really informatie
#   split(.["host"]) %>%
#   lapply(function(x) {
#     # append host to vars to avoid duplicate cols
#     x <- mutate(x, dateOrder = as.numeric(paste0(year, month))) # for reordering dates for display
#     host <- abbrHost(unique(x$host))
#     names(x) <- paste(names(x), host, sep = "_")
#     x <- dplyr::select(x, -contains("year"), -contains("stage"))
#     names(x)[grepl("surveyDates_[A-Z]{2}$", names(x))] <- "surveyDates" # retain common `surveyDate` column. I had to make this hack using grep. Seomthing odd was going on with length of names attribute causing a gsub error. 
#     names(x)[grepl("dateOrder_[A-Z]{2}$", names(x))] <- "dateOrder" # retain common col
#     x
#   }) %>%
#   Reduce(function(x, y) {merge(x, y, by = c("surveyDates", "dateOrder"), 
#                                             all = TRUE)}, .) %>%
#   arrange(dateOrder)
# 
# # reorder date variable for display and convert to trapping period
# lvWide$surveyDates <- reorder(lvWide$surveyDates, lvWide$dateOrder) 
# 
# # format trapping period
# 
# 
# # round columns with complex numbers to 2 digits
# lvIdx <- c(grep("pPres", names(lvWide)), 
#            grep("mnDens", names(lvWide)),
#            grep("sdDens", names(lvWide)))
# 
# lvWide[, lvIdx] <- apply(lvWide[, lvIdx], FUN = round, digits = 2, MARGIN = 2)
# lvWide <- apply(lvWide, FUN = zeroToHyphen, MARGIN = 2) %>%
#   data.frame()
# 
# # add the latex maths env. plus/minus symbol to stdev cols
# sdColIdx <- grep("sdDens", names(lvWide))
# lvWide[, sdColIdx] <- apply(X = lvWide[, sdColIdx],
#                               FUN = function(x){
#                                 # paste latex code only to elements that are not hyphens or NA
#                                 hIdx <- which(x == "-" | is.na(x))
#                                 x <- as.character(x)
#                                 if (length(hIdx) > 0 & length(hIdx) < length(x)){
#                                   x[-hIdx] <- paste("$\\pm$", x[-hIdx])
#                                   x
#                                 } else {
#                                   if (length(hIdx) == length(x)){
#                                     x
#                                     } else {# if length hIdx == 0
#                                       paste("$\\pm$", x)
#                                     }
#                                   }
#                                 }, MARGIN = 2)
# 
# lvLtx <- lvWide %>%
#   mutate(nPres_BW = pasteNonHyphens(nPres_BW, paste0("(", pPres_BW, ")")),
#          nPres_FB = pasteNonHyphens(nPres_FB, paste0("(", pPres_FB, ")")),
#          nPres_LW = pasteNonHyphens(nPres_LW, paste0("(", pPres_LW, ")")),
#          nPres_SR = pasteNonHyphens(nPres_SR, paste0("(", pPres_SR, ")")),
#          nPres_VC = pasteNonHyphens(nPres_VC, paste0("(", pPres_VC, ")")),
#          mnDens_BW = pasteNonHyphens(mnDens_BW, sdDens_BW),
#          mnDens_FB = pasteNonHyphens(mnDens_FB, sdDens_FB),
#          mnDens_LW = pasteNonHyphens(mnDens_LW, sdDens_LW),
#          mnDens_SR = pasteNonHyphens(mnDens_SR, sdDens_SR),
#          mnDens_VC = pasteNonHyphens(mnDens_VC, sdDens_VC)) %>%
#   dplyr::select(surveyDates,
#                 n_LW, nPres_LW, mnDens_LW,
#                 n_SR, nPres_SR, mnDens_SR,
#                 n_VC, nPres_VC, mnDens_VC,
#                 n_FB, nPres_FB, mnDens_FB,
#                 n_BW, nPres_BW, mnDens_BW) %>%
#   apply(FUN = fixNA, MARGIN = 2) %>%
#   data.frame()
# 
# lvColNames <- names(lvLtx) %>%
#   gsub("_[A-Z]{2}$", "", .) %>%
#   gsub("surveyDates", "Trapping period", .) %>%
#   gsub("n$", "N", .) %>%
#   gsub("nPres", "NP", .) %>%
#   gsub("mnDens", "MD", .) %>%
#   paste(collapse = " & ") %>%
#   paste("\\\\")
# 
# # xtable::add.to.row() argument... the hosts in multicols need to be in order in mAltLtx object above
# lvHeader <- list()
# lvHeader$pos <- list(0, 0)
# lvHeader$command <- c(
#   "& \\multicolumn{3}{c}{\\textit{Diplotaxis} spp.} & \\multicolumn{3}{c}{\\textit{C. maritima}} & \\multicolumn{3}{c}{vol. Canola} & \\multicolumn{3}{c}{Forage brassica} & \\multicolumn{3}{c}{\\textit{H. incana}} 
#   \\\\\n \\cmidrule(lr){2-4}\n \\cmidrule(lr){5-7}\n \\cmidrule(lr){8-10}\n \\cmidrule(lr){11-13}\n \\cmidrule(lr){14-16}\n",
#   lvColNames)
# 
# lvCaption <- "Summary data for sites in South Australia where brassicaceous plants were sampled for \\textit{Plutella} larvae during autumn surveys over three seasons. At each location, brassica plants were sampled using a sweep net and the number of \\textit{Plutella} larvae recorded. 
# Sea rocket was sampled for larvae by hand beating the equivalent of 20 plants into a tray. At each location sampled were taken from a single brassica species, with the exception that \\textit{Diplotaxis tenufolia} and \\textit{D. muralis} commonly occurred in mixed populations
# and are grouped together. N = the number of sites sampled, NP = the number of sites where larvae were recorded with the proportion of positive sites in parentheses,
# MD = mean density of larvae per 20 sweeps (20 beat samples for sea rocket) for positive sites only $\\pm$ standard deviation."
# 
# xtable(lvLtx,
#        caption = lvCaption,
#        label = "Summary data for autumn sampling of brassica plants for larvae") %>%
#   print(include.colnames = FALSE,
#         include.rownames = FALSE,
#         table.placement = "p",
#         caption.placement = "top",
#         caption.width = "24cm", # a manual hack.
#         booktabs = TRUE,
#         add.to.row = lAltHeader,
#         scalebox = 1,
#         sanitize.text.function = function(x){x})
# 
# # ############################################
# # make latex table of MOTH autumn survey data  
# # ############################################
# 
# mWide <- byYearHostDf %>%
#   filter(stage == "moths",
#          host != "sea rocket") %>% # remove pheromone trap data from sea rocket. Traps were placed way back in dune and data not really informative really informatie
#   split(.["host"]) %>%
#   lapply(function(x) {
#     # append host to vars to avoid duplicate cols
#     x <- mutate(x, dateOrder = as.numeric(paste0(year, month))) # for reordering dates for display
#     host <- abbrHost(unique(x$host))
#     names(x) <- paste(names(x), host, sep = "_")
#     x <- dplyr::select(x, -contains("year"), -contains("stage"))
#     names(x)[grepl("surveyDates_[A-Z]{2}$", names(x))] <- "surveyDates" # retain common `surveyDate` column. I had to make this hack using grep. Seomthing odd was going on with length of names attribute causing a gsub error. 
#     names(x)[grepl("dateOrder_[A-Z]{2}$", names(x))] <- "dateOrder" # retain common col
#     x
#   }) %>%
#   Reduce(function(x, y) {merge(x, y, by = c("surveyDates", "dateOrder"), 
#                                all = TRUE)}, .) %>%
#   mutate(surveyDates = gsub("20144", "10-14 Mar -- 7-10 Apr 2014 (28--29)", as.character(dateOrder)) %>%
#            gsub("20154", "6-18 Mar -- 6-16 Apr 2015 (30$\\\\pm$1)",. ) %>%
#            gsub("20164", "7-10 Mar -- 4-7 Apr 2016 (28)",. )) %>%
#   arrange(dateOrder)
# 
# # reorder date variable for display and convert to `trapping period`
# mWide$surveyDates <- reorder(mWide$surveyDates, mWide$dateOrder)
# 
# # round columns with complex numbers to 2 digits
# mIdx <- c(grep("pPres", names(mWide)), 
#            grep("mnDens", names(mWide)),
#            grep("sdDens", names(mWide)))
# 
# mWide[, mIdx] <- apply(mWide[, mIdx], 
#                        FUN = round, digits = 2, MARGIN = 2)
# mWide <- apply(mWide, FUN = zeroToHyphen, MARGIN = 2) %>%
#   data.frame()
# 
# # add the latex maths env. plus/minus symbol to stdev cols
# sdColIdx <- grep("sdDens", names(mWide))
# mWide[, sdColIdx] <- apply(X = mWide[, sdColIdx],
#                             FUN = function(x){
#                               # paste latex code only to elements that are not hyphens or NA
#                               hIdx <- which(x == "-" | is.na(x))
#                               x <- as.character(x)
#                               if (length(hIdx) > 0 & length(hIdx) < length(x)){
#                                 x[-hIdx] <- paste("$\\pm$", x[-hIdx])
#                                 x
#                               } else {
#                                 if (length(hIdx) == length(x)){
#                                   x
#                                 } else {# if length hIdx == 0
#                                   paste("$\\pm$", x)
#                                 }
#                               }
#                             }, MARGIN = 2)
# 
# mLtx <- mWide %>%
#   mutate(nPres_BW = pasteNonHyphens(nPres_BW, paste0("(", pPres_BW, ")")),
#          nPres_FB = pasteNonHyphens(nPres_FB, paste0("(", pPres_FB, ")")),
#          nPres_LW = pasteNonHyphens(nPres_LW, paste0("(", pPres_LW, ")")),
#          nPres_NL = pasteNonHyphens(nPres_NL, paste0("(", pPres_NL, ")")),
#          nPres_VC = pasteNonHyphens(nPres_VC, paste0("(", pPres_VC, ")")),
#          mnDens_BW = pasteNonHyphens(mnDens_BW, sdDens_BW),
#          mnDens_FB = pasteNonHyphens(mnDens_FB, sdDens_FB),
#          mnDens_LW = pasteNonHyphens(mnDens_LW, sdDens_LW),
#          mnDens_NL = pasteNonHyphens(mnDens_NL, sdDens_NL),
#          mnDens_VC = pasteNonHyphens(mnDens_VC, sdDens_VC)) %>%
#   dplyr::select(surveyDates,
#                 n_LW, nPres_LW, mnDens_LW,
#                 n_VC, nPres_VC, mnDens_VC,
#                 n_FB, nPres_FB, mnDens_FB,
#                 n_BW, nPres_BW, mnDens_BW,
#                 n_NL, nPres_NL, mnDens_NL) %>%
#   apply(FUN = fixNA, MARGIN = 2) %>%
#   data.frame()
# 
# mColNames <- names(mLtx) %>%
#   gsub("_[A-Z]{2}$", "", .) %>%
#   gsub("surveyDates", "Trapping period (days)", .) %>%
#   gsub("n$", "N", .) %>%
#   gsub("nPres", "NP", .) %>%
#   gsub("mnDens", "MD", .) %>%
#   paste(collapse = " & ") %>%
#   paste("\\\\")
# 
# # xtable::add.to.row() argument... the hosts in multicols need to be in order in mAltLtx object above
# mHeader <- list()
# mHeader$pos <- list(0, 0)
# mHeader$command <- c(
#   "& \\multicolumn{3}{c}{\\textit{Diplotaxis} spp.} & \\multicolumn{3}{c}{vol. Canola} & \\multicolumn{3}{c}{Forage brassica} & \\multicolumn{3}{c}{\\textit{H. incana}} & \\multicolumn{3}{c}{Nil brassicas} 
#   \\\\\n \\cmidrule(lr){2-4}\n \\cmidrule(lr){5-7}\n \\cmidrule(lr){8-10}\n \\cmidrule(lr){11-13}\n \\cmidrule(lr){14-16}\n",
#   mColNames)
# 
# mCaption <- "Summary of pheromone trapping data for sites in South Australia where \\textit{Plutella} moths were trapped during a four week period in autumn over three seasons. At each location, a single pheromone trap was placed either directly adjacent stands of brassica plants or remote from brassicas during March,
# and the number of \\textit{Plutella} moths recorded approximately four weeks later in April. N = the number of sites, NP = the number of sites where moths were recorded with the proportion of positive sites in parentheses,
# MD = mean density of moths per trap for positive sites only $\\pm$ standard deviation."
# 
# xtable(mLtx,
#        caption = mCaption,
#        label = "Summary data for autumn pheromone trapping for moths in autumn") %>%
#   print(include.colnames = FALSE,
#         include.rownames = FALSE,
#         table.placement = "p",
#         caption.placement = "top",
#         caption.width = "24cm", # a manual hack.
#         booktabs = TRUE,
#         add.to.row = mHeader,
#         scalebox = 1,
#         sanitize.text.function = function(x){x})
# 

# # ###########################################################################################
# # Code for alternative (time cols) table format below
# 
# 
# # ############################################
# # make latex table for MOTH autumn survey data
# # ############################################
# 
# # reformat to wide, for presenting table as left-to-right time series
# mothsWide <- byYearHostDf %>%
#   filter(stage == "moths",
#          host != "sea rocket") %>% # remove pheromone trap data from sea rocket. Traps were placed way back in dune and data not really informative really informatie
#   split(.["year"]) %>%
#   lapply(function(x) {
#     # append year to vars to avoid duplicate cols
#     year <- unique(x$year)
#     names(x) <- paste(names(x), year, sep = "_")
#     x <- dplyr::select(x, -contains("year"), -contains("stage"))
#     names(x)[grepl("host_\\d{4}", names(x))] <- "host" # retain common host column. I had to make this hack using grep. Seomthing odd was going on with length of names attribute causing a gsub error. 
#     x
#   }) %>%
#   Reduce(function(x, y) {merge(x, y, by = "host", all = TRUE)}, .) %>%
#   # manually reorder for display
#   mutate(hostOrder = gsub("lincoln weed", 1, host) %>%
#            gsub("forage brassica", 3, .) %>%
#            gsub("vol. canola", 2, .) %>%
#            gsub("buchan weed", 4, .) %>%
#            gsub("wards weed", 5, .) %>%
#            gsub("nil", 6, .)) %>%
#   arrange(hostOrder)
# 
# # round cols with complex numbers to 2 digits
# mIdx <- c(grep("pPres", names(mothsWide)),
#           grep("mnDens", names(mothsWide)),
#           grep("sdDens", names(mothsWide))) %>% sort()
# 
# mothsWide[, mIdx] <- apply(mothsWide[, mIdx],
#                            FUN = round, digits = 2, MARGIN = 2)
# 
# # add the latex maths env. plus/minus symbol to sd cols
# sdCols <- c("sdDens_2014", "sdDens_2015", "sdDens_2016") 
# mothsWide[, sdCols] <- apply(X = mothsWide[, sdCols], 
#                              FUN = function(x){
#                                x[!is.na(x)] <- paste("$\\pm$", x[!is.na(x)])
#                                x
#                              }, 
#                              MARGIN = 2)
# 
# # format latex table (comine some cols)
# mothsLtx <- mothsWide %>%
#   mutate(nPres_2014 = fixNA(paste(nPres_2014, paste0("(", pPres_2014, ")"))),
#          nPres_2015 = fixNA(paste(nPres_2015, paste0("(", pPres_2015, ")"))),
#          nPres_2016 = fixNA(paste(nPres_2016, paste0("(", pPres_2016, ")"))),
#          mnDens_2014 = fixNA(paste(mnDens_2014, sdDens_2014)),
#          mnDens_2015 = fixNA(paste(mnDens_2015, sdDens_2015)),
#          mnDens_2016 = fixNA(paste(mnDens_2016, sdDens_2016))) %>%
#   dplyr::select(host, 
#                 surveyDates_2014, n_2014, nPres_2014, mnDens_2014,
#                 surveyDates_2015, n_2015, nPres_2015, mnDens_2015,
#                 surveyDates_2016, n_2016, nPres_2016, mnDens_2016)
# 
# sDateIdx <- grep("surveyDates", names(mothsLtx))
# mothsLtx[, sDateIdx] <- apply(mothsLtx[, sDateIdx],
#                               FUN = fillDates, MARGIN = 2)
# 
# # for moth data, convert survey dates to `trapping period`
# sDates  <- mothsLtx[1, sDateIdx] %>% unlist() # survey dates
# tPer14  <- sDates[grep("2014", sDates)] %>%   # trapping period, 2014
#   gsub(sDates[grep("2014", sDates)], "10-14 Mar to 7-10 Apr (28 - 29 days)")
# tPer15  <- sDates[grep("2015", sDates)] %>%   
#   gsub(sDates[grep("2015", sDates)], "6-18 Mar to 6-16 Apr (30 $\\pm$ 1 days)")
# tPer16  <- sDates[grep("2016", sDates)] %>% 
#   gsub(sDates[grep("2016", sDates)], "7-10 Mar to 4-7 Apr (28 days)") 
# 
# mothsLtx[, "surveyDates_2014"] <- tPer14
# mothsLtx[, "surveyDates_2015"] <- tPer15
# mothsLtx[, "surveyDates_2016"] <- tPer16
# 
# # generate \multirow
# tpColWidth <- "2.5cm" # trapping period column width
# pbox14 <- paste0("\\parbox{", tpColWidth, "}{10-14 Mar to\\\\7-10 Apr\\\\(28 - 29 days)}")
# pbox15 <- paste0("\\parbox{", tpColWidth, "}{6-18 Mar to\\\\6-16 Apr\\\\(30 $\\pm$ 1 days)}")
# pbox16 <- paste0("\\parbox{", tpColWidth, "}{7-10 Mar to\\\\4-7 Apr\\\\(28 days)}")
# 
# mothsLtx[, "surveyDates_2014"] <- make_multirow(mothsLtx[, "surveyDates_2014"],
#                                                 pbox14)
# mothsLtx[, "surveyDates_2015"] <- make_multirow(mothsLtx[, "surveyDates_2015"],
#                                                 pbox15)
# mothsLtx[, "surveyDates_2016"] <- make_multirow(mothsLtx[, "surveyDates_2016"],
#                                                 pbox16)
# 
# # generate latex table
# colNames <- names(mothsLtx) %>%
#   gsub("_\\d{4}$", "", .) %>%
#   gsub("host", "Plant sp.", .) %>%
#   gsub("surveyDates", "Trap period", .) %>%
#   gsub("n$", "N", .) %>%
#   gsub("nPres", "NP", .) %>%
#   gsub("mnDens", "MD", .) %>%
#   paste(collapse = " & ") %>%
#   paste("\\\\")
# 
# header <- list() # for multirow header line using add.to.row()
# header$pos <- list(0, 0)
# header$command <- c(
#   "& \\multicolumn{4}{c}{2014} & \\multicolumn{4}{c}{2015} & \\multicolumn{4}{c}{2016} \\\\\n \\cmidrule(lr){2-5}\n \\cmidrule(lr){6-9}\n \\cmidrule(lr){10-13}\n",
#   colNames)
# 
# mAlign <- c("r", "l", 
#             paste0("p{", tpColWidth, "}"), "l", "r", "r", # 2014. `p` for paragraph wrapping
#             paste0("p{", tpColWidth, "}"), "l", "r", "r",
#             paste0("p{", tpColWidth, "}"), "l", "r", "r")
# 
# mCaption <- "Summary data for \\textit{P. xylostella} moths trapped in single pheromone traps placed at sites in South Australia, during a four week period in autumn over three years, 
# grouped according to the predominant brassica species at each site. Two brassica species were sampled in all three years.
# N = the number of sites where trapping was conducted, NP = the number of sites with moths recorded with the proportion of positive sites in parentheses,
# MD = mean density of moths per trap for positive sites $\\pm$ standard deviation."
# 
# mCaption <- "Summary data for sites where \\textit{P. xylostella} moths were collected in pheromone traps associated with various wild brassica plants in South Australia during a four week period in autumn across three seasons. 
# Only a single brassica species was sampled at each site. N = the number of sites sampled, NP = the number of sites where larvae were recorded with the proportion of positive sites in parentheses,
# MD = mean density of larvae per 20 sweeps (or per 20 beat samples for sea rocket) for positive sites $\\pm$ standard deviation."
# 
# xtable(mothsLtx,
#        align = mAlign,
#        caption = mCaption,
#        label = "autumn surveys moths") %>%
#   print(include.colnames = FALSE,
#         include.rownames = FALSE,
#         caption.placement = "top",
#         booktabs = TRUE,
#         add.to.row = header,
#         scalebox = 0.7,
#         sanitize.text.function = function(x){x})
# 
# 
# # End moths
# 
# # ##############################################
# # make latex table for LARVAL autumn survey data
# # ##############################################
# 
# # reformat to wide, for presenting table as left-to-right time series
# lvWide <- byYearHostDf %>%
#   filter(stage == "larvae",
#          host != "nil") %>% # remove pheromone trap data from sea rocket. Traps were placed way back in dune and data not really informative really informatie
#   split(.["year"]) %>%
#   lapply(function(x) {
#     # append year to vars to avoid duplicate cols
#     year <- unique(x$year)
#     names(x) <- paste(names(x), year, sep = "_")
#     x <- dplyr::select(x, -contains("year"), -contains("stage"))
#     names(x)[grepl("host_\\d{4}", names(x))] <- "host" # retain common host column. I had to make this hack using grep. Seomthing odd was going on with length of names attribute causing a gsub error. 
#     names(x)[grepl("month_\\d{4}", names(x))] <- "month"
#     mutate(x, HostMon = paste0(host, month)) # unique host x month variable for correct merging
#   }) %>%
#   Reduce(function(x, y) {merge(x, y, by = "HostMon", all = TRUE)}, .) %>%
#   dplyr::select(-contains("host", ignore.case = FALSE), -contains("month")) %>%
#   mutate(month = gsub("[a-zA-Z. ]", "", HostMon),
#          host  = gsub("\\d{1}$", "", HostMon), 
#          # manually reorder for display
#          hostOrder = gsub("lincoln weed", 1, host) %>%
#            gsub("sea rocket", 2, .) %>%
#            gsub("forage brassica", 3, .) %>%
#            gsub("vol. canola", 4, .) %>%
#            gsub("buchan weed", 5, .) %>%
#            gsub("wards weed", 6, .)) %>%
#   arrange(hostOrder)
# 
# # find columns with complex numbers and round to 2 digits
# lvIdx <- c(grep("pPres", names(lvWide)),
#            grep("mnDens", names(lvWide)),
#            grep("sdDens", names(lvWide))) %>% sort()
# 
# lvWide[, lvIdx] <- apply(lvWide[, lvIdx],
#                          FUN = round, digits = 2, MARGIN = 2)
# 
# # add the latex maths env. plus/minus symbol to sd cols
# sdCols <- c("sdDens_2014", "sdDens_2015", "sdDens_2016") 
# lvWide[, sdCols] <- apply(X = lvWide[, sdCols],
#                           FUN = function(x){
#                             x[!is.na(x)] <- paste("$\\pm$", x[!is.na(x)])
#                             x
#                           },
#                           MARGIN = 2)
# 
# 
# 
# # fill the dates
# #sDateIdx <- grep("surveyDates", names(lvWide))
# mar <- grep("3", lvWide$month)
# apr <- grep("4", lvWide$month)  
# lvWide[mar, "surveyDates_2014"] <- fillDates(lvWide[mar, "surveyDates_2014"])
# lvWide[apr, "surveyDates_2014"] <- fillDates(lvWide[apr, "surveyDates_2014"])
# lvWide[mar, "surveyDates_2015"] <- fillDates(lvWide[mar, "surveyDates_2015"])
# lvWide[apr, "surveyDates_2015"] <- fillDates(lvWide[apr, "surveyDates_2015"])
# lvWide[mar, "surveyDates_2016"] <- fillDates(lvWide[mar, "surveyDates_2016"])
# lvWide[apr, "surveyDates_2016"] <- fillDates(lvWide[apr, "surveyDates_2016"])
# 
# # format latex table
# lvLtx <- lvWide %>%
#   mutate(nPres_2014 = fixNA(paste(nPres_2014, paste0("(", pPres_2014, ")"))),
#          nPres_2015 = fixNA(paste(nPres_2015, paste0("(", pPres_2015, ")"))),
#          nPres_2016 = fixNA(paste(nPres_2016, paste0("(", pPres_2016, ")"))),
#          mnDens_2014 = fixNA(paste(mnDens_2014, sdDens_2014)),
#          mnDens_2015 = fixNA(paste(mnDens_2015, sdDens_2015)),
#          mnDens_2016 = fixNA(paste(mnDens_2016, sdDens_2016))) %>%
#   dplyr::select(host, 
#                 surveyDates_2014, n_2014, nPres_2014, mnDens_2014,
#                 surveyDates_2015, n_2015, nPres_2015, mnDens_2015,
#                 surveyDates_2016, n_2016, nPres_2016, mnDens_2016)
# 
# 
# # remove duplicate hosts in host column (top-aligns, instead if mid aligns if multirow used)
# lvLtx$host[duplicated(lvLtx$host)] <- ""
# 
# # remove year from dates
# sDates <- grep("surveyDates", names(lvLtx))
# lvLtx[, sDates] <- apply(lvLtx[, sDates], 
#                          FUN = function(x) {
#                            gsub(" \\d{4}$", "", x)
#                          }, MARGIN = 2)
# 
# # Now generate latex code
# lvColNames <- names(lvLtx) %>%
#   gsub("_\\d{4}$", "", .) %>%
#   gsub("host", "Plant sp.", .) %>%
#   gsub("surveyDates", "Survey dates", .) %>%
#   gsub("n$", "N", .) %>%
#   gsub("nPres", "NP", .) %>%
#   gsub("mnDens", "MD", .) %>%
#   paste(collapse = " & ") %>%
#   paste("\\\\")
# 
# lvHeader <- list()
# lvHeader$pos <- list(0, 0)
# lvHeader$command <- c(
#   "& \\multicolumn{4}{c}{2014} & \\multicolumn{4}{c}{2015} & \\multicolumn{4}{c}{2016} \\\\\n \\cmidrule(lr){2-5}\n \\cmidrule(lr){6-9}\n \\cmidrule(lr){10-13}\n",
#   lvColNames)
# 
# lvAlign <- c("r", "l", 
#              "l", "l", "r", "r", # 2014. `p` for paragraph wrapping
#              "l", "l", "r", "r",
#              "l", "l", "r", "r")
# 
# lvCaption <- "Summary data for sites where \\textit{P. xylostella} larvae were collected during sampling of wild brassica plants in South Australia during March and April over three seasons, 
# grouped according to brassica species sampled. Only a single brassica species was sampled at each site. N = the number of sites sampled, NP = the number of sites where larvae were recorded with the proportion of positive sites in parentheses,
# MD = mean density of larvae per 20 sweeps (or per 20 beat samples for sea rocket) for positive sites $\\pm$ standard deviation."
# xtable(lvLtx,
#        align = lvAlign,
#        caption = lvCaption,
#        label = "autumn larvae") %>%
#   print(include.colnames = FALSE,
#         include.rownames = FALSE,
#         caption.placement = "top",
#         booktabs = TRUE,
#         add.to.row = lvHeader,
#         scalebox = 1,
#         sanitize.text.function = function(x){x})
# 

# ###############
# OLD CODE


# # read in the autumn survey data for 2014, 2015 and 2016
# dat2014 <- read_excel("C:/UserData/Kym/PhD/Data/autumnSurveys/autumnSurveys2014/autumnSurveyData2014.xlsx") %>%
#   filter(!is.na(lat), !is.na(long)) %>% # remove sites with no data
#   # average the count data for sites with two traps, round to integer (for Poission modelling) 
#   group_by(siteNum) %>%
#   summarise(year = unique(year),
#             consensusSiteID = unique(consensusSiteID),
#             region = unique(region),
#             #nTraps = length(site),
#             lat = mean(lat),
#             long = mean(long),
#             host = unique(host),
#             date0 = unique(date0),
#             larvae0 = round(mean(larvae0[!is.na(larvae0)])),
#             date1 = date1[1],
#             larvae1 = round(mean(larvae1[!is.na(larvae1)])),
#             moths1  = round(mean(moths1[!is.na(moths1)])))
# 
# dat2015 <- read_excel("C:/UserData/Kym/PhD/Data/autumnSurveys/autumnSurveys2015/autumnSurveyData2015.xlsx") %>%
#   filter(!is.na(lat), !is.na(long)) %>%
#   dplyr::select(siteNum, year, consensusSiteID, region, lat, long, host, 
#                 date0, larvae0, date1, larvae1, moths1)
# 
# dat2016 <- read_excel("C:/UserData/Kym/PhD/Data/autumnSurveys/autumnSurveys2016/autumnSurveyData2016.xlsx") %>%
#   filter(!is.na(lat), !is.na(long)) %>%
#   dplyr::select(siteNum, year, consensusSiteID, region, lat, long, host, 
#                 date0, larvae0, date1, larvae1, moths1)
# 
# # combine data, split march and april, re-combine, convert to long format
# autumnData <- bind_rows(dat2014, dat2015, dat2016) %>%
#   arrange(year, siteNum) 
# datM <- dplyr::select(autumnData, -date1, -larvae1, -moths1) %>%
#   rename(date = date0, larvae = larvae0)
# datA <- dplyr::select(autumnData, -date0, -larvae0) %>%
#   rename(date = date1, larvae = larvae1, moths = moths1)
# 
# autumnData <- bind_rows(datM, datA) %>%
#   melt(measure.vars = c("larvae", "moths"),
#        variable.name = "stage",
#        value.name = "count") %>%
#   mutate(month = lubridate::month(date),
#          # add siteCode for merging with other dfs
#          siteCode = paste(siteNum, year, sep = "_"),
#          # add a sampling timepoint factor
#          monYear = paste(str_pad(month, width = 2, side = "left", pad = 0),
#                          year, sep = "-") %>%
#            as.factor()) %>%
#   filter(!is.na(count)) %>% # no need to keep NAs
#   dplyr::select(siteNum, siteCode, consensusSiteID, region, year, month, date, monYear, 
#                 lat, long, host, stage, count)
# 
# # write the workspace object to a file for use in other scripts (`autumnClimate.R`)
# dir <- "C:/UserData/Kym/PhD/Data/autumnSurveys/"
# write.table(autumnData, paste0(dir, "allYearsAutumnSurveyData.txt"),
#             row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# 
# 
