# #####################################################################################
# Summarise Plutella species and Wolbachia infection genotype data
# K. Perry, 31/3/2017
#
# This script does:
# 1. make a summary table of full collection data by population 
# 2. make a summary table of P. australiana occurrence by host and state
# 3. export the latex code for both tables
# 4. make a geographic map figure displaying the distribution of genotypes
# #####################################################################################

library(tidyverse)
library(magrittr)
library(readxl)
library(reshape2)
library(sp)
library(rgdal)
library(maps)
library(mapdata)
library(maptools)
library(mapplots)
library(mapproj) # for ggplot2 coord_map()
library(lattice)
library(ggplot2)
library(RColorBrewer)
library(xtable)


# #########################################################################
# Make a summary table of full collection and genotype data by population #
# #########################################################################

colTypes <- c("n", "t", "t", "n", "n", "n", "n", "n", "n", "t", "t", "n", "t", "t", "t",
           "t", "t", "t", "n", "t", "n", "t", "t", "t", "t", "t") %>%
  gsub("^n$", "numeric", .) %>%
  gsub("^t$", "text", .)

# read in genotype data 
# include only the samples assayed for Plutella species (including unsuccessful assays)
genotypes <- read_excel("C:/UserData/Kym/PhD/RAD2/genotypesMasterRAD2.xlsx",
                        sheet = "genotypes-master", col_types = colTypes) %>%
  filter(!is.na(popCode), # exclude empty plate/positions
         genotypedForPlutellaSp == 1) %>%
  # to avoid any confusion, remove the raw data column with wolbachia genotypes coded with faint bands in parentheses
  dplyr::select(-origWolbachia) %>% 
  mutate(wolbachia = paste0(`wol_16S-2`, `wol_16S-6`) %>%
             gsub("00", "0", . ) %>%
             gsub("01", "0", .) %>%
             gsub("10", "0", .) %>%
             gsub("11", "1", .) %>%
             as.numeric()) # presence of two clear bands is encoded as `infected` with wolbachia
  
# generate summary table of collection details and genotypes by population sample 
genotypes$species[is.na(genotypes$species)] <- "unknown"  # replace null genotypes
genotypes$popCode[genotypes$popCode == 127] <- 128        # combine Baird's Bay with Calca population
popsDf <- genotypes %>%
  group_by(popCode, species, wolbachia) %>%
  summarise(location = unique(location),
            state = unique(state),
            year = unique(year),
            season = unique(season),
            hostType = unique(hostType) %>%
              gsub("c$", "canola", .) %>%
              gsub("f$", "forage", .) %>% # combine forage with canola
              gsub("v$", "vegetables", .) %>%
              gsub("w$", "wild", .), 
            nInd = length(popCode),
            nPxyl = length(which(species == "Px")),
            nPxyl_uninfected = length(intersect(which(species == "Px"), which(wolbachia == 0))),
            nPxyl_wol = length(intersect(which(species == "Px"), which(wolbachia == 1))),
            nPaus = length(which(species == "Paus")),
            nPsp = length(which(species == "unknown"))) %>%
  # now combine the genotypes for each population sample 
  group_by(popCode) %>%
  summarise(location = unique(location),
            state = unique(state),
            year = unique(year),
            season = unique(season),
            hostType = unique(hostType) %>%
              gsub("c$", "canola", .) %>%
              gsub("f$", "canola", .) %>% # combine forage with canola
              gsub("v$", "vegetables", .) %>%
              gsub("w$", "wild", .),
            nIndS = sum(nInd), # num inds screened
            nIndG = nIndS - sum(nPsp), # num inds successfully genotyped for Plutella species (all were genotyped for wol)
            nPxyl = sum(nPxyl),
            pPxyl = nPxyl / nIndG,
            nPxyl_uninfected = sum(nPxyl_uninfected),
            nPxyl_wol = sum(nPxyl_wol),
            pPxyl_wol = nPxyl_wol / nPxyl,
            nPaus = sum(nPaus),
            pPaus = nPaus / nIndG) %>% # proportion Paus among genotyped inds
  # create a unique popName (RAD2 study pops have one assigned already, but many of these pops don't)
  mutate(popName = paste(location, state, paste0(season, year), hostType),
         season = gsub("s$", "spring", season) %>%
           gsub("a$", "autumn", .)) %>%
  # get lats/longs
  merge(read_excel("C:/UserData/Kym/PhD/RAD2/genotypesMasterRAD2.xlsx", 
                   sheet = "pops-master") %>%
          filter(popCode != 9.1) %>%
          mutate(popCode = as.numeric(popCode),
                 collectionDatePOSIX = collectionDate, # retain for arranging Df.
                 collectionDate = format(as.Date(collectionDate), "%b-%Y")) %>%
          dplyr::select(popCode, latitude, longitude, host, 
                        stageCollected, collectionDate, collectionDatePOSIX), by = "popCode")


pasteDegrees <- function(x, dig = 2){
  format(round(x, dig), nsmall = dig) %>%
    paste0("\\SI{", ., "}{\\degree}")
}

wrapParenth <- function(x){
  paste0("{(}", x, "{)}")
}

# calculate grand totals for the bottom table row
grandTots <- popsDf %>%
  summarise(nIndG = sum(nIndG),
            nPaus = sum(nPaus),
            nPxyl = sum(nPxyl),
            nPxyl_wol = sum(nPxyl_wol)) %>%
  mutate(pop = NA,
         collectionDate = NA, 
         latitude = NA,
         longitude = NA,
         host = NA,
         pPaus = wrapParenth(sprintf("%.2f", nPaus / nIndG)),
         pPxyl = wrapParenth(sprintf("%.2f", nPxyl / nIndG)),
         pPxyl_wol = wrapParenth(sprintf("%.2f", nPxyl_wol / nPxyl))) %>%
  dplyr::select(pop, collectionDate, latitude, longitude, host, # the blank col
                nIndG, nPaus, pPaus, nPxyl, pPxyl, nPxyl_wol, pPxyl_wol)

# format for latex code
genotypesLatex <- popsDf %>%
  mutate(pop = paste(location, state),
         host = gsub("Lincoln weed", "Sand rocket", host) %>%
           gsub("Dog weed", "Wall rocket", .) %>%
           gsub("vol. Canola", "Wild canola", .),
         latitude  = pasteDegrees(latitude),
         longitude = pasteDegrees(longitude),
         nPxyl = nPxyl,
         pPxyl = wrapParenth(format(round(pPxyl, 2), nsmall = 2)),
         nPaus = nPaus,
         pPaus = wrapParenth(format(round(pPaus, 2), nsmall = 2)),
         nPxyl_wol = nPxyl_wol,
         pPxyl_wol = wrapParenth(format(round(pPxyl_wol, 2), nsmall = 2))) %>%
  dplyr::select(state, location, # purely for arranging
                pop, 
                collectionDate, collectionDatePOSIX,
                latitude, longitude,
                host,
                nIndG,
                nPaus, pPaus,
                nPxyl, 
                pPxyl, 
                nPxyl_wol, pPxyl_wol) %>%
  arrange(state, location, collectionDatePOSIX) %>%
  dplyr::select(-state, -location, -collectionDatePOSIX) 

# replace a single cell value
genotypesLatex$pPxyl_wol[genotypesLatex$pPxyl_wol == "{(} NaN{)}"] <- "{\\textendash}"

tCmd   <- unlist(grandTots) %>% as.character() 
totCmd <- c("\\multicolumn{5}{l}{\\textbf{\\textit{Total}}}", tCmd[6:12])

# make a multirow header using the add.to.row() argument
gHeader <- list()
gHeader$pos <- list(0, 0, 75)
gHeader$command = c(
  paste(c("& & & & & & ",
          "\\multicolumn{2}{c}{\\textit{P. australiana}} & ",
          "\\multicolumn{4}{c}{\\textit{P. xylostella}} \\\\\n",
          "\\cmidrule(lr){7-8}\n", "\\cmidrule(lr){9-12}\n"),
        collapse = ""), 
  paste(c("Location & ", "\\makecell{Collection \\\\ date} & ",
          "{Latitude} & ", "{Longitude} & ", "Host & ",
          "{\\makecell{No. \\\\ genotyped}} & ",
          "\\multicolumn{2}{c}{No. ($f$)} & ", 
          "\\multicolumn{2}{c}{No. ($f$)} & ", 
          "\\multicolumn{2}{c}{\\makecell{No. ($f$) \\\\ \\textit{wol}-infected}} \\\\"),
        collapse = ""),
  paste0(paste(totCmd, collapse = " & "), "\\\\\n"))

algn <- c("l", 
          "l", # pop 
          "c", # collection date
          "S", # latitude
          "S", # longitude
          "l", # host
          "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]", # nIndG
          "S[table-number-alignment=right,table-figures-decimal=0,table-figures-integer=4]", # nPaus,
          "@{}S[table-align-text-pre=false,table-figures-integer=1,round-mode=places,round-precision=2]@{\\space\\space\\space}", # pPaus
          "S[table-number-alignment=right,table-figures-decimal=0,table-figures-integer=4]", # nPxyl
          "@{}S[table-align-text-pre=false,table-figures-integer=1,round-mode=places,round-precision=2]", # pPxyl
          "S[table-number-alignment=right,table-figures-decimal=0,table-figures-integer=3]", # nPxyl_wol 
          "@{}S[table-align-text-pre=false,table-figures-integer=1,round-mode=places,round-precision=2]@{\\space}") # pPxyl_wol


# genoCap <- "Summary of \\textit{Plutella} collections from Australia. 
# For each location, the numbers and frequency ($f$) of each species and the \\textit{Wolbachia} infection status of \\textit{P. xylostella} 
# are presented. All \\textit{P. australiana} individuals were infected with \\textit{Wolbachia}."
# xtable(genotypesLatex, 
#        align = algn,
#        caption = genoCap,
#        lab = "tab:genotypes") %>% 
#   print.xtable(include.rownames = FALSE,
#                include.colnames = FALSE,
#                table.placement = "p",
#                caption.placement = "top",
#                add.to.row = gHeader,
#                NA.string = "",
#                scalebox = 0.65,
#                booktabs = TRUE,
#                sanitize.text.function = function(x){x})

# A version for submission to BMC Evolutionary Biology 
# Note, BMC requires legends under tables.
# I implement this using package threeparttable commands in Latex.
# I don't hard code them in here.
# I also use package adjustbox rather than \scalebox

# The only change here is to the caption
bmcCap <- "Summary of \\textit{Plutella} collections from Australia. 
Presented are the number and frequency ($f$) of each species and the number and frequency of \\textit{Wolbachia}-infected individuals among the \\textit{P. xylostella} individuals."
xtable(genotypesLatex,
       align = algn,
       caption = bmcCap,
       lab = "tab:genotypes") %>%
  print.xtable(include.rownames = FALSE,
               include.colnames = FALSE,
               table.placement = "p",
               caption.placement = "top",
               add.to.row = gHeader,
               NA.string = "",
               #scalebox = 0.65,
               booktabs = TRUE,
               sanitize.text.function = function(x){x})

# ###############################################################
# Make a bar plot of the Plutella genotypes for a presentation #
# ###############################################################

moltenSpecies <- popsDf %>%
  melt(id.vars = c("location", "state", "year", "hostType"),
       measure.vars  = c("nPxyl", "nPaus"),
       variable.name = "species",
       value.name    = "nInd") %>%
  mutate(location2 = paste(location, state, year),
         hostType = gsub("forage", "Canola", hostType) %>% # combine forage with canola
           gsub("canola", "Canola", .) %>%
           gsub("vegetables", "Vegetable Brassicas", .) %>%
           gsub("wild", "Wild Brassicas", .),
         species = gsub("nPxyl", "P. xylostella (light shade)", species) %>%
           gsub("nPaus", "P. australiana (dark shade)", .)) %>%
  # reorder variables
  arrange(hostType, state, location, year, desc(species)) %>%
  mutate(order = 1: nrow(.)) %>%
  dplyr::rename(`Host type` = hostType, Species = species)

# hack: Add a space to Narrogin WA 2015 for one of the hosts
# ... to remove duplicate location/year
narroginWildIdx <- row.names(moltenSpecies[moltenSpecies$location2 == "Narrogin WA 2015" &
                                             moltenSpecies$`Host type` == "Wild Brassicas", ]) %>% as.numeric()
moltenSpecies$location2[narroginWildIdx] <- "Narrogin WA  2015"
moltenSpecies$location2 %<>% reorder(moltenSpecies$order)
moltenSpecies$Species   %<>% reorder(moltenSpecies$order) 

myPalette = brewer.pal(8, "Accent")#[c(3, 4, 6)]
#jpeg("PausHostCollectionsBarpot.jpg", height = 480, width = 720, quality = 150)
pdf("PausHostCollectionsBarpot.pdf", height = 6, width = 11)
ggplot(moltenSpecies, aes(x = location2, y = nInd, 
                          fill = `Host type`, alpha = Species)) +
  geom_bar(stat = "identity", colour = 'lightgrey') +
  labs(x = "Collection locations", y = "No. individuals genotyped") +
  theme_classic() +
  theme(axis.title   = element_text(face = "bold", size = 12),
        axis.text.x  = element_text(angle = 60, hjust = 1, size = 8),
        legend.title = element_text(face = "bold", size = 12),
        legend.text  = element_text(size = 12)) +
  scale_alpha_manual(values = c(0.35, 1)) +
  scale_fill_manual(values = myPalette)
dev.off()

# End make plot
# ###########################################################

# ##############################################################
# Summarise sex ratio data for infected and uninfected genotypes
# ##############################################################

naToUnresolved <- function(x){
  x[is.na(x)] <- "unresolved"
  x
}

sexRatio <- genotypes %>%
  group_by(species, wolbachia) %>%
  summarise(females = length(which(gender == "f")),
            males   = length(which(gender == "m")),
            sexUnresolved = length(which(gender == "u"))) %>%
  ungroup() %>%
  mutate(species   = naToUnresolved(species),
         wolbachia = naToUnresolved(wolbachia) %>%
           gsub(0, "uninfected", .) %>%
           gsub(1, "INFECTED", .),
         sexRatioFM = round(females / males, 2))

write.csv(sexRatio, "sexRatio.csv", row.names = FALSE)

# ########################################################################################
# Extract the collectors of each population for acknowledgement in the P. aus manuscript
# ########################################################################################

collPath <- file.path("C:", "UserData", "Kym", "PhD", "genotypes", "genotypesMasterRAD2.xlsx")
collections <- read_excel("C:/UserData/Kym/PhD/RAD2/genotypesMasterRAD2.xlsx",
                        sheet = "pops-master") %>%
  filter(popCode != 9.1) %>%
  mutate(popCode  = as.numeric(popCode)) %>%
  dplyr::select(popCode, collector) %>%
  merge(popsDf, by = "popCode")

collectorNames <- unique(collections$collector) %>% sort()
collectorNames
# acknowledge  <- collectorNames[collectorNames != "Kym Perry"] # remove myself



# ######################################################################
# Make summary table of P. australiana frequency across hosts and states 
# ######################################################################

# to summarise the overall numbers of indvs genotyped ...:
data_frame(nInd = length(genotypes$location),
           Pxyl = length(which(genotypes$species == "Px")),
           Paus = length(which(genotypes$species == "Paus")),
           Unresolved  = length(which(genotypes$species == "unknown")),
           Pxyl_pr = Pxyl / nInd,
           Paus_pr = Paus / nInd,
           pUnresolved  = Unresolved  / nInd)

# for future reference, a df of the 30 unresolved samples (not genotyped for Plutella sp. at COI)
# unresolved <- genotypes %>%
#   # to avoid any confusion, remove the raw data column with wolbachia genotypes coded with faint bands in parentheses 
#   dplyr::select(-origWolbachia) %>% 
#   mutate(wolbachia = paste0(`wol_16S-2`, `wol_16S-6`) %>%
#            gsub("00", "0", . ) %>%
#            gsub("01", "0", .) %>%
#            gsub("10", "0", .) %>%
#            gsub("11", "1", .) %>% # presence of two clear bands is encoded as `infected` with wolbachia
#            as.numeric()) %>%
#   filter(species == "unknown")

# A tibble: 1 × 7
# nInd  Pxyl  Paus   Psp   Pxyl_pr    Paus_pr     Psp_pr
# <int> <int> <int> <int>     <dbl>      <dbl>      <dbl>
#   1  1477  1300   147    30 0.8801625 0.09952607 0.02031144

# ... by host type (across all states) 
byHostDf <- popsDf %>%
  group_by(hostType) %>%
  summarise(nPops = length(hostType), # capital N = num sites
            nPopsPa = length(which(nPaus > 0)),
            pPopsPa = format(round(nPopsPa/nPops, 2), nsmall = 2),
            nIndS = sum(nIndS), # num ind screened
            nIndG = sum(nIndG), # num ind genotyped
            nIndPa  = sum(nPaus),
            pIndPa  = format(round(nIndPa/nIndG, 2), nsmall = 2)) %>%
  arrange(desc(pIndPa))

# ... by state (across all host types)
byStateDf <- popsDf %>%
  group_by(state) %>%
  summarise(nPops = length(state), # capital prefixes = `sampling locations`
            nPopsPa = length(which(nPaus > 0)),
            pPopsPa = format(round(nPopsPa/nPops, 2), nsmall = 2),
            nIndS = sum(nIndS), # lower case prefixes = `individuals`
            nIndG = sum(nIndG),
            nIndPa = sum(nPaus),
            pIndPa = format(round(nIndPa/nIndG, 2), nsmall = 2)) %>%
  arrange(desc(pIndPa))

# ... format for latex
genoSumLatex <- bind_rows(rename(byStateDf, group = state),
                          rename(byHostDf,  group = hostType))  %>%
  mutate(group = gsub("canola", "Canola", group) %>%
           gsub("vegetables", "Vegetables", .) %>%
           gsub("wild", "Wild brassicas", .) %>%
           gsub("NSW", "New South Wales (NSW)", .) %>%
           gsub("QLD", "Queenland (QLD)", .) %>%
           gsub("SA", "South Australia (SA)", .) %>%
           gsub("VIC", "Victoria (VIC)", .) %>%
           gsub("TAS", "Tasmania (TAS)", .) %>%
           gsub("WA", "Western Australia (WA)", .),
         pPopsPa = wrapParenth(pPopsPa),
         pIndPa   = wrapParenth(pIndPa)) %>%
  dplyr::select(group, nPops, nPopsPa, pPopsPa, nIndG, nIndPa, pIndPa)
  
# create new rows for `add.to.rows` argument for genotype summary
gsNms <- list()
gsNms$pos <- list(0, 0, 6)
gsNms$command <- c(
  # header row
  paste("\\multicolumn{1}{l}{Group} & ", # multicolumn for 1 column is a good trick to decouple header spacing from \quad command!  
        "{\\makecell{No. \\\\ collection \\\\ locations}} & ",
        "\\multicolumn{2}{c}{\\makecell{No. ($f$) \\\\ locations \\\\ \\textit{P.aus}}} & ",
        "{\\makecell{No. \\\\ individuals \\\\ genotyped}} & ",
        "\\multicolumn{2}{c}{\\makecell{No. ($f$) \\\\ \\textit{P.aus}}} \\\\\n", collapse = ""),
  "\\multicolumn{7}{l}{\\textbf{\\textit{Australian state}}} \\\\\n",
  "\\multicolumn{7}{l}{\\textbf{\\textit{Brassica host type}}} \\\\\n"
  )

gsAlgn <- c("l", 
            ">{\\quad}l", 
            "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]",
            "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]",
            "@{}S[table-align-text-pre=false,table-figures-integer=1,round-mode=places,round-precision=2]",
            "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]",
            "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]",
            "@{}S[table-align-text-pre=false,table-figures-integer=1,round-mode=places,round-precision=2]")
gsDig <- c(0, 
           0, 0, 0, 2, 0, 0, 2)
gsCap <- "Number and frequency ($f$) of \\textit{P. australiana} in \\textit{Plutella} collections from different Australian states and brassica host types.
"

xtable(genoSumLatex,
       caption = gsCap,
       align   = gsAlgn,
       digits  = gsDig,
       lab = "tab:genotypesHostState") %>% 
  print.xtable(floating = TRUE,
               include.rownames = FALSE,
               include.colnames = FALSE,
               add.to.row = gsNms,
               table.placement = "p",
               caption.placement = "top",
               NA.string = "",
               scalebox = 0.8,
               booktabs = TRUE,
               sanitize.text.function = function(x){x})


# #################################################################################################
# 21/10/2017: A revised table summarising by host only (for submission to BMC Evolutionary Biology 
# ... (and possibly thesis)
# #################################################################################################

hostSumLatex <- byHostDf %>%
  mutate(hostType = gsub("vegetables", "Vegetable crops", hostType) %>%
           gsub("canola", "Canola crops", .) %>%
           gsub("forage", "Forage brassicas", .) %>%
           gsub("wild", "Wild brassicas", .),
         pPopsPa = wrapParenth(pPopsPa),
         pIndPa   = wrapParenth(pIndPa)) %>%
  dplyr::select(hostType, nPops, nPopsPa, pPopsPa, nIndG, nIndPa, pIndPa)

hSumNms <- list()
hSumNms$pos <- list(0)
hSumNms$command <- paste(
  "\\multicolumn{1}{l}{Host} & ", # multicolumn for 1 column is a good trick to decouple header spacing from \quad command!  
  "{\\makecell{No.\\\\locations}} & ",
  "\\multicolumn{2}{c}{\\makecell{No. \\textit{P.aus}\\\\locations}} & ",
  "{\\makecell{No.\\\\genotyped}} & ",
  "\\multicolumn{2}{c}{\\makecell{No. \\textit{P.aus}}}\\\\\n", collapse = "")


hSumAlgn <- c("l", 
            #">{\\quad}l",
            "l",
            "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]",
            "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]",
            "@{}S[table-align-text-pre=false,table-figures-integer=1,round-mode=places,round-precision=2]",
            "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]",
            "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]",
            "@{}S[table-align-text-pre=false,table-figures-integer=1,round-mode=places,round-precision=2]")

hSumDig <- c(0, 0, 0, 0, 2, 0, 0, 2)
hSumCap <- "Frequency of \\textit{P. australiana} in \\textit{Plutella} collected from different locations 
and brassica host types."
xtable(hostSumLatex,
       caption = hSumCap,
       align   = hSumAlgn,
       digits  = hSumDig,
       lab = "tab:genotypesByHost") %>% 
  print.xtable(floating = TRUE,
               include.rownames = FALSE,
               include.colnames = FALSE,
               add.to.row = hSumNms,
               table.placement = "h",
               caption.placement = "top",
               NA.string = "",
               #scalebox = 1.0,
               booktabs = TRUE,
               sanitize.text.function = function(x){x})


# #####################################################################################
# plot an unprojected geographic map of Plutella genotypes among 75 population samples 
# ... plot P. xylostella as small white circles
# ... plot populations with P. australiana present as large pies
# ... defaults to a longlat projection
# #####################################################################################

# read in a shape file with the Australian coastline and state borders 
# ... downloaded online from `Fertility of Australian Soils` page
old <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "genotypes")
setwd("C:/UserData/Kym/PhD/Data/geoSpatialData/nsaasr9nnd_02211a04ec_alb132")
ausPolygon <- readOGR(dsn = ".", layer = "aust_cd66states")
setwd(old)

# create pies
paMelt <- popsDf %>% 
  filter(nPaus > 0) %>% 
  melt(id.vars = c("longitude", "latitude"),
       measure.vars = c("nPxyl", "nPaus"),
       variable.name = "species",
       value.name = "numInd")

pxMelt <- popsDf %>% 
  filter(nPaus == 0) %>% 
  melt(id.vars = c("longitude", "latitude"),
       measure.vars = c("nPxyl", "nPaus"),
       variable.name = "species",
       value.name = "numInd")

xyzPx <- mapplots::make.xyz(x = pxMelt[,"longitude"],
                            y = pxMelt[, "latitude"],
                            z = pxMelt$numInd, group = pxMelt$species)

xyzPa <- mapplots::make.xyz(x = paMelt[,"longitude"],
                            y = paMelt[, "latitude"],
                            z = paMelt$numInd, group = paMelt$species)

# grab the five RADseq populations
radPops <- popsDf %>%
  filter(location %in% c("Boomi", "Calca", "Esperance", "Gilgandra", "Goulburn"), 
         nPaus > 0)
coordinates(radPops) <- ~longitude + latitude

# generate gridlines and gridline labels
gridLines <- gridlines(ausPolygon, 
                       easts = seq( 120, 160, by = 10),
                       norths = seq(-40, -10, by = 10), ndisc = 200)
proj4string(gridLines) <- CRS("+proj=longlat") # set proj prints pretty labels
gridLabels  <- labels(gridLines, side = 2:3)

# make a scale bar
# function from here: http://www.flutterbys.com.au/stats/tut/tut5.4.html
scalebar <- function(loc,length,unit="km",division.cex=0.8,bg="white",border="black",...) {
  lablength <-length
  length<-length/111
  x <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1]
  y <- c(0,length/(10*3:1))+loc[2]
  rect(x[1]-diff(x)[1]/4,y[1]-(y[2]-y[1]),x[5]+strwidth(paste(round(x[5]*111-loc[1]*111,0),unit))/2+diff(x)[1]/4,y[4]+(y[4]-y[1])/2, col=bg,border=border)
  cols <- rep(c("black","white"),2)
  for (i in 1:4) rect(x[i],y[1],x[i+1],y[2],col=cols[i])
  for (i in 1:5) segments(x[i],y[2],x[i],y[3])
  labels <- round(x[c(1,3)]*111-loc[1]*111,0)
  labels <- append(labels, paste(round(x[5]*111-loc[1]*111,0),unit))
  text(x[c(1,3,5)],y[4],labels,adj=0.5,cex=division.cex)
}

greyPalette <- brewer.pal(9, "Greys")
pieColours  <- c(greyPalette[3], greyPalette[9])

# plot, write svg (the points load correctly into inkscape with svg, but slightly off with pdf)

svg("genotypesMap.svg", width = 15, height = 12)
plot(aus, col = 'white', border = 'black', lwd = 0.5)
draw.pie(xyzPx$x, xyzPx$y, xyzPx$z, 
         scale = FALSE, # if true, scales the pie size by sample size
         radius = 0.65, col = pieColours)
# overplot the Pa pies so they are easier to move in inkscape
draw.pie(xyzPa$x, xyzPa$y, xyzPa$z, 
         scale = FALSE, # if true, scales the pie size by sample size
         radius = 0.65, col = pieColours)
lines(gridLines, col = 'grey', lty = 3, lwd = 1)
text(gridLabels, col = 'grey', adj = 1, cex = 1)
# overplot points
points(popsDf[, c("longitude", "latitude")],
       pch = 21, col = 'black', bg = 'black', cex = 1.2)
points(radPops, pch = 22, bg = alpha('green', 0.5), cex = 2)
scalebar(c(115,-43),1000, border = NA, bg = NA, division.cex = 0.9)
dev.off()


# also write a pdf for compiling into the manuscript draft
setwd("C:/UserData/Kym/PhD/thesis/images")
pdf("genotypesMapLocationsOnly.pdf", width = 15, height = 12)
plot(aus, col = 'white', border = 'black', lwd = 0.5)
draw.pie(xyzPx$x, xyzPx$y, xyzPx$z, 
         scale = FALSE, # if true, scales the pie size by sample size
         radius = 0.65, col = pieColours)
# overplot the Pa pies so they are easier to move in inkscape
# draw.pie(xyzPa$x, xyzPa$y, xyzPa$z, 
#          scale = FALSE, # if true, scales the pie size by sample size
#          radius = 0.65, col = pieColours)
lines(gridLines, col = 'grey', lty = 3, lwd = 1.5)
text(gridLabels, col = 'grey', adj = 1, cex = 1)
text(x = paPops[, "longitude"], y = paPops[, "latitude"],
           labels = paPops[, "location"])
points(x = paPops[, "longitude"], y = paPops[, "latitude"],
       pch = 21, bg = 'pink', cex = 1.5)
# overplot points
# points(popsDf[, c("longitude", "latitude")],
#        pch = 21, col = 'black', bg = 'black', cex = 1.2)
points(radPops, pch = 22, bg = alpha('green', 0.5), cex = 2)
scalebar(c(115,-43),1000, border = NA, bg = NA, division.cex = 0.9)
dev.off()
setwd(old)

# also additionally, write a map with pies scaled by sample size


# also write a pdf for compiling into the manuscript draft
#setwd("C:/UserData/Kym/PhD/thesis/images")
pdf("genotypesMapScaledPies.pdf", width = 15, height = 12)
plot(aus, col = 'white', border = 'black', lwd = 0.5)
draw.pie(xyzPx$x, xyzPx$y, xyzPx$z, 
         scale = TRUE, # if true, scales the pie size by sample size
         radius = 0.65, col = pieColours)
# overplot the Pa pies so they are easier to move in inkscape
draw.pie(xyzPa$x, xyzPa$y, xyzPa$z, 
         scale = TRUE, # if true, scales the pie size by sample size
         radius = 0.65, col = pieColours)
lines(gridLines, col = 'grey', lty = 3, lwd = 1.5)
text(gridLabels, col = 'grey', adj = 1, cex = 1)
# overplot points
# points(popsDf[, c("longitude", "latitude")],
#        pch = 21, col = 'black', bg = 'black', cex = 1.2)
points(radPops, pch = 22, bg = alpha('green', 0.5), cex = 2)
scalebar(c(115,-43),1000, border = NA, bg = NA, division.cex = 0.9)
dev.off()


# End script 

# # #######################################################
# # Make a reprojected version of the above genotypes map
# # ... reprojection is simply a data transformation to 
# # ... represent a sphere on a flat surface 
# # ######################################################
# 
# # read in a shape file with the Australian coastline and state borders 
# # ... downloaded online from `Fertility of Australian Soils` page
# old <- getwd()
# setwd("C:/UserData/Kym/PhD/Data/geoSpatialData/nsaasr9nnd_02211a04ec_alb132")
# ausPoly <- readOGR(dsn = ".", layer = "aust_cd66states")
# setwd(old)
# 
# # check projection
# proj4string(ausPoly)
# # [1] NA
# # you can't reproject from NA. You need to first set one.
# # set projection as `longlat` with a datum (WGS84 is common)
# origProj <- CRS("+proj=longlat +ellps=WGS84")
# proj4string(ausPoly) <- origProj
# 
# # now, we need to assign a projection (there are many) using a valid CRS string
# # you can get projection information from the `rgdal` library
# # ... make a dataframe of EPSG projections
# EPSG <- rgdal::make_EPSG()
# # we'll select the `GDA_1994` projection:
# # ... extract the proj4string info 
# gdaPrj4 <- subset(EPSG, code == 3107)[, "prj4"] # GDA_1994_South_Australia_Lambert WKID: 3107 Authority: EPSG
# gdaPrj4
# # [1] "+proj=lcc +lat_1=-28 +lat_2=-36 +lat_0=-32 +lon_0=135 +x_0=1000000 +y_0=2000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
# # ... define desired projection
# newProj <- CRS(gdaPrj4) 
# 
# # reproject Australia polygon
# ausPolyProj  <- spTransform(ausPoly, newProj)
# #plot(ausCoast)
# #plot(ausCoastSP)
# 
# # reproject data
# popsSP <- SpatialPointsDataFrame(
#   coords = popsDf[,  c("longitude", "latitude")],
#   data   = popsDf,
#   proj4string = origProj
# )
# popsProj <- spTransform(popsSP, newProj)
# 
# # generate reprojected pies
# meltPops <- popsDf %>%
#   melt(id.vars = c("longitude", "latitude"),
#        measure.vars = c("nPxyl", "nPaus"),
#        variable.name = "species",
#        value.name = "numInd")
# 
# meltPopsSP <- SpatialPointsDataFrame(
#   coords = meltPops[, c("longitude", "latitude")],
#   data   = meltPops,
#   proj4string = origProj
# )
# meltPopsProj <- spTransform(meltPopsSP, newProj)      
# 
# xyzProj <- make.xyz(x = coordinates(meltPopsProj)[,"longitude"],
#                     y = coordinates(meltPopsProj)[, "latitude"],
#                     z = meltPopsProj$numInd, group = meltPopsProj$species)
#        
# 
# # generate gridlines, reproject
# gl <- gridlines(ausPoly,
#                 easts = seq(110, 160, by = 10),
#                 norths = seq(-40, -10, by = 10), ndisc = 200)
# glProj <- spTransform(gl, newProj)
# glLab <- labels(glProj, 
#                 labelCRS = origProj, # use untranformed lat/long labels
#                 side = 2:3)
# 
# # plot
# plot(ausPolyProj, lwd = 1, col = "white")
# lines(glProj, col = alpha("black", 0.5), lty = 3, lwd = 1.5)
# text(glLab, col = "grey", adj = 1, cex = 0.8)
# #points(meltPopsProj, pch = 19, cex = 1)
# draw.pie(xyzProj$x, xyzProj$y, xyzProj$z, scale = FALSE,
#          radius = 100000,
#          col = c(alpha("white", 0.5), gryPal[9]))

# End script
# #######################################################################

# Old code
# # read in a shape file with the Australian coastline and state borders 
# # ... downloaded online from `Fertility of Australian Soils` page
# old <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "genotypes")
# setwd("C:/UserData/Kym/PhD/Data/geoSpatialData/nsaasr9nnd_02211a04ec_alb132")
# aus <- readOGR(dsn = ".", layer = "aust_cd66states")
# setwd(old)
# 
# # plot the frequency of Plutella genotypes as pies on a map
# # ... convert data to long format, make `xyz` object for the pies
# paPops <- popsDf %>% 
#   filter(nPaus > 0)
# 
# paPopsMelt <- paPops %>% # pops with P. australiana in long format
#   melt(id.vars = c("longitude", "latitude"),
#        measure.vars = c("nPxyl", "nPaus"),
#        variable.name = "species",
#        value.name = "numInd")
# 
# xyz <- mapplots::make.xyz(x = paPopsMelt[,"longitude"],
#                           y = paPopsMelt[, "latitude"],
#                           z = paPopsMelt$numInd, group = paPopsMelt$species)
# 
# pxPops <- filter(popsDf, nPaus == 0)
# coordinates(pxPops) <- pxPops[, c("longitude", "latitude")] 
# # ... converted the above to `sp` class using the coordinates method 
# # ... alt, I could have piped the above cmd to SpatialPointsDataFrame()
# 
# # generate gridlines and labels
# gridLines <- gridlines(aus, 
#                        easts = seq( 120, 160, by = 10),
#                        norths = seq(-40, -10, by = 10), ndisc = 200)
# # ... set a projection for gridlines (allows printing of `pretty` labels)
# proj4string(glines) <- CRS("+proj=longlat")
# gridLabels  <- labels(gridLines, side = 2:3)
# 
# # grab the 5 sympatric pops used for a pop gen study 
# # ... (for highlighting by overplotting green circles)
# sym5Pops <- popsDf %>%
#   filter(location %in% c("Boomi", "Calca", "Esperance", "Gilgandra", "Goulburn"),
#          nPaus > 0)
# coordinates(sym5Pops) <- sym5Pops[, c("longitude", "latitude")] # convert to `sp`
# 
# # plot, export pdf and svg plot versions for importing into inkscape
# # ... (I use inkscape to move the overplotted pies, which is why I overplot the point locations)
# # ... set a colour palette for the pies
# gryPal <- brewer.pal(9, "Greys")
# pieCol  <- c(alpha("white", 0.5), gryPal[9])
# 
# pdf("genotypesMap.pdf", width = 10, height = 8)
# plot(aus)
# # plot the P.x pops as small white semi-transparent circles
# points(pxPops, pch = 21, col = 'black', bg = pieCol[1], cex = 2)
# draw.pie(xyz$x, xyz$y, xyz$z, 
#          scale = TRUE, # if true, scales the pie size by sample size
#          radius = 0.7, col = pieCol)
# # plot all pops as black points
# points(popsDf[, c("longitude", "latitude")], # I didn't convert popsDf to an `sp` object because I need it for other things.
#        pch = 21, col = 'black', bg = 'black')
# # overplot green circles to highlight 5 populations
# points(sym5Pops, pch = 21, bg = alpha('green', 0.5), cex = 1.5)
# text(x = paPops[, "longitude"], y = paPops[, "latitude"],
#      labels = paPops[, "location"])
# lines(glines, col = 'grey', lty = 3, lwd = 1.5)
# text(glabs, col = 'grey', adj = 1, cex = 0.8)
# #genotypesPlot <- recordPlot() # record the plot as graphics object
# dev.off()
# 
# # export an svg as well
# svg("genotypesMap.svg", width = 10, height = 8)
# #genotypesPlot
# plot(aus)
# # plot the P.x pops as small white semi-transparent circles
# points(pxPops, pch = 21, col = 'black', bg = pieCol[1], cex = 2)
# draw.pie(xyz$x, xyz$y, xyz$z, 
#          scale = FALSE, # if true, scales the pie size by sample size
#          radius = 1, col = pieCol)
# # plot all pops as black points
# points(popsDf[, c("longitude", "latitude")], # I didn't convert popsDf to an `sp` object because I need it for other things.
#        pch = 21, col = 'black', bg = 'black')
# # overplot green circles to highlight 5 populations
# points(sym5Pops, pch = 21, bg = alpha('green', 0.5), cex = 1.5)
# text(x = paPops[, "longitude"], y = paPops[, "latitude"],
#      labels = paPops[, "location"])
# lines(glines, col = 'grey', lty = 3, lwd = 1.5)
# text(glabs, col = 'grey', adj = 1, cex = 0.8)
# dev.off()

