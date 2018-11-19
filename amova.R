# ######################################################################################
# AMOVA analysis for the P. xylostella dataset, 59 pops, 842 ind, 1008 SNPs
# K. Perry
# Last updated: 2/5/2018
#
# run hierachical AMOVA models using the `ade4` implementation in poppr::poppr.amova()

# Notes:
# We will run several hierarchical model structures:
# 1) 
# 1) Test effects of year for those locations that were resampled in 2014 and 2015 
# ######################################################################################

library(tidyverse)
library(readxl)
library(geosphere)
library(adegenet)
library(poppr)
library(xtable)


# ##################################################################################
# define helper functions to extract population metadata from the population names #
# ##################################################################################

getPop <- function(x){
  x %>% str_sub(1, 9)
}
getHost <- function(x){
  x %>% 
    str_sub(8, 8) %>%
    str_replace_all(
      c("^v$" = "Vegetables", "^c$" = "Canola", "^w$" = "Wild", "^f$" = "Forage")
      )
}
getYear <- function(x){
  x %>% 
    str_sub(5, 6) %>%
    str_replace_all(c("14" = "2014", "15" = "2015"))
}
getSeason <- function(x){
  x %>%
    str_sub(7, 7) %>%
    str_replace_all(c("a$" = "Autumn", "s$" = "Spring"))
}

# this function turns 18 resampled populations into their 9 sites sampled in 2014 and 2015)
getSite <- function(x, df = popsMaster){
  popStrToSite <- function(str){
    df[df$popString == str, "Location"]
  }
  patt <- paste(unique(x), collapse = "|")
  x %>%
    str_sub(1, 9) %>%
    str_replace("dedT14svx", "newT15svx") %>% # this site has different location names. deduplicating here. 
    str_replace_all(patt, popStrToSite)
}

# function to extract locations less than a threshold distance from each other from a matrix with distances
nearSites <- function(distm, thres = 1) {
  if (is.null(colnames(distm))) stop("Distance matrix must have column names")
  lapply(colnames(distm), function(x) {
    y <- distm[, x]
    val <- y[y < thres]
    if(length(val) > 0) {
      data_frame(siteA = names(val),
                 siteB = x,
                 distance = val)
    }
  }) %>%
    bind_rows()
}

# ###############
# Organise data #
# ###############

# read in population metadata
popsMaster <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "popsMasterRAD2.csv") %>%
  read_csv() %>%
  filter(Species == "Pxyl", Location != "Gilgandra") # drop two Gilgandra samples
  
# ============================================
# Data for analysis of all heirarchical levels  
# ============================================

# create a genind object for all P xylostella individuals
gen833 <- file.path(
  "C:", "UserData", "Kym", "PhD", "RAD2", "pgdSpider", "Px.59pops.833ind.1032SNPs.gen"
  ) %>%
  read.genepop(ncode = 3)

# replace with correct population names (genepop files tend to use the first sample name)
pop(gen833) <- pop(gen833) %>% getPop()

# set the strata slot of the genind object
strataDf <- data_frame(
  pop    = factor(pop(gen833)),
  # for the purposes of AMOVA, count forage crops as canola crops
  host   = getHost(pop), #%>% 
    #str_replace("Forage", "Canola") %>% factor(),
  season = factor(getSeason(pop)),
  year   = factor(getYear(pop))
  )
strata(gen833) <- strataDf


# =====================================================================================
# Data for a refined temporal analysis for locations resampled across seasons and years
# =====================================================================================

# First Identify the populations (locations) that were sampled in both 2014 & 2015
# Construct a pairwise geographic distance matrix between 2014 and 2015 populations
pops2014 <- popsMaster %>% filter(Year == 2014) 
pops2015 <- popsMaster %>% filter(Year == 2015)
geodist <- distm(x = pops2014[, c("Longitude", "Latitude")], # must be x y order
                 y = pops2015[, c("Longitude", "Latitude")])
geodist <- geodist / 1000 # convert to km
rownames(geodist) <- pops2014$popString # can't pipe these
colnames(geodist) <- pops2015$popString

# extract names of populations from sites resampled across 2014 and 2015 (< 1km threshold) 
(pops1km <- nearSites(geodist, thres = 1)) # locations < 1 kilometre apart across years
# A tibble: 8 x 3
# siteA     siteB distance
# <chr>     <chr>      <dbl>
#   1 bunQ14scx bunQ15scx  0.1242016
# 2 gatQ14svx gatQ15svx  0.0000000
# 3 mnjW14svx mnjW15svx  0.0000000
# 4 dedT14svx newT15svx  0.0000000 # note that these have different location names: correct prior to publication
# 5 vrgS14svx vrgS15svx  0.1158270
# 6 wbeV14svx wbeV15svx  0.0000000
# 7 picS14swx wlkS15awx  0.0000000 # spring 14 to autumn 15. note that picS14swx is Walkers Beach - see KPRAD2 notes 25/2/2017
# 8 picS14swx wlkS15swx  0.0000000 # spring 14 to spring 15. note that picS14swx is Walkers Beach - see KPRAD2 notes 25/2/2017 
#> we could exclude # 7 above (wlkS15awx), the additional time slice, as # 8 has the same location spring14-spring15

# raising the threshold to 10 km identifies an additional 3 sites
(pops10km <- nearSites(geodist, thres = 10)) # locations < 10 kilometres apart from years 
# A tibble: 11 x 3
# siteA     siteB  distance
# <chr>     <chr>     <dbl>
#   1 bunQ14scx bunQ15scx 0.1242016
# 2 gatQ14svx gatQ15svx 0.0000000
# 3 mnjW14svx mnjW15svx 0.0000000
# 4 mthS14scx mthS15scx 6.7957436
# 5 picS14awx mthS15scx 6.4806991  # remove this one from analysis: already have a mount hope comparison, and this is diff crop type.
# 6 dedT14svx newT15svx 0.0000000  # note that each year has a different location name: correct prior to publication
# 7 vrgS14svx vrgS15svx 0.1158270
# 8 wbeV14svx wbeV15svx 0.0000000
# 9 werN14svx werN15svx 7.5791641  # the werombi site was 7.6 km away. Include to represent NSW, but note in manuscript 
# 10 picS14swx wlkS15awx 0.0000000 # autumn & spring: note that picS14swx is Walkers Beach
# 11 picS14swx wlkS15swx 0.0000000 # autumn & spring: note that picS14swx is Walkers Beach

# NOTE a couple of additional locations were sampled at two time points within season
# ... wlkS15awx to wlkS15swx  (autumn 2015 to spring 2015)
# ... colebatch and tintinara (autumn 2015 to spring 2015)

# conduct analysis for locations sampled in both 2014 and 2015
# create a genind object for 14 populations (7 locations resampled across years)
# first create a vector of the pop names
resampledPopsYears <- c(pops1km$siteA, pops1km$siteB) %>% 
  unique() %>%
  .[which(. != "wlkS15awx")] # exclude the walkers beach autumn population

# now create a new genind object for the 14 pops (7 per year) \
# \ by subsetting the main genind file 
df833 <- genind2df(gen833)
df195 <- df833 %>%
  mutate(ind = row.names(.)) %>% # keep the individual names
  filter(pop %in% resampledPopsYears)

gen195 <- df2genind(df195 %>% dplyr::select(-ind, -pop), ncode = 3) 
indNames(gen195)   <- df195$ind
pop(gen195) <- getPop(df195$ind) 
View(genind2df(gen195)) # check the ind/pop names look correct

# also, write a list of these sample names for other programs, and a populations file
rad2Path <- file.path("C:", "UserData", "Kym", "PhD", "RAD2")
df195$ind %>%
  write.table(file.path(rad2Path, "Px14pops195ind.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
df195$ind %>%
  data_frame(ind = ., pop = getPop(ind)) %>%
  write.table(file.path(rad2Path, "Px14pops195indPopulations.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)

# now, set the hierchary in the strata slot of the genind object
# add an additional "location" factor for the hierarchy
strataDf195 <- data_frame(
  pop    = factor(pop(gen195)),
  # for the purposes of AMOVA, count forage crops as canola crops
  host   = getHost(pop), #%>% 
  #str_replace("Forage", "Canola") %>% factor(),
  season = factor(getSeason(pop)),
  year   = factor(getYear(pop))
) %>%
  left_join(popsMaster %>% 
              dplyr::select(pop = popString, location = Location, state = State),
            by = "pop") %>%
  # Rename Deddincgton to Newstead, because both pops are from the same location 
  mutate(location = str_replace(location, "Deddington", "Newstead"))
strata(gen195) <- strataDf195


# ###########################################################################
# Run AMOVA model using pegas implementation: spatial and temporal analysis #
# ###########################################################################

# Do not permute to test significance. Most of the signal is noise and 

# Model 1: ~year/host/pop
# P. xylotella, 59 populations, 833 individuals, 1032 SNPs
amovaModel1YHP <- poppr.amova(
  gen833, hier = ~year/host/pop, cutoff = 0.2, method = "pegas", within = FALSE
)
amovaModel1YHP

# Model 8: ~year/location (test whether spatial or temporal variance is greater)
amovaModel8ResampledLocsYL <- poppr.amova(
  gen195, hier = ~year/location, cutoff = 0.3, method = "pegas", within = FALSE
)
amovaModel8ResampledLocsYL

#> Note: other exploratory analyses conducted are retained at the bottom of this script

# ######################################################
# Make an AMOVA table using Latex code for publication #
# ######################################################

# The model structures to report are:
# -- Model 1: P. xylostella 59 pops, 833 ind, 1032 SNPs, ~year/host/pop
# -- Model 8: Temporal anaysis of 7 reampled locations: P. xylostella, 7 locations, 14 pops, 195 individuals. ~ year/location

# We don't need the total Mean Square Dev in AMOVA output
total2NA <- function(x){
  x[length(x)] <- NA
  x
}

# Make a dataframe of Model 1
varCompTotal <- sum(amovaModel1YHP$varcomp)
modelA <- amovaModel1YHP$tab %>% 
  mutate(
    factor  = row.names(.) %>%
      str_replace_all(c("year" = "Year", "host" = "Host", "pop" = "Population")),
    MSD = total2NA(MSD),
    varComp = c(amovaModel1YHP$varcomp, varCompTotal) %>% round(4), # added the total
    `%`     = ((100 * varComp) / varCompTotal) %>% round(4)
    ) %>%
  dplyr::select(factor, df, SS = SSD, MS = MSD, `Var. comp.` = varComp, `%`)

# Make a dataframe of Model 8
varCompTotalModel8 <- sum(amovaModel8ResampledLocsYL$varcomp)
modelB <- amovaModel8ResampledLocsYL$tab %>% 
  mutate(
    factor  = row.names(.) %>%
      str_replace_all(c("year" = "Year", "location" = "Location")),
    MSD = total2NA(MSD),
    varComp = c(amovaModel8ResampledLocsYL$varcomp, varCompTotalModel8) %>% round(4), # added the total
    `%`     = ((100 * varComp) / varCompTotalModel8) %>% round(4)
  ) %>%
  dplyr::select(factor, df, SS = SSD, MS = MSD, `Var. comp.` = varComp, `%`)


# output Latex code with both models in a single table
amovaCap <- c("Analysis of molecular variance (AMOVA) under two hierachical model strutures. 
In Model A, all 59 populations collected from four \\textit{Brassica} host types in 2014 and 2015 were analyzed 
and variance was partitioned among years, among host within years and among populations within host.
In Model B, populations from 7 locations sampled in both 2014 and 2015 were analyzed and variance was partitioned among years 
              and among locations within years.")
amovaDig <- c(0, 0, 0, 3, 3, 3, 3)
amovaAln <- c(
  "l", 
  "l", # use ">{\\quad}l" to indent the first column under sub-headers
  "S[table-number-alignment=center,table-figures-integer=3,table-figures-decimal=0]\n",
  "S[round-mode=places,round-precision=3,table-number-alignment=center,table-figures-integer=5,table-figures-decimal=3]\n",
  "S[round-mode=places,round-precision=3,table-number-alignment=center,table-figures-integer=3,table-figures-decimal=3]\n",
  "S[round-mode=places,round-precision=3,table-number-alignment=center,table-figures-integer=3,table-figures-decimal=3]\n",
  "S[round-mode=places,round-precision=2,table-number-alignment=center,table-figures-integer=3,table-figures-decimal=2]\n"
  )

# set up the header names, and additional rows with AMOVA model names
amovaAddRows <- list()
amovaAddRows$pos = list(0, 0,0,  5, 5)
amovaAddRows$command = c(
  "\\multicolumn{1}{l}{AMOVA summary}\\\\\n",
  "\\multicolumn{6}{l}{\\textsc{Model A}}\\\\\n",
  "Source & {df} & {SS} & {MS} & {Est. var.} & {\\si{\\percent}}\\\\\n",
  "\\multicolumn{6}{l}{\\textsc{Model B}}\\\\\n",
  "Source & {df} & {SS} & {MS} & {Est. var.} & {\\si{\\percent}}\\\\\n"
  )

# LaTeX table code  
modelA %>%
  bind_rows(modelB) %>%
  xtable(caption = amovaCap,
         digits  = amovaDig, 
         align   = amovaAln,
		 lab     = "tab:amova") %>%
  print(include.rownames = FALSE,
        include.colnames = FALSE,
        caption.placement = "top",
        table.placement   = "p",
        add.to.row = amovaAddRows,
        booktabs = TRUE,
        NA.string = "",
        sanitize.text.function = function(x){x})

# End script
# ############################################################################################################################# 
# 
# # ################################
# # EXploratory AMOVA models below #
# # ################################
# 
# # adding `within = TRUE` option (model 1)
# amovaYearHPInd <- poppr.amova(
#   gen833, hier = ~year/host/pop, cutoff = 0.2, method = "pegas", within = TRUE
# )
# amovaYearHPInd
# 
# # Model 2: ~host/year/pop
# amovaHostYP <- poppr.amova(
#   gen833, hier = ~host/year/pop, cutoff = 0.2, method = "pegas", within = FALSE
# )
# amovaHostYP
# 
# # Model 3: ~host/pop (ignoring year)
# amovaHostPop <- poppr.amova(
#   gen833, hier = ~host/pop, cutoff = 0.2, method = "pegas", within = FALSE
# )
# amovaHostPop
# 
# # Model 4: ~year/pop (ignoring host)
# amovaYearPop <- poppr.amova(
#   gen833, hier = ~year/pop, cutoff = 0.2, method = "pegas", within = FALSE
# )
# amovaYearPop
# 
# # Add some models incorporating season information (autumn versus spring)
# # Model 5: ~year/season/pop (ignoring host)
# amovaYearSP <- poppr.amova(
#   gen833, hier = ~year/season/pop, cutoff = 0.2, method = "pegas", within = FALSE
# )
# amovaYearSP
# 
# # Model 6: ~year/season/host/pop
# amovaYearSHP <- poppr.amova(
#   gen833, hier = ~year/season/host/pop, cutoff = 0.2, method = "pegas", within = FALSE
# )
# amovaYearSHP
# 
# # Model 7: ~year/host/season/pop
# amovaYearSHP <- poppr.amova(
#   gen833, hier = ~year/season/host/pop, cutoff = 0.2, method = "pegas", within = FALSE
# )
# amovaYearSHP
# 
# # Now Model 8, within = TRUE
# amovaResampledLocationsYearLocWithin <- poppr.amova(
#   gen195, hier = ~year/location, cutoff = 0.2, method = "pegas", within = TRUE
# )
# amovaResampledLocationsYearLocWithin
# 
# # Model 9 (reverse): ~location/year
# amovaResampledLocationsLocYearPop <- poppr.amova(
#   gen195, hier = ~location/year, cutoff = 0.2, method = "pegas", within = FALSE
# )
# amovaResampledLocationsLocYearPop
# 
# # Now Model with within = TRUE
# # Model 9 : ~location/year
# amovaResampledLocYearWithin <- poppr.amova(
#   gen195, hier = ~location/year, cutoff = 0.2, method = "pegas", within = TRUE
# )
# amovaResampledLocYearWithin
# 
# # =================================================================================================================
# # What happens if you analysis the data for 2014 and 2015 separately? looking for spatial & host-related structure)
# # =================================================================================================================
# 
# # Analyse the 2014 dataset: create new genind object for 2014
# df833 <- genind2df(gen833)
# df434 <- df833 %>%
#   mutate(ind = row.names(.)) %>% # keep the individual names
#   filter(getYear(pop) == "2014")
# gen434 <- df2genind(df434 %>% dplyr::select(-ind, -pop), ncode = 3) 
# indNames(gen434)   <- df434$ind
# pop(gen434) <- getPop(df434$ind) 
# View(genind2df(gen434)) # check the ind/pop names look correct
# 
# strataDf434 <- data_frame(
#   pop    = factor(pop(gen434)),
#   # for the purposes of AMOVA, count forage crops as canola crops
#   host   = getHost(pop) %>% str_replace("Forage", "Canola") %>% factor(),
#   season = factor(getSeason(pop)),
#   year   = factor(getYear(pop))
# ) %>%
#   left_join(popsMaster %>% 
#               dplyr::select(pop = popString, location = Location, state = State),
#             by = "pop") %>%
#   # Rename Deddington to Newstead, because both pops are from the same location 
#   mutate(location = str_replace(location, "Deddington", "Newstead"))
# strata(gen434) <- strataDf434
# 
# # Model 10: ~host/pop
# amova2014HostPop <- poppr.amova(
#   gen434, hier = ~host/pop, cutoff = 0.3, method = "pegas", within = TRUE
# )
# amova2014HostPop
# #> Running just 2014 and running just the year hierarchy does not change the result
# 
# # Analyse the 2015 dataset: create new genind object for 2015
# df833 <- genind2df(gen833)
# df399 <- df833 %>%
#   mutate(ind = row.names(.)) %>% # keep the individual names
#   filter(getYear(pop) == "2015")
# gen399 <- df2genind(df399 %>% dplyr::select(-ind, -pop), ncode = 3) 
# indNames(gen399)   <- df399$ind
# pop(gen399) <- getPop(df399$ind) 
# View(genind2df(gen399)) # check the ind/pop names look correct
# 
# strataDf399 <- data_frame(
#   pop    = factor(pop(gen399)),
#   # for the purposes of AMOVA, count forage crops as canola crops
#   host   = getHost(pop) %>% str_replace("Forage", "Canola") %>% factor(),
#   season = factor(getSeason(pop)),
#   year   = factor(getYear(pop))
# ) %>%
#   left_join(popsMaster %>% 
#               dplyr::select(pop = popString, location = Location, state = State),
#             by = "pop") %>%
#   # Rename Deddington to Newstead, because both pops are from the same location 
#   mutate(location = str_replace(location, "Deddington", "Newstead"))
# strata(gen399) <- strataDf399
# 
# # Model 11: ~host/pop
# amova2015HostPop <- poppr.amova(
#   gen399, hier = ~host/pop, cutoff = 0.3, method = "pegas", within = TRUE
# )
# amova2015HostPop
# #> Running just 2015 and running just the year hierarchy does not change the result
# 
# # ==================================================================
# # Finally, try the seasonal AMOVA comparison (autumn versus spring)
# # ==================================================================
# 
# # There were two locations resampled in autumn and spring in 2015
# #> Colebatch, colS15afx. Tintinara, tinS15sfx
# 
# # Make a genind object
# df28 <- df833 %>%
#   mutate(ind = row.names(.)) %>% # keep the individual names
#   filter(pop %in% c("colS15afx", "tinS15sfx"))
# gen28 <- df2genind(df28 %>% dplyr::select(-ind, -pop), ncode = 3) 
# indNames(gen28)   <- df28$ind
# pop(gen28) <- getPop(df28$ind) 
# View(genind2df(gen28)) # check the ind/pop names look correct
# 
# # Add the strata
# strataDf28 <- data_frame(
#   pop    = factor(pop(gen28)),
#   # for the purposes of AMOVA, count forage crops as canola crops
#   host   = getHost(pop) %>% str_replace("Forage", "Canola") %>% factor(),
#   season = factor(getSeason(pop)),
#   year   = factor(getYear(pop))
# ) %>%
#   left_join(popsMaster %>% 
#               dplyr::select(pop = popString, location = Location, state = State),
#             by = "pop") %>%
#   # Rename Deddington to Newstead, because both pops are from the same location 
#   mutate(location = str_replace(location, "Deddington", "Newstead"))
# strata(gen28) <- strataDf28
# 
# 
# # Model 12: ~season (ignoring host)
# amovaSeason <- poppr.amova(
#   gen28, hier = ~season, cutoff = 0.3, method = "pegas", within = FALSE
# )
# amovaSeason
# #> We can't include location as not enough degrees of freedom in the model.
# #> model table shows season only makes up 0.3% of variance anyway.
# 
# 
# # End script
# # #############################################
# # LEGACY CODE BELOW
# # #############################################
# 
# # # create genind object for all Px samples, 59 pops, 842 inds, 1008 SNPs
# # gen842 <- file.path(
# #   "C:", "UserData", "Kym", "PhD", "RAD2", "pgdSpider", "Px.59pops.842ind.1008SNPs.gen"
# # ) %>%
# #   read.genepop(ncode = 3)
# # 
# # # genepop uses the first sample for the pop names - replace with the correct pop strings 
# # pop(gen842) <- pop(gen842) %>% getPop()
# # 
# # # set the strata slot of the genind object
# # strataDf <- data_frame(
# #   pop  = pop(gen842),
# #   host = getHost(pop),
# #   season = getSeason(pop),
# #   year = getYear(pop)
# # )
# # strata(gen842) <- strataDf
# # 
# # # run AMOVA ~host/pop
# # amovaHP <- poppr.amova(gen842, hier = ~host/pop, cutoff = 0.3, method = "ade4", within = FALSE)
# # amovaHP
# # amovaHPTest <- randtest(amovaHP, nrepet = 999)
# # 
# # # output as table
# # # note: poppr.amova and randtest output use different labels - homogenise for merging.
# # # note: unhelpfully, different labels are used between slots - homogenise for merging. 
# # amovaHPTestDf <- data_frame(
# #   varSource = amovaHPTest$names %>%
# #     str_replace("Variations between host", "Between host") %>%
# #     str_replace("Variations between samples", "Between samples Within host") %>%
# #     str_replace("Variations within samples", "Within samples"),
# #   pVal = amovaHPTest$pvalue
# # )
# #   
# # (amovaHPSummary <- amovaHP$results %>%
# #   mutate(varSource = rownames(.),
# #          varSource = str_trim(varSource)) %>%
# #   filter(varSource != "Total") %>%
# #   # add variance components
# #   full_join(amovaHP$componentsofcovariance %>%
# #               mutate(varSource = rownames(.) %>%
# #                        str_replace("Variations  ", "") %>%
# #                        str_trim()) %>%
# #               filter(varSource != "Total variations")) %>%
# #   # add fixation indices
# #   full_join(amovaHP$statphi %>%
# #               mutate(varSource = rownames(.) %>%
# #                        str_replace("Phi-samples-total", "Within samples") %>%
# #                        str_replace("Phi-samples-host", "Between samples Within host") %>%
# #                        str_replace("Phi-host-total", "Between host"))) %>%
# #   # add significance values
# #   full_join(amovaHPTestDf) %>%
# #   dplyr::select(varSource, Df, `Sum Sq`, varComp = Sigma, varPerc = `%`, 
# #                 Phi, # Phi is the fixation indices
# #                 pVal))
# # 
# # # function to add the amova subscripts
# # fixationSubscripts <- function(x,
# #                                collections_total, # these vars are values in variables `source of variation` values.
# #                                samples_collections,
# #                                samples_total) {
# #   anchor <- function(y) {
# #     paste0("^", y, "$")
# #   }
# #   x %>%
# #     str_replace(anchor(collections_total), "$F_\\\\textsc\\{ct\\}$") %>%
# #     str_replace(anchor(samples_collections), "$F_\\\\textsc\\{sc\\}$") %>%
# #     str_replace(anchor(samples_total), "$F_\\\\textsc\\{st\\}$")
# # }
# # 
# # # format LateX table
# # (amovaHPLatex <- amovaHPSummary %>%
# #   mutate(
# #     fixationSymbol = fixationSubscripts(
# #       varSource, 
# #       collections_total   = "Between host",
# #       samples_collections = "Between samples Within host",
# #       samples_total       = "Within samples"
# #       ),
# #     varSource = varSource %>% 
# #       str_replace("^Between host$", "Among host types") %>%
# #       str_replace("^Between samples Within host$", "Among populations within host types") %>%
# #       str_replace("^Within samples$", "Within populations"), # or "Among individuals within host types
# #     varComp = round(varComp, 4),
# #     varPerc = round(varPerc, 4)) %>%
# #   dplyr::select(varSource, Df, `Sum Sq`, varComp, varPerc, fixationSymbol, Phi, pVal)) # all p values are non-significant (# remember to code in asterisks for significance levels in other tables)
# # 
# # # 
# # # amovaHPAlign <- c("l",
# # #                   "l",
# # #                   "S[round-mode=places,round-precision=0,table-number-alignment=center,table-figures-integer=3,table-figures-decimal=1]",
# # #                   "S[table-figures-integer=1,table-figures-decimal=4]",
# # #                   "S[table-figures-integer=1,table-figures-decimal=4]",
# # #                   "l",
# # #                   "S[round-mode=places,round-precision=4,table-figures-integer=1,table-figures-decimal=4]")
# # 
# # # set up header row for add.to.rows (the same headers can be used for both tables)
# # colNms <- list()
# # colNms$pos <- list(0)
# # colNms$command <- paste(
# #   "Source of variation & ", "{df} & ", "{\\makecell{Sum of\\\\squares}} & ", "{\\makecell{Variance\\\\components}} & ", 
# #   "{\\makecell{Percentage\\\\of variance}} & ", "\\multicolumn{2}{c}{\\makecell{Fixation\\\\index}} & ", "{$p$-value}\\\\\n"
# # )
# # 
# # amovaHPAlign <- c("l",
# #                   "l",
# #                   "S[round-mode=places,round-precision=0,table-number-alignment=center,table-figures-integer=4,table-figures-decimal=1]",
# #                   "S[round-mode=places,round-precision=3,table-number-alignment=center,table-figures-integer=5,table-figures-decimal=4]",
# #                   "S[round-mode=places,round-precision=3,table-number-alignment=center,table-figures-integer=3,table-figures-decimal=4]",
# #                   "S[round-mode=places,round-precision=3,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=4]",
# #                   "S[table-text-alignment=left]@{}",
# #                   "@{\\;{=}\\;}S[round-mode=places,round-precision=4,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=4]", # note column separated by `=` sign with whitespace
# #                   "S[round-mode=places,round-precision=4,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=4]")
# # 
# # amovaHPDigits <- c(0, 0, 0, 3, 3, 3, 0, 4, 3)
# # capHostPop <- "\\textsc{AMOVA} to compare genetic variation in nuclear SNPs among \\textit{P. xylostella} populations ($n$= 842 individuals) 
# # in four groups. Each group is a \\textit{Brassica} host type sampled. In the model, genetic variation is partitioned among host types, 
# # among populations within host types and within populations." 
# # xtable(amovaHPLatex, 
# #        caption = capHostPop, 
# #        align = amovaHPAlign,
# #        digits = amovaHPDigits,
# #        lab = "tab:fst") %>% 
# #   print.xtable(include.rownames = FALSE,
# #                include.colnames = FALSE,
# #                caption.placement = "top",
# #                table.placement = "p",
# #                add.to.row = colNms,
# #                #NA.string = "--",
# #                scalebox = 0.9,
# #                booktabs = TRUE, 
# #                sanitize.text.function = function(x){x})
# # 
# # # ######################################################################################
# # # Model 2A: Temporal comparisons for sites that were resampled in 2014 and 2015
# # # Model structure: Among years, Among sites within years, Within sites (populations)
# # # 9 resampled sites, 248 ind: ~Year/Site 
# # # ######################################################################################
# # 
# # # create a genind object for all Px samples, 59 pops, 842 inds, 1008 SNPs
# # gen248 <- file.path(
# #   "C:", "UserData", "Kym", "PhD", "RAD2", "pgdSpider", "Px.18pops.248ind.1008SNPs.gen"
# # ) %>%
# #   read.genepop(ncode = 3)
# # 
# # # genepop uses the first sample for the pop names - replace with the correct pop strings 
# # pop(gen248) <- pop(gen248) %>% getPop()
# # 
# # # set the strata slot of the genind object
# # strata248 <- data_frame(
# #   pop  = pop(gen248),
# #   site = getSite(pop),
# #   host = getHost(pop),
# #   season = getSeason(pop),
# #   year = getYear(pop)
# # )
# # strata(gen248) <- strata248
# # 
# # # run AMOVA ~ year/site
# # amovaYearS <- poppr.amova(gen248, hier = ~year/site, cutoff = 0.3, method = "ade4", within = FALSE)
# # # it drops 79 loci at cutoof of 30%
# # amovaYearS
# # amovaYearSTest <- randtest(amovaYearS, nrepet = 999)
# # 
# # # output as table
# # # note: poppr.amova and randtest output use different labels - homogenise for merging.
# # # note: unhelpfully, different labels are used between slots - homogenise for merging. 
# # amovaYearSTestDf <- data_frame(
# #   varSource = amovaYearSTest$names %>%
# #     str_replace("Variations between year", "Between year") %>%
# #     str_replace("Variations between samples", "Between samples Within year") %>%
# #     str_replace("Variations within samples", "Within samples"),
# #   pVal = amovaYearSTest$pvalue
# # )
# # 
# # (amovaYearSSummary <- amovaYearS$results %>%
# #   mutate(varSource = rownames(.),
# #          varSource = str_trim(varSource)) %>%
# #   filter(varSource != "Total") %>%
# #   # add variance components
# #   full_join(amovaYearS$componentsofcovariance %>%
# #               mutate(varSource = rownames(.) %>%
# #                        str_replace("Variations  ", "") %>%
# #                        str_trim()) %>%
# #               filter(varSource != "Total variations")) %>%
# #   # add fixation indices
# #   full_join(amovaYearS$statphi %>%
# #               mutate(varSource = rownames(.) %>%
# #                        str_replace("Phi-samples-total", "Within samples") %>%
# #                        str_replace("Phi-samples-year", "Between samples Within year") %>%
# #                        str_replace("Phi-year-total", "Between year"))) %>%
# #   # add significance values
# #   full_join(amovaYearSTestDf) %>%
# #   dplyr::select(varSource, Df, `Sum Sq`, varComp = Sigma, varPerc = `%`, 
# #                 Phi, # Phi is the fixation indices
# #                 pVal))
# # 
# # # format LateX table for Model 2A (~year/site)
# # amovaYearSLatex <- amovaYearSSummary %>%
# #   mutate(
# #     fixationSymbol = fixationSubscripts(
# #       varSource,
# #       collections_total   = "Between year",
# #       samples_collections = "Between samples Within year",
# #       samples_total       = "Within samples"
# #     ),
# #     varSource = varSource %>% 
# #       str_replace("^Between year$", "Among years") %>%
# #       str_replace("^Between samples Within year$", "Among populations within years") %>%
# #       str_replace("^Within samples$", "Within populations"), # or "Among individuals within populations
# #     varComp = round(varComp, 4),
# #     varPerc = round(varPerc, 4)) %>%
# #   dplyr::select(varSource, Df, `Sum Sq`, varComp, varPerc, fixationSymbol, Phi, pVal) # all p values are non-significant (# remember to code in asterisks for significance levels in other tables)
# # 
# # 
# # # set up header row for add.to.rows (the same headers can be used for both tables)
# # colNms <- list()
# # colNms$pos <- list(0)
# # colNms$command <- paste(
# #   "Source of variation & ", "{df} & ", "{\\makecell{Sum of\\\\squares}} & ", "{\\makecell{Variance\\\\components}} & ", 
# #   "{\\makecell{Percentage\\\\of variance}} & ", "\\multicolumn{2}{c}{\\makecell{Fixation\\\\index}} & ", "{$p$-value}\\\\\n"
# #   )
# # 
# # 
# # # align xolumn on equals sign
# # # https://tex.stackexchange.com/questions/4964/align-equals-sign-in-table
# # amovaYearSAlign <- c("l",
# #                   "l",
# #                   "S[round-mode=places,round-precision=0,table-number-alignment=center,table-figures-integer=4,table-figures-decimal=1]",
# #                   "S[round-mode=places,round-precision=3,table-number-alignment=center,table-figures-integer=5,table-figures-decimal=4]",
# #                   "S[round-mode=places,round-precision=3,table-number-alignment=center,table-figures-integer=3,table-figures-decimal=4]",
# #                   "S[round-mode=places,round-precision=3,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=4]",
# #                   "S[table-text-alignment=left]@{}",
# #                   "@{\\;{=}\\;}S[round-mode=places,round-precision=4,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=4]", # note column separated by `=` sign with whitespace
# #                   "S[round-mode=places,round-precision=4,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=4]")
# # 
# # amovaYearSDigits <- c(0, 0, 0, 3, 3, 3, 0, 4, 4)
# # capYearSite <- "\\textsc{AMOVA} to compare genetic variation in nuclear SNPs among \\textit{P. xylostella} populations ($n$= 248 individuals) 
# # in nine groups. Each group is a location in Australia sampled in 2014 and 2015. In the model, genetic variation is partitioned among years, 
# # among populations within years and within populations." 
# # xtable(amovaYearSLatex, 
# #        caption = capYearSite, 
# #        align = amovaYearSAlign,
# #        digits = amovaYearSDigits,
# #        lab = "tab:amovaYearSite") %>% 
# #   print.xtable(include.rownames = FALSE,
# #                include.colnames = FALSE,
# #                caption.placement = "top",
# #                table.placement = "p",
# #                add.to.row = colNms,
# #                #NA.string = "--",
# #                scalebox = 0.9,
# #                booktabs = TRUE, 
# #                sanitize.text.function = function(x){x})
# # 
# # 
# # 
# # # ######################################################################################
# # # Model 2B: Temporal comparisons for sites that were resampled in 2014 and 2015
# # # Model structure (reverse 2A): 
# # # ... Among sites, Among years within sites, Within years (populations)
# # # 9 resampled sites, 248 ind: ~Year/Site 
# # # ######################################################################################
# # 
# # # run AMOVA ~ site/year
# # amovaSiteY <- poppr.amova(gen248, hier = ~site/year, cutoff = 0.3, method = "ade4", within = FALSE)
# # # it drops 79 loci at cutoof of 30%
# # amovaSiteY
# # amovaSiteYTest <- randtest(amovaSiteY, nrepet = 999)
# # 
# # 
# # # output as table
# # # note: poppr.amova and randtest output use different labels - homogenise for merging.
# # # note: unhelpfully, different labels are used between slots - homogenise for merging. 
# # amovaSiteYTestDf <- data_frame(
# #   varSource = amovaSiteYTest$names %>%
# #     str_replace("Variations between site", "Between site") %>%
# #     str_replace("Variations between samples", "Between samples Within site") %>%
# #     str_replace("Variations within samples", "Within samples"),
# #   pVal = amovaSiteYTest$pvalue
# # )
# # 
# # (amovaSiteYSummary <- amovaSiteY$results %>%
# #   mutate(varSource = rownames(.),
# #          varSource = str_trim(varSource)) %>%
# #   filter(varSource != "Total") %>%
# #   # add variance components
# #   full_join(amovaSiteY$componentsofcovariance %>%
# #               mutate(varSource = rownames(.) %>%
# #                        str_replace("Variations  ", "") %>%
# #                        str_trim()) %>%
# #               filter(varSource != "Total variations")) %>%
# #   # add fixation indices
# #   full_join(amovaSiteY$statphi %>%
# #               mutate(varSource = rownames(.) %>%
# #                        str_replace("Phi-samples-total", "Within samples") %>%
# #                        str_replace("Phi-samples-site", "Between samples Within site") %>%
# #                        str_replace("Phi-site-total", "Between site"))) %>%
# #   # add significance values
# #   full_join(amovaSiteYTestDf) %>%
# #   dplyr::select(varSource, Df, `Sum Sq`, varComp = Sigma, varPerc = `%`, 
# #                 Phi, # Phi is the fixation indices
# #                 pVal))
# # 
# # 
# # 
# # # ####################################################################
# # # digression: write lists of populations to run a STRUCTURE analysis 
# # # ... to compare resampled populations across time slices
# # # October 2017 
# # # ####################################################################
# # 
# # # # resampled vegetable locations only (sep by less than 1 km, + werombi at 7.6km)
# # # old <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "amova")
# # # setwd(file.path("C:", "UserData", "Kym", "PhD", "RAD2"))
# # # pops10km %>%
# # #   mutate(crop = str_sub(pop2014, 8, 8)) %>%
# # #   filter(crop == "v") %>%
# # #   dplyr::select(pop2014, pop2015) %>%
# # #   c() %>% unlist() %>% as.character() %>%
# # #   write.table("resampled6Locations12PopsVeges.txt", 
# # #               row.names = FALSE, col.names = FALSE, quote = FALSE)
# # # setwd(old)
# # # 
# # # # resampled locations, all, ACROSS YEARS only (sep by less than 10km)
# # # # ... here we exclude additional autumn time slices to restrict to cross-year comparisons
# # # setwd(file.path("C:", "UserData", "Kym", "PhD", "RAD2"))
# # # pops10km %>%
# # #   mutate(crop = str_sub(pop2014, 8, 8)) %>%
# # #   filter(pop2014 != "picS14awx",      # exclude this comparison: have mount hope already
# # #          pop2015 != "wlkS15awx") %>%  # exclude this comparison: cross/year comparisons only here
# # #   dplyr::select(pop2014, pop2015) %>%
# # #   c() %>% unlist() %>% as.character() %>%
# # #   write.table("resampled9Locations18PopsYears.txt",
# # #               row.names = FALSE, col.names = FALSE, quote = FALSE)
# # # setwd(old)
# # # 
# # # # resampled locations, all (sep by less than 10km), including the additional autumn time slices
# # # setwd(file.path("C:", "UserData", "Kym", "PhD", "RAD2"))
# # # pops10km %>%
# # #   mutate(crop = str_sub(pop2014, 8, 8)) %>%
# # #   filter(pop2014 != "picS14awx") %>%      # exclude this comparison: have mount hope already
# # #          dplyr::select(pop2014, pop2015) %>%
# # #   c() %>% unlist() %>% as.character() %>%
# # #   c("colS15afx", "tinS15sfx") %>% # include colebatch & tintinara: these are the same location
# # #   write.table("resampled11Locations22Pops.txt",
# # #               row.names = FALSE, col.names = FALSE, quote = FALSE)
# # # setwd(old)
# # 
# # # ############
# # # ... OLD BELOW ...
# # 
# # # # now read in a pairwise Fst matrix to use for genetic distances and copy to upper tri
# # # copyLowerToUpperTriMatrix <- function(m){
# # #   m[upper.tri(m)] <- t(m)[upper.tri(m)]
# # #   m
# # # }
# # # 
# # # pxFst59pops <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "genepop", "Px842.txt.ST2.xlsx") %>%
# # #   read_excel(range = "A4:BH63") %>%
# # #   dplyr::select(-pop) %>%
# # #   as.matrix() %>%
# # #   copyLowerToUpperTriMatrix()
# # # 
# # # # use the original genepop input file to assign the names, to ensure correct order
# # # # note, a copy of the gen file in gen842 (but renamed to Px842.txt) was used to run genepop
# # # popsPwFst <- unique(pop(gen842)) %>% getPop()
# # # rownames(pxFst59pops) <- popsPwFst
# # # colnames(pxFst59pops) <- popsPwFst
# # 
# # 
# # # 15/11/2017: Next: 
# # # ... supply pairwise Fst as the distance matrix
# # # ... create genind objects for 2014 and 2015 datasets, to run separately
# # # create genind object for the 18 resampled populations
# # 
# # 
# # # # BELOW: Developing a function to automatically extract and join the results of AMOVa and k randtest. 
# # # 19/11/2017: status: a few str_replaces need trouble shooting to join correctly. Nearly there, stopping now though.
# # # makeAmovaTableAde4 <- function(amova = amovaSiteY, 
# # #                                krandtest = amovaSiteYTest, 
# # #                                groups  = "site",  # name of hierarchy 1 (highest level)
# # #                                samples = "year")  # name of hierarchy 2 
# # # {
# # #   varGroup <- paste("Variations between", groups)  
# # #   varSampB <- paste("Between", samples, "Within", groups) # I know about the capital W - there for joining dfs.
# # #   varSampW <- paste("Within", samples)
# # #   
# # #   testDf <- data_frame(
# # #     varSource = krandtest$names %>%
# # #       str_replace(varGroup, paste("Between", groups)) %>%
# # #       str_replace("Variations between samples", varSampB) %>%  
# # #       str_replace("Variations within samples",  varSampW),
# # #     pVal = krandtest$pvalue
# # #   )
# # #   
# # #   phiGroup <- paste("Phi", groups, "total", sep = "-")
# # #   phiSampB <- paste("Phi-samples", groups, sep = "-")
# # #   phiSampW <- "Phi-samples-total" # always a std string for this line
# # #   
# # #   amova$results %>%
# # #     mutate(varSource = rownames(.),
# # #            varSource = str_trim(varSource)) %>%
# # #     filter(varSource != "Total") %>%
# # #     # add variance components
# # #     full_join(amova$componentsofcovariance %>%
# # #                 mutate(varSource = rownames(.) %>%
# # #                          str_replace("Variations  ", "") %>%
# # #                          str_trim()) %>%
# # #                 filter(varSource != "Total variations")) %>%
# # #     # add fixation indices
# # #     full_join(amova$statphi %>%
# # #                 mutate(varSource = rownames(.) %>%
# # #                          str_replace(phiSampW, varSampW) %>%
# # #                          str_replace(phiSampB, varSampB) %>%
# # #                          str_replace(phiGroup, paste("Between", groups)))) %>%
# # #     # add significance values
# # #     full_join(testDf) %>%
# # #     dplyr::select(varSource, Df, varComp = Sigma, varPerc = `%`,
# # #                   Phi, # Phi is the fixation indices
# # #                   pVal)
# # #   
# # # }
# # 
# # 
