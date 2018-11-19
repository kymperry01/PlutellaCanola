# #######################################################################
# Pairwise Fst matrices and heat maps for P. xylostella for publication
# K Perry, 12/4/2018
# 
# This script analyses the following datasets:
# (i)   P. australiana, 5 pops, 52 individuals
# (ii)  P. xylostella, 31 pops, 440 individuals
# (iii) P. xylostella, 28 pops, 402 individuals
#
# For each dataset, this script does:
# (i)   create a pairwise genetic (pw Fst) and geographic (km) distance matrix
# (ii)  P. australiana dataset: output latex table code
# (iii) P. xylostella datasets: output pairwise Fst heatmap figures
# (iv)  Assesses isolation by distance for the P. xylostella datasets using a mantel test

# Notes:
# This script uses the workspace in `genepop` project.
# adegenet expects a `.gen` file extension and text in the header line
# adegenet uses the last sample name in each pop as the pop name, whereas genepop uses the first
#
# ########################################################################

library(tidyverse)
library(adegenet)
library(ade4)
library(readxl)
library(geosphere)
library(xtable)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

# ###########################
# set up functions and data #
# ###########################

# function to rename 13 samples from Px to Paus (the genepop file has old sample names)
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

# function to rename pops for publication
renamePausPops <- function(x) {
  gsub("booN14sca", "Boomi", x) %>%
    gsub("calS14awa", "Calca", .) %>%
    gsub("espW14sca", "Esperance", .) %>%
    gsub("gilN14swa", "Gilgandra", .) %>%
    gsub("glbN15sca", "Goulburn", .)
}

extractPop <- function(x){
  str_sub(x, 1, 9)
}

# function to replace missing elements of a vector with values
fillValues <- function(x){
  if(length(unique(x[!is.na(x)])) > 1) warning(
    "More than one unique value to replace missing elements. Recycling replacement values."
  )
  if(all(is.na(x))) {
    x
  } else {
    x[is.na(x)] <- x[!is.na(x)]
    x 
  }
} 

dupValuesToNA <- function(x){
  x[duplicated(x)] <- NA
  x
}

# function to sort the order of a pop pair string, to output "popA & popB" regardless of pop order in the string
unifyPopPairOrder <- function(x) {
  str_split(x, " & ") %>%
    lapply(sort) %>%
    vapply(paste, collapse = " & ", FUN.VALUE = character(1))
}

popNames <- function(x, popsDf = popsMaster){
  vapply(
    X = x,
    FUN = function(x) {
      popsDf$locState[which(popsDf$popString == x)]
    },
    FUN.VALUE = character(1)
  )
}

formatLabels <- function(x) {
  format(x, nsmall = 4, scientific = FALSE) %>%
    gsub("NA", "", .)
}

# standard error function
se <- function(x, na.rm = TRUE){
  if(na.rm) x <- x[!is.na(x)]
  sd(x) / sqrt(length(x))
}


# read in geographic plotting order for P. xylostella
plotOrderPx <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "popsMasterRAD2PlotOrderPxyl.xlsx") %>%
  read_excel() %>%
  dplyr::select(popString, orderGeo, orderYearGeo)

# read in pops master to extract names and coordinates, and add plot order
popsMaster <- read_csv("C:/UserData/Kym/PhD/RAD2/popsMasterRAD2.csv") %>%
  mutate(locState = paste(Location, State)) %>%
  merge(plotOrderPx, by = "popString", all = TRUE) %>%
  arrange(Year, Host, orderGeo) %>%
  mutate(orderYearHostGeo = 1:nrow(.))


# ###########################################################################
# Plutella australiana dataset: 5 populations, 52 individuals, 974 SNP loci #
# ###########################################################################

# read in the original genepop input file to ensure correct order of population names
genindPa52 <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "genepop", "Pa52.gen") %>%
  read.genepop(ncode = 3)

pops52 <- data_frame(
  popString = unique(pop(genindPa52)) %>%
    renamePxToPaus() %>%
    extractPop(),
  genepopOrder = 1:length(popString)
  ) %>% # retain the original order of pops in the genepop file prior to merge 
  merge(popsMaster, by = "popString", all = FALSE) %>%
  mutate(popString = renamePausPops(popString)) %>%
  arrange(genepopOrder)
stopifnot(nrow(pops52) == length(unique(pop(genindPa52))))

# create genetic distance half matrix
fst52 <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "genepop", "Pa52.txt.ST2.xlsx") %>%
  read_excel(range = "A5:F10", col_names = TRUE) %>%
  dplyr::select(-pop) %>%
  as.matrix() %>%
  magrittr::set_rownames(pops52$popString) %>%
  magrittr::set_colnames(pops52$popString)

# create pairwise geographic distance half-matrix
geo52 <- distm(x = pops52[, c("Longitude", "Latitude")])# defaults to the Haversine method
geo52 <- round(geo52 / 1000, 0) # convert to km 
rownames(geo52) <- renamePausPops(pops52$popString)
colnames(geo52) <- renamePausPops(pops52$popString)
stopifnot(identical(rownames(fst52), rownames(geo52)))

# create single pairwise matrix with Fst below the diagonal and geographic distance in upper diagonal
fst52[upper.tri(fst52)] <- geo52[upper.tri(geo52)]

# test formatting for latex
wrapBrace <- function(x){
  paste0("{", x, "}")
}

# format a siunitx right-aligned table for Latex
wrapMultiCol <- function(x) {
  paste0("\\multicolumn{1}{c}{", x, "}")
}

fst52Latex <- format(fst52, digits = 4, nsmall = 0)
fst52Upper <- fst52Latex[upper.tri(fst52Latex)] %>% as.integer() %>% wrapBrace()
  as.integer() %>% wrapBrace() 
fst52Latex[upper.tri(fst52Latex)] <- fst52Upper
fst52Latex[grep("NA$", fst52Latex)] <- "{--}"
colnames(fst52Latex) <- wrapMultiCol(colnames(fst52Latex)) # multicolumn{1} decouples names from column alignment.
colnames(fst52Latex)

fst52Align <- c("l",
              "S[table-column-width=0.5cm,table-number-alignment=right,table-text-alignment=right,table-figures-integer=1,table-figures-decimal=4]",
              "S[table-column-width=0.5cm,table-number-alignment=right,table-text-alignment=right,table-figures-integer=1,table-figures-decimal=4]",
              "S[table-column-width=0.5cm,table-number-alignment=right,table-text-alignment=right,table-figures-integer=1,table-figures-decimal=4]",
              "S[table-column-width=0.5cm,table-number-alignment=right,table-text-alignment=right,table-figures-integer=1,table-figures-decimal=4]",
              "S[table-column-width=0.5cm,table-number-alignment=right,table-text-alignment=right,table-figures-integer=1,table-figures-decimal=4]")


xtable(fst52Latex, 
       caption = "Pairwise comparisons of Weir and Cockerham's \\citep{weir1984} $F_{\\textsc{st}}$ (below diagonal) and 
       geographic distance in kilometres (above diagonal) among populations of \\textit{P. australiana} from five locations.", 
       align = fst52Align,
       lab = "tab:fst") %>% 
  print.xtable(include.rownames = TRUE,
               caption.placement = "top",
               table.placement = "p",
               NA.string = "--",
               booktabs = TRUE, 
               sanitize.text.function = function(x){x})


# ##############################################################
# P. xylostella, 59 populations, 833 individuals, 1032 SNPs    #
# Plot heatmaps of the pairwise Fst for 2014 and 2015 datasets #
# ##############################################################
# 15/4/2018

# ===========================================
# check the significance of the exactG tests
# ===========================================

# all the 2014 pairwise tests had p.value = 1.000
# read in the 2015 exact G test
exactG399 <- file.path(
  "C:", "UserData", "Kym", "PhD", "RAD2", "genepop", "px399ind1032SNPs",
  "px399.txt.2G2.xlsx"
  ) %>%
  read_excel(range = "A8:E385", col_names = FALSE) %>%
  set_names(c("popA", "&", "popB", "chi2", "pValue")) %>%
  mutate(popPair = paste(popA, "&", popB),
         popUnified = unifyPopPairOrder(popPair))

# now test the different pVal corrections
temp <- exactG399$pValue
pAdj <- sapply(p.adjust.methods, function(meth) p.adjust(temp, meth)) 
View(pAdj)
#> no adjusted pValues are significant using any correction.
#> all non-significant


# extract population names from the original input genepop file (to ensure correct label order of pops in matrix)
genindPx434 <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "genepop", "px434.gen") %>%
  read.genepop(ncode = 3)

pops434 <- data_frame(
  popString = unique(pop(genindPx434)) %>%
    renamePxToPaus() %>%
    extractPop(),
  genepopOrder = 1:length(popString)
  ) %>% # retain the original order of pops in the genepop file prior to merge
  left_join(popsMaster, by = "popString") %>%
  arrange(genepopOrder) # left_join keeps order but this makes explicit.
stopifnot(nrow(pops434) == length(unique(pop(genindPx434))))

fst434 <- file.path(
  "C:", "UserData", "Kym", "PhD", "RAD2", "genepop","px434.txt.ST2.xlsx"
) %>%
  read_excel(range = ("A5:AF36"), col_names = TRUE) %>%
  dplyr::select(-pop) %>%
  as.matrix() %>%
  magrittr::set_rownames(pops434$popString) %>%
  magrittr::set_colnames(pops434$popString) %>%
  # convert matrix to long df for ggplot
  melt(varnames = c("popA", "popB"), 
       value.name = "fst") %>%
  mutate(Year = 2014)

# extract population names from the original input genepop file (to ensure correct label order of pops in matrix)
genindPx399 <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "genepop", "px399.gen") %>%
  read.genepop(ncode = 3)

pops399 <- data_frame(
  popString = unique(pop(genindPx399)) %>%
    renamePxToPaus() %>%
    extractPop(),
  genepopOrder = 1:length(popString)
  ) %>% # retain the original order of pops in the genepop file prior to merge
  left_join(popsMaster, by = "popString") %>%
  arrange(genepopOrder) # left_join keeps order but this makes explicit.
stopifnot(nrow(pops399) == length(unique(pop(genindPx399))))

fst399 <- file.path(
  "C:", "UserData", "Kym", "PhD", "RAD2", "genepop","px399.txt.ST2.xlsx"
) %>%
  read_excel(range = ("A5:AC33"), col_names = TRUE) %>%
  dplyr::select(-pop) %>%
  as.matrix() %>%
  magrittr::set_rownames(pops399$popString) %>%
  magrittr::set_colnames(pops399$popString) %>%
  # convert matrix to long df for ggplot
  melt(varnames = c("popA", "popB"), 
       value.name = "fst") %>%
  mutate(Year = 2015)

# combine for drawing a faceted plot with 2014 & 2015 datasets
fst833 <- fst434 %>%
  bind_rows(fst399) %>%
  mutate(popPair = paste(
    popA, popB, sep = " & "),
    popPairUnified = unifyPopPairOrder(popPair)
    ) %>%
  # add geographic plottng order for population "A"
  left_join(
    popsMaster %>%
      filter(Species == "Pxyl") %>%
      dplyr::select(popA = popString, orderPopA = orderYearGeo),
    by = "popA"
  ) %>%
  arrange(Year, orderPopA, popPairUnified) %>%
  mutate(orderPopA = 1:nrow(.)) %>%
  # add geographic plotting order for population "B"
  left_join(
    popsMaster %>%
      filter(Species == "Pxyl") %>%
      dplyr::select(popB = popString, orderPopB = orderYearGeo),
    by = "popB"
  ) %>%
  # ====================================================================
  # # to remove the `self-pairwise` in the plot, run this filtering step
  # filter(!popA %in% c("bunQ14scx", "bunQ15scx"), # popA is the x axis, popB the y-axis
  #        !popB %in% c("nhpW14scx", "kalW15scx")) %>%
  # ====================================================================
  # Population A or B needs to be in reverse order to plot the triangular matrix correctly
  arrange(Year, desc(orderPopB), popPairUnified) %>%
  mutate(orderPopB = 1:nrow(.)) %>%
  # reorder the population variables for plotting
  mutate(popA = reorder(popA, .$orderPopA),
         popB = reorder(popB, .$orderPopB))

# now fill (duplicate) the Fst values for all population pairs, because \ 
# \ reordering variables means the triangular half matrix is no longer symmetrical
fst833Filled <- fst833 %>%
  split(.$popPairUnified) %>%
  map(function(x){
    stopifnot(nrow(x) <= 2) # should be 1 or 2 rows ("popA & popA" pairs are not duplicated)
    x %>% mutate(fst = fillValues(fst))
  }) %>%
  bind_rows() %>%
  # now remove values for the duplicate pairs to plot a nice triangular half matrix
  arrange(popPairUnified) %>% 
  split(.$popPairUnified) %>%
  map(function(x){
    stopifnot(nrow(x) <= 2) # should be one or two rows
    x %>% mutate(fst = dupValuesToNA(fst))
  }) %>%
  bind_rows()

# plot datasets on separate panels in facted plot
# (how to adjust colour_bar: https://ggplot2.tidyverse.org/reference/guide_colourbar.html)
myPal <- brewer.pal(9, "RdBu")
(plotFst <- fst833Filled %>%
  ggplot(aes(x = popA, popB, fill = fst)) +
  geom_tile(colour = "white") + 
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 13),
        axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 12),
        aspect.ratio = 1) +  # for square tiles
  guides(fill = guide_colorbar(
    barwidth = 0.75, barheight = 10, ticks = FALSE,
    title.position = "top",
    title.hjust = 0, title.vjust = -10,
    draw.ulim = TRUE, draw.llim = TRUE
    )) +
  ggtitle(expression(paste("Pairwise ", italic(F)[ST], " comparisons for ", italic(`P. xylostella`), " populations for 1032SNPs"))) +
  labs(x = NULL, y = NULL,
       fill = expression(italic(F)[ST])) +
  scale_x_discrete(position = "bottom",
                   labels = popNames(levels(fst833$popA))) +
  scale_y_discrete(labels = popNames(levels(fst833$popB))) +
  # manually get a brewer palette by setting high and how values from the brewer red blue diverging palette.
  scale_fill_gradient2(low  = myPal[3], # start at 3 in pallete because the lowest value below zero is only half the magnitude of high above zero 
                       high = myPal[9], na.value = "white",
                       limits = c(-11, 37) / 1000,
                       breaks = seq(-10, 35, by = 5) / 1000,
                       labels = seq(-10, 35, by = 5) / 1000) +
  facet_wrap(~ Year, scales = "free", nrow = 1))
#thesisImagesPath <- file.path("C:", "UserData", "Kym", "PhD", "thesis", "images")
file.path("heatmapPxPairwiseFstPx2014-15Draft.pdf") %>%
  ggsave(width = 15, height = 7)

# =================================================
# calculate the range of Fst values for manuscript text
# =================================================

# overall
tibble(minFst = min(fst833$fst, na.rm = TRUE),
       maxFst = max(fst833$fst, na.rm = TRUE),
       mnFst = mean(fst833$fst, na.rm = TRUE),
       sdFst = sd(fst833$fst, na.rm = TRUE),
       medFst = median(fst833$fst, na.rm = TRUE))

# 2014
tibble(minFst = min(fst434$fst, na.rm = TRUE),
       maxFst = max(fst434$fst, na.rm = TRUE),
       mnFst = mean(fst434$fst, na.rm = TRUE),
       sdFst = sd(fst434$fst, na.rm = TRUE),
       medFst = median(fst434$fst, na.rm = TRUE))

# 2015
tibble(minFst = min(fst399$fst, na.rm = TRUE),
       maxFst = max(fst399$fst, na.rm = TRUE),
       mnFst = mean(fst399$fst, na.rm = TRUE),
       sdFst = sd(fst399$fst, na.rm = TRUE),
       medFst = median(fst399$fst, na.rm = TRUE))

# 2015 - excluding southend
fst399Nosouthend <- fst399 %>%
  filter(popA != "sthS15awx", popB != "sthS15awx")
tibble(minFst = min(fst399Nosouthend$fst, na.rm = TRUE),
       maxFst = max(fst399Nosouthend$fst, na.rm = TRUE),
       mnFst = mean(fst399Nosouthend$fst, na.rm = TRUE),
       sdFst = sd(fst399Nosouthend$fst, na.rm = TRUE),
       medFst = median(fst399Nosouthend$fst, na.rm = TRUE))

# 2015 - southend only
southend <- fst399 %>%
  filter(popA == "sthS15awx" | popB == "sthS15awx")
tibble(minFst = min(southend$fst, na.rm = TRUE),
       maxFst = max(southend$fst, na.rm = TRUE),
       mnFst = mean(southend$fst, na.rm = TRUE),
       sdFst = sd(southend$fst, na.rm = TRUE),
       medFst = median(southend$fst, na.rm = TRUE))

# The number of pairwise comparisons (duplicate comparisons have NAs.)
# total comparisons
# > length(fst399$fst[!is.na(fst399$fst)])
# [1] 378
# total comparisons excluding the 27 (n-1) southend comparisons
# > length(fst399$fst[!is.na(fst399$fst)]) - 27
# [1] 351

# ################################################################
# Now, plot a heatmap of geographic distance between populations #
# ################################################################

# ====
# 2014
# ====
# create a geographic distance matrix
geodist434 <- distm(x = pops434[, c("Longitude", "Latitude")]) # defaults to the Haversine method
geodist434 <- round(geodist434 / 1000, 0) %>% # convert to km
  magrittr::set_rownames(pops434$popString) %>%
  magrittr::set_colnames(pops434$popString) %>%
  # convert matrix to long df for ggplot
  melt(varnames = c("popA", "popB"),
       value.name = "distanceKm") %>%
  mutate(Year = 2014)

# ====
# 2015
# ====
# create a geographic distance matrix
geodist399 <- distm(x = pops399[, c("Longitude", "Latitude")])# defaults to the Haversine method
geodist399 <- round(geodist399 / 1000, 0) %>% # convert to km 
  magrittr::set_rownames(pops399$popString) %>%
  magrittr::set_colnames(pops399$popString) %>%
  # convert matrix to long df for ggplot
  melt(varnames = c("popA", "popB"),
       value.name = "distanceKm") %>%
  mutate(Year = 2015)

# ====================================================
# combine for drawing a faceted plot for both datasets
# ====================================================

geodist833 <- geodist434 %>%
  bind_rows(geodist399) %>%
  mutate(popPair = paste(
    popA, popB, sep = " & "),
    popPairUnified = unifyPopPairOrder(popPair)
  ) %>%
  # add geographic plottng order for population "A"
  left_join(
    popsMaster %>%
      filter(Species == "Pxyl") %>%
      dplyr::select(popA = popString, orderPopA = orderYearGeo),
    by = "popA"
  ) %>%
  arrange(Year, orderPopA, popPairUnified) %>%
  mutate(orderPopA = 1:nrow(.)) %>%
  # add geographic plotting order for population "B"
  left_join(
    popsMaster %>%
      filter(Species == "Pxyl") %>%
      dplyr::select(popB = popString, orderPopB = orderYearGeo),
    by = "popB"
  ) %>%
  # ====================================================================
# # to remove the `self-pairwise` in the plot, run this filtering step
# filter(!popA %in% c("bunQ14scx", "bunQ15scx"), # popA is the x axis, popB the y-axis
#        !popB %in% c("nhpW14scx", "kalW15scx")) %>%
# ====================================================================
# Population A or B needs to be in reverse order to plot the triangular matrix correctly
arrange(Year, desc(orderPopB), popPairUnified) %>%
  mutate(orderPopB = 1:nrow(.)) %>%
  # reorder the population variables for plotting
  mutate(popA = reorder(popA, .$orderPopA),
         popB = reorder(popB, .$orderPopB))


# remove the duplicate values to create a triangular half matrix
geodist833Filled <- geodist833 %>%
  split(.$popPairUnified) %>%
  map(function(x){
    stopifnot(nrow(x) <= 2) # should be 1 or 2 rows ("popA & popA" pairs are not duplicated)
    x %>% mutate(distanceKm = fillValues(distanceKm))
  }) %>%
  bind_rows() %>%
  # now remove values for the duplicate pairs to plot a nice triangular half matrix
  arrange(popPairUnified) %>% 
  split(.$popPairUnified) %>%
  map(function(x){
    stopifnot(nrow(x) <= 2) # should be one or two rows
    x %>% mutate(distanceKm = dupValuesToNA(distanceKm))
  }) %>%
  bind_rows()

# now draw a faceted plit with grey scale palette
greyPal <- brewer.pal(9, "Greys")
(plotGeo <- geodist833Filled %>%
  ggplot(aes(x = popA, popB, fill = distanceKm)) +
  geom_tile(colour = "white") + 
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 12),
        aspect.ratio = 1) +  # for square tiles
  guides(fill = guide_colorbar(
    barwidth = 0.75, barheight = 10, ticks = FALSE,
    title.position = "top",
    title.hjust = 0, title.vjust = -10,
    draw.ulim = TRUE, draw.llim = TRUE
  )) +
  ggtitle(expression(paste("Geographic distance in Km among populations of ", italic(`P. xylostella`)))) +
  labs(x = NULL, y = NULL,
       fill = "Distance (Km)") +
  scale_x_discrete(position = "bottom",
                   labels = popNames(levels(geodist833$popA))) +
  scale_y_discrete(labels = popNames(levels(geodist833$popB))) +
  # manually get a brewer palette by setting high and how values from the brewer red blue diverging palette.
  scale_fill_gradient(low  = greyPal[1], # start at 3 in pallete because the lowest value below zero is only half the magnitude of high above zero 
                      high = greyPal[8], na.value = "white"
                      #limits = c(-11, 37) / 1000,
                      #breaks = seq(-10, 35, by = 5) / 1000,
                       #labels = seq(-10, 35, by = 5) / 1000
                      ) +
  facet_wrap(~ Year, scales = "free", nrow = 1))
# thesisImagesPath <- file.path("C:", "UserData", "Kym", "PhD", "thesis", "images")
file.path("heatmapPxPairwiseGeodist2014-15Draft.pdf") %>%
  ggsave(width = 15, height = 7)

# ============================================================================
# now arrange the fst and geodist plots on a single panel to allow comparison
# ============================================================================

library(grid)
library(gridExtra)

# save the legend from each plot separately.
saveLegend <- function(gplot){
  tmp <- ggplot_gtable(ggplot_build(gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
legendFst <- saveLegend(plotFst)
legendGeo <- saveLegend(plotGeo)

# now make plots without legends (so that panels are all the same size) \
# and without the x-axis labels for the fst plot (the geo plot axis will suffice)
plotFstGrid <- plotFst +
  theme(legend.position = "none",
        #axis.text.x = element_blank(),
        plot.title  = element_blank())
plotGeoGrid <- plotGeo +
  theme(legend.position = "none",
        #axis.text.x = element_blank(),
        plot.title  = element_blank())

# now make nested grobs for each plot with a legend
legendWidth <- unit(0.05, "npc")
grobPlotFst <- arrangeGrob(
  plotFstGrid, legendFst, nrow = 1,
  widths  = unit.c(unit(1, "npc") - legendWidth, legendWidth),
  heights = unit.c(unit(1, "npc")) 
  )
grobPlotGeo <- arrangeGrob(
  plotGeoGrid, legendGeo, nrow = 1,
  widths = unit.c(unit(1, "npc") - legendWidth, legendWidth),
  heights = unit.c(unit(1, "npc"))
  )

pdf("heatmapPxPairwiseFstGeodist.pdf", width = 14, height = 12)
grid.arrange(grobPlotFst, grobPlotGeo, nrow = 2)
dev.off()
#> the legends are not saving correctly in the rendered pdf, but look fine on screen.
#> might need t adjust in inkscape.

# ==========================================================
# Calculate the range of geographic distances in the dataset
# ==========================================================

# 2014
# max dist:
geodist434 %>%
  filter(distanceKm == max(distanceKm)) %>%
  mutate(popNamesA = popNames(popA),
         popNamesB = popNames(popB))


tibble(minKm= min(geodist434$distanceKm, na.rm = TRUE),
       maxKm = max(geodist434$distanceKm, na.rm = TRUE),
       mnKm  = mean(geodist434$distanceKm, na.rm = TRUE),
       seKm  = se(geodist434$distanceKm, na.rm = TRUE),
       sdKm  = sd(geodist434$distanceKm, na.rm = TRUE),
       medKm = median(geodist434$distanceKm, na.rm = TRUE))

# 2015
# max dist
geodist399 %>%
  filter(distanceKm == max(distanceKm)) %>%
  mutate(popNamesA = popNames(popA),
         popNamesB = popNames(popB))

tibble(minKm = min(geodist399$distanceKm, na.rm = TRUE),
       maxKm = max(geodist399$distanceKm, na.rm = TRUE),
       mnKm  = mean(geodist399$distanceKm, na.rm = TRUE),
       seKm  = se(geodist399$distanceKm, na.rm = TRUE),
       sdKm  = sd(geodist399$distanceKm, na.rm = TRUE),
       medKm = median(geodist399$distanceKm, na.rm = TRUE))


# ======================
# Isolation by distance
# ======================
# 14/5/2018

# Assess evidence for isolation by distance for the P. xylostella datasets
# Analyse the 2014 and 2015 datasets separately 
# Method: Create pairwise distance matrices (as follows), then run a mantel test:
# -- geodist: natural log of geographic distance
# -- gendist: slatkin's linearised Fst (== Fst/(1-Fst)) 

# ============
# 2014 dataset:
# create genetic distance matrix
genDist434 <- file.path(
  "C:", "UserData", "Kym", "PhD", "RAD2", "genepop","px434.txt.ST2.xlsx"
  ) %>% # note: this code is duplicated earlier in script, but intermediate object not saved!
  read_excel(range = ("A5:AF36"), col_names = TRUE) %>%
  dplyr::select(-pop) %>%
  as.matrix() %>%
  magrittr::set_rownames(pops434$popString) %>%
  magrittr::set_colnames(pops434$popString)
genDist434 <- genDist434 / (1 - genDist434) # apply slatkin's linearisation transformation
genDist434 <- as.dist(genDist434)

# create geographic distance matrix (some code duplicted earlier in script but no intermediate objects saved)
geoDist434 <- geosphere::distm(x = pops434[, c("Longitude", "Latitude")]) # defaults to the Haversine method
geoDist434 <- log(geoDist434 / 1000)  %>% # natural log of distance in km
  magrittr::set_rownames(pops434$popString) %>%
  magrittr::set_colnames(pops434$popString) %>%
  as.dist()
stopifnot(attr(genDist434, "Labels") == attr(geoDist434, "Labels")) # check the population labels match!

# test relationship isolation by distance using a mantel test
mantel.rtest(genDist434, geoDist434, nrepet = 10000)
#> mantel's r = 0.11355, simulated p-value 0.1315868
lm434 <- lm(genDist434 ~ geoDist434) 
summary(lm434)
plot(genDist434 ~ geoDist434, xlab = "Ln Km", ylab = "Fst/(1-Fst)")
abline(coef(lm434), lty = 2, lwd = 1)

# =============  
# 2015 dataset
# create genetic distance matrix
genDist399 <- file.path(
  "C:", "UserData", "Kym", "PhD", "RAD2", "genepop","px399.txt.ST2.xlsx"
) %>%
  read_excel(range = ("A5:AC33"), col_names = TRUE) %>%
  dplyr::select(-pop) %>%
  as.matrix() %>%
  magrittr::set_rownames(pops399$popString) %>%
  magrittr::set_colnames(pops399$popString)
genDist399 <- genDist399 / (1 - genDist399) # apply slatkin's linearisation transformation
genDist399 <- as.dist(genDist399)

# create geographic distance matrix (some code duplicted earlier in script but no intermediate objects saved)
geoDist399 <- geosphere::distm(x = pops399[, c("Longitude", "Latitude")]) # defaults to the Haversine method
geoDist399[geoDist399 == 0] <- 1 # offset zero distance to 1 metre to allow log.
geoDist399 <- log(geoDist399 / 1000)  %>% # natural log of distance in km
  magrittr::set_rownames(pops399$popString) %>%
  magrittr::set_colnames(pops399$popString) %>%
  as.dist()
stopifnot(attr(genDist399, "Labels") == attr(geoDist399, "Labels")) # check the population labels match!

# test relationship isolation by distance using a mantel test
mantel.rtest(genDist399, geoDist399, nrepet = 10000)
#> mantel's r = -0.09012718, simulated p-value 0.8222178
lm399 <- lm(genDist399 ~ geoDist399) 
summary(lm399)
plot(genDist399 ~ geoDist399, xlab = "Ln Km", ylab = "Fst/(1-Fst)")
abline(coef(lm399), lty = 2, lwd = 1)

#> No evidence of isolation by distance in either year.


# End script
# ######################################################

