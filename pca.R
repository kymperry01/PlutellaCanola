# ###############################################################
# Principal components analysis for P. xylostella RADseq 2 study
# K. Perry, 18/5/2018

# Notes:
# This script conducts PCA for 833 P. xylostella individuals from 59 populations
# RADseq was performed and 1032 nuclear SNPs genotyped for population genetics.
# The dataset comprises samples collected in 2014 and 2015.
# Each year is analyzed separately
# The SNPs have are genotyped in at least 80% of individuals
# Missing SNPs will need to be imputed for PCA.

# This script does:
# Performs PCA analysis using R package adegenet 
# Creates plots for publication.

# To do:
# Run analysis for different loci missingness and using different groups
# #######################################################################

library(ade4)
library(adegenet)
library(tidyverse)
library(ggplot2)

# ###############
# Organise data #
# ###############

# read in populations master for labelling the pops.
popsMaster <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "popsMasterRAD2.csv") %>%
  read_csv() %>%
  mutate(locState = paste(Location, State))

# function to replace the popString with the popName
popNames <- function(x, popsDf = popsMaster){
  vapply(
    X = x,
    FUN = function(x) {
      popsDf$locState[which(popsDf$popString == x)]
    },
    FUN.VALUE = character(1)
  )
}

# some helper functions to extract population metadata out of sample name strings
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

############################################################
# 2014 dataset: 31 populations, 434 individuals, 1032 SNPs #
############################################################

gen434 <- file.path(
  "C:", "UserData", "Kym", "PhD", "RAD2", "pgdSpider", "Px.31pops.434ind.1032SNPs.gen"
  ) %>% 
  read.genepop(ncode = 3)

# set the `population` grouping factor (genepop uses first sample name - convert to the pop name instead)
pop(gen434) <- getPop(pop(gen434))

# replace missing values using the mean (note: this reduces the variance)
allFreq434 <- scaleGen(gen434, NA.method = "mean")

# run basic pca without scaling
pca434 <- dudi.pca(allFreq434, cent = FALSE, scale = FALSE, scannf = FALSE, nf = 3)

# check the bar plot of eigenvalues
barplot(pca2014$eig[1:50], main = "PCA 2014 eigenvalues", col = heat.colors(50))

# plot the pca
myCols <- funky(nPop(gen434))
file.path("C:", "UserData", "Kym", "PhD", "RAD2", "pca",
          "pcaPx31pops434indDraft.pdf") %>%
  pdf(height = 7, width = 7)
s.class(
  pca434$li, fac = pop(gen434), col = alpha(myCols, 0.5), 
  label = popNames(levels(pop(gen434))),
  clabel = 0.8,   # cex for pop labels 
  cellipse = 1,   # cex for the ellipse # 1 takes to 50% of the way between points (1.5 is default)
  cstar  = 1,     # cex for the wheel spokes (1 takes it to the points)
  cpoint = 2,     # cex for points,
  grid = FALSE    # remove grid
)
#title("PCA of 31 P. xylostella populations from 2014")
add.scatter.eig(
  pca434$eig[1:20], nf = 3, xax = 1, yax = 2, posi = "bottomright"
)
dev.off()


############################################################
# 2015 dataset: 28 populations, 399 individuals, 1032 SNPs #
############################################################

gen399 <- file.path(
  "C:", "UserData", "Kym", "PhD", "RAD2", "pgdSpider", "Px.28pops.399ind.1032SNPs.gen"
) %>% 
  read.genepop(ncode = 3)

# set the `population` grouping factor (genepop uses first sample name - convert to the pop name instead)
pop(gen399) <- getPop(pop(gen399))

# replace missing values using the mean (note: this reduces the variance)
allFreq399 <- scaleGen(gen399, NA.method = "mean")

# run basic pca without scaling
pca399 <- dudi.pca(allFreq399, cent = FALSE, scale = FALSE, scannf = FALSE, nf = 3)

# check the bar plot of eigenvalues
barplot(pca399$eig[1:50], main = "PCA 2015 eigenvalues", col = heat.colors(50))

# plot the pca
myCols399 <- funky(nPop(gen399))
file.path("C:", "UserData", "Kym", "PhD", "RAD2", "pca",
          "pcaPx28pops399indDraft.pdf") %>%
  pdf(height = 7, width = 7)
s.class(
  pca399$li, fac = pop(gen399), col = alpha(myCols399, 0.5), 
  label = popNames(levels(pop(gen399))),
  clabel = 0.8,   # cex for pop labels 
  cellipse = 1,   # cex for the ellipse # 1 takes to 50% of the way between points (1.5 is default)
  cstar  = 1,     # cex for the wheel spokes (1 takes it to the points)
  cpoint = 2,     # cex for points,
  grid = FALSE    # remove grid
)
#title("PCA of 28 P. xylostella populations from 2015")
add.scatter.eig(
  pca399$eig[1:20], nf = 3, xax = 1, yax = 2, posi = "bottomright"
)
dev.off()

# ##############################################################################
# Try running a PCA on just the 7 locations sampled in 2014 and 2015 (14 pops) #
# ##############################################################################

# Is there greater spatial variation within or across years

# Create a new genind object with the 14 pops.
# now create a new genind object for the 14 pops (7 per year) \
# \ by subsetting the genind file for all samples 
gen833 <- file.path(
  "C:", "UserData", "Kym", "PhD", "RAD2", "pgdSpider", "Px.59pops.833ind.1032SNPs.gen"
) %>% 
  read.genepop(ncode = 3)
pop(gen833) <- getPop(pop(gen833))

# read in a list of the 14 populations from the 7 resampled locations
resampledPopsYears <- file.path(
  "C:", "UserData", "Kym", "PhD", "RAD2", "Px14pops195ind.txt"
) %>%
  readLines() %>%
  getPop() %>%
  unique()

df833 <- genind2df(gen833)
df195 <- df833 %>%
  mutate(ind = row.names(.)) %>% # keep the individual names
  filter(pop %in% resampledPopsYears)

gen195 <- df2genind(df195 %>% dplyr::select(-ind, -pop), ncode = 3) 
indNames(gen195)   <- df195$ind
pop(gen195) <- getPop(df195$ind) 
View(genind2df(gen195)) # check the ind/pop names look correct

# # now, set the hierchary in the strata slot of the genind object
# # add an additional "location" factor for the hierarchy
strataDf195 <- data_frame(
  pop    = factor(pop(gen195)),
  host   = getHost(pop), #%>%
  season = factor(getSeason(pop)),
  year   = factor(getYear(pop))
) %>%
  left_join(popsMaster %>%
              dplyr::select(pop = popString, location = Location, state = State),
            by = "pop") %>%
  # Rename Deddincgton to Newstead, because both pops are from the same location
  mutate(location = str_replace(location, "Deddington", "Newstead"))
strata(gen195) <- strataDf195

# replace missing values using the mean (note: this reduces the variance)
allFreq195 <- scaleGen(gen195, NA.method = "mean")

# run basic pca without scaling
pca195 <- dudi.pca(allFreq195, cent = FALSE, scale = FALSE, scannf = FALSE, nf = 3)

# check the bar plot of eigenvalues
barplot(pca195$eig[1:50], main = "PCA 195 eigenvalues", col = heat.colors(50))

# plot the pca
myCols195 <- funky(nPop(gen195))
file.path("C:", "UserData", "Kym", "PhD", "RAD2", "pca",
          "PxPCA14pops195indDraft.pdf") %>%
  pdf(height = 7, width = 7)
# you can try altering `fac` arg below to `strata(gen195)$year` or `strata(gen195)$location`
# \ but it doesn;t change the result.   
s.class(
  pca195$li, fac = strata(gen195)$location, col = alpha(myCols195, 0.5), 
  label  = levels(strata(gen195)$location),
  clabel = 0.8,   # cex for pop labels 
  cellipse = 1,   # cex for the ellipse # 1 takes to 50% of the way between points (1.5 is default)
  cstar  = 1,     # cex for the wheel spokes (1 takes it to the points)
  cpoint = 2,     # cex for points,
  grid = FALSE    # remove grid
)
#title("PCA of 195 individuals from 7 locations sampled in both 2014 and 2015")
add.scatter.eig(
  pca195$eig[1:20], nf = 3, xax = 1, yax = 2, posi = "bottomright"
)
dev.off()






# UP TO HERE:
# To do:
# Could try imputing the values ased on allele frequencies in each pop (put a few imputations in a supplementary figure?)
# Read to docuentation - understand the arguments I've used there
# Write the results up.
# Look into bootstrapping the PCA - multivariate techniques book.
# finish revising the code in this script.




# read in genepop file
genind2014 <- file.path(
  "C:", "UserData", "Kym", "PhD", "RAD2", "pgdSpider","Px.1008SNPs.440ind.2014.gen"
  ) %>%
  read.genepop(ncode = 3)

# set grouping factor
pop(genind2014) <- pop(genind2014) %>%
  str_sub(1, 9)

# replace missing values using the mean (note: this reduces the variance)
x <- scaleGen(genind2014, NA.method = "mean")
class(x)

# run pca (without scaling)
pca2014 <- dudi.pca(x, cent = FALSE, scale = FALSE, scannf = FALSE, nf = 3)

# check the barplot of eigenvalues
barplot(pca2014$eig[1:50], main = "PCA 2014 eigenvalues",
        col = heat.colors(50))
#> There are no outstanding PCs.

# plot the pca
myCols <- funky(nPop(genind2014))
file.path("C:", "UserData", "Kym", "PhD", "RAD2", "pca",
          "pcaPx2014DRAFT.pdf") %>%
  pdf(height = 7, width = 7)
s.class(
  pca2014$li, fac = pop(genind2014), col = alpha(myCols, 0.5), 
  label = popNames(levels(pop(genind2014))),
  clabel = 0.5,     # cex for pop labels 
  cellipse = 1, # cex for the ellipse # 1 takes to 50% of the way between points (1.5 is default)
  cstar  = 1,     # cex for the wheel spokes (1 takes it to the points)
  cpoint = 2,     # cex for points,
  grid = FALSE    # remove grid
  )
title("PCA of 31 P. xylostella populations from 2014")
add.scatter.eig(
  pca2014$eig[1:20], nf = 3, xax = 1, yax = 2, posi = "bottomright"
  )
dev.off()

# #########################################################
# 2015 dataset: 28 populations, 402 individuals, 1008 SNPs
# #########################################################

# read in genepop file
genind2015 <- file.path(
  "C:", "UserData", "Kym", "PhD", "RAD2", "pgdSpider","Px.1008SNPs.402ind.2015.gen"
) %>%
  read.genepop(ncode = 3)

# set grouping factor
pop(genind2015) <- pop(genind2015) %>%
  str_sub(1, 9)

# replace missing values using the mean (note: this reduces the variance)
y <- scaleGen(genind2015, NA.method = "mean")

# run pca (without scaling)
pca2015 <- dudi.pca(y, cent = FALSE, scale = FALSE, scannf = FALSE, nf = 3)

# check the barplot of eigenvalues
barplot(pca2015$eig[1:50], main = "PCA 2015 eigenvalues",
        col = heat.colors(50))
#> PCs not prominent.

# plot the pca
myCols15 <- funky(nPop(genind2015))
file.path("C:", "UserData", "Kym", "PhD", "RAD2", "pca",
          "pcaPx2015DRAFT.pdf") %>%
  pdf(height = 7, width = 7)
s.class(
  pca2015$li, fac = pop(genind2015), col = alpha(myCols15, 0.75), 
  label = popNames(levels(pop(genind2015))),
  clabel = 0.3,     # cex for pop labels 
  cellipse = 1, # cex for the ellipse # 1 takes to 50% of the way between points (1.5 is default)
  cstar  = 1,     # cex for the wheel spokes (1 takes it to the points)
  cpoint = 1.5,     # cex for points,
  grid = FALSE    # remove grid
)
title("PCA of 28 P. xylostella populations from 2015")
add.scatter.eig(
  pca2014$eig[1:20], nf = 3, xax = 1, yax = 2, posi = "bottomleft"
)
dev.off()

# Next:
# find way to permute the pca. If can, then can use the method of imputing the missing alleles rather thn use the mean
# try a pca on 2015 data without Southend.



# ####################################
# CODE BELOW for imputing values rather than using the mean (to avoid reducing variance)


# #######################
# pca func
# plotPCA <- function(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.440ind.2014.gen",
#                     code = 3,
#                     name = "Px2014", # a name for this analysis
#                     group = "pop",
#                     na.threshold = 1, # 0-1, proportion of NAs allowed. Loci that equal or exceed this prop in ANY pop are dropped (1 = loci is dropped if ANY pop has all NAs)
#                     scale_labs = 1, 
#                     write_pdf = FALSE, 
#                     write_loadplot = FALSE) {
#   
#   # define function to blacklist SNP loci that equal or exceed a threshold proportion of NAs in ANY pop
#   missLoci <- function(x, na.threshold = 1) {
#     nInd <- nrow(x) * (na.threshold)
#     miss <- apply(x, FUN = function(x) {sum(is.na(x))}, MARGIN = 2)
#     names(miss[which(miss >= nInd)])
#   }
#   # define function to impute missing allele counts without reducing the pop variance (randomly sample from a probability distribution of allele frequencies for each population)
#   imputeVals <- function(x){
#     miss <- is.na(x)
#     p <- c(sum(x[!miss] == 0), 
#            sum(x[!miss] == 1), 
#            sum(x[!miss] == 2)) / sum(!miss)
#     x[miss] <- sample(c(0, 1, 2), size = sum(miss), prob = p, replace = TRUE)
#     x
#   }
#   
#   # read in path to a genepop file or a `genind` object
#   if(class(gp) == "genind") {
#     genind <- gp
#   } else if (is.character(gp)) {
#     genind <- read.genepop(gp, code)
#     }
#   
#   # ###############################################
#   # CHRIS ... you'll need to block out the following code chunk .. just grouping my data into Factors for PCA 
#   
#   # set the grouping factor (correct the pops using gsub first ..genepop infile has first ind name as pop name)
#   if (group == "pop") {
#     pop <- gsub("gilN14swx-10m.40", "gilN14swa-10m.40", pop(genind2014)) %>% # this line a hack to rename one mislabelled Px pop (added 12/2/2017)
#       gsub("-\\d{2}.\\.\\d{2}$", "", .)
#   } else if (group == "host") {
#     pop <- gsub("-\\d{2}.\\.\\d{2}$", "", pop(genind)) %>%
#       str_sub(8, 8) %>%
#       gsub("^v$", "Vegetables", .) %>%
#       gsub("^w$", "Wild hosts", .) %>%
#       gsub("^c$", "Canola", .) %>%
#       gsub("^f$", "Forage", .)
#   } else if (group == "season") {
#     pop <- gsub("-\\d{2}.\\.\\d{2}$", "", pop(genind)) %>%
#       str_sub(7, 7) %>%
#       gsub("^a$", "Autumn", .) %>%
#       gsub("^s$", "Spring", .) 
#   } else if (group == "state") {
#     pop <- gsub("-\\d{2}.\\.\\d{2}$", "", pop(genind)) %>%
#       str_sub(4, 4) %>%
#       gsub("^S$", "SA", .) %>%
#       gsub("^N$", "NSW", .) %>%
#       gsub("^W$", "WA", .) %>%
#       gsub("^V$", "VIC", .) %>%
#       gsub("^T$", "TAS", .) %>%
#       gsub("^Q$", "QLD", .)
#   } else if (group == "year") {
#     pop <- gsub("-\\d{2}.\\.\\d{2}$", "", pop(genind)) %>%
#       str_sub(5, 6) %>%
#       paste0("20", .)
#   } else if (group == "species") {
#     pop <- gsub("-\\d{2}.\\.\\d{2}$", "", pop(genind)) %>%
#       str_sub(9, 9) %>%
#       gsub("^x$", "Pxyl", .) %>%
#       gsub("^a$", "Paus", .) 
#   }
#   
#   # CHRIS, end comment out code chunk
#   # ###################################
#   
#   # assign the pop grouping factor to genind object
#   pop(genind) <- pop
#   npop <- nPop(genind)
#   
#   # drop SNP loci that exceed a threshold proportion of NAs in ANY pop
#   blLoci <- tab(genind) %>% 
#     as.data.frame() %>%
#     split(pop(genind)) %>% # split by population
#     lapply(FUN = missLoci, na.threshold = na.threshold) %>%
#     unlist () %>%
#     gsub("\\.\\d{3}$", "", .) %>% # remove the 3-digit allele suffix, keep SNP name
#     unique()
#   
#   # blacklist any loci with one allele (don't know why there are some in Px/Pa file!)
#   blLoci1All <- lapply(seppop(genind, res.type = "genind"), function(x) {
#     which(nAll(x) < 2) %>%
#       names()
#     }) %>%
#     unlist() %>%
#     as.character() %>%
#     unique()
#   
#   blIdx <- which(locNames(genind) %in% unique(c(blLoci, blLoci1All))) 
#   if (length(blIdx) > 0) {  
#    genind <- genind[, loc = -blIdx]
#   }
#   
#   # for each loci, impute missing allele counts from a probability distribution of allele counts in each pop
#   # (for bi-allelic data, impute one allele, assign the other allele count as (2 - imputedVal))
#   spop <- seppop(genind, res.type = "genind")
#   ptab <- lapply(spop, function(x) {
#     
#     sloc <- seploc(x, res.type = "matrix")
#     lapply(sloc, function(y) {
#       
#       y[, 1] <- imputeVals(y[, 1])
#       y[, 2] <- 2 - y[, 1]
#       y
#       }) %>%
#       Reduce(function(x, y) {cbind(x, y)}, .) # h stack loci matrices
#   }) %>%
#     Reduce(function(x, y) {rbind(x, y)}, .) # v stack pop matrices
#   
# 
#   # recontruct a new genind object from the allele table with imputed missing values
#   genind <- new("genind", tab = ptab)
#   pop(genind) <- pop
#   
#   # perform the PCA on scaled and centred frequencies
#   pca <- dudi.pca(genind, cent = TRUE, scale = TRUE, scannf = FALSE, nf = 3)
#   barplot(pca$eig[1:50], 
#           main = "PCA eigenvalues", col = heat.colors(50))
#   
#   col <- funky(npop)
#   s.class(pca$li, fac = as.factor(pop), col = col, clabel = scale_labs)
#   title(paste0("PCA of ", paste(name, paste0(nPop(genind), "pops"),
#                                 paste0(nInd(genind), "inds"), 
#                                 paste0(nLoc(genind), "loci"), sep = ", ")))
#   add.scatter.eig(pca$eig[1:20], 3, 1, 2)
#   
#   #loadingplot(pca$c1^2)
#   
#   if (write_pdf) {
#     pdf(paste0("pca.", name, ".pdf"))
#     col <- funky(npop)
#     s.class(pca$li, fac = as.factor(pop), col = col, clabel = scale_labs)
#     title(paste0("PCA of ", paste(name, paste0(nPop(genind), "pops"),
#                                   paste0(nInd(genind), "inds"), 
#                                   paste0(nLoc(genind), "loci"), sep = ", ")))
#     add.scatter.eig(pca$eig[1:20], 3, 1, 2)
#     dev.off()
#     
#   }
#   if (write_loadplot) {
#     pdf(paste0("pca.loadplot.", name, ".pdf"))
#     loadingplot(pca$c1^2)
#     dev.off()
#   }
#   
# } # End function
# 
# 
# # ************HERE ******************************************
# # define a function to check 'missingness' of loci in each pop
# popAlleleMiss <- function(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.440ind.2014.gen",
#                           code = 3,
#                           group = "pop",
#                           na.threshold = 1) {
#   
#   # define function to blacklist SNP loci that equal or exceed a threshold proportion of NAs in ANY pop
#   missLoci <- function(x, na.threshold = 1) {
#     nInd <- nrow(x) * (na.threshold)
#     miss <- apply(x, FUN = function(x) {sum(is.na(x))}, MARGIN = 2)
#     names(miss[which(miss >= nInd)])
#   }
# 
#   # read in genepop file
#   genind <- read.genepop(gp, code)
#   
#   # set the grouping factor (correct the pops using gsub first ..genepop infile has first ind name as pop name)
#   if (group == "pop") {
#     pop <- gsub("-\\d{2}.\\.\\d{2}$", "", pop(genind))
#   } else if (group == "host") {
#     pop <- gsub("-\\d{2}.\\.\\d{2}$", "", pop(genind)) %>%
#       str_sub(8, 8) %>%
#       gsub("^v$", "Vegetables", .) %>%
#       gsub("^w$", "Wild hosts", .) %>%
#       gsub("^c$", "Canola", .) %>%
#       gsub("^f$", "Forage", .)
#   } else if (group == "season") {
#     pop <- gsub("-\\d{2}.\\.\\d{2}$", "", pop(genind)) %>%
#       str_sub(7, 7) %>%
#       gsub("^a$", "Autumn", .) %>%
#       gsub("^s$", "Spring", .) 
#   } else if (group == "state") {
#     pop <- gsub("-\\d{2}.\\.\\d{2}$", "", pop(genind)) %>%
#       str_sub(4, 4) %>%
#       gsub("^S$", "SA", .) %>%
#       gsub("^N$", "NSW", .) %>%
#       gsub("^W$", "WA", .) %>%
#       gsub("^V$", "VIC", .) %>%
#       gsub("^T$", "TAS", .) %>%
#       gsub("^Q$", "QLD", .)
#   } else if (group == "year") {
#     pop <- gsub("-\\d{2}.\\.\\d{2}$", "", pop(genind)) %>%
#       str_sub(5, 6) %>%
#       paste0("20", .)
#   } else if (group == "species") {
#     pop <- gsub("-\\d{2}.\\.\\d{2}$", "", pop(genind)) %>%
#       str_sub(9, 9) %>%
#       gsub("^x$", "Pxyl", .) %>%
#       gsub("^a$", "Paus", .)  
#     }
# 
#   # assign the pop grouping factor to genind object
#   pop(genind) <- pop
#   npop <- nPop(genind)
#   
#   # drop SNP loci that exceed a threshold proportion of NAs in ANY pop
#   blLoci <- tab(genind) %>% 
#     as.data.frame() %>%
#     split(pop(genind)) %>% # split by population
#     lapply(FUN = missLoci, na.threshold = na.threshold) %>%
#     unlist () %>%
#     gsub("\\.\\d{3}$", "", .) %>% # remove the 3-digit allele suffix, keep SNP name
#     unique()
#   
#   blIdx <- which(locNames(genind) %in% blLoci) 
#   genind <- genind[, loc = -blIdx]
#   
#   splPop <- tab(genind) %>% 
#     as.data.frame() %>%
#     split(pop(genind))
#   
#   lociMiss <- lapply(splPop, function(df) {
#     res <- apply(df, FUN = function(x) {
#       miss <- sum(is.na(x))
#       miss / length(x)
#     }, MARGIN = 2)
#     data_frame(mean = mean(res),
#                sd = sd(res))
#   }) %>%
#     bind_rows() %>%
#     cbind(data.frame(pop = names(splPop)), .)
#   
#   lociMiss
#   ggplot(lociMiss, aes(x = pop, y = mean)) +
#     geom_bar(stat = "identity") +
#     geom_errorbar(aes(ymin = mean - sd,
#                       ymax = mean + sd)) +
#     ggtitle(paste0("Proportion loci missing per pop: Overall mean missing = ", 
#                    mean(lociMiss$mean))) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# } # End function
# 
# # #####################################################
# # Run analyses
# 
# # *********************
# # analysis 1: only drop Loci where ANY pop has all NAs
# # (keep the same na.threshold for the two function calls below)
# plotPCA(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.440ind.2014.gen",
#         group = "pop",
#         na.threshold = 1,
#         name = "Px2014",
#         scale_labs = 0.6, 
#         write_pdf = TRUE)
# popAlleleMiss(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.440ind.2014.gen",
#               group = "pop",
#               na.threshold = 1)
# 
# # *********************
# # analysis 2: only keep loci where all pops have less than 0.5 NAs
# plotPCA(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.440ind.2014.gen",
#         na.threshold = 0.5,
#         group = "pop",
#         name = "Px2014.na0.5",
#         scale_labs = 0.8, 
#         write_pdf = TRUE)
# popAlleleMiss(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.440ind.2014.gen",
#               na.threshold = 0.5)
# 
# 
# # *********************
# # analysis 3: only keep loci where all pops have less than one third NAs
# plotPCA(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.440ind.2014.gen",
#         na.threshold = 0.33,
#         name = "Px2014.na0.33",
#         scale_labs = 0.8, 
#         write_pdf = TRUE)
# popAlleleMiss(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.440ind.2014.gen",
#               na.threshold = 0.33)
# 
# # ********************
# # group by host, only drop Loci where ANY pop has all NAs
# plotPCA(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.440ind.2014.gen",
#         group = "host",
#         name = "Px2014.host",
#         na.threshold = 1,
#         scale_labs = 0.8, 
#         write_pdf = TRUE)
# popAlleleMiss(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.440ind.2014.gen",
#               na.threshold = 1,
#               group = "host")
# 
# # *********************
# # group by state, only drop Loci where ANY pop has all NAs
# plotPCA(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.440ind.2014.gen",
#         group = "state",
#         name = "Px2014.state",
#         na.threshold = 1,
#         scale_labs = 0.3, 
#         write_pdf = TRUE)
# # Pop Miss Fail ..
# popAlleleMiss(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.440ind.2014.gen",
#               na.threshold = 1,
#               group = "state")
# 
# # *********************
# # group by season, only drop Loci where ANY pop has all NAs
# plotPCA(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.440ind.2014.gen",
#         group = "season",
#         name = "Px2014.season",
#         na.threshold = 1,
#         scale_labs = 0.5, 
#         write_pdf = TRUE)
# popAlleleMiss(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.440ind.2014.gen",
#               na.threshold = 1,
#               group = "season")
# 
# # *******************************
# # 2015 dataset
# 
# # analysis 1: only drop Loci where ANY pop has all NAs
# # (keep the same na.threshold for the two function calls below)
# plotPCA(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.402ind.2015.gen",
#         group = "pop",
#         na.threshold = 1,
#         name = "Px2015",
#         scale_labs = 0.6, 
#         write_pdf = TRUE,
#         write_loadplot = TRUE)
# # FAIL pop ALLELE MISS
# popAlleleMiss(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.402ind.2015.gen",
#               na.threshold = 1)
# # ******************************
# # analysis 2: keep only loci with less than 0.5 NAs in all pops
# # (keep the same na.threshold for the two function calls below)
# plotPCA(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.402ind.2015.gen",
#         na.threshold = 0.5,
#         name = "Px2015.na0.5",
#         scale_labs = 0.8, 
#         write_pdf = TRUE)
# # FAIL pop ALLELE MISS
# popAlleleMiss(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008SNPs.402ind.2015.gen",
#               na.threshold = 0.5)
# 
# # *******************************
# # Make a Px / Pa PCA plot
# 
# # Analysis 1: drop any loci with all NAs in ANY pop
# plotPCA(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/pxpa4pops78ind.gen",
#         group = "pop",
#         na.threshold = 1,
#         name = "PxPa4pops78ind",
#         scale_labs = 0.3, 
#         write_pdf = TRUE)
# popAlleleMiss(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/pxpa4pops78ind.gen",
#               na.threshold = 1)
# 
# # Analysis 2: keep only loci with less than 0.5 missing NAs in all pops
# plotPCA(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/pxpa4pops78ind.gen",
#         group = "pop",
#         na.threshold = 0.5,
#         name = "PxPa4pops78ind.na0.5",
#         scale_labs = 0.5, 
#         write_pdf = TRUE)
# popAlleleMiss(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/pxpa4pops78ind.gen",
#               na.threshold = 0.5,
#               group = "pop")
# 
# # Seperate PCA using only canola and vegetable populations
# plotPCA(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1048SNPs.CanVegPops.gen",
#         group = "pop",
#         na.threshold = 1,
#         name = "PxCanVeg",
#         scale_labs = 0.4, 
#         write_pdf = TRUE)
# 
# plotPCA(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1048SNPs.CanVegPops.gen",
#         group = "host",
#         na.threshold = 1,
#         name = "PxCanVeg.Host",
#         scale_labs = 0.5, 
#         write_pdf = TRUE)
# 
# # ################
# # Plot all Px together
# plotPCA(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.1008.SNPs.842.indv.gen",
#         group = "pop",
#         na.threshold = 1,
#         name = "allPx",
#         scale_labs = 0.5, 
#         write_pdf = TRUE)
# 
# # ################
# # 12/2/2017
# # Plot Esperance and Southend pops agaist other Px and Paus samples
# # 10 pops, 128 ind, 656 SNPs, maxMiss 0.9
# plotPCA(gp = "C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.checkEspSth.656SNPs.maxMiss0.9.gen",
#         group = "pop",
#         na.threshold = 1,
#         name = "Px.checkEspSth",
#         scale_labs = 0.5, 
#         write_pdf = TRUE)
# 
# # subset1 ... drop the 8 Paus samples first
# genind <- read.genepop("C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.checkEspSth.656SNPs.maxMiss0.9.gen",
#                        ncode = 3)
# pop(genind) <- gsub("gilN14swx-10m.40", "gilN14swa-10m.40", pop(genind)) %>% # this line a hack to rename one mislabelled Px pop (added 12/2/2017)
#   gsub("-\\d{2}.\\.\\d{2}$", "", .)
# 
# PausIdx <- indNames(genind) %in% c("gilN14swx-01f.26", "gilN14swx-02f.33", "gilN14swx-03f.33",
#                       "gilN14swx-04f.40", "gilN14swx-06m.26", "gilN14swx-07m.26",
#                       "gilN14swx-08m.33", "gilN14swx-10m.40") %>%
#   which()
# 
# noPaus <- genind[-PausIdx, ] 
# plotPCA(gp = noPaus,
#         group = "pop",
#         na.threshold = 1,
#         name = "Px.checkEspSth.dropPaus",
#         scale_labs = 0.5, 
#         write_pdf = TRUE)
# 
# # subset2 ... drop southend
# indNames(genind) # southend are 113:128
# sthIdx <- 113:128
# noSth <- genind[-sthIdx, ]
# plotPCA(gp = noSth,
#         group = "pop",
#         na.threshold = 1,
#         name = "Px.checkEspSth.dropSth",
#         scale_labs = 0.5, 
#         write_pdf = TRUE)
# 
# # subset3 ... drop esperance
# indNames(genind) # esp are 40:51
# espIdx <- 40:51
# noEsp <- genind[-espIdx, ]
# plotPCA(gp = noEsp,
#         group = "pop",
#         na.threshold = 1,
#         name = "Px.checkEspSth.dropEsp",
#         scale_labs = 0.5, 
#         write_pdf = TRUE)
# 
# 
# 
# # **************************************************
# # Make a neighbour joining tree using the ape package
# library(ape)
# genind <- read.genepop("C:/UserData/Kym/PhD/RAD2/pgdSpider/Px.checkEspSth.656SNPs.maxMiss0.9.gen",
#                        ncode = 3)
# d <- dist(tab(genind))
# tre <- nj(d)
# par(xpd = TRUE)
# # vector of tip labels for the outgroup (rooted tree)
# og <- c("gilN14swx-01f.26", "gilN14swx-02f.33",
#         "gilN14swx-03f.33", "gilN14swx-04f.40",
#         "gilN14swx-06m.26", "gilN14swx-07m.26",
#         "gilN14swx-08m.33", "gilN14swx-10m.40")
# pdf("checkPx.NJ.656SNPs.pdf", height = 15, width = 15)
# plot(root(tre, outgroup = og, node = NULL, resolve.root = TRUE), 
#      type = "phylogram", edge.w = 2, cex = 0.5)
# dev.off()
# #edgelabels(tex = round(tre$edge.length, 1), bg = rgb(0.8, 0.8, 1, 0.8))
# 
# 
# 
# 