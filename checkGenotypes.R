# ####################################################################################
# Check the genotypes for 15 individuals RAD sequenced in duplicate for the RAD2 study
# KD Perry, 4/6/2018
# 
# Notes: 
# 15 individuals were sequenced in duplicate.
# This script will run a PCA to check the same indviduals group together.
# We might expect some differences between samples due to different in read depth leading to missing data.
# ########################################################################################################

library(tidyverse)
library(adegenet)
library(readxl)

# =======================================================================================================
# Create a list of individuals sequenced in duplicate to check genotypes and write a file with the names 
# =======================================================================================================

# Read in the list of the duplicate individuals.
# These samples should have the names but different library numbers( `.nn` suffix) and a `d` suffix.
# These are the correct smples names (same as `old` names for these samples - no incorrect genotypes).
dupInd <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "dupInd15.txt") %>%
  readLines()

# Now pull out 15 original sample names (this code copied from `writeSampleMetaDataForSRASUbmission.R`)
# There are 960 samples (945 unique samples) in the data. We duplicate sequenced 15 individuals to check genotypes.
samplesMasterRAD2 <- file.path(RAD2dir, "genotypesMasterRAD2.xlsx") %>%
  read_excel(sheet = "RAD2-samples-master", range = "B1:F961") %>%
  dplyr::rename(sampleNamesExcel = `sampleNames(RADpools)`) %>% # original names made in Excel (not fixed width, and 13 wrong species genotypes)
  # add the DNA plate and position coordinates
  mutate(plate = str_replace(Plate, "JF_1", "JF01"), # Jacinta Frater's plates
         plate = str_replace(plate, "JF_2", "JF02"),
         plate = str_pad(plate, width = 2, side = "left", pad = 0),
         plate = paste0("KP", plate), # The rest are Kym Perry's plates
         plate = str_replace(plate, "^KPJF", "JF"),
         well  = str_pad(Well, width = 3, side = "left", pad = 0),
         DNA_plate_position = paste(plate, well, sep = "_"),
         # now add a column with indNames in the new format, for joining
         sampleNamesOld = renameSamplesRAD2(sampleNamesExcel),    # add sample names in new fixed width format but 13 wrong genotypes
         sampleNamesCorrect = renamePxToPaus(sampleNamesOld)) %>% # the final sample names with 13 corrected species genotypes
  dplyr::select(sampleNamesExcel, sampleNamesOld, sampleNamesCorrect, DNA_plate_position)
# Note: there are 960 samples but 945 unique plate_well combos. \

removeLibSuffix <- function(x){
  x %>% str_sub(1, 13) # whole name except suffix
}
dupSamples30 <- samplesMasterRAD2 %>%
  filter(removeLibSuffix(sampleNamesOld) %in% removeLibSuffix(dupInd))
stopifnot(identical(dupSamples30$sampleNamesOld, dupSamples30$sampleNamesCorrect)) 
#> `old` and `corrected` names are the same
#> use these names to extract the samples from a vcf file

dupIndNames30 <- dupSamples30 %>%
  dplyr::select(sampleNamesCorrect) %>%
  unlist() %>% 
  as.character()

# also write a file with the samples list for use in other programs (e.g. filtering vcf)
dupIndNames30 %>%
  write.table("dupInd30.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

# also write a populations file for vcf to genepop format conversion in pgdspider
# for this, we'll define the `populations` as sample pairs (15 pops with 2 samples in each).
data_frame(ind = sort(dupIndNames30)) %>%
  mutate(pop = str_sub(ind, 1, 13)) %>% # Use the sample name as the pop name 
  write.table("dupInd30PopsFile.txt",
              quote = FALSE, row.names = FALSE, col.names = FALSE) 

# ==================================
# Now run a PCA with the duplicates
# ==================================

# Read in genepop file with 30 samples and 1473 SNPs (created from hard-filtered vcf and converted in pgdspiders - see RAD2 notes line 7474)
gen30 <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "pgdSpider", "Px.30dupInd.1473SNPs.gen") %>%
  adegenet::read.genepop(ncode = 3)

# Set the hierchary in the strata slot of the genind object (use sample pair name as pop name)
strataDf30 <- data_frame(
  pop = pop(gen30) %>% removeLibSuffix() %>% factor()
  )
strata(gen30) <- strataDf30

# replace missing values using the mean (note: this reduces the variance)
allFreq30 <- scaleGen(gen30, NA.method = "mean")

# run basic pca without scaling
pca30 <- dudi.pca(allFreq30, cent = FALSE, scale = FALSE, scannf = FALSE, nf = 3)

# check the bar plot of eigenvalues
barplot(pca30$eig[1:50], main = "PCA 30 ind eigenvalues", col = heat.colors(50))

# plot the pca
myCols30 <- funky(nPop(gen30))
file.path("C:", "UserData", "Kym", "PhD", "RAD2", "Px30DupIndPCA.pdf") %>%
  pdf(height = 7, width = 7)
s.class(
  pca30$li, fac = strata(gen30)$pop, col = alpha(myCols30, 0.5), 
  label  = levels(strata(gen30)$pop),
  clabel = 0.2,   # cex for pop labels 
  cellipse = 1,   # cex for the ellipse # 1 takes to 50% of the way between points (1.5 is default)
  cstar  = 1,     # cex for the wheel spokes (1 takes it to the points)
  cpoint = 2,     # cex for points,
  grid = FALSE    # remove grid
)
title("PCA of 15 P. xylostella individuals RAD sequenced in duplicate (1473 SNPs)",
      cex.main = 1)
add.scatter.eig(
  pca30$eig[1:20], nf = 3, xax = 1, yax = 2, posi = "bottomright", ratio = 0.15
)
dev.off()

# Genotypes match nicely!

# End script.
# ############################################################################







  