# #######################################################################################
# Write sample meta-data for submitting sequences to the NCBI Sequence Read Archive (SRA)
# KD Perry, 27/5/2018

# =================
# This script does:
# =================

# UPDATE:
# This script:
# (1) Writes metadata for BioSamples (ie one row per individual) using fields specified 
#  \  in the NCBI invertebrate biosample template
# (2) Writes SRA metadata with details of each sequencing run (experiment) for each sample
# (3) Writes a master metadata file for our own purposes showing 
#  \  the naming history for each sample (see below).

# ===============================
# Notes on sample naming history:
# ===============================
# The sample naming scheme was designed to enable sample meta data to be extracted from each name string.
# The naming convention is fully explained in the footer of this script. ******* STILL TO ADD *********
# Across the entire RADseq pipeline, there are 3 `versions` of sample names (don't panic - everything is documented)
# In fact there are only 2 versions, but genotypes for 13 samples were later corrected, necessitating changes to the name strings for 13 samples.
# All version use the same metadata encoding scheme and sample numbers, but the differences are as follows:
# Version 1: `sampleNamesV1Excel`: sample names were originally assigned by joining column data in Excel. Hence, numbers and strings were not fixed width 
# Version 2: `sampleNamesV2OriginalGenotypes`: to enable convenient metadata extraction and sub-sampling from sample names (eg using grep etc) \
# \ names were converted to a new convention with fixed width format for all metadata characters (see footer).
# Version 3: `sampleNamesV3Final`: after COI sequencing, 13 genotypes were corrected from P. xylostella to P. australiana.
# \ This meant the character encoding Plutella species needed changing from `x` to `a` in the name strings for 13 samples.

# Inclusion of additional samples from a prevous sequencing effort (RAD1):
# An extra 16 samples from two populations (Calca SA, Picnic Beach SA) were included in the RAD2 study from previous a sequencing effort at ACRF Frome Road in 2014 (`RAD1` samples). \
# \ These 16 samples were originally named under a different scheme (eg Calca-1, Picnic-1 etc).
# \ Samples were first renamed to conform to the `sampleNamesV1EXcel` scheme (see commands in Excel file in local dir \preProcessing\`renameSamplesRAD1.xlsx`)
# \ Then they were renamed to the RAD2 scheme `sampleNamesV2OriginalGenotypes`, and then `sampleNamesV3Final` (though none of these genotypes needed correction)
# \ In their name strings, a library number was assigned starting from 70, to differentiate them from the RAD2 library numbers (which went up to 60)

# ========================================
# When to use which sample name `version`:
# ========================================
# The v1Excel names were used in the genotypesMasterRAD2 spreadsheet contain the 960 samples sequenced at AGRF (Ade/Mel) in 2017 (`RAD2` samples).
# All the raw fastq files (_1.fastq, _2.fastq) had these original Excel sample names.
# After demultiplexing, all raw fastq files were renamed to the fixed width `sampleNamesV2OriginalGenotypes` format \
# \ using bash commands made with the function `renameSamplesRAD2()` (the function definition is included in this script below)
# Then, downstream alignment to the reference genome, genotyping and hard-filtering steps were performed with these V2 names.
# Finally, names with corrected genotypes were used for downstream population genetic analyses, including grouping samples and interpreting results. 
# Input files for use in other programs (eg GENEPOP, STRUCTURE) population were often converted from VCF format using PGDSpider (or sometimes scripts),
# \ which means these input files will sometimes have the `sampleNamesV2OriginalGenotypes` names.
# \ However, great care was taken to ensure populations were correctly grouped according to species using custom R functions to update the names as required.
# The custom R functions `renamePxToPaus()` and `renamePausBackToPx()` (included in this script below) \
# \ were used extensively during pop gen analyses to ensure forward and backward compatibility when working with lists of sample names.
# \ Many lists of sample names were written in local directory `RAD2` using R script `writeSampleNamesRAD2.R` and used to subset files.
# \ These lists contain the `sampleNamesV2OriginalGenotypes` names (because lists were often used to filter and subset VCF files).
# \ In many cases, conversion of sample names to `sampleNamesV3Final` format was the final step before reporting results.

# In summary: 
# (i)   There should never be a need to use the V1Excel names again, other than to match samples with their metadata \
#     \ in the original file `genotypeaMasterRAD2.xlsx`, sheet = "RAD2-samples-master". 
# (ii)  When retrieving samples or their names from fastq, BAM or VCF files, `sampleNamesV2OriginalGenotypes` names must be used.
# (iii) Input file for many downstream population genetic analysis program will also have the version 2 names, but the sample/population grouping 
#     \ according to Plutella species will be correct. The function renamePxToPaus() was used prior to creating the groups.  
# (iv)  Conversion of sample names in sample name lists can be performed using helper functions renamePxToPaus() and renamePausBackToPx() as required.  

# ============================
# For SRA submission, we will:
# ============================
# Use the final V3 sample names with corrected genotypes in the sample metadata.
# Prior to upload, we will rename the raw fastq files with their correct species names

# ==========================
# For our purposes, we will: 
# ==========================
# Retain all names versions on our master metadata file to ensure backward compatibility.

# #######################################################################################

library(tidyverse)
library(readxl)

# ##################
# Helper functions #
# ##################

# function used to rename filenames to the new fixed width format
# function copied from R script `~preProcessing/renameSamplesRAD2.R`:
# see that script for a full explanation of the naming convention.
renameSamplesRAD2 <- function(filename) {
  basename <- str_extract(filename, "^.*-\\d{1,2}(\\.99)?[mfu]{1}")
  suffix <- mapply(str_split, string = filename, pattern = basename) %>%
    lapply(magrittr::extract, 2) %>% 
    unlist #%>%
  #gsub("(\\.fastq$)$", ".fastq.gz", .) # if needed append the .gz to names
  popName <- str_sub(basename, 1, 9)
  smNum <- str_extract(basename, "-\\d{1,2}") %>%
    gsub("[-]", "", .) %>% 
    str_pad(width = 2, side = "left", pad = "0")
  gender <- str_extract(basename, "[mfu]{1}$")
  lib <- str_extract(basename, "\\.\\d{1,2}-") %>%
    gsub("[.-]", "", .) %>% 
    str_pad(width = 2, side = "left", pad = 0)
  dup <- gsub("\\.99[mfu]$", "d", basename) %>% str_extract(., "d$")
  renamed <- paste0(popName, "-", smNum, gender, ".", lib, 
                    dup, suffix) %>% gsub("NA", "", .)
  renamed
}

# function to rename 13 samples from Px to Paus (rename samples for correct population grouping)
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

# function to convert correct (Paus) back to incorrect (Px) names
# This is for backward compatility, where names need to be extracted from fastqs/vcfs
renamePausBackToPx <- function(x) {
  gsub("gilN14swa-01f.26", "gilN14swx-01f.26", x) %>%
    gsub("gilN14swa-06m.26", "gilN14swx-06m.26", .) %>%
    gsub("gilN14swa-07m.26", "gilN14swx-07m.26", .) %>%
    gsub("gilN14swa-08m.33", "gilN14swx-08m.33", .) %>%
    gsub("gilN14swa-02f.33", "gilN14swx-02f.33", .) %>%
    gsub("gilN14swa-03f.33", "gilN14swx-03f.33", .) %>%
    gsub("gilN14swa-04f.40", "gilN14swx-04f.40", .) %>%
    gsub("gilN14swa-10m.40", "gilN14swx-10m.40", .) %>%
    gsub("boyW14sca-05m.03", "boyW14scx-05m.03", .) %>%
    gsub("espW14sca-01f.04", "espW14scx-01f.04", .) %>%
    gsub("espW14sca-06m.11", "espW14scx-06m.11", .) %>%
    gsub("espW14sca-07m.11", "espW14scx-07m.11", .) %>%
    gsub("espW14sca-04m.04", "espW14scx-04m.04", .)
}

getPop <- function(x){
  x %>% str_sub(1, 9)
} 
getSpecies <- function(x){
  species <- c("^x$" = "Plutella xylostella", "^a$" = "Plutella australiana")
  x %>% 
    str_sub(9, 9) %>% 
    str_replace_all(species)
}
getState <- function(x){
  states <- c("^S$" = "South Australia", "^W$" = "Western Australia",
              "^V$" = "Victoria", "^N$" = "New South Wales", "^T$" = "Tasmania",
              "^Q$" = "Queensland")
  x %>% 
    str_sub(4, 4) %>% 
    str_replace_all(states)
}
getSex <- function(x){
  sexes <- c("^m$" = "Male", "^f$" = "Female", "^u$" = "Unknown")
  x %>%
    str_sub(13, 13) %>%
    str_replace_all(sexes)
}
getHost <- function(x) {
  hosts <- c("^v$" = "Vegetables", "^w$" = "Wild", "^c$" = "Canola", "^f$" = "Forage")
  x %>%
    str_sub(8, 8) %>%
    str_replace_all(hostType)
}
getRADseqLibrary <- function(x) {
  x %>%
    str_extract("\\d{2}d?$") %>%
    str_replace("d$", "") # strip the `d` suffix from 15 technical duplicate sample names
}

# #################################################################################################################
# Summarise meta-data for the 887 P. xylostella and P. australana samples included in the RAD2 sequencing studies # 
# #################################################################################################################

# Create a meta-data master file of all samples RAD sequenced: 960 sequenced in RAD2 (AGRF, 2017) and 16 samples in RAD1 (ACRF, 2014)
# The sequenced samples include:
# (1) 887 P. xylostella and P. australiana samples that were included in two final published manuscripts, comprising:
#   \ (a) 833 P. xylostella individuals included in the RAD2 P. xylostella manauscript
#   \ (b) 99  P. australiana (52) and P. xylostella (47) individuals included in the RAD2 cryptic species manuscript (Perry et al 2018, BMC Evol. Biol.)
# note: The P.x samples in both studies overlap, but there are 2 extra Gilgandra samples included in the Cryptic species study (hence n = 887, not 885).
# (2) The xxx samples that were sequenced but later excluded from the studies.
#   \ **** STILL TO ADD ***()

# ==============================================================================
# read in list of samples included in the final analyses for 2 published studies
# ==============================================================================

# These lists were written in separate R script `writeSampleNames.R`
# note: lists use version 2 sample names with some incorrect species genotypes 
# This is correct as sample name lists required backward compatilibility to extract samples from VCFs etc)
RAD2path   <- file.path("C:", "UserData", "Kym", "PhD", "RAD2") 
Px833indV2names   <- file.path(RAD2path, "Px833ind.txt") %>% readLines()   # 833 P. xylostella samples
PxPa99indV2names  <- file.path(RAD2path, "pxpa.5pops.99ind.samples.txt") %>% readLines() # 47 Px and 52 P. australiana samples from the cryptic species study. (45 of the 47 Px samples overlap with the Px833ind)
PxPa887indV2names <- c(Px833indV2names, PxPa99indV2names) %>% unique() # 833 Px samples + 2 Px samples from Gilgandra + 52 Pa samples
PxPa16indRAD1


# ===============================================================================================
# create a dataframe with the metadata for all sequenced samples (included and excluded samples)
# ===============================================================================================

# read in metadata for population collections (which has new sample names)
# Note: the file "popsMasterRAD2.csv" was written in a separate script: `writeSampleNames.R`.
popsMasterRAD2 <- file.path(RAD2dir, "popsMasterRAD2.csv") %>%
  read_csv()

# 16 samples from RAD1 study were included in the RAD2 study (but only 12 kept after filtering)
# Get the original RAD1 sample names, and the new sample names in fixed width format 
samplesMasterRAD1 <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "preProcessing", "renameSamplesRAD1.xlsx") %>%
  read_excel() %>%
  dplyr::select(name, renamed) %>%
  mutate(name    = str_replace(name, "_\\d\\.fastq", ""),         # names as per a scheme from a previous sequencing study (RAD1)
         renamed = str_replace(renamed, "_\\d\\.fastq", "")) %>%  # renamed to the RAD2 scheme, as per `sampleNamesV1Excel` format
  dplyr::select(sampleNamesRAD1Study = name, sampleNamesV1Excel = renamed) %>%
  distinct(sampleNamesRAD1Study, .keep_all = TRUE) %>%
  mutate(sampleNamesV2OriginalGenotypes = renameSamplesRAD2(sampleNamesV1Excel), # converted to the new fixed width format 
         sampleNamesV3Final = renamePxToPaus(sampleNamesV2OriginalGenotypes),    # note: V2 and V3 names are the same as no genotypes needed correcting
         DNA_plate_position = "single tubes",
         RADseqLibrary = getRADseqLibrary(sampleNamesV3Final)) # these samples were prepared in single tubes, but add values for `DNA_plate_position` var. for joining 

# read in the RAD2 genotypes data for each individual (has the old sample names in the original format)
# There are 960 samples (945 unique samples) in the data. We duplicate sequenced 15 individuals to check genotypes.
samplesMasterRAD2 <- file.path(RAD2dir, "genotypesMasterRAD2.xlsx") %>%
  read_excel(sheet = "RAD2-samples-master", range = "B1:F961") %>%
  dplyr::rename(sampleNamesV1Excel = `sampleNames(RADpools)`) %>% # original names made in Excel (not fixed width, and 13 wrong species genotypes)
  # add the DNA plate and position coordinates
  mutate(plate = str_replace(Plate, "JF_1", "JF01"), # Jacinta Frater's plates
         plate = str_replace(plate, "JF_2", "JF02"),
         plate = str_pad(plate, width = 2, side = "left", pad = 0),
         plate = paste0("KP", plate), # The rest are Kym Perry's plates
         plate = str_replace(plate, "^KPJF", "JF"),
         well  = str_pad(Well, width = 3, side = "left", pad = 0),
         DNA_plate_position = paste(plate, well, sep = "_"),
         # now add a column with indNames in the new format, for joining
         sampleNamesV2OriginalGenotypes = renameSamplesRAD2(sampleNamesV1Excel),  # add sample names in new fixed width format but 13 wrong genotypes
         sampleNamesV3Final = renamePxToPaus(sampleNamesV2OriginalGenotypes),
         RADseqLibrary = getRADseqLibrary(sampleNamesV3Final)) %>% # the final sample names with 13 corrected species genotypes
  dplyr::select(sampleNamesV1Excel, sampleNamesV2OriginalGenotypes, sampleNamesV3Final, DNA_plate_position, RADseqLibrary)


# =============================================================================================================
# create a data_frame showing which studies each `included` sample was used in (for joining to the metadata df)
# =============================================================================================================

# The codes for each study are:
# RAD2PxPa99ind: A population genetic study of 47 Px and 52 P. australiana samples: `Cryptic Plutella species` (Perry et al 2018 BMC Evol Biol) 
# RAD2Px833ind:  A population genetic study of 833 P. xylostella samples from Australia 
# RAD1PxPa16ind: RAD1 study which had 69 Px and 3 P. australiana samples. 16 were added to the RAD2 study, and 12 retained for the analysis. \
# \ RAD1 study was published in our paper for the 7th International DBM workshop: `Genome-wide SNP discovery`  (Perry et al 2017, pre-print online) 
PxPa99indDf <- data_frame(sampleNamesV2OriginalGenotypes =  PxPa99indV2names,
                       RAD2PxPa99ind = "RAD2PxPa99ind")
Px833indDf  <- data_frame(sampleNamesV2OriginalGenotypes =  Px833indV2names,
                       RAD2Px833ind  = "RAD2Px833ind")
# get the 12 (of 16) RAD1 samples that were included in the RAD2 studies
PxPa12ind  <- samplesMasterRAD1$sampleNamesV2OriginalGenotypes[samplesMasterRAD1$sampleNamesV2OriginalGenotypes %in% PxPa887indV2names]
PxPa12indDf <- data_frame(sampleNamesV2OriginalGenotypes = PxPa12ind, 
                          RAD1PxPa12ind = "RAD1PxPa12ind") 

# create a variable summarising all studies that each samples was included in 
includedSamplesStudiesDf <- Px833indDf %>%
  full_join(PxPa99indDf, by = "sampleNamesV2OriginalGenotypes") %>%
  full_join(PxPa12indDf, by = "sampleNamesV2OriginalGenotypes") %>%
  mutate(studiesIncludedIn = paste(RAD2Px833ind, RAD2PxPa99ind, RAD1PxPa12ind, sep = ", ") %>%
           str_replace_all(" NA,?", "") %>% 
           str_replace_all(",$", "") %>%
           str_replace_all("NA, ", "")) %>%
  dplyr::select(sampleNamesV2OriginalGenotypes, studiesIncludedIn)
  
  
# ====================================================================================
# finally, create and write a master file with the metadata for all sequenced samples  
# ====================================================================================

samplesMaster976ind <- samplesMasterRAD2 %>%
  full_join(samplesMasterRAD1) %>%
  dplyr::select(sampleNamesRAD1Study, everything()) %>% # put the RAD1 column first (chrono order)
  mutate(popString = getPop(sampleNamesV3Final)) %>%
  left_join(popsMasterRAD2, by = "popString") %>%
  left_join(includedSamplesStudiesDf)


samplesMaster976ind %>%
  write.table("test.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
# write a fixed width file
colWidths <- sapply(samplesMaster976ind, nchar) %>% 
  data.frame() %>% 
  sapply(max, na.rm = TRUE)

fixedWidth <- function(x, maxWidth){
  if(is.numeric(x)){
    x <- str_pad(x, width = maxWidth, side = "left", pad = 0)
  }
  x %>%
    format(width = maxWidth + 2, scientific = FALSE)
}
samplesMaster976ind %>%
  purrr::map2_dfr(.x = ., .y = colWidths, .f = fixedWidth) %>%
  write.table("samplesRAD2MetadataMaster976ind.txt", sep = "\t",
              quote = FALSE, col.names = TRUE, row.names = FALSE)

# UP TO HERE!!!

# create the metadata master file for the 887 samples
# POSSIBLY: Ad the original EXcel names
metaData887Ind <- samplesMasterRAD2 %>%
  full_join(samplesMasterRAD1) %>%
  filter(sampleNamesOld %in% sampleNamesOldPxPa887ind) # just keep the samples used in RAD2 studies


metaDataMaster887 <- data_frame(
  sampleNamesOld = sampleNamesOldPxPa887ind
  ) %>% # need all samples, not just the RAD2 samplesMaster (that's why can't pipe from sampleMaster object)
  mutate(sampleNamesCorrect = renamePxToPaus(sampleNamesOld),
         popString = getPop(sampleNamesCorrect)) %>%
  # add the population metadata
  left_join(samplesMasterRAD2) %>% # doesn't have the RAD1 sample
  left_join(popsMaster, by = "popString") %>% 
  ### THERE ARE BLOODY NAs in teh data here .... WHY???!!!
  mutate(checkNames = renameSamplesRAD2(sampleNamesExcel))


# ==================================================
# List sample names with sequences for upload to SRA
# ==================================================

RAD2dir <- file.path("C:", "UserData", "Kym", "PhD", "RAD2") 
Px833ind  <- file.path(RAD2dir, "Px833ind.txt") %>% readLines() # 833 P. xylostella samples from RAD2 study
PxPa99ind <- file.path(RAD2dir, "pxpa.5pops.99ind.samples.txt") %>% readLines() # 99 Px and Pa from the RAD2 cryptic species study. (overlaps with the Px833ind)
# sample names have some incorrect genotypes (but the lists were used for backward compatilibility to extract early files)



# combine the sample names and update them to the correct names
sampleNamesOldPxPa887ind <- c(Px833ind, PxPa99ind) %>% unique() 


# ##############################################################################################
# Write BioSample metadata for submission of sequences to the NCBI Sequence Read Archive (SRA) #
# ##############################################################################################


# ======================================================================================
# Create a master data_frame with the fields (and order) matching the Biosample Template
# ======================================================================================

SRAdir <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "SRAsubmission")
biosampleTemplate <- file.path(SRAdir, "Invertebrate.1.0.biosampleTemplate.xlsx") %>%
  read_excel(range = "A13:Y13", col_names = TRUE, col_types = "text") %>%
  set_names(str_replace(names(.), "\\*", ""))
#> need the names just to match the order of columns in my output file

metaDataMaster887 %>%
  filter(is.na(sampleNamesExcel))
  
  # these samples are in the samples list but not in the samplesMaster
  #> > oldSampleNamesPxPa887ind[!oldSampleNamesPxPa887ind %in% samplesMaster$sampleNamesOld]
  [1] "calS14awx-14u.67" "calS14awx-15u.68" "calS14awx-16u.69" "calS14awx-17u.71" "picS14awx-03u.67" "picS14awx-04u.68"
  [7] "picS14awx-05u.69" "picS14awx-06u.70" "picS14awx-07u.71" "picS14awx-08u.72" "calS14awa-19u.70" "calS14awa-20u.72"

# Create the metadata master file for all the samples  
metadata887Df <- data_frame(indNameWrongGenotypes = PxPa887ind) %>%
  mutate(indNameCorrected = renamePxToPaus(indNameWrongGenotypes),# use the new correct species genotypes
         popString = getPop(indNameCorrected)) %>%
  left_join(popsMaster, by = "popString") %>%
  mutate(species = getSpecies(indNameCorrected),
         # Locations: need to add Town names somewhere
         state   = getState(indNameCorrected),
         bioproject_accession = "PRJNA471964",
         isolate = "Not applicable",
         breed   = "Not applicable",
         geo_loc_name = paste("Australia", state, sep = ":"),
         tissue = "Whole body",
         biomaterial_provider = paste(
           "Kym Perry, School of Agriculture, Food and Wine, University of Adelaide, Adelaide, South Australia, 5005"
           ),
         env_biome = "Agricultural production landscape",
         lat_lon   = paste(paste(round(Latitude,  2), "N"),
                           paste(round(Longitude, 2), "W")),
         sex = getSex(indNameCorrected),
         collected_by = paste(collector, "Kym Perry", sep = ", ") %>% # Add KP as a collector for all pops
                    str_replace("Kym Perry, Kym Perry", "Kym Perry")
         ) %>%
  # Add origcollector, collectionDate, townNames, hostSpecies, stageCollected
  left_join(genotypesMaster)

# now format the metadata with columns in the same order as the bioSampleTemplate
# for easy cuttig and pasting the data into the template 

# UP TO HERE: FIX YOUR NAs!

biosample887Df <- metadata887Df %>%
  mutate(
    #sample_title = paste(paste(species, collapse = "_"), "RADseq", "Sample", sep = "_"), FIX!
    isolation_source = paste(host, "plants")
  ) %>%
  dplyr::select(
    popString,
    sample_name  = indNameCorrected,
    sample_title, # ADDED
    bioproject_accession, 
    organism = species,
    isolate, breed,
    host = hostSpecies,
    isolation_source, # ADDED
    collection_date = collectionDate,
    #geo_loc_name = ADD paste country:Town State,
    tissue, 
    # age, altitude, ADD blanks
    biomaterial_provider, collected_by, 
    # depth, ADD blanks
    dev_stage = stageCollected, env_biome, 
    # identified_by ADD blanks,
    lat_lon, sex,
    # specimen_voucher, temp ADD blanks,
    # description ADD a description
    DNA_plate_position) %>%
  arrange(str_sub(popString, 1, 8), # the location without the species! (species is chr 9 in the string)
          desc(organism), sample_name) %>% # to get the P. xylostella samples in order
  dplyr::select(-popString) 

View(biosample887Df)

# UP TO HERE
# TO DO:
# Add sample title? Plutella_xylostella_RADseq_Sample
# Check all the fields.
# write csv with teh fields in the correct order as per the biosample template.
  










