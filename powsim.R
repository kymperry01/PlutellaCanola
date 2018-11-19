# #########################################################
# Produce a LaTex table of powsim analysis for publication
# KD Perry, 20/5/2018
# 
# POWSIM Ryman (2006) is a program that uses simulations to calculate 
# the statistical power of a set of loci to detect true population structure 
# for pre-defined Fst values. 
# 
# This script does:
# 1. Generates input files for POWSIM (Powsim is run separately, outside this script)
# 2. Read in a summary of results and output code for a LaTeX table for publication
# #########################################################

library(tidyverse)
library(readxl)
library(adegenet)
library(reshape2)
library(xtable)


# ###########################
# set up powsim input files #
# ###########################

# =========================================================
# define function to convert genepop to powsim input format
# =========================================================
genepopToPowsim <- function(
  genepopFilePath = file.path("C:", "UserData", "Kym", "PhD", "RAD2", "pgdSpider"),
  nCode = 3, # whether alleles in genepop file are encoded with 2 or 3 integers
  outDir = file.path("C:", "UserData", "Kym", "PhD", "RAD2", "powsim"),
  outName = "addName",
  dataDescription = "P. xylostella SNP data", # short (~3 word) description for the data for powsim header line 2."
  numLoci, #  
  numPops, # number of pops used for the powsim simulation (max 30 currently)
  effectivePopSize = 500,
  gensDrift = 10,
  sampleSizes, # an integer vector of length numPops, with the sample sizes for each pop.
  numSimRuns,
  fisherParams = c(1000, 100, 1000), # burnins, demems, per run
  runChiFisher = c(1, 1),
  eraseDetailedOut = 0
){
  # checks
  stopifnot(file.exists(genepopFilePath))
  stopifnot(file.exists(outDir))
  stopifnot(length(sampleSizes) == numPops)
  stopifnot(length(fisherParams) == 3)
  
  # read in the genepop file, convert to allele frequencies
  alleleFreq <- adegenet::read.genepop(genepopFilePath, ncode = nCode) %>%
    as_tibble() %>%
    summarise_all(sum, na.rm = TRUE) %>%
    t() %>%
    as.data.frame() %>% # as_tibble drops the row.names()
    set_names("alleleCount") %>%
    mutate(locusAllele = row.names(.),
           locus  = str_split(locusAllele, "\\.") %>%
             map_chr(magrittr::extract2, 1),
           allele = str_split(locusAllele, "\\.") %>%
             map_chr(magrittr::extract2, 2)) %>%
    split(.$locus) %>%
    lapply(function(x) {
      x %>%
        mutate(alleleFreq = alleleCount / sum(alleleCount),
               alleleNum  = 1:length(alleleCount),
               allelesPerLocus = length(alleleCount))
    }) %>%
    bind_rows() %>%
    dplyr::arrange(as.numeric(str_replace(locus, "SNP_", ""))) %>%
    dplyr::select(locus, allele, alleleNum, allelesPerLocus, alleleCount, alleleFreq)
  
  # convert to powsim format
  powFreq <- alleleFreq %>%
    reshape2::dcast(locus + allelesPerLocus ~ alleleNum, value.var = "alleleFreq") %>%
    dplyr::rename(numAlleles = allelesPerLocus) %>%
    dplyr::arrange(as.numeric(str_replace(locus, "SNP_", "")))
  
  # write powsim input file to directory
  params <- list(# named for human readability only
    l1  = paste("POWSIM input file generated on", Sys.time(), 
                "using R script `powsimToGenepop.R` by Kym D. Perry 2018"),
    l2  = paste(c(dataDescription,  "with", numLoci, "loci and", numPops, "populations"), collapse = " "),
    l3  = paste(c(fisherParams, "Iteration/permutation factors for Fisher's exact test"), collapse = " "),
    l4  = paste(c(runChiFisher, "Run chi-square and/or Fisher's exact test (1/0)"), collapse = " "),
    l5  = paste(c(eraseDetailedOut, "Erase detailed output (1/0)"), collapse = " "),
    l6  = "",
    l7  = paste(c(numLoci, "Loci"), collapse = " "),
    l8  = paste(c(numPops, "Populations"), collapse = " "),
    l9  = paste(c(effectivePopSize, "Effective population size when drifting apart"), collapse = " "),
    l10 = paste(c(gensDrift, "Generations of drift (t)"), collapse = " "),
    l11 = paste(c(sampleSizes, "Sample sizes after drift"), collapse = " "), # an integer vector of sample sizes for each sub-population
    l12 = paste(c(numSimRuns, "Number of simulation runs"), collapse = " "),
    l13 = "",
    l14 = "Number of allele frequencies for each locus in base population",
    powFreq %>% dplyr::select(-locus)
  )
  
  params %>% 
    lapply(write.table, file.path(outDir, outName), # outSuffix designed to identify param sets ... manually delete suffix prior to running powsim
           append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
} # end function

# ============================
# Set up parameters for powsim
# ============================

globPath <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "pgdSpider")
outPath  <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "powsim")

# get sample sizes for the datasets for use in the powsim analysis.
px434SampleSizes <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "px434SampleSize.txt") %>%
  read.table() %>%
  set_names(c("pop", "n"))

px399SampleSizes <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "px399SampleSize.txt") %>%
  read.table() %>%
  set_names(c("pop", "n"))

# function to calculate Fst for given effective pop size and generations drift (Nei 1987)
# use this to find Ne, t to give target Fst value
fst <- function(Ne, t){
  x <- 1 - (1 - 1 / (2 * Ne)) ^ t
  round(x, 5)
}

# P. xylostella is a species with large effective population sizes. 
# 5000 is probably a reasonable assumption. Other invert papers have used 2000, 4000 fr Ne.
fst(5000, 10)  #> fst = 0.001. Use 10 gens
fst(5000, 50)  #> fst = 0.005. Use 50 gens

# =======================================================
# Write powsim input files (- and then run powim locally)
# =======================================================

# px 31 pops, 434ind, 1032SNPs, target fst = 0.001
fst(5000, 10)  #run 10 gens drift with Ne 5000 
genepopToPowsim(genepopFilePath = file.path(globPath, "Px.31pops.434ind.1032SNPs.gen"),
                outDir  = file.path(outPath, "px434.1032SNPs"),
                outName = "powsim.in.fst001.px434",
                dataDescription = "P.xylostella SNP data",
                numLoci = 1032, 
                numPops = 30, 
                effectivePopSize = 5000,
                gensDrift = 10,
                sampleSizes = px434SampleSizes$n[-31],
                numSimRuns = 1000,
                eraseDetailedOut = 1)

# px 31 pops, 434ind, 1032SNPs, target fst = 0.0027. The global Fst of our dataset
fst(5000, 27)  #run 27 gens drift with Ne 5000
genepopToPowsim(genepopFilePath = file.path(globPath, "Px.31pops.434ind.1032SNPs.gen"),
                outDir  = file.path(outPath, "px434.1032SNPs"),
                outName = "powsim.in.fst0027.px434",
                dataDescription = "P.xylostella SNP data",
                numLoci = 1032,
                numPops = 30,
                effectivePopSize = 5000,
                gensDrift = 27,
                sampleSizes = px434SampleSizes$n[-31],
                numSimRuns = 1000)

# px 31 pops, 434ind, 1032SNPs, target fst = 0.005
fst(5000, 50)  #run 50 gens drift with Ne 5000 
genepopToPowsim(genepopFilePath = file.path(globPath, "Px.31pops.434ind.1032SNPs.gen"),
                outDir  = file.path(outPath, "px434.1032SNPs"),
                outName = "powsim.in.fst005.px434",
                dataDescription = "P.xylostella SNP data",
                numLoci = 1032, 
                numPops = 30, 
                effectivePopSize = 5000,
                gensDrift = 50,
                sampleSizes = px434SampleSizes$n[-31],
                numSimRuns = 1000,
                eraseDetailedOut = 1)

# px 31 pops, 434ind, 1032SNPs, target fst = 0.010
fst(5000, 100)  #run 253 gens drift with Ne 5000 
genepopToPowsim(genepopFilePath = file.path(globPath, "Px.31pops.434ind.1032SNPs.gen"),
                outDir  = file.path(outPath, "px434.1032SNPs"),
                outName = "powsim.in.fst01.px434",
                dataDescription = "P.xylostella SNP data",
                numLoci = 1032, 
                numPops = 30, 
                effectivePopSize = 5000,
                gensDrift = 100,
                sampleSizes = px434SampleSizes$n[-31],
                numSimRuns = 1000,
                eraseDetailedOut = 1)

# px 31 pops, 434ind, 1032SNPs, target fst = 0.025
fst(5000, 253)  #run 253 gens drift with Ne 5000 
genepopToPowsim(genepopFilePath = file.path(globPath, "Px.31pops.434ind.1032SNPs.gen"),
                outDir  = file.path(outPath, "px434.1032SNPs"),
                outName = "powsim.in.fst025.px434",
                dataDescription = "P.xylostella SNP data",
                numLoci = 1032, 
                numPops = 30, 
                effectivePopSize = 5000,
                gensDrift = 253,
                sampleSizes = px434SampleSizes$n[-31],
                numSimRuns = 1000,
                eraseDetailedOut = 1)

# Px 31 pops, 434 ind, 1032SNPs, target Fst = 0.0027 (th global Fst for the 2014 dataset)
fst(5000, 27) #> 0.0027 21/5/2018
genepopToPowsim(genepopFilePath = file.path(globPath, "Px.31pops.434ind.1032SNPs.gen"),
                outDir  = file.path(outPath, "px434.1032SNPs"),
                outName = "powsim.in.fst0027.px434",
                dataDescription = "P.xylostella SNP data",
                numLoci = 1032, 
                numPops = 30, 
                effectivePopSize = 5000,
                gensDrift = 27,
                sampleSizes = px434SampleSizes$n[-31],
                numSimRuns = 1000,
                eraseDetailedOut = 1)


# px 28 pops, 399ind, 1032SNPs, target fst = 0.001
fst(5000, 10)  #run 10 gens drift with Ne 5000 
genepopToPowsim(genepopFilePath = file.path(globPath, "Px.28pops.399ind.1032SNPs.gen"),
                outDir  = file.path(outPath, "px399.1032SNPs"),
                outName = "powsim.in.fst001.px399",
                dataDescription = "P.xylostella SNP data",
                numLoci = 1032, 
                numPops = 28, 
                effectivePopSize = 5000,
                gensDrift = 10,
                sampleSizes = px399SampleSizes$n,
                numSimRuns = 1000,
                eraseDetailedOut = 1)

fst(5000, 50)  #run 50 gens drift with Ne 5000 
genepopToPowsim(genepopFilePath = file.path(globPath, "Px.28pops.399ind.1032SNPs.gen"),
                outDir  = file.path(outPath, "px399.1032SNPs"),
                outName = "powsim.in.fst005.px399",
                dataDescription = "P.xylostella SNP data",
                numLoci = 1032, 
                numPops = 28, 
                effectivePopSize = 5000,
                gensDrift = 50,
                sampleSizes = px399SampleSizes$n,
                numSimRuns = 1000,
                eraseDetailedOut = 1)

fst(5000, 100)  #run 100 gens drift with Ne 5000 
genepopToPowsim(genepopFilePath = file.path(globPath, "Px.28pops.399ind.1032SNPs.gen"),
                outDir  = file.path(outPath, "px399.1032SNPs"),
                outName = "powsim.in.fst01.px399",
                dataDescription = "P.xylostella SNP data",
                numLoci = 1032, 
                numPops = 28, 
                effectivePopSize = 5000,
                gensDrift = 100,
                sampleSizes = px399SampleSizes$n,
                numSimRuns = 1000,
                eraseDetailedOut = 1)

fst(5000, 253)  #run 50 gens drift with Ne 5000 
genepopToPowsim(genepopFilePath = file.path(globPath, "Px.28pops.399ind.1032SNPs.gen"),
                outDir  = file.path(outPath, "px399.1032SNPs"),
                outName = "powsim.in.fst025.px399",
                dataDescription = "P.xylostella SNP data",
                numLoci = 1032, 
                numPops = 28, 
                effectivePopSize = 5000,
                gensDrift = 253,
                sampleSizes = px399SampleSizes$n,
                numSimRuns = 1000,
                eraseDetailedOut = 1)

# Px 28 pops, 399 ind, 1032SNPs, target Fst = 0.0056 (th global Fst for the 2015 dataset)
fst(5000, 56) #> 0.00558 21/5/2018
genepopToPowsim(genepopFilePath = file.path(globPath, "Px.28pops.399ind.1032SNPs.gen"),
                outDir  = file.path(outPath, "px399.1032SNPs"),
                outName = "powsim.in.fst0056.px399",
                dataDescription = "P.xylostella SNP data",
                numLoci = 1032, 
                numPops = 28, 
                effectivePopSize = 5000,
                gensDrift = 56,
                sampleSizes = px399SampleSizes$n,
                numSimRuns = 1000,
                eraseDetailedOut = 1)



# #########################################
# Output results summary in a LaTeX table #
# #########################################

# read in the summarised powsim results
powsimRes <- file.path(
  "C:", "UserData", "Kym", "PhD", "RAD2", "powsim", "Px.1032SNPs.results", "PxPowsimResultsSummary.xlsx"
) %>%
  read_excel(range = "A5:E9", col_names = FALSE) %>%
  set_names(c("Fst", "fisher2014", "chi2014", "fisher2015", "chi2015"))


# output LaTex table code
powCaption <- "Power analysis for 1032 SNP marker loci identified in Australian populations of \\textit{P. xylostella}.
The probability of SNP loci detecting true population differentiation at predefined $F_\\textsc{st}$ values 
according to Fisher’s Exact and Chi-Squared tests. 
Analyses were conducted in POWSIM assuming an effective population size ($N_e$) of 5000.
Simulations for $F_\\textsc{st}$ = 0.0027 were conducted for 2014 data only 
and $F_\\textsc{st}$ = 0.0056 for 2015 data only, corresponding to the global $F_\\textsc{st}$ values for these years." 
powDigits  <- c(0, 4, 4, 4, 4, 4)
powAlign   <- c("l", # just the row.names which aren't printed
  "S[table-number-alignment=center,table-figures-integer=1,table-figures-decimal=4]\n",
  "S[table-number-alignment=center,table-figures-integer=1,table-figures-decimal=4]\n",
  "S[table-number-alignment=center,table-figures-integer=1,table-figures-decimal=4]\n",
  "S[table-number-alignment=center,table-figures-integer=1,table-figures-decimal=4]\n",
  "S[table-number-alignment=center,table-figures-integer=1,table-figures-decimal=4]\n"
)

# set up header rows
powAddRows <- list()
powAddRows$pos <- list(0, 0)
powAddRows$command <- c(
  paste("& \\multicolumn{2}{c}{2014} & \\multicolumn{2}{c}{2015}\\\\\n",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}\n", collapse = ""),
  "$F_\\textsc{st}$ & {Fisher's Exact} & {Chi-Square} & {Fisher's Exact} & {Chi-Square}\\\\\n"
)

# Latex code
powsimRes %>%
  xtable(caption = powCaption,
         digits  = powDigits,
         align   = powAlign,
		 lab   = "tab:powsim") %>%
  print(include.rownames = FALSE, 
        include.colnames = FALSE,
        table.placement = "p",
        caption.placement = "top",
        booktabs = TRUE,
        NA.string = "\\textendash",
        add.to.row = powAddRows,
        sanitize.text.function = function(x){x})
  
# End script
# ###################################################################
