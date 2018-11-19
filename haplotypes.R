# #########################################################################################
# plot haplotype network for Plutella australiana and P. xylostella
# K. Perry, 8/5/2017
#
# Notes:
# This analysis plots haplotype networks for 27 Plutella australiana individuals 
#  ... from study `KPRAD2`, which was sequenced at the COI gene by S. Baxter, and
#  ... 57 individuals from landry and herberts study downloaded as .tsv file.
# Simon downloaded trace files, checked chromatograms, dropped low qual samples,
# ... aligned together with kp samples in geneious v10.0, trimmed to common length, 
#  ... and exported a new fasta file.
# 
# This script does:
# ... reads in the fasta file with trimmed/aligned COI sequecnes for 89 Paus individuals,
# ... ands new sample names, writes new fasta files and plots hapotype networks in R pegas.
# 
# #########################################################################################

library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(reshape2)
library(pegas)
library(RColorBrewer)
library(xtable)


# ###############################
# Write new formatted fasta files 
# ###############################

# read in the cleaned Landry and Perry COI sequences for each species
# ... (done separately by species, because the species cannot easily be determined from the combined fasta)
paDNA <- read.dna("C:/UserData/Kym/PhD/RAD2/haplotypes/Landry_Perry_Pa_Extract.fasta",
               format = "fasta")
pxDNA <- read.dna("C:/UserData/Kym/PhD/RAD2/haplotypes/Landry_Perry_Px_Extract (modified) (modified).fasta",
                  format = "fasta")

# make a dataframe containing the sequences and names
seqDf <- data.frame(oldNames = labels(paDNA),
                    seq = sapply(paDNA, paste, collapse = "")) %>%
  # filter out the landry Pa consensus sequence
  filter(oldNames != "Landry_Pa_Consensus_extraction") %>%
  mutate(species = "P.aus") %>%
  # combine with the px sequences
  bind_rows(data.frame(oldNames = labels(pxDNA),
                       seq = sapply(pxDNA, paste, collapse = "")) %>%
              mutate(species = "P.x")) %>%
  magrittr::set_rownames(NULL)

# rename the samples
# ... first, write the old names to file
write.csv(seqDf[, c("oldNames", "species")], "seqNames.csv", row.names = FALSE) 
# ... manually add new names to this file in a separate column, 
# ... then change filename (to avoid overwrite), read back in
# ... note: to find the sample locations, I cross-referenced the old names with the sample names in `COI summary sheet``
seqNames <- read.csv("C:/UserData/Kym/PhD/RAD2/haplotypes/seqNamesNew.csv")
seqDf %<>% merge(seqNames, by = c("oldNames", "species")) %>%
  # add a column identifying landry's and perry's samples
  mutate(study = NA)
seqDf$study[grep("-\\d{2}$", seqDf$newNames)] <- "landry"
seqDf$study[grep("-\\d{2}$", seqDf$newNames, invert = TRUE)] <- "perry"

# create additional seqName variables to facilitate plotting haplotype by `study` and `species` factor in R pegas 
# ... allows us to determine how many haplotypes discovered by landry, perry
seqDf %<>% mutate(ID = paste(newNames, species, study, sep = "_")) %>%
  arrange(species, study, newNames) %>%
  dplyr::select(ID, seq, species, study)

# write new fasta files for all samples, P.x samples, and P.aus samples
pxpaFA <- c(rbind(paste0("> ", as.character(seqDf$ID)), as.character(seqDf$seq)))
write.table(pxpaFA, "pxpaCOI.fasta", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

pxDf <- filter(seqDf, species == "P.x")
pxFA <- c(rbind(paste0("> ", as.character(pxDf$ID)), as.character(pxDf$seq)))
write.table(pxFA, "pxCOI.fasta", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

paDf <- filter(seqDf, species == "P.aus")
paFA <- c(rbind(paste0("> ", as.character(paDf$ID)), as.character(paDf$seq)))
write.table(paFA, "paCOI.fasta", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)


# ##################################
# Plot haplotype networks in R pegas
# 
# ... assign labels
# ... plot network with sample sizes
# ... plot network with the haplotype labels
# ... plot outer circles with smaller sizes (singleton)
# ... import to inkscape
# ##################################


# *********************
# P.x haplotype network
# *********************
pxCOI <- read.dna("C:/UserData/Kym/PhD/RAD2/haplotypes/pxCOI.fasta",
                  format = "fasta")
xh <- haplotype(pxCOI, labels = NULL)
plot(xh)
xnet <- haploNet(xh)

# set the new labels for xh and xnet
# ... create vector of hap names in the order to assign names to haps 
# ... below, print the names/freq attributes here for transparency 
pxHapNames <- c("PxCOI04", "PxCOI01", "PxCOI05", "PxCOI03", "PxCOI02")
attr(xh, "dimnames")[[1]] 
# [1] "I"   "II"  "III" "IV"  "V"
attr(xh, "dimnames")[[1]] <- pxHapNames
attr(xh, "dimnames")[[1]]
# [1] "PxCOI4" "PxCOI1" "PxCOI5" "PxCOI3" "PxCOI2"

attr(xnet, "freq")
# [1] 23 76  1  1  1
attr(xnet, "labels")
# [1] "I"   "II"  "III" "IV"  "V"
attr(xnet, "labels") <- pxHapNames
attr(xnet, "labels")
# [1] "PxCOI4" "PxCOI1" "PxCOI5" "PxCOI3" "PxCOI3"

# plot a network with the hap labels, and sample sizes
plot(xnet, threshold = 0)
# plot a network labelled with the sample size
xnet2 <- xnet
attr(xnet2, "labels") <- attr(xnet, "freq")
plot(xnet2, threshold = 0)

# export some pdfs
pdf("PxHapNet.pdf")
plot(xnet, threshold = 0)
dev.off()

pdf("PxHapNetSize.pdf")
plot(xnet2, threshold = 0)
dev.off()

# export some svgs
svg("PxHapNet.svg")
plot(xnet, threshold = 0)
dev.off()

svg("PxHapNetSize.svg")
plot(xnet2, threshold = 0)
dev.off()

# ********************************
# P. australiana haplotype network
# ********************************
paCOI <- read.dna("C:/UserData/Kym/PhD/RAD2/haplotypes/paCOI.fasta",
                  format = "fasta")
ah <- haplotype(paCOI, labels = NULL)
plot(ah)
anet <- haploNet(ah)


# set the new haplotype names for ah and anet
# . label the Pa haplotypes according to the nomenclature in draft HapNetWorkTable.xlsx
# ... view the draft haplotype network table in this file:
# . hapnetTab <- "C:/UserData/Kym/PhD/RAD2/haplotypes/HapNetworkTable.xlsx"
# ... find the index of individuals with each `haplotype` 
# ... (to id the samples, use the names in the dna.bin object, from which haplotype object was derived)
# ... then, match the samples with genbaknk accession numbers, and \
# ... \ cross reference these by vewing the genbank accession numbers for each sample in this file:
# . gbPath <-  "C:/UserData/Kym/PhD/RAD2/haplotypes/genbankAccessionsPerryPlutellaCOI.txt"

# index the paCOI object to find the accession numbers for haplotypes [[1]] to [[9]]
# ... note that we're indexing the default order of haplotypes
labels(paCOI)[attr(ah, "index")[[1]]] # n = 74 samples, Pa COI 01
labels(paCOI)[attr(ah, "index")[[2]]] # [1] "LNSWA731-05_P.aus_landry", Pa COI 09
labels(paCOI)[attr(ah, "index")[[3]]] # [1] "LSM1299-11_P.aus_landry" , Pa COI 08
labels(paCOI)[attr(ah, "index")[[4]]] # [1] "MCCAA2949-12_P.aus_landry", pa COI 04
labels(paCOI)[attr(ah, "index")[[5]]] # n = 6 samples, Pa COI 02
labels(paCOI)[attr(ah, "index")[[6]]] # [1] "PHLCA920-11_P.aus_landry", pa COI 05
labels(paCOI)[attr(ah, "index")[[7]]] # [1] "dypWAU1011Ca_P.aus_perry", matches accession MF151883, Pa COI 07
labels(paCOI)[attr(ah, "index")[[8]]] # [1] "dypWAU1012Aa_P.aus_perry", matches accession MF151885, Pa COI 03
labels(paCOI)[attr(ah, "index")[[9]]] # [1] "wagNSW0201Ga_P.aus_perry", matches accession MF151836, pa COI 06

# ... create vector of hap names in the order above to assign names to haps 
# ... below, also print the names/freq attributes here for transparency 
paHapNames <- c("PaCOI01", "PaCOI09", "PaCOI08", "PaCOI04", "PaCOI02",
                "PaCOI05", "PaCOI07", "PaCOI03", "PaCOI06")
attr(ah, "dimnames")[[1]] 
# [1] "I"    "II"   "III"  "IV"   "V"    "VI"   "VII"  "VIII" "IX" 
attr(ah, "dimnames")[[1]] <- paHapNames
attr(ah, "dimnames")[[1]]
# [1] "PaCOI01" "PaCOI09" "PaCOI08" "PaCOI04" "PaCOI02" "PaCOI05" "PaCOI07" "PaCOI03" "PaCOI06"

attr(anet, "freq")
# [1] 74  1  1  1  6  1  1  1  1
attr(anet, "labels")
# [1] "I"    "II"   "III"  "IV"   "V"    "VI"   "VII"  "VIII" "IX"  
attr(anet, "labels") <- paHapNames
attr(anet, "labels")
# [1] "PaCOI01" "PaCOI09" "PaCOI08" "PaCOI04" "PaCOI02" "PaCOI05" "PaCOI07" "PaCOI03" "PaCOI06"

# plot a network with the hap labels, and sample sizes
plot(anet, threshold = 0)
# plot a network labelled with the sample size
anet2 <- anet
attr(anet2, "labels") <- attr(anet, "freq")
plot(anet2, threshold = 0)

# export some pdfs
pdf("PaHapNet.pdf")
plot(anet, threshold = 0)
dev.off()

pdf("PaHapNetSize.pdf")
plot(anet2, threshold = 0)
dev.off()

# export some svgs
svg("PaHapNet.svg")
plot(anet, threshold = 0)
dev.off()

svg("PaHapNetSize.svg")
plot(anet2, threshold = 0)
dev.off()

# calculate sample sizes, frequency by study (landry vs perry)
(PxFreqStudy <- haploFreq(pxCOI, split = "_", what = 3))
(PaFreqStudy <- haploFreq(paCOI, split = "_", what = 3))

data.frame(landryPx = PxFreqStudy[, "landry"],
           perryPx  = PxFreqStudy[, "perry"])
# landryPx perryPx
# 1       14       9
# 2       41      35
# 3        1       0
# 4        1       0
# 5        1       0

data.frame(N_landryPx = sum(PxFreqStudy[, "landry"]),
           N_perryPx  = sum(PxFreqStudy[, "perry"]))
# N_landryPx N_perryPx
# 1         58        44

data.frame(landryPaus = PaFreqStudy[, "landry"],
           perryPaus  = PaFreqStudy[, "perry"])
# landryPaus perryPaus
# 1         42        32
# 2          1         0
# 3          1         0
# 4          1         0
# 5          4         2
# 6          1         0
# 7          0         1
# 8          0         1
# 9          0         1

data.frame(N_landryPaus = sum(PaFreqStudy[, "landry"]),
           N_perryPaus  = sum(PaFreqStudy[, "perry"]))
# N_landryPaus N_perryPaus
# 1           50          37

# #####################################################
# extract the meta data for perry P.x and P.aus samples
# ... for submission to genbank
# #####################################################

# first read in the metadata for each sample and population
# ... variable `pos` contains the unique plate and well coordinates for each individual 
genotypesPath <- "C:/UserData/Kym/PhD/RAD2/genotypesMasterRAD2.xlsx"
s <- "skip"
t <- "text"
colT <- c(s,t,t,s,s,s,s,s,"numeric",t,t,s,s,t,t,t,t,s,s,s,s,s,s,s,s,s)  
metaData <- read_excel(genotypesPath, sheet = "genotypes-master",
                              col_types = colT) %>%
  mutate(plate = gsub("JF_2", "JF02", plate) %>%
           gsub("JF_1", "JF01", .) %>%
           str_pad(side = "left", width = 2, pad = 0),
         position = str_pad(position, side = "left", width = 3, pad = 0),
         pos = paste0(plate, position)) %>%
  dplyr::select(pos, popCode, location, state, host, hostType, gender) %>%
  merge(read_excel(genotypesPath, sheet = "pops-master") %>%
          filter(popCode != 9.1) %>%
          mutate(popCode = as.numeric(popCode),
                 collectionDate = format(as.Date(collectionDate), "%d-%b-%Y")) %>%
          dplyr::select(popCode, latitude, longitude,  
                        stageCollected, collectionDate, collector), by = c("popCode"))

# get the perry seqNames, extract well and plate positions, merge with metaData
renameHost <- function(x) {
  gsub("c", "Canola", x) %>%
    gsub("wt", "Wild turnip", .) %>%
    gsub("v", "Brassica vegetables", .) %>%
    gsub("lw", "Diplotaxis sp.", .) %>%
    gsub("f", "Forage brassica", .) %>%
    gsub("wr", "Wild radish", .)
}

seqPerry <- seqDf %>%
  filter(study == "perry") %>%
  mutate(ID = gsub("_(.)+", "", ID),
         pos = str_extract(ID,"(JF)?[0-9]{4}[A-H]{1}"),
         plate = str_extract(pos, "(JF)?\\d{2}"),
         well = str_extract(pos, "\\d{2}[A-H]$")) %>%
  dplyr::select(pos, ID, species, seq) %>%
  merge(metaData, by = "pos") %>%
  mutate(Country = "Australia",
         Lat_Lon = paste(round(latitude, 2), "N", 
                         round(longitude, 2), "E"),
         host = renameHost(host),
         Isolation_source = paste(location, state),
         species = gsub("P.x", "Plutella xylostella", species) %>%
           gsub("P.aus", "Plutella australiana", .),
         gender = gsub("f$", "Female", gender) %>%
           gsub("m$", "Male", .) %>%
           gsub("u$", "Unknown", .)) %>%
  dplyr::select(Sequence_ID = ID, Collection_date = collectionDate, 
                Collected_by = collector, Country, Lat_Lon, 
                Host = host, Isolate = pos, Isolation_source, Species = species, Sex = gender, Sequence = seq)

write.csv(seqPerry, "perry.81PxPa.612bpCOI.metadata.csv")
# write a fasta file with just the perry samples (for submission to genbank)
perryFA <- c(rbind(paste0("> ", as.character(seqPerry$Sequence_ID)), as.character(seqPerry$Sequence)))
write.table(perryFA, "perry.81PxPa.612bpCOI.fasta", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)


# #######################################################################
# Make a Latex table for the GenBank accession numbers for the Haplotypes
# ... supplementary tables
# #######################################################################

# ====================
# Plutella xylostella 
# ====================

# Read in the summary table (from Simon Baxter) containing an example accession 
# for each haplotype
hapAccFile <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "haplotypes", "HapNetworkTable.xlsx")
PxAccessions <- read_excel(hapAccFile, sheet = "Sheet1", skip = 2)
PxAccessions <- PxAccessions[1:5, 1:8]

# format latex table
PxAccLatex <- PxAccessions %>%
  dplyr::select(-Study) %>%
  rename(`{No. individuals}` = n,
         `Sequence reference` = sequenceRefNo)

colWidth <- "0.75cm"
PausAccAlign <- c("l",
                  "l",
                  paste0(">{\\centering\\arraybackslash}p{", colWidth, "}"),
                  paste0("@{}>{\\centering\\arraybackslash}p{", colWidth, "}"),
                  paste0("@{}>{\\centering\\arraybackslash}p{", colWidth, "}"),
                  paste0("@{}>{\\centering\\arraybackslash}p{", colWidth, "}"),
                  "S[round-mode=places,round-precision=0,table-number-alignment=center, table-figures-decimal=0,table-figures-integer=2]", 
                  "l")

# Jess Saw Haplotypes:
# PxMt01 = DQ394347
# PxMt02 = DQ394348
# PxMt06 = DQ394352

# make a multi-row header using the add to rows argument
# make a multirow header using the add.to.row() argument
PxHeader <- list()
PxHeader$pos <- list(0, 0)
PxHeader$command = c(
  "& \\multicolumn{4}{c}{Nucleotide position} & \\\\\n \\cmidrule(lr){2-5}\n",
  paste(paste(names(PxAccLatex), collapse = " & "), "\\\\\n", collapse = "")
)

PxAccCap <- "The four variable nucleotide sites among the five \\textit{P. xylostella} 
613 bp COI haplotypes identified in 102 individuals from Australia.
Shown are sequences from this study and re-analysed sequences from Landry and Hebert (2013) downloaded
from dx.doi.org//10.5883//DS-PLUT1.
Three haplotypes correspond to those reported by Saw et al. (2006): 
PxCOI01/PxMt01, GenBank accession: DQ394347; PxCOI02/PxMt06, GenBank accession: DQ394352;
PxCOI04/PxMt02, GenBank accession: DQ394348. 
Nucleotide positions were determined from sequence MF151841.
Only positions that differ from haplotype PxCOI01 are shown."

xtable(PxAccLatex, 
       caption = PxAccCap, 
       align = PxAccAlign,
       lab = "tab:PxAccessions") %>% 
  print.xtable(include.rownames = FALSE,
               include.colnames = FALSE,
               add.to.row = PxHeader,
               caption.placement = "top",
               table.placement = "h",
               booktabs = TRUE, 
               sanitize.text.function = function(x){x})

  
# ######################
# Plutella australiana #
# ######################

# format latex table
PausAccLatex <- read_excel(hapAccFile, sheet = "Sheet1", skip = 11) %>%
  filter(Haplotype != "Total") %>%
  dplyr::select(-Study) %>%
  rename(`{No. individuals}` = n,
         `Sequence reference` = sequenceRefNo)

names(PausAccLatex)

# make a multi-row header using the add to rows argument
# make a multirow header using the add.to.row() argument
PaHeader <- list()
PaHeader$pos <- list(0, 0)
PaHeader$command = c(
  "& \\multicolumn{8}{c}{Nucleotide position} & \\\\\n \\cmidrule(lr){2-9}\n",
  paste(paste(names(PausAccLatex), collapse = " & "), "\\\\\n", collapse = "")
  )

colWidth <- "0.7cm"
# I add the @{} space, otherwise cols too squashed to centre exactly.
PausAccAlign <- c("l",
                  "l",
                  paste0(">{\\centering\\arraybackslash}p{", colWidth, "}"),
                  paste0("@{}>{\\centering\\arraybackslash}p{", colWidth, "}"),
                  paste0("@{}>{\\centering\\arraybackslash}p{", colWidth, "}"),
                  paste0("@{}>{\\centering\\arraybackslash}p{", colWidth, "}"),
                  paste0("@{}>{\\centering\\arraybackslash}p{", colWidth, "}"),
                  paste0("@{}>{\\centering\\arraybackslash}p{", colWidth, "}"),
                  paste0("@{}>{\\centering\\arraybackslash}p{", colWidth, "}"),
                  paste0("@{}>{\\centering\\arraybackslash}p{", colWidth, "}"),
                  "S[round-mode=places,round-precision=0,table-number-alignment=center, table-figures-decimal=0,table-figures-integer=2]", 
                  "l")

PausAccCap <- "The eight variable nucleotide sites among the nine \\textit{P. australiana} 
613 bp COI haplotypes identified in 87 individuals from Australia.
Haplotypes PaCOI01 and PaCOI02 were identified among sequences from this study and Landry and Hebert (2013), 
and PaCOI04, PaCOI05, PaCOI08 and PaCOI09 were identified from Landry and Hebert (2013).
Nucleotide positions were determined from sequence MF151865.
Only the positions that differ from haplotype PaCOI01 are shown."

xtable(PausAccLatex, 
       caption = PausAccCap, 
       align = PausAccAlign,
       lab = "tab:PausAccessions") %>% 
  print.xtable(include.rownames = FALSE,
               include.colnames = FALSE,
               add.to.row = PaHeader,
               caption.placement = "top",
               table.placement = "h",
               booktabs = TRUE, 
               sanitize.text.function = function(x){x})


# # End script
# ##################################################################

