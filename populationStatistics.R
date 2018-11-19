# ###########################################################################
# Create population summary statistics tables for populations genetic datasets
# K. Perry, 17/4/2018
# 
# This script analyses the following datasets:
# -- Paus.5pop5.53ind,  ~allSites(410780) / ~SNPs(978)
# -- pxpa.5pops.101ind, ~allSites(305531) / ~SNPs(708)
# -- 

# For each dataset, this script does:
# -- reads in population level summary parameters (ar, size, Ho, He, Fis) generated using R diveRsity:basicstats() function
# -- reads in additional summary parameters generated using vcfstats
# -- combine stats into a dataframe, formats
# -- outputs latex table code, which has separate upper/lower panels for the SNP/allSites statistics 
#
# Notes:
# -- this script uses the workspace in local directory `\RAD2\populationStatistics`
# -- the allSites files include indels
# -- the pxpa5pops datasets includes samples from 5 locations, but is strictly 10 pops/table rows (2 species).
# -- reading in files with list.files may fai if any files are open (creates duplicated name with ~ in dir)
# ########################################################################################

library(tidyverse)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(xtable)

# ###########################
# define functions and data #
# ###########################

# function to rename the 13 mislabelled Paus samples
# ... (initially labelled as Px, after sequencing, identified as Paus)
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

# utility functions to extract population meta-data from the individuals names

extractPop <- function(ind) {
  stringr::str_sub(ind, 1, 9) %>%
    gsub("boyW14sca", "espW14sca", .) # ensures one boyup sample is grouped with the esperance pop (P. aus only)
}
extractHost <- function(x) {
  stringr::str_sub(x, 8, 8) %>%
    gsub("^v$", "Vegetables", .) %>%
    gsub("^c$", "Canola", .) %>%
    gsub("^f$", "Forage", .) %>%
    gsub("^w$", "Wild brassicas", .)
}
extractYear <- function(x) {
  stringr::str_sub(x, 5, 6) %>%
    paste0("20", .) %>% 
    as.numeric()
}
extractSeason <- function(x) {
  stringr::str_sub(x, 7, 7) %>%
    gsub("^s$", "Spring", .) %>%
    gsub("^a$", "Autumn", .)
}
extractSpecies <- function(x) {
  str_sub(x, 9, 9) %>%
  gsub("^a$", "P. aus", .) %>%
  gsub("^x$", "P. x", .)
}
extractAusState <- function(x) {
  str_sub(x, 4, 4) %>%
    gsub("^N$", "NSW", .) %>%
    gsub("^S$", "SA",  .) %>%
    gsub("^V$", "VIC", .) %>%
    gsub("^W$", "WA",  .) %>%
    gsub("^T$", "TAS", .) %>%
    gsub("^Q$", "QLD", .)
}
extractSampleNum <- function(x) {
  str_sub(x, 11, 12) %>%
    as.numeric()
}
# pretty names for publication
popNames <- function(x, popsDf = popsMaster){
  vapply(
    X = x,
    FUN = function(x) {
      popsDf$locationState[which(popsDf$popString == x)]
    },
    FUN.VALUE = character(1)
  ) %>% as.character()
}

# ###############
# organise data #
# ###############

# read in population metadata and format location name for displaying in a table
popsMaster <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "popsMasterRAD2.csv") %>%
  read_csv() %>%
  mutate(locationState = paste(Location, State))

# define local directories containing the summary output from various programs
vcfPath      <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "vcftools")   # vcftools/vcfstats summaries
hrfPath      <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "hierfstats") # hierfstat summaries
popStatsPath <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "populationStatistics") 

# ##################################################################################
# Generate a population statistics table for dataset: PxPa, 5pops, 99ind, allSites #
# ##################################################################################

# read in the population summary table generated using R package `hierfstat`
PxPaHierfstatAllsites <- file.path(popStatsPath, "PxPa.5pops.99ind.allSites.hierfstat.basicstats.pop.txt") %>%
  read.table(header = TRUE) %>%
  mutate(pop = extractPop(pop))

# read in the vcftools depth summary
PxPaVcfDepthAllsites <- file.path(vcfPath, "PxPa.5pops.99ind.genotypeGVCFs.allSites.hardFilt.maxMiss0.7.idepth") %>%
  read.table(header = TRUE) %>%
  set_names(c("ind", "sites", "idepth"))

# read in the vcfstats summary
PxPaVcfStatsAllsites  <- file.path(vcfPath, "PxPa.5pops.99ind.genotypeGVCFs.allSites.hardFilt.maxMiss0.7.vcfstats.grepSummary.txt") %>%
  read.table(header = TRUE)

# merge results into a data_frame summarising the population level genetic diversity statistics
PxPaStatisticsAllsites <- PxPaVcfStatsAllsites %>%
  full_join(PxPaVcfDepthAllsites, by = "ind") %>%
  mutate(ind = renamePxToPaus(ind),
         pop = extractPop(ind),
         species = extractSpecies(pop)) %>%
  group_by(pop) %>%
  summarise(species = unique(species),
            nInd   = length(ind),
            nSites = round(mean(sites), 1),
            pSites = round(nSites / 305136, 2),
            idepth = round(mean(idepth), 1),
            homAA  = round(mean(homAA), 1),
            hetRA  = round(mean(hetRA), 1),
            homRR  = round(mean(homRR), 1),
            ref_count =   round(mean(ref_count), 1),
            snp_count =   round(mean(snp_count), 1),
            indel_count = round(mean(indel_count), 1),
            private = round(mean(private), 1)) %>% 
  full_join(PxPaHierfstatAllsites, by = "pop") %>%
  mutate(popName = popNames(pop))#%>%
  #left_join(popLabs, by = "pop")

# output LaTex table code
PxPaStatsAllsitesLatex <- PxPaStatisticsAllsites %>%
  dplyr::select(Population = popName,
                `\\multicolumn{1}{l}{Species}` = species,
                `{\\textit{n}}` = indPerLoc, 
                `{Sites}` = nSites, # wrap S column headers in curly braces, so it doesn't try to treat as a number
                `{\\makecell{Site\\\\depth}}` = idepth,
                `{SNPs}`= snp_count,
                `{Indels}` = indel_count,
                `{\\makecell{Private\\\\sites}}` = private,
                `$H_{\\textsc{o}}$` = Ho,
                `$H_{\\textsc{s}}$` = Hs,
                `$F_{\\textsc{is}}$` = Fis)

PxPaStatsAllsitesLatex$Population[duplicated(PxPaStatsAllsitesLatex$Population)] <- NA
PxPaAllsitesAlign <- c(
  "l", "l", ">{\\itshape}l",
  "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=3]", # center align S columns without decimals by pasting this option onto relevant column aligns, 
  "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=6]", 
  "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]", 
  "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=4]", 
  "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=4]", 
  "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=3]",
  "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]", 
  "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]", 
  "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]"
  )
PxPaAllsitesDigits <- c(0, 0, 0, 1, 0, 0, 0, 0, 0, 3, 3, 3)

PxPaAllsitesCaption <- "Population statistics for variant and invariant sites for sympatric populations of 
\\textit{P. australiana} (\\textit{P. aus}) and \\textit{P. xylostella} (\\textit{P. x}) from five locations."

xtable(PxPaStatsAllsitesLatex,
       caption = PxPaAllsitesCaption,
       align   = PxPaAllsitesAlign,
       digits  = PxPaAllsitesDigits,
       lab = "tab:pxpaPopstatsAS") %>% 
  print.xtable(floating = TRUE,
               include.rownames = FALSE,
               include.colnames = TRUE,
               table.placement = "h",
               caption.placement = "top",
               NA.string = "",
               booktabs = TRUE,
               sanitize.text.function = function(x){x})  

# ################################################################################
# Calculate mean depth per individual across genotyped sites for three datasets: #
# (for reporting in the P. australiana manuscript)
# ################################################################################
# 30/3/2018

# -- PxPa, 5 pops, 99 ind, allSites
mean(PxPaVcfDepthAllsites$idepth)
#> [1] 36.71231

# -- PxPa, 5 pops, 99 ind, variants
PxPaVcfDepthSNPs <- file.path(vcfPath, "PxPa.5pops.99ind.genotypeGVCFs.variants.hardFilt.maxMiss0.7.idepth") %>%
  read.table(header = TRUE) %>%
  set_names(c("ind", "sites", "idepth"))
mean(PxPaVcfDepthSNPs$idepth)
#> [1] 36.19314

# -- Pa, 5 pops, 52 ind, 974 SNPs
Pa52indVcfDepth974SNPs <- file.path(vcfPath, "Pa.5pops.52ind.genotypeGVCFs.variants.hardFilt.maxMiss0.7.idepth") %>%
  read.table(header = TRUE) %>%
  set_names(c("ind", "sites", "idepth"))
mean(Pa52indVcfDepth974SNPs$idepth)
#> [1] 33.74253

# #################################################
# Boxplot of the proportion of heterozygous sites #
# #################################################

# Dataset: PxPa, 5 pops, 99 individuals, allSites with indels removed
# Uses the output of python script (Simon Martin, UK): parseVCF.py, evaluateCalls.py, 

# read in data
PxPaPyStats <- file.path(
  popStatsPath, "PxPa.5pops.99ind.genotypeGVCFs.allSites.hardFilt.maxMiss0.7.rmIndels.recode.out"
  ) %>%
  read.table(skip = 123, header = TRUE) %>%
  set_names(c("ind", "calls", "percent", "hom", "het")) %>%
  mutate(ind = renamePxToPaus(ind),
         pop = extractPop(ind),
         species = extractSpecies(ind),
         propHet = het / calls) %>%
  mutate(popName = popNames(pop)) %>%
  # strip the state abbreviation from pop name
  mutate(pop2 = str_split(popName, " ") %>%
           map_chr(magrittr::extract, 1))
  
greypal <-brewer.pal(9, "Greys")
dw <- 0.7 # dodge width
yBreaks <- seq(0, 15, by = 1) / 1000
lineWd <- 0.5 # line width for boxplot and axes  
ggplot(PxPaPyStats, aes(x = pop2, y = propHet, fill = species)) + 
  geom_boxplot(position = position_dodge(width = dw),
               outlier.shape = NA,
               width = 0.7, lwd = lineWd) + 
  geom_point(position = position_jitterdodge(dodge.width = dw,
                                             jitter.width = 0.1),
             size = 1.75) +
  theme_classic(base_size = 16) +
  theme(axis.title.x = element_text(#face = "bold", 
                                    colour = 'black',
                                    margin = margin(t = 20)), # move text away from axis
        axis.title.y = element_text(#face = "bold", 
                                    colour = 'black',
                                    margin = margin(r = 20)),
        axis.text.x  = element_text(colour = 'black'),
        axis.text.y  = element_text(colour = 'black'),
        axis.line    = element_line(size = lineWd),
        axis.ticks   = element_line(size = lineWd), 
        panel.grid.major.x = element_line(colour = alpha(greypal[3], 0.5), linetype = 1),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        #axis.line = element_line(size = 0.6, color = "black"),
        legend.text  = element_text(face = "italic"),
        legend.position = c(0.85, 0.95),
        aspect.ratio = 1) +
  labs(x = "Location", y = "Proportion of heterozygous sites", fill = NULL) +
  scale_y_continuous(breaks = yBreaks) +
  scale_fill_manual(values = greypal[c(6,3)],
                    labels = c("P. australiana", "P. xylostella")) +
  # reorder the x variable geographically West to East
  scale_x_discrete(limits = c("Esperance", "Calca", "Goulburn", "Gilgandra", "Boomi"))
ggsave("heterozygosityBoxplotFinalProof.pdf", width = 6, height = 6) 


# #################################################################################################
# generate a population statistics table for the P. xylostella RADseq datasets from 2014 and 2015 #
# #################################################################################################

# list individuals names used in the P. xylostella study 
PxSamples <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "px833Ind.txt") %>%
  read.table() %>% unlist() %>% as.character()

# read in a file with the plot geographic plotting order for the Pxyl pops
popsOrderPx <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "popsMasterRAD2PlotOrderPxyl.xlsx") %>%
  read_excel() %>%
  dplyr::select(popString, Year, orderGeo, orderYearGeo)

# =============================================================
# 2014 SNPs dataset: 31 populations, 434 individuals, 1032 SNPs
# =============================================================

# list the 2014 samples
PxSamples2014 <- PxSamples[extractYear(PxSamples) == 2014] 

# read in the population statistics output from other programs
hrfPath <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "hierfstat") 
Px2014HierfstatSNPs <- file.path(hrfPath, "Px.31pops.434ind.1032SNPs.hierfstat.basicstats.popsummary.txt") %>% 
  read.table(header = TRUE) %>%
  mutate(pop = extractPop(pop))

# read in vcftools indv site depth summary
Px2014VcftoolsDepthSNPs <- file.path(vcfPath, "Px.59pops.833ind.genotypeGVCFs.variants.hardFilt.maxMiss0.8.maf0.01.idepth") %>%
  read.table(header = TRUE) %>%
  set_names(c("ind", "sites", "idepth")) %>%
  filter(ind %in% PxSamples2014) %>%
  mutate(ind = as.character(ind)) # to avoid warnings about factor levels.

# read in vcfstats summary
Px2014VcfstatsSNPs <- file.path(vcfPath, "Px.31pops.434ind.2014.genotypeGVCFs.variants.hardFilt.maxMiss0.8.maf0.01.vcfstats.grepSummary.txt") %>%
  read.table(header = TRUE) %>%
  mutate(ind = as.character(ind))

# Now merge the stastics into a population summary
Px2014SummarySNPs <- Px2014VcfstatsSNPs %>%
  full_join(Px2014VcftoolsDepthSNPs, by = "ind") %>%
  mutate(pop = extractPop(ind)) %>%
  group_by(pop) %>%
  summarise(nInd   = length(ind),
            nSites = round(mean(sites), 1),
            pSites = round(nSites / 1032, 2),
            idepth = round(mean(idepth), 1),
            homAA  = round(mean(homAA), 1),
            hetRA  = round(mean(hetRA), 1),
            homRR  = round(mean(homRR), 1),
            ref_count = round(mean(ref_count), 1),
            snp_count = round(mean(snp_count), 1)) %>%
  full_join(Px2014HierfstatSNPs, by = "pop") %>%
  mutate(popName = popNames(pop)) #%>%
  # merge(popsOrderPx %>%
  #         rename(pop = popString)) %>%
  # arrange(orderYearGeo)

# ========================================================================================
# 2014 All sites dataset: 31 populations, 434 individuals, 590086 confidently-called sites
# ========================================================================================

# list the 2014 samples
PxSamples2014 <- PxSamples[extractYear(PxSamples) == 2014] 

# read in the population statistics output from other programs
hrfPath <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "hierfstat") 
Px2014HierfstatAllsites <- file.path(hrfPath, "Px.31pops.434ind.2014.allSites.hierfstat.basicstats.pop.txt") %>% 
  read.table(header = TRUE) %>%
  mutate(pop = extractPop(pop))

# read in vcftools indv site depth summary
vcfPath <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "vcftools")
Px2014VcftoolsDepthAllsites <- file.path(vcfPath, "Px.59pops.833ind.genotypeGVCFs.allSites.hardFilt.maxMiss0.8.idepth") %>%
  read.table(header = TRUE) %>%
  set_names(c("ind", "sites", "idepth")) %>%
  filter(ind %in% PxSamples2014) %>%
  mutate(ind = as.character(ind)) # to avoid warnings about factor levels.

# read in vcfstats summary
Px2014VcfstatsAllsites <- file.path(vcfPath, "Px.31pops.434ind.2014.genotypeGVCFs.allSites.hardFilt.maxMiss0.8.vcfstats.grepSummary.txt") %>%
  read.table(header = TRUE) %>%
  mutate(ind = as.character(ind))

# Now merge the stastics into a population summary
Px2014SummaryAllsites <- Px2014VcfstatsAllsites %>%
  full_join(Px2014VcftoolsDepthAllsites, by = "ind") %>%
  mutate(pop = extractPop(ind)) %>%
  group_by(pop) %>%
  summarise(nInd   = length(ind),
            nSites = round(mean(sites), 1),
            pSites = round(nSites / 590086, 2),
            idepth = round(mean(idepth), 1),
            homAA  = round(mean(homAA), 1),
            hetRA  = round(mean(hetRA), 1),
            homRR  = round(mean(homRR), 1),
            ref_count = round(mean(ref_count), 1),
            snp_count = round(mean(snp_count), 1),
            indel_count = round(mean(indel_count), 1),
            private = round(mean(private), 1)) %>%
  full_join(Px2014HierfstatAllsites, by = "pop") %>%
  mutate(popName = popNames(pop)) #%>%
# merge(popsOrderPx %>%
#         rename(pop = popString)) %>%
# arrange(orderYearGeo)

# ==============================================================
# 2015 SNPs dataset: 31 populations, 399 individuals, 1032 SNPs
# ==============================================================

# list the 2015 samples
PxSamples2015 <- PxSamples[extractYear(PxSamples) == 2015] 

# read in the population statistics output from other programs
hrfPath <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "hierfstat") 
Px2015HierfstatSNPs <- file.path(hrfPath, "Px.28pops.399ind.1032SNPs.hierfstat.basicstats.popsummary.txt") %>% 
  read.table(header = TRUE) %>%
  mutate(pop = extractPop(pop))

# read in vcftools indv site depth summary
vcfPath <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "vcftools")
Px2015VcftoolsDepthSNPs <- file.path(vcfPath, "Px.59pops.833ind.genotypeGVCFs.variants.hardFilt.maxMiss0.8.maf0.01.idepth") %>%
  read.table(header = TRUE) %>%
  set_names(c("ind", "sites", "idepth")) %>%
  filter(ind %in% PxSamples2015) %>%
  mutate(ind = as.character(ind)) # to avoid warnings about factor levels.

# read in vcfstats summary
Px2015VcfstatsSNPs <- file.path(vcfPath, "Px.28pops.399ind.2015.genotypeGVCFs.variants.hardFilt.maxMiss0.8.maf0.01.vcfstats.grepSummary.txt") %>%
  read.table(header = TRUE) %>%
  mutate(ind = as.character(ind))

# Now merge the stastics into a population summary
Px2015SummarySNPs <- Px2015VcfstatsSNPs %>%
  full_join(Px2015VcftoolsDepthSNPs, by = "ind") %>%
  mutate(pop = extractPop(ind)) %>%
  group_by(pop) %>%
  summarise(nInd   = length(ind),
            nSites = round(mean(sites), 1),
            pSites = round(nSites / 1032, 2),
            idepth = round(mean(idepth), 1),
            homAA  = round(mean(homAA), 1),
            hetRA  = round(mean(hetRA), 1),
            homRR  = round(mean(homRR), 1),
            ref_count = round(mean(ref_count), 1),
            snp_count = round(mean(snp_count), 1)) %>%
  full_join(Px2015HierfstatSNPs, by = c("pop")) %>%
  mutate(popName = popNames(pop)) #%>%
# merge(popsOrderPx %>%
#         rename(pop = popString)) %>%
# arrange(orderYearGeo)

# sea rocket pops: (picS14swx, picS14swx, sthS15awx, wlkS15awx, wlkS15swx)
# ========================================================================================
# 2015 All sites dataset: 28 populations, 399 individuals, 590086 confidently-called sites
# ========================================================================================

# list the 2015 samples
PxSamples2015 <- PxSamples[extractYear(PxSamples) == 2015] 

# read in the population statistics output from other programs
hrfPath <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "hierfstat") 
Px2015HierfstatAllsites <- file.path(hrfPath, "Px.28pops.399ind.2015.allSites.hierfstat.basicstats.pop.txt")  %>% 
  read.table(header = TRUE) %>%
  mutate(pop = extractPop(pop))

# read in vcftools indv site depth summary
vcfPath <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "vcftools")
Px2015VcftoolsDepthAllsites <- file.path(vcfPath, "Px.59pops.833ind.genotypeGVCFs.allSites.hardFilt.maxMiss0.8.idepth") %>%
  read.table(header = TRUE) %>%
  set_names(c("ind", "sites", "idepth")) %>%
  filter(ind %in% PxSamples2015) %>%
  mutate(ind = as.character(ind)) # to avoid warnings about factor levels.

# read in vcfstats summary
Px2015VcfstatsAllsites <- file.path(vcfPath, "Px.28pops.399ind.2015.genotypeGVCFs.allSites.hardFilt.maxMiss0.8.vcfstats.grepSummary.txt") %>%
  read.table(header = TRUE) %>%
  mutate(ind = as.character(ind))

# Now merge the stastics into a population summary
Px2015SummaryAllsites <- Px2015VcfstatsAllsites %>%
  full_join(Px2015VcftoolsDepthAllsites, by = "ind") %>%
  mutate(pop = extractPop(ind)) %>%
  group_by(pop) %>%
  summarise(nInd   = length(ind),
            nSites = round(mean(sites), 1),
            pSites = round(nSites / 590086, 2),
            idepth = round(mean(idepth), 1),
            homAA  = round(mean(homAA), 1),
            hetRA  = round(mean(hetRA), 1),
            homRR  = round(mean(homRR), 1),
            ref_count = round(mean(ref_count), 1),
            snp_count = round(mean(snp_count), 1),
            indel_count = round(mean(indel_count), 1),
            private = round(mean(private), 1)) %>%
  full_join(Px2015HierfstatAllsites, by = "pop") %>%
  mutate(popName = popNames(pop)) #%>%
# merge(popsOrderPx %>%


# ======================================================
# Summarise values (mean, range) for the manuscript text
# ======================================================
summariseStats <- function(df) {
  df %>%
    summarise(
      nPops   = nrow(.),
      meanHo  = mean(Ho), sdHo  = sd(Ho), minHo = min(Ho), maxHo = max(Ho),
      meanHs  = mean(Hs), sdHs  = sd(Hs), minHs = min(Hs), maxHs = max(Hs),
      meanFis = mean(Fis),sdFis = sd(Fis), minFis = min(Fis), maxFis = max(Fis)
      ) %>%
    map_dfr(round, 4) %>%
    map_dfr(format, scientific = FALSE)
}

# ALL SITES: average stats for populations across years, hosts, seasons
tmpDfAllsites <- Px2014SummaryAllsites %>%
  bind_rows(Px2015SummaryAllsites) %>%
  mutate(host = extractHost(pop),
         year = extractYear(pop),
         season = extractSeason(pop))
# overall
tmpDfAllsites %>% summariseStats()
# A tibble: 1 x 13
#        nPops meanHo   sdHo  minHo  maxHo meanHs   sdHs  minHs  maxHs meanFis  sdFis  minFis maxFis
#        <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>   <chr>  <chr>   <chr>  <chr>
#   1    59   0.0092  0.0002 0.0088 0.0097 0.0096 0.0001 0.0091 0.0099  0.0408 0.0152 -0.0131  0.071
tmpDfAllsites %>%
  split(.$year) %>%
  map_dfr(summariseStats, .id = "year")
# # A tibble: 2 x 14
# year nPops meanHo   sdHo  minHo  maxHo meanHs   sdHs  minHs  maxHs meanFis  sdFis  minFis maxFis
# <chr> <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>   <chr>  <chr>   <chr>  <chr>
# 1  2014    31 0.0093 0.0002 0.0089 0.0097 0.0097 0.0001 0.0095 0.0099  0.0392 0.0158  0.0104  0.071
# 2  2015    28 0.0092 0.0002 0.0088 0.0096 0.0096 0.0001 0.0091 0.0098  0.0426 0.0146 -0.0131 0.0634
tmpDfAllSites %>%
  split(.$host) %>%
  map_dfr(summariseStats, .id = "host")
# # A tibble: 4 x 14
#             host nPops meanHo   sdHo  minHo  maxHo meanHs   sdHs  minHs  maxHs meanFis  sdFis  minFis maxFis
# <chr> <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>   <chr>  <chr>   <chr>  <chr>
# 1         Canola    29 0.0092 0.0002  0.009 0.0097 0.0096 0.0001 0.0095 0.0098  0.0426 0.0115  0.0104 0.0619
# 2         Forage     3 0.0092 0.0003 0.0089 0.0096 0.0096 0.0001 0.0096 0.0097  0.0443 0.0164  0.0278 0.0607
# 3     Vegetables    15 0.0092 0.0003 0.0088 0.0097 0.0096 0.0001 0.0094 0.0098  0.0413 0.0183  0.0107  0.071
# 4 Wild brassicas    12 0.0093 0.0002  0.009 0.0097 0.0096 0.0002 0.0091 0.0099  0.0349 0.0189 -0.0131 0.0579
tmpDfAllsites %>%
  split(.$season) %>%
  map_dfr(summariseStats, .id = "season")
# # A tibble: 2 x 14
#   season  nPops meanHo   sdHo  minHo  maxHo meanHs   sdHs  minHs  maxHs meanFis  sdFis  minFis maxFis
# <chr> <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>   <chr>  <chr>   <chr>  <chr>
# 1 Autumn     7 0.0092 0.0003 0.0089 0.0097 0.0096 0.0003 0.0091 0.0099   0.032 0.0244 -0.0131  0.0607
# 2 Spring    52 0.0092 0.0002 0.0088 0.0097 0.0096 0.0001 0.0094 0.0098   0.042 0.0134  0.0104  0.071

# SNPs: average stats for populations across years, hosts, seasons
tmpDfSNPs <- Px2014SummarySNPs %>%
  bind_rows(Px2015SummarySNPs) %>%
  mutate(host = extractHost(pop),
         year = extractYear(pop),
         season = extractSeason(pop))

tmpDfSNPs %>%
  split(.$year) %>%
  map_dfr(summariseStats, .id = "year")
# # A tibble: 2 x 14
# year nPops meanHo   sdHo  minHo  maxHo meanHs   sdHs  minHs  maxHs meanFis  sdFis  minFis maxFis
# <chr> <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>   <chr>  <chr>   <chr>  <chr>
#   1  2014    31 0.2074 0.0075 0.1949 0.2233 0.2058 0.0042 0.2007 0.2205 -0.0079 0.0179 -0.0459 0.0207
# 2  2015    28 0.2057 0.0056 0.1948 0.2164 0.2045 0.0041 0.1953 0.2116 -0.0069 0.0149 -0.0598 0.0153

tmpDfSNPs %>%
  split(.$host) %>%
  map_dfr(summariseStats, .id = "host")
# # A tibble: 4 x 14
# host nPops meanHo   sdHo  minHo  maxHo meanHs   sdHs  minHs  maxHs meanFis  sdFis  minFis maxFis
# <chr> <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>   <chr>  <chr>   <chr>  <chr>
# 1         Canola    29 0.2061 0.0062 0.1948 0.2172 0.2051 0.0029 0.1985 0.2104 -0.0053 0.0151 -0.0385 0.0207
# 2         Forage     3 0.2046 0.0041  0.202 0.2093 0.2032 0.0005 0.2026 0.2036 -0.0057 0.0125 -0.0202 0.0019
# 3     Vegetables    15 0.2068 0.0064 0.1959 0.2182 0.2053 0.0036 0.1994 0.2116 -0.0083 0.0164 -0.0459  0.016
# 4 Wild brassicas    12 0.2081 0.0089 0.1968 0.2233 0.2057 0.0072 0.1953 0.2205 -0.0119 0.0207 -0.0598 0.0125

tmpDfSNPs %>%
  split(.$season) %>%
  map_dfr(summariseStats, .id = "season")
# # A tibble: 2 x 14
# season nPops meanHo   sdHo  minHo  maxHo meanHs   sdHs  minHs  maxHs meanFis  sdFis  minFis maxFis
# <chr> <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>   <chr>  <chr>   <chr>  <chr>
#   1 Autumn     7 0.2079 0.0105 0.1968 0.2233 0.2051 0.0093 0.1953 0.2205 -0.0142 0.0258 -0.0598 0.0117
# 2 Spring    52 0.2064 0.0061 0.1948 0.2182 0.2052 0.0031 0.1985 0.2116 -0.0065 0.0148 -0.0459 0.0207



# LEGACY (But correct numbers)
# summariseStats(Px2014SummaryAllsites)
# #> meanHo  minHo  maxHo meanFis minFis maxFis
# #> 0.0093 0.0089 0.0097  0.0392 0.0104  0.071
# summariseStats(Px2014SummarySNPs)
# #> meanHo  minHo  maxHo meanFis  minFis maxFis
# #> 0.2074 0.1949 0.2233 -0.0079 -0.0459 0.0207
# summariseStats(Px2015SummaryAllsites)
# #> meanHo  minHo  maxHo meanFis  minFis maxFis
# #> 0.0092 0.0088 0.0096  0.0426 -0.0131 0.0634
# summariseStats(Px2015SummarySNPs)
# #> meanHo  minHo  maxHo meanFis  minFis maxFis
# #> 0.2057 0.1948 0.2164 -0.0069 -0.0598 0.0153


# ===========================================================================================
# Summarise mean depth per individual across genotyped sites for the P. xylostella manuscript
# ===========================================================================================
# 6/5/2018
PxVcftoolsDepthAllsites <- file.path(
  vcfPath, "Px.59pops.833ind.genotypeGVCFs.allSites.hardFilt.maxMiss0.8.idepth"
  ) %>%
  read.table(header = TRUE) %>%
  set_names(c("ind", "sites", "idepth")) %>%
  mutate(ind = as.character(ind)) # to avoid warnings about factor levels.

PxVcftoolsDepthSNPs <- file.path(
  vcfPath, "Px.59pops.833ind.genotypeGVCFs.variants.hardFilt.maxMiss0.8.maf0.01.idepth"
) %>%
  read.table(header = TRUE) %>%
  set_names(c("ind", "sites", "idepth")) %>%
  mutate(ind = as.character(ind)) # to avoid warnings about factor levels.

# -- Px, 59 pops, 833 ind, ALL sites
mean(PxVcftoolsDepthAllsites$idepth) # [1] 33.37822
# -- Px, 59 pops, 833 ind, SNPs
mean(PxVcftoolsDepthSNPs$idepth)     # [1] 33.95487

# ==============================================
# Plot gene diversity vs observed heterozygosity
# ==============================================

# a function to highlight the sea rocket populations
extractSeaRocketPops <- function(x){
  srPops <- c("picS14awx" = "Sea rocket", 
              "picS14swx" = "Sea rocket", "sthS15awx" = "Sea rocket", 
              "wlkS15awx" = "Sea rocket", "wlkS15swx" = "Sea rocket")
  x <- str_replace_all(x, srPops)
  y <- rep("Other", length(x))
  y[which(x == "Sea rocket")] <- "Sea rocket"
  y
  }
  

PxAllJoined <- Px2014SummaryAllsites %>%
  mutate(dataset = "allSites", Year = 2014) %>%
  bind_rows(Px2014SummarySNPs %>% mutate(dataset = "SNPs", Year = 2014)) %>%
  bind_rows(Px2015SummarySNPs %>% mutate(dataset = "SNPs", Year = 2015)) %>%
  bind_rows(Px2015SummaryAllsites %>% mutate(dataset = "allSites", Year = 2015)) %>%
  mutate(datasetYear = paste(dataset, Year),
         season = extractSeason(pop),
         hostType  = extractHost(pop),
         seaRocket = extractSeaRocketPops(pop))


# to get equal coords, plot datasets separately first, then use grid.arrange
(plotSNPs <- PxAllJoined %>%
  mutate(Year = factor(Year)) %>%
  filter(dataset == "SNPs") %>%
  ggplot(aes(x = Hs, y = Ho#, 
             #colour = seaRocket
             )) +
    geom_point(colour = alpha("red", 0.8)) +
    geom_abline(slope = 1, lty = 2) +
    theme_grey(base_size = 14) +
    theme(aspect.ratio = 1) +
    coord_fixed(ratio = 1) +
    labs(x = "Gene diversity (Expected heterozygosity)", 
         y = "Observed heterozygosity") +
    scale_x_continuous(limits = c(0.1925, 0.225),
                       breaks = seq(195, 225, by = 5) / 1000,
                       labels = seq(195, 225, by = 5) / 1000) +
    scale_y_continuous(limits = c(0.1925, 0.225),
                       breaks = seq(195, 225, by = 5) / 1000,
                       labels = seq(195, 225, by = 5) / 1000) +
    facet_wrap(~datasetYear, scales = "fixed", nrow = 1))


(plotAllsites <- PxAllJoined %>%
  mutate(Year = factor(Year)) %>%
  filter(dataset == "allSites") %>%
  ggplot(aes(x = Hs, y = Ho#, 
             #colour = seaRocket
             )) +
    geom_point(colour = "blue") +
    geom_abline(slope = 1, lty = 2) +
    theme_grey(base_size = 14) +
    theme(aspect.ratio = 1) +
    coord_fixed(ratio = 1) +
    labs(x = "Gene diversity (Expected heterozygosity)", 
         y = "Observed heterozygosity") +
    scale_x_continuous(limits = c(0.00875, 0.01000),
                       breaks = seq(875, 1000, by = 25) / 1e5,
                       labels = seq(875, 1000, by = 25) / 1e5) +
    scale_y_continuous(limits = c(0.00875, 0.01000),
                       breaks = seq(875, 1000, by = 25) / 1e5,
                       labels = seq(875, 1000, by = 25) / 1e5) +
    facet_wrap(~datasetYear, scales = "fixed", nrow = 1))

# # save a plot legend
# saveLegend <- function(gplot){
#   tmp <- ggplot_gtable(ggplot_build(gplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   legend
# }
# legendSNPs <- saveLegend(plotSNPs)

# draw plots together
pdf("PxObservedVExpectedHetDraft.pdf", width = 15, height = 5)
gridExtra::grid.arrange(plotAllsites, plotSNPs, nrow = 1)
dev.off()

# ===================================================
# Now produce Latex tables for 2014 and 2015 datasets
# ===================================================

# ====
# 2014
# ====

# create distinct variable names for joining SNP and allSites datasets 
# add prefix _SNPs to the variables from the SNP dataset
Px2014Latex <- Px2014SummarySNPs %>%
  set_names(paste(names(.), "SNPs", sep = "_")) %>%
  dplyr::rename(pop = pop_SNPs) %>%
  full_join(Px2014SummaryAllsites, by = "pop") %>%
  # arrange the table by state then location alphabetically
  left_join(popsMaster %>% 
              dplyr::select(pop = popString, Location, State), 
            by = "pop") %>%
  arrange(State, Location) %>%
  dplyr::select(pop, # as a check only
                Population = popName,
                # variables from the allsites dataset
                indPerLoc, nSites, idepth, snp_count, indel_count, private, Ho, Hs, Fis, 
                # variables from the SNP dataset
                indPerLoc_SNPs, idepth_SNPs, Ho_SNPs, Hs_SNPs, Fis_SNPs)

Px2014Digits <- c(0, # rownames, whether printed or not
                  0, # popnames
                  1, 0, 0, 0, 0, 0, # numeric variables start here 
                  4, 4, 4, # Ho, Hs, Fis
                  1, 0, 
                  4, 4, 4) # Ho, Hs, Fis

Px2014Align <- c("l", "l\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=1]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=6,table-figures-decimal=0]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=0]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=4,table-figures-decimal=0]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=4,table-figures-decimal=0]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=0]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=0,table-figures-decimal=4]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=0,table-figures-decimal=4]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=0,table-figures-decimal=4]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=1]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=0]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=0,table-figures-decimal=4]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=0,table-figures-decimal=4]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=0,table-figures-decimal=4]\n")

Px2014Caption <- "Welcome to 2014."
# Px2014Caption <- "Genetic diversity statistics for 31 populations of \\textit{P. xylostella} collected from Australia in 2014 
# based on 1008 SNP sites. 
# Statistics include population means for the number of individuals genotyped per locus (\\textit{n}), 
# site depth, observed heterozygosity ($H_{\\textsc{o}}$), gene diversity ($H_{\\textsc{s}}$) 
# and Nei's inbreeding coefficient, $F_{\\textsc{is}}$."

# make a multirow header using the add.to.row() argument
# wrap S column headers in curly braces, so it doesn't try to treat as a number
headerNames <- c(
  # All sites data headings
  "Population", "{$N$}", "{Sites}", "{\\makecell{Site\\\\depth}}", "{SNPs}", "{Indels}", "{\\makecell{Private\\\\sites}}",
  "{$H_{\\textsc{o}}$}", "{$H_{\\textsc{s}}$}", "{$F_{\\textsc{is}}$}",
  # SNP data headings
  "{$N$}", "{\\makecell{Site\\\\depth}}", "{$H_{\\textsc{o}}$}", "{$H_{\\textsc{s}}$}", "{$F_{\\textsc{is}}$}"
  ) %>%
  paste(collapse = " & ")

Px2014Header <- list()
Px2014Header$pos <- list(0, 0)
Px2014Header$command <- c(
  " & \\multicolumn{9}{c}{All variant and invariant sites}",
  " & \\multicolumn{5}{c}{1032 SNP variants}\\\\\n",
  "\\cmidrule(lr){2-10} ", "\\cmidrule(lr){11-15}\n"
  ) %>% paste(collapse = "") %>%
  c(headerNames %>% paste0("\\\\\n"))

# export Latex table code
Px2014Latex %>%
  dplyr::select(-pop) %>%
  xtable(caption = Px2014Caption,
         align   = Px2014Align,
         digits  = Px2014Digits,
         lab = "tab:Px2014PopulationStatistics") %>% 
  print.xtable(floating = TRUE,
               include.rownames = FALSE,
               include.colnames = FALSE,
               table.placement = "p",
               caption.placement = "top",
               add.to.row = Px2014Header,
               NA.string = "",
               scalebox = 0.5,
               booktabs = TRUE,
               sanitize.text.function = function(x){x})



# ====
# 2015
# ====

# create distinct variable names for joining SNP and allSites datasets 
# add prefix _SNPs to the variables from the SNP dataset
Px2015Latex <- Px2015SummarySNPs %>%
  set_names(paste(names(.), "SNPs", sep = "_")) %>%
  dplyr::rename(pop = pop_SNPs) %>%
  full_join(Px2015SummaryAllsites, by = "pop") %>%
  # arrange the table by state then location alphabetically
  left_join(popsMaster %>% 
              dplyr::select(pop = popString, Location, State), 
            by = "pop") %>%
  arrange(State, Location) %>%
  dplyr::select(pop, # as a check only
                Population = popName,
                # variables from the allsites dataset
                indPerLoc, nSites, idepth, snp_count, indel_count, private, Ho, Hs, Fis, 
                # variables from the SNP dataset
                indPerLoc_SNPs, idepth_SNPs, Ho_SNPs, Hs_SNPs, Fis_SNPs)

Px2015Digits <- c(0, # rownames, whether printed or not
                  0, # popnames
                  1, 0, 0, 0, 0, 0, # numeric variables start here 
                  4, 4, 4, # Ho, Hs, Fis
                  1, 0, 
                  4, 4, 4) # Ho, Hs, Fis

Px2015Align <- c("l", "l\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=1]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=6,table-figures-decimal=0]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=0]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=4,table-figures-decimal=0]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=4,table-figures-decimal=0]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=0]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=0,table-figures-decimal=4]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=0,table-figures-decimal=4]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=0,table-figures-decimal=4]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=1]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=0]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=0,table-figures-decimal=4]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=0,table-figures-decimal=4]\n",
                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-integer=0,table-figures-decimal=4]\n")

Px2015Caption <- "Welcome to 2015."
# Px2015Caption <- "Genetic diversity statistics for 31 populations of \\textit{P. xylostella} collected from Australia in 2014 
# based on 1008 SNP sites. 
# Statistics include population means for the number of individuals genotyped per locus (\\textit{n}), 
# site depth, observed heterozygosity ($H_{\\textsc{o}}$), gene diversity ($H_{\\textsc{s}}$) 
# and Nei's inbreeding coefficient, $F_{\\textsc{is}}$."

# make a multirow header using the add.to.row() argument
# wrap S column headers in curly braces, so it doesn't try to treat as a number
headerNames <- c(
  # All sites data headings
  "Population", "{$N$}", "{Sites}", "{\\makecell{Site\\\\depth}}", "{SNPs}", "{Indels}", "{\\makecell{Private\\\\sites}}",
  "{$H_{\\textsc{o}}$}", "{$H_{\\textsc{s}}$}", "{$F_{\\textsc{is}}$}",
  # SNP data headings
  "{$N$}", "{\\makecell{Site\\\\depth}}", "{$H_{\\textsc{o}}$}", "{$H_{\\textsc{s}}$}", "{$F_{\\textsc{is}}$}"
) %>%
  paste(collapse = " & ")

Px2015Header <- list()
Px2015Header$pos <- list(0, 0)
Px2015Header$command <- c(
  " & \\multicolumn{9}{c}{All variant and invariant sites}",
  " & \\multicolumn{5}{c}{1032 SNP variants}\\\\\n",
  "\\cmidrule(lr){2-10} ", "\\cmidrule(lr){11-15}\n"
) %>% paste(collapse = "") %>%
  c(headerNames %>% paste0("\\\\\n"))

# export Latex table code
Px2015Latex %>%
  dplyr::select(-pop) %>%
  xtable(caption = Px2015Caption,
         align   = Px2015Align,
         digits  = Px2015Digits,
         lab = "tab:Px2015PopulationStatistics") %>% 
  print.xtable(floating = TRUE,
               include.rownames = FALSE,
               include.colnames = FALSE,
               table.placement = "p",
               caption.placement = "top",
               add.to.row = Px2015Header,
               NA.string = "",
               scalebox = 0.5,
               booktabs = TRUE,
               sanitize.text.function = function(x){x})


# =======================================================
# Box plot for gene diversity (= expected heterozygosity)
# =======================================================
# 28/4/2018
# some indviduals aong our samples appear to be related.
# Observed het among pops is similar, however \
# some pops (e.g. southend) show reduced gene diversity.
# Therefore, we'll visualise gene diversity using boxplots

# read in the 2014 data
Px2014VcftoolsHet <- file.path(
  popStatsPath, "Px.31pops.434ind.2014.genotypeGVCFs.allSites.hardFilt.maxMiss0.8.het") %>%
  read.table(header = TRUE) %>%
  set_names(c("ind", "homO", "homE", "nSites", "Fis")) %>%
  mutate(pop = extractPop(ind))

Px2014VcftoolsHet %>%
  ggplot(aes(x = pop, y = homO)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# read in the 2015 data
Px2015VcftoolsHet <- file.path(
  popStatsPath, "Px.28pops.399ind.2015.genotypeGVCFs.allSites.hardFilt.maxMiss0.8.het") %>%
  read.table(header = TRUE) %>%
  set_names(c("ind", "homO", "homE", "nSites", "Fis")) %>%
  mutate(pop = extractPop(ind))

Px2015VcftoolsHet %>%
  ggplot(aes(x = pop, y = homE)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))  
#> no differences on expected or observed homozygosity on a per-individual basis.
# ###############################################################################

















# # OLD BELOW
# # #############
# # SNPs
# # #############
# 
# # read in hierfstat derived population summary (nei's Fis)
# px31HfSNP <- file.path(hrfPath, "px31pops2014.SNPs.hierfstat.basicstats.pop.txt") %>%
#   read.table(header = TRUE)
# 
# # read in genepop derived gene diversities
# px31GpSNP <- file.path(popStatsPath, "genepopGeneDiversitiesPops.DIV.xlsx") %>%
#   read_excel(sheet = "px31popsSNPs", skip = 5)
# 
# # read in vcftools indv site depth summary
# px31VdSNP <- file.path(vfPath, "Px.59pops.842ind.genotypeGVCFs.variants.hardFilt.maxMiss0.7.idepth") %>%
#   read.table(header = TRUE) %>%
#   filter(INDV %in% pxInd2014)
# 
# # read in vcfstats summary
# px31VsSNP <- file.path(vfPath, "Px.31pops.440ind.2014.genotypeGVCFs.variants.hardFilt.maxMiss0.7.vcfstats.grepSummary.txt") %>%
#   read.table(header = TRUE)
# 
# # merge results into a dataframe summarising population statistics
# px31StatsSNP <- px31VsSNP %>%
#   merge(px31VdSNP %>% setNames(c("ind", "sites", "idepth")), by = "ind") %>%
#   mutate(ind = renamePxToPaus(ind),
#          pop = extractPop(ind)) %>%
#   group_by(pop) %>%
#   summarise(nInd   = length(ind),
#             nSites = round(mean(sites), 1),
#             pSites = round(nSites / 1008, 2),
#             idepth = round(mean(idepth), 1),
#             homAA  = round(mean(homAA), 1),
#             hetRA  = round(mean(hetRA), 1),
#             homRR  = round(mean(homRR), 1),
#             ref_count = round(mean(ref_count), 1),
#             snp_count = round(mean(snp_count), 1)) %>% 
#   merge(px31HfSNP %>% mutate(pop = extractPop(pop)), by = "pop") %>%
#   merge(px31GpSNP %>% mutate(pop = extractPop(Sample)) %>% 
#           rename(gpFis = Fis), by = "pop") %>%
#   merge(popLabs) %>%
#   merge(popsOrderPx %>%
#           rename(pop = popString)) %>%
#   arrange(orderYearGeo)
# 
# px31StatsSNP %>%
#   ggplot(aes(x = Ho, Hs)) +
#   geom_point() +
#   geom_abline(intercept = 0, slope = 1)
# 
# # format for latex
# # ... # wrap S column headers in curly braces so siuntix treats as text
# px31SNPLatex <- px31StatsSNP %>%
#   dplyr::select(Population = popS,
#                 `{\\textit{n}}` = indPerLoc, 
#                 `{\\makecell{No. \\\\ sites}}` = nSites,
#                 `{\\makecell{Site \\\\ depth}}` = idepth,
#                 `$H_{\\textsc{o}}$` = Ho,
#                 `$H_{\\textsc{s}}$` = Hs,
#                 `$F_{\\textsc{is}}$` = Fis)
# 
# px31SNPAlg <- c("l", "l", 
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=1,table-figures-integer=2]",
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=3]",
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]",
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]", 
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]", 
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]")
# px31SNPDig <- c(0, 0, 1, 0, 0, 4, 4, 4)
# px31SNPCap <- "Genetic diversity statistics for 31 populations of \\textit{P. xylostella} collected from Australia in 2014 
# based on 1008 SNP sites. 
# Statistics include population means for the number of individuals genotyped per locus (\\textit{n}), 
# site depth, observed heterozygosity ($H_{\\textsc{o}}$), gene diversity ($H_{\\textsc{s}}$) 
# and Nei's inbreeding coefficient, $F_{\\textsc{is}}$."
# 
# xtable(px31SNPLatex,
#        caption = px31SNPCap,
#        align   = px31SNPAlg,
#        digits  = px31SNPDig,
#        lab = "tab:px31PopstatsSNP") %>% 
#   print.xtable(floating = TRUE,
#                include.rownames = FALSE,
#                include.colnames = TRUE,
#                table.placement = "p",
#                caption.placement = "top",
#                NA.string = "",
#                scalebox = 0.8,
#                booktabs = TRUE,
#                sanitize.text.function = function(x){x})
# 
# 
# 
# # #############
# # All sites
# # #############
# 
# # read in hierfstat derived population summary (nei's Fis)
# px31HfAS <- file.path(hfPath, "px31pops2014.allSites.hierfstat.basicstats.pop.txt") %>%
#   read.table(header = TRUE) # ... DONE
# 
# # read in genepop derived gene diversities
# px31GpAS <- file.path(popStatsPath, "genepopGeneDiversitiesPops.DIV.xlsx") %>%
#   read_excel(sheet = "px31popsAllsites", skip = 5) # ... DONE
# 
# # read in vcftools indv site depth summary
# px31VdAS <- file.path(vfPath, "Px.31pops.440ind.2014.genotypeGVCFs.allSites.hardFilt.maxMiss0.7.idepth") %>%
#   read.table(header = TRUE) # ... DONE
# 
# # read in vcfstats summary
# px31VsAS <- file.path(vfPath, "Px.31pops.440ind.2014.genotypeGVCFs.allSites.hardFilt.maxMiss0.7.vcfstats.grepSummary.txt") %>%
#   read.table(header = TRUE)
# 
# # merge results into a dataframe summarising population statistics
# px31StatsAS <- px31VsAS %>%
#   merge(px31VdAS %>% setNames(c("ind", "sites", "idepth")), by = "ind") %>%
#   mutate(ind = renamePxToPaus(ind), # not needed, but no harm
#          pop = extractPop(ind)) %>%
#   group_by(pop) %>%
#   summarise(nInd   = length(ind),
#             nSites = round(mean(sites), 1),
#             pSites = round(nSites / 680919, 2),
#             idepth = round(mean(idepth), 1),
#             homAA  = round(mean(homAA), 1),
#             hetRA  = round(mean(hetRA), 1),
#             homRR  = round(mean(homRR), 1),
#             ref_count =   round(mean(ref_count), 1),
#             snp_count =   round(mean(snp_count), 1),
#             indel_count = round(mean(indel_count), 1),
#             private = round(mean(private), 1)) %>% 
#   merge(px31HfAS %>% mutate(pop = extractPop(pop)), by = "pop") %>%
#   merge(px31GpAS %>% mutate(pop = extractPop(Sample)) %>% 
#           rename(gpFis = Fis), by = "pop") %>%
#   merge(popLabs) %>%
#   merge(popsOrderPx %>%
#           rename(pop = popString)) %>%
#   arrange(orderYearGeo)
# 
# px31StatsAS %>%
#   ggplot(aes(x = Hs, y = Ho)) +
#   geom_point() +
#   geom_abline(intercept = 0, slope = 1)
# # observed het mostly lower than expected
# 
# # calculate site depth per individual for manuscript
# px31VdAS %>%
#   summarise(n = nrow(.),
#             mnDepthPerInd = mean(MEAN_DEPTH))
# # n mnDepthPerInd
# # 1 440      32.10415
# 
# 
# 
# # format for latex
# px31ASLatex <- px31StatsAS %>%
#   dplyr::select(Population = popS,
#                 #`\\multicolumn{1}{l}{Species}` = species,
#                 `{\\textit{n}}` = indPerLoc, 
#                 `{\\makecell{No.\\\\sites}}` = nSites, # wrap S column headers in curly braces, so it doesn't try to treat as a number
#                 `{\\makecell{Site\\\\depth}}` = idepth,
#                 `{SNPs}`= snp_count,
#                 `{Indels}` = indel_count,
#                 `{\\makecell{Private\\\\sites}}` = private,
#                 `$H_{\\textsc{o}}$` = Ho,
#                 `$H_{\\textsc{s}}$` = Hs,
#                 `$F_{\\textsc{is}}$` = Fis)
# 
# px31ASAlg <- c("l", "l",
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=3]", # center align S columns without decimals by pasting this option onto relevant column aligns, 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=6]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=4]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=4]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=3]",
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]")
# px31ASDig <- c(0, 0, 1, 0, 0, 0, 0, 0, 4, 4, 4)
# px31ASCap <- "Genetic diversity statistics for 31 populations of \\textit{P. xylostella} collected from Australia in 2014 
# based on \\num{667243} confidently called variant and invariant sites. 
# Statistics include population means for the number of individuals genotyped per locus (\\textit{n}), 
# site depth, observed heterozygosity ($H_{\\textsc{o}}$), gene diversity ($H_{\\textsc{s}}$), the number of 
# SNP and indel and private sites, and Nei's inbreeding coefficient, $F_{\\textsc{is}}$."
# 
# xtable(px31ASLatex,
#        caption = px31ASCap,
#        align   = px31ASAlg,
#        digits  = px31ASDig,
#        lab = "tab:px31PopstatsAS") %>% 
#   print.xtable(floating = TRUE,
#                include.rownames = FALSE,
#                include.colnames = TRUE,
#                table.placement = "p",
#                caption.placement = "top",
#                NA.string = "",
#                scalebox = 0.8,
#                booktabs = TRUE,
#                sanitize.text.function = function(x){x})  
# 
# 
# 
# # ##################################################################################
# # generate a population statistics table for dataset: px, 28pops, 402 ind, 2015 
# # ##################################################################################
# 
# # a list of individuals from 2015
# pxInd2015  <- pxSamples[str_sub(pxSamples, 5, 6) == "15"]
# 
# # ###############################
# # px 28pops, 402 ind, 1008  SNPs
# # ###############################
# 
# # read in hierfstat derived population summary (nei's Fis)
# px28HfSNP <- file.path(hfPath, "px28pops2015.SNPs.hierfstat.basicstats.pop.txt") %>%
#   read.table(header = TRUE)
# 
# # read in genepop derived gene diversities
# px28GpSNP <- file.path(popStatsPath, "genepopGeneDiversitiesPops.DIV.xlsx") %>%
#   read_excel(sheet = "px28popsSNPs", skip = 5)
# 
# # read in vcftools indv site depth summary
# px28VdSNP <- file.path(vfPath, "Px.59pops.842ind.genotypeGVCFs.variants.hardFilt.maxMiss0.7.idepth") %>%
#   read.table(header = TRUE) %>% 
#   filter(INDV %in% pxInd2015)
# 
# # read in vcfstats summary
# px28VsSNP <- file.path(vfPath, "Px.28pops.402ind.2015.genotypeGVCFs.variants.hardFilt.maxMiss0.7.vcfstats.grepSummary.txt") %>%
#   read.table(header = TRUE)
# 
# # merge results into a dataframe summarising population statistics
# px28StatsSNP <- px28VsSNP %>%
#   merge(px28VdSNP %>% setNames(c("ind", "sites", "idepth")), by = "ind") %>%
#   mutate(ind = renamePxToPaus(ind),
#          pop = extractPop(ind)) %>%
#   group_by(pop) %>%
#   summarise(nInd   = length(ind),
#             nSites = round(mean(sites), 1),
#             pSites = round(nSites / 1008, 2),
#             idepth = round(mean(idepth), 1),
#             homAA  = round(mean(homAA), 1),
#             hetRA  = round(mean(hetRA), 1),
#             homRR  = round(mean(homRR), 1),
#             ref_count = round(mean(ref_count), 1),
#             snp_count = round(mean(snp_count), 1)) %>% 
#   merge(px28HfSNP %>% 
#           mutate(pop = extractPop(pop)), 
#         by = "pop") %>%
#   merge(px28GpSNP %>% 
#           mutate(pop = extractPop(Sample)) %>% 
#           rename(gpFis = Fis), 
#         by = "pop") %>%
#   merge(popLabs) %>%
#   merge(popsOrderPx %>%
#           rename(pop = popString)) %>%
#   arrange(orderYearGeo)
# 
# px28StatsSNP %>%
#   ggplot(aes(x = Hs, y = Ho)) +
#   geom_point() +
#   geom_abline(intercept = 0, slope = 1)
# 
# 
# 
# # format for latex
# # ... # wrap S column headers in curly braces so siuntix treats as text
# px28SNPLatex <- px28StatsSNP %>%
#   dplyr::select(Population = popS,
#                 `{\\textit{n}}` = indPerLoc, 
#                 `{\\makecell{No. \\\\ sites}}` = nSites,
#                 `{\\makecell{Site \\\\ depth}}` = idepth,
#                 `$H_{\\textsc{o}}$` = Ho,
#                 `$H_{\\textsc{s}}$` = Hs,
#                 `$F_{\\textsc{is}}$` = Fis)
# 
# px28SNPAlg <- c("l", "l", 
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=1,table-figures-integer=2]",
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=3]",
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]",
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]", 
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]", 
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]")
# px28SNPDig <- c(0, 0, 1, 0, 0, 4, 4, 4)
# px28SNPCap <- "Genetic diversity statistics for 28 populations of \\textit{P. xylostella} collected from Australia in 2015 
# based on 1008 SNP sites. 
# Statistics include population means for the number of individuals genotyped per locus (\\textit{n}), 
# site depth, observed heterozygosity ($H_{\\textsc{o}}$), gene diversity ($H_{\\textsc{s}}$) 
# and Nei's inbreeding coefficient, $F_{\\textsc{is}}$."
# 
# xtable(px28SNPLatex,
#        caption = px28SNPCap,
#        align   = px28SNPAlg,
#        digits  = px28SNPDig,
#        lab = "tab:px28PopstatsSNP") %>% 
#   print.xtable(floating = TRUE,
#                include.rownames = FALSE,
#                include.colnames = TRUE,
#                table.placement = "p",
#                caption.placement = "top",
#                NA.string = "",
#                scalebox = 0.8,
#                booktabs = TRUE,
#                sanitize.text.function = function(x){x})
# 
# 
# 
# # ###############################
# # px 28pops, 402 ind, allSites
# # ###############################
# 
# # read in hierfstat derived population summary (nei's Fis)
# 
# px28HfAS <- file.path(hfPath, "px28pops2015.allSites.hierfstat.basicstats.pop.txt") %>%
#   read.table(header = TRUE) # ... DONE
# 
# # read in genepop derived gene diversities
# px28GpAS <- file.path(popStatsPath, "genepopGeneDiversitiesPops.DIV.xlsx") %>%
#   read_excel(sheet = "px28popsAllsites", skip = 5) # ... DONE
# 
# # read in vcftools indv site depth summary
# px28VdAS <- file.path(vfPath, "Px.28pops.402ind.2015.genotypeGVCFs.allSites.hardFilt.maxMiss0.7.idepth") %>%
#   read.table(header = TRUE) # ... DONE
# 
# # HERE.
# # read in vcfstats summary
# px28VsAS <- file.path(vfPath, "Px.28pops.402ind.2015.genotypeGVCFs.allSites.hardFilt.maxMiss0.7.vcfstats.grepSummary.txt") %>%
#   read.table(header = TRUE)
# 
# # merge results into a dataframe summarising population statistics
# px28StatsAS <- px28VsAS %>%
#   merge(px28VdAS %>% setNames(c("ind", "sites", "idepth")), by = "ind") %>%
#   mutate(ind = renamePxToPaus(ind), # not needed, but no harm
#          pop = extractPop(ind)) %>%
#   group_by(pop) %>%
#   summarise(nInd   = length(ind),
#             nSites = round(mean(sites), 1),
#             pSites = round(nSites / 653166, 2),
#             idepth = round(mean(idepth), 1),
#             homAA  = round(mean(homAA), 1),
#             hetRA  = round(mean(hetRA), 1),
#             homRR  = round(mean(homRR), 1),
#             ref_count =   round(mean(ref_count), 1),
#             snp_count =   round(mean(snp_count), 1),
#             indel_count = round(mean(indel_count), 1),
#             private = round(mean(private), 1)) %>% 
#   merge(px28HfAS %>% mutate(pop = extractPop(pop)), by = "pop") %>%
#   merge(px28GpAS %>% mutate(pop = extractPop(Sample)) %>% 
#           rename(gpFis = Fis), by = "pop") %>%
#   merge(popLabs) %>%
#   merge(popsOrderPx %>%
#           rename(pop = popString)) %>%
#   arrange(orderYearGeo)
# 
# px28StatsAS %>%
#   ggplot(aes(x = Hs, y = Ho)) +
#   geom_point() +
#   geom_abline(intercept = 0, slope = 1)
# #> observed lower than expected
# 
# 
# # format for latex
# px28ASLatex <- px28StatsAS %>%
#   dplyr::select(Population = popS,
#                 `{\\textit{n}}` = indPerLoc, 
#                 `{\\makecell{No.\\\\sites}}` = nSites, # wrap S column headers in curly braces, so it doesn't try to treat as a number
#                 `{\\makecell{Site\\\\depth}}` = idepth,
#                 `{SNPs}`= snp_count,
#                 `{Indels}` = indel_count,
#                 `{\\makecell{Private\\\\sites}}` = private,
#                 `$H_{\\textsc{o}}$` = Ho,
#                 `$H_{\\textsc{s}}$` = Hs,
#                 `$F_{\\textsc{is}}$` = Fis)
# 
# px28ASAlg <- c("l", "l",
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=3]", # center align S columns without decimals by pasting this option onto relevant column aligns, 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=6]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=4]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=4]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=3]",
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]")
# px28ASDig <- c(0, 0, 1, 0, 0, 0, 0, 0, 4, 4, 4)
# px28ASCap <- "Genetic diversity statistics for 28 populations of \\textit{P. xylostella} collected from Australia in 2015 
# based on \\num{667243} confidently called variant and invariant sites. 
# Statistics include population means for the number of individuals genotyped per locus (\\textit{n}), 
# site depth, observed heterozygosity ($H_{\\textsc{o}}$), gene diversity ($H_{\\textsc{s}}$), the number of 
# SNP and indel and private sites, and Nei's inbreeding coefficient, $F_{\\textsc{is}}$."
# 
# xtable(px28ASLatex,
#        caption = px28ASCap,
#        align   = px28ASAlg,
#        digits  = px28ASDig,
#        lab = "tab:px28PopstatsAS") %>% 
#   print.xtable(floating = TRUE,
#                include.rownames = FALSE,
#                include.colnames = TRUE,
#                table.placement = "p",
#                caption.placement = "top",
#                NA.string = "",
#                scalebox = 0.8,
#                booktabs = TRUE,
#                sanitize.text.function = function(x){x})  
# 
# # #########################################################################
# # Calculate individual site depth for P. xylostella samples, for manuscript
# # 4/11/2017
# # 1008 SNPS
# # read in the idepth file for 1008 SNPs.
# px59VdVR <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "populationStatistics",
#                       "Px.59pops.842ind.genotypeGVCFs.variants.hardFilt.maxMiss0.7.idepth") %>%
#   read.table(header = TRUE)
# # calculate mean depth per individual
# px59VdVR %>%
#   summarise(mean(MEAN_DEPTH))
# 
# # merge the .idepth files for 2014 and 2015 (these were derived from the hard filtered vcf with 842 ind anyway)
# # NOTE: this is interim. I need to re-generate the idepth files using consistent sites across years.
# px31VdAS %>% 
#   bind_rows(px28VdAS) %>%
#   summarise(depthPerInd = mean(MEAN_DEPTH))
# # 1    32.94992
# 
# 
# # ########################################################################
# # plot observed vs expected heterozygosity for the allSites and SNP files
# # ... this will be used to visualise the statistics tables
# # ... if informative, a publication figure (prob not)
# # ########################################################################
# 
# fmtLabels <- function(x) sprintf("%.4f", x)
# 
# # allSites dataset
# pxStatsAS <- px31StatsAS %>%
#   bind_rows(px28StatsAS) 
# minHet <- min(c(pxStatsAS$Ho, pxStatsAS$Hs))
# maxHet <- max(c(pxStatsAS$Ho, pxStatsAS$Hs))
# 
# pxStatsAS %>%
#   mutate(Year = factor(Year)) %>%
#   ggplot(aes(x = Hs, y = Ho, colour = Year)) +
#   geom_point() +
#   facet_wrap(~Year) +
#   scale_x_continuous(limits = c(minHet, maxHet),
#                      labels = fmtLabels) +
#   scale_y_continuous(limits = c(minHet, maxHet),
#                      labels = fmtLabels) +
#   geom_abline(slope = 1, lty = 5) +
#   coord_equal()
# 
# # SNP dataset
# pxStatsSNP <- px31StatsSNP %>%
#   bind_rows(px28StatsSNP)
# 
# minHet2 <- min(c(pxStatsSNP$Ho, pxStatsSNP$Hs))
# maxHet2 <- max(c(pxStatsSNP$Ho, pxStatsSNP$Hs))
# pxStatsSNP %>%
#   mutate(Year = factor(Year)) %>%
#   ggplot(aes(x = Hs, y = Ho, colour = Year)) +
#   geom_point() +
#   facet_wrap(~Year) +
#   scale_x_continuous(limits = c(minHet2, maxHet2),
#                      labels = fmtLabels) +
#   scale_y_continuous(limits = c(minHet2, maxHet2),
#                      labels = fmtLabels) +
#   geom_abline(slope = 1, lty = 5) +
#   coord_equal()
# 
# # now look at the Px Pa dataset
# minHet3 <- min(c(pxpaStatsAS$Ho, pxpaStatsAS$Hs))
# maxHet3 <- max(c(pxpaStatsAS$Ho, pxpaStatsAS$Hs))
# pxpaStatsAS %>%
#   ggplot(aes(x = Hs, y = Ho, colour = species)) +
#   geom_point() +
#   #facet_wrap(~species) +
#   scale_x_continuous(limits = c(minHet3, maxHet3),
#                      labels = fmtLabels) +
#   scale_y_continuous(limits = c(minHet3, maxHet3),
#                      labels = fmtLabels) +
#   geom_abline(slope = 1, lty = 5) +
#   coord_equal()
# 
# minHet4 <- min(c(pxpaStatsSNP$Ho, pxpaStatsSNP$Hs))
# maxHet4 <- max(c(pxpaStatsSNP$Ho, pxpaStatsSNP$Hs))
# pxpaStatsSNP %>%
#   ggplot(aes(x = Hs, y = Ho, colour = species)) +
#   geom_point() +
#   #facet_wrap(~species) +
#   scale_x_continuous(limits = c(minHet4, maxHet4),
#                      labels = fmtLabels) +
#   scale_y_continuous(limits = c(minHet4, maxHet4),
#                      labels = fmtLabels) +
#   geom_abline(slope = 1, lty = 5) +
#   coord_equal()
# 
# # summarise range of Fis values
# pxStatsAS %>% 
#   group_by(Year) %>%
#   summarise(minFis = round(min(Fis), 4),
#             maxFis = round(max(Fis), 4))
# 
# pxStatsSNP %>% 
#   group_by(Year) %>%
#   summarise(minFis = round(min(Fis), 4),
#             maxFis = round(max(Fis), 4))
# 
# # summarise the range of heterozygosities in the allsites datasets
# pxStatsAS %>%
#   group_by(Year) %>%
#   summarise(minHetO = round(min(Ho), 4),
#             maxHetO = round(max(Ho), 4)) 
# 
# 
# # summarise the range of private sites (excluding the two outliers, for purposes of discussing this point in the manuscript) 
# pxStatsAS %>%
#   filter(!pop %in% c("werN15svx", "picS14awx")) %>%
#   group_by(Year) %>%
#   summarise(minPrivate = min(private),
#             maxPrivate = max(private)) 
# 
# # ############################################################################################################
# # Create a box plot of the proportion of heterozygous sites: Px, 2014 440 ind, and Px 2015 402 ind, all sites
# # ... derived from vcftools, vcftool --het and output of python script
# # ... this is made from the allsites excluding indels, and the output of python script
# # ... check for differences in observed het by site
# # ############################################################################################################
# 
# old <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "populationStatistics")
# outThesis <- file.path("C:", "UserData", "Kym", "PhD", "thesis", "images") 
# 
# pyStatsPx28 <- file.path(
#   vfPath, 
#   "Px.28pops.402ind.2015.genotypeGVCFs.allSites.hardFilt.maxMiss0.7.rmIndels.recode.pyStats.txt"
#   ) %>% 
#   read.table(skip = 423, header = FALSE, stringsAsFactors = FALSE) %>%
#   setNames(c("ind", "nCalls", "percent", "hom", "het"))
# 
# pyStatsPx31 <- file.path(
#   vfPath, 
#   "Px.31pops.440ind.2014.genotypeGVCFs.allSites.hardFilt.maxMiss0.7.rmIndels.recode.pyStats.txt"
# ) %>% 
#   read.table(skip = 461, header = FALSE, stringsAsFactors = FALSE) %>%
#   setNames(c("ind", "nCalls", "percent", "hom", "het"))
# 
# file.path(old, "hetBoxPxStudyDRAFT.pdf") %>%
#   pdf(height = 8, width = 12)
# pyStatsPx31 %>%
#   bind_rows(pyStatsPx28) %>%
#   mutate(pop = extractPop(ind),
#          host = extractHost(ind),
#          prHet = het / nCalls) %>%
#   merge(popsOrderPx %>%
#           rename(pop = popString)) %>%
#   # add host order
#   arrange(host, prHet) %>%
#   mutate(orderHost = 1:nrow(.)) %>%
#   mutate(pop = reorder(pop, .$orderHost)) %>% # it won't reorder by orderHost directly within the same call to mutate()
#   ggplot(aes(x = pop, y = prHet, fill = host)) +
#   geom_boxplot() +
#     #outlier.shape = NA) + # to avoid duplicating outlier points 
#   #geom_jitter(width = 0.1) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#   labs(x = NULL, y = "Proportion of heterozygous bi-allelic sites")
# dev.off()
#   
# # ####################################################
# # # find the 9 outlier `high het samples for analysis#
# # 10/12/2017
# highHet9 <- pyStatsPx31 %>%
#   bind_rows(pyStatsPx28) %>%
#   mutate(pop = extractPop(ind),
#          host = extractHost(ind),
#          prHet = het / nCalls) %>%
#   arrange(desc(prHet)) %>%
#   filter(prHet > 0.01)
# 
# # find 9 random samples with `normal` levels of het
# normHet <- pyStatsPx31 %>%
#   bind_rows(pyStatsPx28) %>%
#   mutate(pop = extractPop(ind),
#          host = extractHost(ind),
#          prHet = het / nCalls) %>%
#   filter(prHet < 0.01)
# 
# # grab 9 random samples with normal Het.
# normHet9 <- sample(normHet$ind, size = 9, replace = FALSE)
# normHet9Df <- normHet[grepl(paste(normHet9, collapse = "|"), normHet$ind), ]
# normHet9Df
# # End grab high het samples
# ####################################################################


# plot ordered by year and then geographically from west to east

pdf("hetBoxPx_OrderedByYearGeo_DRAFT.pdf", width = 12, height = 8)
# #greypal <-brewer.pal(9, "Greys")
ggplot(indStats842, aes(x = reorder(popSY, indStats842$orderYearGeo), 
                        y = propHet, fill = host)) + 
  #geom_boxplot(position = position_dodge(width = dw)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  labs(x = NULL, y = "Proportion of heterozygous sites") +
  theme_bw(base_size = 10) +
  theme(axis.title.x = element_text(face = "bold", 
                                    colour = 'black', 
                                    margin = margin(t = 20)), # move text away from axis
        axis.title.y = element_text(face = "bold", colour = 'black', margin = margin(r = 20)),
        axis.text.x  = element_text(colour = 'black', angle = 60, hjust = 1),
        axis.text.y  = element_text(colour = 'black'),
        #panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        #axis.line = element_line(size = 0.6, color = "black"),
        legend.text  = element_text(face = "italic")) +
        #legend.position = c(0.85, 0.94)) +
  scale_fill_brewer(palette = "Accent")
dev.off()

# plot ordered by host, then geographically from west to east
pdf("hetBoxPx_OrderedByHostGeo_DRAFT.pdf", width = 12, height = 8)
# #greypal <-brewer.pal(9, "Greys")
ggplot(indStats842, aes(x = reorder(popSY, indStats842$orderHostGeo) , 
                        y = propHet, fill = host)) + 
  #geom_boxplot(position = position_dodge(width = dw)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  labs(x = NULL, y = "Proportion of heterozygous sites") +
  theme_bw(base_size = 10) +
  theme(axis.title.x = element_text(face = "bold", 
                                    colour = 'black', 
                                    margin = margin(t = 20)), # move text away from axis
        axis.title.y = element_text(face = "bold", colour = 'black', margin = margin(r = 20)),
        axis.text.x  = element_text(colour = 'black', angle = 60, hjust = 1),
        axis.text.y  = element_text(colour = 'black'),
        #panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        #axis.line = element_line(size = 0.6, color = "black"),
        legend.text  = element_text(face = "italic")) +
  #legend.position = c(0.85, 0.94)) +
  scale_fill_brewer(palette = "Accent")
dev.off()

# a plot ordered by host then prop het from smallest to largest
pdf("hetBoxPx_OrderedByHostHet_DRAFT.pdf", width = 12, height = 8)
# #greypal <-brewer.pal(9, "Greys")
ggplot(indStats842, aes(x = reorder(popSY, indStats842$orderHostHet) , 
                        y = propHet, fill = host)) + 
  #geom_boxplot(position = position_dodge(width = dw)) + 
  geom_boxplot() +
  #geom_jitter() +
  labs(x = NULL, y = "Proportion of heterozygous sites") +
  theme_bw(base_size = 10) +
  theme(axis.title.x = element_text(face = "bold", 
                                    colour = 'black', 
                                    margin = margin(t = 20)), # move text away from axis
        axis.title.y = element_text(face = "bold", colour = 'black', margin = margin(r = 20)),
        axis.text.x  = element_text(colour = 'black', angle = 60, hjust = 1),
        axis.text.y  = element_text(colour = 'black'),
        #panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        #axis.line = element_line(size = 0.6, color = "black"),
        legend.text  = element_text(face = "italic")) +
  #legend.position = c(0.85, 0.94)) +
  scale_fill_brewer(palette = "Accent")
dev.off()

# a plot grouped by host
pdf("hetBoxPx_GroupedByHost_DRAFT.pdf", width = 12, height = 8)
# #greypal <-brewer.pal(9, "Greys")
ggplot(indStats842, aes(x = host, y = propHet, fill = host)) + 
  #geom_boxplot(position = position_dodge(width = dw)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  labs(x = NULL, y = "Proportion of heterozygous sites") +
  theme_bw(base_size = 10) +
  theme(axis.title.x = element_text(face = "bold", 
                                    colour = 'black', 
                                    margin = margin(t = 20)), # move text away from axis
        axis.title.y = element_text(face = "bold", colour = 'black', margin = margin(r = 20)),
        axis.text.x  = element_text(colour = 'black', angle = 60, hjust = 1),
        axis.text.y  = element_text(colour = 'black'),
        #panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        #axis.line = element_line(size = 0.6, color = "black"),
        legend.text  = element_text(face = "italic")) +
  #legend.position = c(0.85, 0.94)) +
  scale_fill_brewer(palette = "Accent")
dev.off()

# merge with the COI/ px fixed SNPs data.

# # End script
# # #######################################################################
# 
# ################
# OLD CODE BELOW
# ###############

# # #########################################################################
# # Generate a population statistics table for dataset: PxPa, 5pops, 100ind #
# # #########################################################################
# 
# # #############
# # SNPs
# # #############
# 
# # read in hierfstat derived population summary (nei's Fis)
# pxpaHfSNP <- file.path(hfPath, "pxpa100ind707.hierfstat.basicstats.pop.txt") %>%
#   read.table(header = TRUE)
# 
# # read in genepop derived gene diversities
# pxpaGpSNP <- file.path(popStatsPath, "genepopGeneDiversitiesPops.DIV.xlsx") %>%
#   read_excel(sheet = "pxpa100ind707SNPs", skip = 5)
# 
# # read in vcftools indv site depth summary
# pxpaVdSNP <- file.path(vfPath, "pxpa.5pops.100ind.genotypeGVCFs.variants.hardFilt.maxMiss0.7.idepth") %>%
#   read.table(header = TRUE)
# 
# # read in vcfstats summary
# pxpaVsSNP <- file.path(vfPath, "pxpa.5pops.100ind.genotypeGVCFs.variants.hardFilt.maxMiss0.7.vcfstats.grepSummary.txt") %>%
#   read.table(header = TRUE)
# 
# # merge results into a dataframe summarising population statistics
# pxpaStatsSNP <- pxpaVsSNP %>%
#   merge(pxpaVdSNP %>% setNames(c("ind", "sites", "idepth")), by = "ind") %>%
#   mutate(ind = renamePxToPaus(ind),
#          pop = extractPop(ind),
#          species = extractSpecies(pop)) %>%
#   group_by(pop) %>%
#   summarise(species = unique(species),
#             nInd   = length(ind),
#             nSites = round(mean(sites), 1),
#             pSites = round(nSites / 707, 2),
#             idepth = round(mean(idepth), 1),
#             homAA  = round(mean(homAA), 1),
#             hetRA  = round(mean(hetRA), 1),
#             homRR  = round(mean(homRR), 1),
#             ref_count = round(mean(ref_count), 1),
#             snp_count = round(mean(snp_count), 1)) %>% 
#   merge(pxpaHfSNP %>% mutate(pop = extractPop(pop)), by = "pop") %>%
#   merge(pxpaGpSNP %>% mutate(pop = extractPop(Sample)) %>% 
#           rename(gpFis = Fis), by = "pop") %>%
#   merge(popLabs)
# 
# 
# # format for latex
# # ... # wrap S column headers in curly braces so siuntix treats as text
# pxpaSNPLatex <- pxpaStatsSNP %>%
#   dplyr::select(Population = popS,
#                 `\\multicolumn{1}{l}{Species}` = species,
#                 `{\\textit{n}}` = indPerLoc, 
#                 `{\\makecell{No. \\\\ sites}}` = nSites,
#                 `{\\makecell{Site \\\\ depth}}` = idepth,
#                 `$H_{\\textsc{o}}$` = Ho,
#                 `$H_{\\textsc{s}}$` = Hs,
#                 `$F_{\\textsc{is}}$` = Fis)
# 
# pxpaSNPLatex$Population[duplicated(pxpaSNPLatex$Population)] <- NA
# pxpaSNPAlg <- c("l", 
#                 "l", ">{\\itshape}l",
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=1,table-figures-integer=2]",
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=3]",
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]",
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]", 
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]", 
#                 "S[table-column-width=1.2cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]")
# pxpaSNPDig <- c(#0, 
#   0, 0, 0, 1, 0, 0, 3, 3, 3)
# pxpaSNPCap <- "Population statistics for 707 shared bi-allelic SNPs for sympatric populations of \\textit{Plutella} species from five locations. 
# Statistics include population means for the number of individuals genotyped per locus (\\textit{n}), 
# site depth, observed heterozygosity ($H_{\\textsc{o}}$), gene diversity ($H_{\\textsc{s}}$), 
# and Nei's inbreeding coefficient, $F_{\\textsc{is}}$."
# 
# xtable(pxpaSNPLatex,
#        caption = pxpaSNPCap,
#        align  = pxpaSNPAlg,
#        digits = pxpaSNPDig,
#        lab = "tab:pxpaPopstatsSNP") %>% 
#   print.xtable(floating = TRUE,
#                include.rownames = FALSE,
#                include.colnames = TRUE,
#                table.placement = "p",
#                caption.placement = "top",
#                NA.string = "",
#                scalebox = 0.8,
#                booktabs = TRUE,
#                sanitize.text.function = function(x){x})
# 
# 
# # #################
# # All Sites
# # #################
# 
# # read in hierfstat derived population summary (nei's Fis)
# pxpaHfAS <- file.path(hfPath, "pxpa.5pops.100ind.allSites.hierfstat.basicstats.pop.txt") %>%
#   read.table(header = TRUE)
# 
# # read in genepop derived gene diversities
# pxpaGpAS <- file.path(popStatsPath, "genepopGeneDiversitiesPops.DIV.xlsx") %>%
#   read_excel(sheet = "pxpa100indAllsites", skip = 5)
# 
# # read in vcftools indv site depth summary
# pxpaVdAS <- file.path(vfPath, "pxpa.5pops.100ind.genotypeGVCFs.allSites.hardFilt.maxMiss0.7.idepth") %>%
#   read.table(header = TRUE)
# 
# # read in vcfstats summary
# pxpaVsAS <- file.path(vfPath, "pxpa.5pops.100ind.genotypeGVCFs.allSites.hardFilt.maxMiss0.7.vcfstats.grepSummary.txt") %>%
#   read.table(header = TRUE)
# 
# # merge results into a dataframe summarising population statistics
# pxpaStatsAS <- pxpaVsAS %>%
#   merge(pxpaVdAS %>% setNames(c("ind", "sites", "idepth")), by = "ind") %>%
#   mutate(ind = renamePxToPaus(ind),
#          pop = extractPop(ind),
#          species = extractSpecies(pop)) %>%
#   group_by(pop) %>%
#   summarise(species = unique(species),
#             nInd   = length(ind),
#             nSites = round(mean(sites), 1),
#             pSites = round(nSites / 305136, 2),
#             idepth = round(mean(idepth), 1),
#             homAA  = round(mean(homAA), 1),
#             hetRA  = round(mean(hetRA), 1),
#             homRR  = round(mean(homRR), 1),
#             ref_count =   round(mean(ref_count), 1),
#             snp_count =   round(mean(snp_count), 1),
#             indel_count = round(mean(indel_count), 1),
#             private = round(mean(private), 1)) %>% 
#   merge(pxpaHfAS %>% mutate(pop = extractPop(pop)), by = "pop") %>%
#   merge(pxpaGpAS %>% mutate(pop = extractPop(Sample)) %>% 
#           rename(gpFis = Fis), by = "pop") %>%
#   merge(popLabs)
# 
# # format for latex
# pxpaASLatex <- pxpaStatsAS %>%
#   dplyr::select(Population = popS,
#                 `\\multicolumn{1}{l}{Species}` = species,
#                 `{\\textit{n}}` = indPerLoc, 
#                 `{Sites}` = nSites, # wrap S column headers in curly braces, so it doesn't try to treat as a number
#                 `{\\makecell{Site\\\\depth}}` = idepth,
#                 `{SNPs}`= snp_count,
#                 `{Indels}` = indel_count,
#                 `{\\makecell{Private\\\\sites}}` = private,
#                 `$H_{\\textsc{o}}$` = Ho,
#                 `$H_{\\textsc{s}}$` = Hs,
#                 `$F_{\\textsc{is}}$` = Fis)
# 
# pxpaASLatex$Population[duplicated(pxpaASLatex$Population)] <- NA
# pxpaASAlg <- c("l", "l", ">{\\itshape}l",
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=3]", # center align S columns without decimals by pasting this option onto relevant column aligns, 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=6]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=4]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=4]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=0,table-figures-integer=3]",
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]", 
#                "S[table-column-width=1.1cm,table-number-alignment=center,table-figures-decimal=3,table-figures-integer=1]")
# pxpaASDig <- c(0, 0, 0, 1, 0, 0, 0, 0, 0, 3, 3, 3)
# # pxpaASCap <- "Population statistics for sympatric populations of \\textit{Plutella} species from five locations in Australia. 
# # Statistics include population means for the number of individuals genotyped per locus (\\textit{n}), 
# # site depth, observed heterozygosity ($H_{\\textsc{o}}$), gene diversity ($H_{\\textsc{s}}$), the number of 
# # SNP and indel and private sites, and Nei's inbreeding coefficient, $F_{\\textsc{is}}$."
# 
# pxpaASCap <- "Population statistics for variant and invariant sites for sympatric populations of 
# \\textit{P. australiana} (\\textit{P. aus}) and \\textit{P. xylostella} (\\textit{P. x}) from five locations. 
# Statistics presented include population means for the number of individuals genotyped per locus (\\textit{n}), 
# observed heterozygosity ($H_{\\textsc{o}}$), gene diversity ($H_{\\textsc{s}}$) and Nei's inbreeding coefficient, $F_{\\textsc{is}}$."
# 
# xtable(pxpaASLatex,
#        caption = pxpaASCap,
#        align  = pxpaASAlg,
#        digits = pxpaASDig,
#        lab = "tab:pxpaPopstatsAS") %>% 
#   print.xtable(floating = TRUE,
#                include.rownames = FALSE,
#                include.colnames = TRUE,
#                table.placement = "p",
#                caption.placement = "top",
#                NA.string = "",
#                scalebox = 0.8,
#                booktabs = TRUE,
#                sanitize.text.function = function(x){x})  





# OLD CODE TO CHECK OUTLIER SAMPLES WITH HIGH HET
# ################################################################################################
# Create a master list of the all individuals and their statistics
# ... there are some outliers. This will be used to identify them.
# .. will use the output of a python script to look at teh proportion of heterozygous individuals
# ################################################################################################

# # read in the output of the python script for biallelic site calls
# pyFile31 <- "Px.31pops.440ind.2014.genotypeGVCFs.allSites.hardFilt.maxMiss0.7.rmIndels.recode.stats.txt"
# pyFile28 <- "Px.28pops.402ind.2015.genotypeGVCFs.allSites.hardFilt.maxMiss0.7.rmIndels.recode.stats.txt"
# 
# pyStats31 <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "populationStatistics", pyFile31) %>%
#   read.table(skip = 461, stringsAsFactors = FALSE) %>%
#   setNames(c("ind", "calls", "percent", "hom", "het")) %>%
#   mutate(ind = renamePxToPaus(ind))
# 
# pyStats28 <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "populationStatistics", pyFile28) %>%
#   read.table(skip = 423, stringsAsFactors = FALSE) %>%
#   setNames(c("ind", "calls", "percent", "hom", "het")) %>%
#   mutate(ind = renamePxToPaus(ind))
# 
# pyStats59 <- bind_rows(pyStats31, pyStats28)
# 
# # orderStates <- data_frame( # order the states geographically with increasing longitude.
# #   state = c("WA", "SA", "VIC", "TAS", "NSW", "QLD"),
# #   stateOrder = gsub("WA", 1, state) %>%
# #     gsub("SA", 2, .) %>%
# #     gsub("VIC", 3, .) %>%
# #     gsub("TAS", 4, .) %>%
# #     gsub("NSW", 5, .) %>%
# #     gsub("QLD", 6, .)
# # )
# 
# # merge with the individual statistics from vcfstats 
# indStats842 <- px31VsAS %>%
#   bind_rows(px28VsAS) %>% 
#   merge(allPyStats, by = "ind") %>%
#   mutate(propHet = het / calls, # from the python script
#          pop  = extractPop(ind),
#          host = extractHost(ind),
#          year = extractYear(ind),
#          state = extractAusState(ind)) %>%
#   merge(popLabs, by = "pop") %>%
#   # merge with the geographic plotting order
#   merge(popsOrderPx %>%
#           rename(pop = popString), by = "pop") %>%
#   # now set plotting order for reordering variables
#   arrange(year, orderYearGeo) %>%
#   mutate(orderYearGeo = 1:nrow(.)) %>%
#   arrange(host, orderYearGeo) %>%
#   mutate(orderHostGeo = 1:nrow(.)) %>%
#   arrange(host, propHet) %>%
#   mutate(orderHostHet = 1: nrow(.))
# 
# # extract `outlier` sample
# # ... 3 esperance samples look different on structure plots
# esp3 <- c("espW14scx-08u.11",  "espW14scx-10u.18", "espW14scx-13u.46")
# # ... the entire southend population (n=16) looks different on structure plots
# # ... 9 samples have much higher het than all other samples.
# indStats842 %>%
#   filter(propHet > 0.01) %>%
#   dplyr::select(ind)
# write.table("outliersHighHet.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
# 
# indStats842 %>%
#   filter(propHet > 0.01 | ind %in% esp3 | popS == "Southend SA") %>%
#   write.table("outliersHighHet3Epsp16Southend.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
# 
# # now, look at the `COI_summary` sheet to see which samples have been sequenced at COI
# # ... first define function to update the sample names to the RAD2 syntax
# 
# renameSamplesRAD2 <- function(filename) {
#   basename <- str_extract(filename, "^.*-\\d{1,2}(\\.99)?[mfu]{1}")
#   suffix <- mapply(str_split, string = filename, pattern = basename) %>%
#     lapply(magrittr::extract, 2) %>% 
#     unlist #%>%
#   #gsub("(\\.fastq$)$", ".fastq.gz", .) # if needed append the .gz to names
#   popName <- str_sub(basename, 1, 9)
#   smNum <- str_extract(basename, "-\\d{1,2}") %>%
#     gsub("[-]", "", .) %>% 
#     str_pad(width = 2, side = "left", pad = "0")
#   gender <- str_extract(basename, "[mfu]{1}$")
#   lib <- str_extract(basename, "\\.\\d{1,2}-") %>%
#     gsub("[.-]", "", .) %>% 
#     str_pad(width = 2, side = "left", pad = 0)
#   dup <- gsub("\\.99[mfu]$", "d", basename) %>% str_extract(., "d$")
#   renamed <- paste0(popName, "-", smNum, gender, ".", lib, 
#                     dup, suffix) %>% gsub("NA", "", .)
#   renamed
# }
# 
# # ... first read in the position/well coords for each sample in RAD2 genotypes master
# sampleCoords <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "genotypesMasterRAD2.xlsx") %>%
#   read_excel(sheet = "RAD2-samples-master") %>%
#   mutate(position = paste(Plate, Well, sep = "."),
#          sampleNames = renameSamplesRAD2(`sampleNames(RADpools)`), # update from original RADpools syntax to the correc format
#          sampleNames = renamePxToPaus(sampleNames)) %>%            # rename the incorrect samples 
#   dplyr::select(position, sampleNames)
# 
# # now match the sample names  COI with coords
# COI <- file.path("C:", "UserData", "Kym", "PhD", "RAD2", "haplotypes", "COI_Summary.xlsx") %>%
#   read_excel() %>%
#   mutate(position = paste(plate, position, sep = ".")) %>%
#   filter(RAD2study == 1) %>% # leaves 51 samples
#   merge(sampleCoords, by = "position", all = FALSE) 
# 
# 
# %>% # drops 1 sample - Calca 7 - so, no matter
#   dplyr::select(-species) %>%
#   # you need to run function `renameSamples.R` function on sample names to update to RAD2 ..
#   filter(!is.na(RAD2study),
#          # remove samples with uncertain species
#          `Likely species` != "FAIL",
#          # remove the five poor qual samples Simon excluded
#          `COI genotype name` %in% gsub("_\\(reversed\\)", "", labels(co1dna))) %>%
#   mutate(species_x_wolb = paste(`Likely species`, 
#                                 `Likely wolbachia`, sep = "_"),
#          ind = `COI genotype name`,
#          location = paste(location, state))

# nonoutliers <- indStats842 %>%
#   filter(propHet < 0.01)
# nonoutliers <- nonoutliers[1:10, ]
# bind_rows(outliers, nonoutliers) %>% View()

# ####################################################################
# compare the population statistics generated by hierfstat and genepop
# pxpaStatsSNP %>%
#   reshape2::melt(id.vars = "pop", # use all non-measured variables
#                  measure.vars  = c("Fis", "gpFis"),
#                  variable.name = "program",
#                  value.name = "Fis") %>%
#   mutate(program = gsub("gpFis", "genepop", program) %>%
#            gsub("Fis", "hierfstat", .)) %>%
#   merge(pxpaStatsSNP %>%
#           reshape2::melt(id.vars = "pop", # use all non-measured variables
#                          measure.vars  = c("Ho", "1-Qintra"),
#                          variable.name = "program",
#                          value.name = "Ho") %>%
#           mutate(program = gsub("1-Qintra", "genepop", program) %>%
#                    gsub("Ho", "hierfstat", .)), by = c("pop", "program")) %>%
#   merge(pxpaStatsSNP %>%
#           reshape2::melt(id.vars = "pop", # use all non-measured variables
#                          measure.vars  = c("Hs", "1-Qinter"),
#                          variable.name = "program",
#                          value.name = "Hs") %>%
#           mutate(program = gsub("1-Qinter", "genepop", program) %>%
#                    gsub("Hs", "hierfstat", .)), by = c("pop", "program")) %>%
#   reshape2::melt(variable.name = "parameter",
#                  value.name = "parValue") %>%
#   ggplot2::ggplot(aes(x = pop, y = parValue, fill = program)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   facet_wrap(~parameter) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1))

# # compare observed and expected H
# pxpaStatsSNP %>%
#   ggplot(aes(x = Ho, y = Hs, shape = species)) +
#   geom_point(size = 2) +
#   geom_line(aes(x = Ho, y = Ho, lty = species))

# End compare programs(pxpa, 5pops, 100ind)
# #########################################################################

# ##############################################################################
# compare the population statistics generated by hierfstat and genepop, allSites
# pxpaStatsAS %>%
#   reshape2::melt(id.vars = "pop", # use all non-measured variables
#                  measure.vars  = c("Fis", "gpFis"),
#                  variable.name = "program",
#                  value.name = "Fis") %>%
#   mutate(program = gsub("gpFis", "genepop", program) %>%
#            gsub("Fis", "hierfstat", .)) %>%
#   merge(pxpaStatsAS %>%
#           reshape2::melt(id.vars = "pop", # use all non-measured variables
#                          measure.vars  = c("Ho", "1-Qintra"),
#                          variable.name = "program",
#                          value.name = "Ho") %>%
#           mutate(program = gsub("1-Qintra", "genepop", program) %>%
#                    gsub("Ho", "hierfstat", .)), by = c("pop", "program")) %>%
#   merge(pxpaStatsAS %>%
#           reshape2::melt(id.vars = "pop", # use all non-measured variables
#                          measure.vars  = c("Hs", "1-Qinter"),
#                          variable.name = "program",
#                          value.name = "Hs") %>%
#           mutate(program = gsub("1-Qinter", "genepop", program) %>%
#                    gsub("Hs", "hierfstat", .)), by = c("pop", "program")) %>%
#   reshape2::melt(variable.name = "parameter",
#                  value.name = "parValue") %>%
#   ggplot2::ggplot(aes(x = pop, y = parValue, fill = program)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   facet_wrap(~parameter) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1))
# 
# # # compare observed and expected H
# pxpaStatsAS %>%
#   ggplot(aes(x = Ho, y = Hs, shape = species)) +
#   geom_point(size = 2) +
#   geom_line(aes(x = Ho, y = Ho, lty = species))

# End compare programs, pxpa AS
# ###########################################################################

