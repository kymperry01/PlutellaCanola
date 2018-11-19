# #########################################################################
# Summarise the experimental cross data for P. australiana x P. xylostella 
# 30/5/2017
# #########################################################################

library(dplyr)
library(magrittr)
library(readxl)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(xtable)

# #################
# define functions
# #################

# function to wrap values in parentheses in latex
wrapParenth <- function(x){
  paste0("{(}", x, "{)}")
} 

# function for standard error of the mean (omitting NAs)
se <- function(x) {sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))}

# ########################
# Summarise F0 cross data
# ########################

# read in raw data for the F0 crosses
wd <- "C:/UserData/Kym/PhD/Data/P.australiana/"
rawF0 <- read_excel(paste0(wd, "matedCrossesv2.5Formatted v2.xlsx"),
                    sheet = "Mated Cross Data F0", skip = 4)

# some columns contain different data types in different rows.
# identify rows with numeric count data (as opposed to presence-only)
presVal <- 999 # value coding `presence` data
countEggsRowsIdx <- intersect(which(rawF0$eggsV != presVal), # index of rows with egg count data in both columns (eggV, eggsD)
                              which(rawF0$eggsD != presVal))

# calculate total eggs in rows with count data in both columns
totalEggs <- rawF0$eggsV[countEggsRowsIdx] + rawF0$eggsD[countEggsRowsIdx]
# calculate total proportion of dessicated eggs (among total eggs)
prDudEggs   <- rawF0$eggsD[countEggsRowsIdx] / totalEggs


# summarise data by cross type (PaPa, PxPx, PaPx, PxPa)
# ... capital `N` = number of crosses that produced eggs/larvae/adults
# ... capital `A` = number of crosses assessed for that variable
# ... capital `P` = proportion of crosses that produced ... var
summaryF0 <- rawF0 %>%
  mutate(eggsTot = NA,
         eggsTot = replace(eggsTot, 
                           list = countEggsRowsIdx, values = totalEggs),
         pEggsD  = NA,
         pEggsD  = replace(pEggsD, 
                           list = countEggsRowsIdx, values = prDudEggs)) %>% # the proportion of desiccated eggs in each cross
  group_by(cross) %>%
  # some of these variables are no longer needed ...
  summarise(N = length(cross), # num crosses set up
            AeggsV = sum(!is.na(eggsV)), # n crosses assessed for eggsV
            NeggsV = sum(eggsV[!is.na(eggsV)] > 0), # n crosses producing eggsV
            PeggsV = round(NeggsV / AeggsV, 3), # prop. crosses producing eggsV
            AeggsD = sum(!is.na(eggsD)),
            NeggsD = sum(eggsD[!is.na(eggsD)] > 0), # num crosses producing eggsD
            PeggsD = round(NeggsD / AeggsD, 3),
            NeggsDCount = sum(!eggsD %in% c(NA, presVal)), # num of crosses with count data for eggsD (for mean calc.)
            mnNrEggsD = mean(eggsD[!eggsD %in% c(NA, presVal)]),
            sdNrEggsD = sd(eggsD[!eggsD %in% c(NA, presVal)]),
            mnNrEggsV = mean(eggsV[!eggsV %in% c(NA, presVal)]),
            sdNrEggsV = sd(eggsV[!eggsV %in% c(NA, presVal)]),
            AeggsTot = sum(!is.na(eggsV) | !is.na(eggsD)),
            NeggsTot = sum(eggsV[!is.na(eggsV) > 0] | eggsD[!is.na(eggsD) > 0]), # numb crosses that produced eggs (V or D)
            PeggsTot = round(NeggsTot / AeggsTot, 3),
            NeggsTotCount = sum(!eggsTot %in% c(NA, presVal)), # number of crosses with egg count data (for mean calc.)
            mnNrEggsTot = round(mean(eggsTot[!is.na(eggsTot)]), 3), # mean eggs produced for those crosses with egg count data (not NA)
            sdNrEggsTot = sd(eggsTot[!is.na(eggsTot)]),
            seNrEggsTot = se(eggsTot),
            Alarv  = sum(!is.na(larvae)),
            Nlarv  = sum(larvae[!is.na(larvae)] > 0),
            Plarv  = round(Nlarv / Alarv, 3),
            # there are no count data for larvae
            Aadul  = sum(!is.na(adults)),
            Nadul  = sum(adults[!is.na(adults)] > 0),
            Padul  = round(Nadul / Aadul, 3),
            NadulCount = sum(!adults %in% c(NA, presVal)), # number of crosses with adult count data (for mean calc.)
            mnNrAdul = round(mean(adults[!adults %in% c(NA, presVal)]), 3), # mean num adults produced for crosses with count data
            sdNrAdul = round(sd(adults[!adults %in% c(NA, presVal)]), 3),
            seNrAdul = se(adults[!adults %in% c(NA, presVal)]),
            totAdultsF = sum(adultsF[!adultsF %in% c(NA, presVal)]),
            totAdultsM = sum(adultsM[!adultsM %in% c(NA, presVal)])) 

# ##################################################################################
# Summarize the number of rows with count data that were used to calculate the means
# ##################################################################################

nrow(summaryF0) # 266
summaryF0 %>% filter(adults != 999)   %>% nrow()
# (225 / 266) * 100 = 84.6%
summaryF0 %>% filter(!is.na(eggsTot)) %>% nrow()
# (213 / 266) * 100 = 80.1%

#write.csv(summaryF0, "summaryF0.csv",
 #         quote = FALSE, row.names = FALSE)

# export latex code for the cross data
f0Latex <- summaryF0 %>%
  mutate(sortOrder = c(1, 3, 4, 2),
         # to align the crosses by the 'x', make separate columns for M,F
         crossF = gsub("PafPam", "\\\\textit{P.aus}\\\\female", cross) %>%
           gsub("PafPxm", "\\\\textit{P.aus}\\\\female", .) %>%
           gsub("PxfPam", "\\\\textit{P.x}\\\\female", .) %>%
           gsub("PxfPxm", "\\\\textit{P.x}\\\\female", .),
         crossM = gsub("PafPam", "x \\\\textit{P.aus}\\\\male", cross) %>%
           gsub("PafPxm", "x \\\\textit{P.x}\\\\male", .) %>%
           gsub("PxfPam", "x \\\\textit{P.aus}\\\\male", .) %>%
           gsub("PxfPxm", "x \\\\textit{P.x}\\\\male", .),
         PeggsTot = wrapParenth(PeggsTot),
         Padul    = wrapParenth(Padul),
         mnNrEggsTot = round(mnNrEggsTot, 2),
         sdNrEggsTot = round(sdNrEggsTot, 2),
         seNrEggsTot = round(seNrEggsTot, 2),
         mnNrAdul   = round(mnNrAdul, 2),
         sdNrAdul   = round(sdNrAdul, 2),
         seNrAdul   = round(seNrAdul, 2),
         mnNrEggsSD = paste(mnNrEggsTot, sdNrEggsTot, sep = " \\pm "),
         mnNrEggsSE = paste(mnNrEggsTot, seNrEggsTot, sep = " \\pm "),
         mnNrAdulSD = paste(mnNrAdul, sdNrAdul, sep = " \\pm "),
         mnNrAdulSE = paste(mnNrAdul, seNrAdul, sep = " \\pm ")) %>%
  arrange(sortOrder) %>%
  dplyr::select(crossF, crossM, N, 
                NeggsTot, PeggsTot,
                Nadul, Padul,
                mnNrEggsSE, mnNrAdulSE)
    
# set up multicolumn headers
# set up add.to.rows for the table sub-headers
f0Nms <- list()
f0Nms$pos <- list(0)
f0Nms$command <- c(
  # paste(c("& & & ", "\\multicolumn{3}{c}{Crosses producing egg progeny} & ", "\\multicolumn{3}{c}{Crosses producing adult progeny} \\\\\n"), collapse = "") %>%
  #   paste0("\\cmidrule(lr){4-6}\n", "\\cmidrule(lr){7-9}\n"),
  paste(c("\\multicolumn{2}{l}{Cross (\\female{ } x \\male)} & ", 
          "{\\makecell{No.\\\\replicates}} & ", 
          "\\multicolumn{2}{c}{\\makecell{No. reps\\\\eggs}} & ", 
          "\\multicolumn{2}{c}{\\makecell{No. reps\\\\adults}} & ", 
          "{\\makecell{Mean $\\pm$ SEM \\\\ no. eggs}} & ",
          "{\\makecell{Mean $\\pm$ SEM \\\\ no. adults}}\\\\\n"), collapse = ""))

  
f0Algn <- c("l", 
            "l", "@{\\space}l", 
            "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]",
            "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]",
            "@{}S[table-align-text-pre=false,table-figures-integer=1,round-mode=places,round-precision=2]",
            "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]",
            "@{}S[table-align-text-pre=false,table-figures-integer=1,round-mode=places,round-precision=2]",
            "S[separate-uncertainty,table-figures-uncertainty=1,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=2]",
            "S[separate-uncertainty,table-figures-uncertainty=1,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=2]")

f0Dig <- c(0, 
           0, 0, 0, 0, 2, 2, 0, 2, 2)
# f0Cap <- "Progeny production of reciprocal single pair crosses of \\textit{P. australiana} and \\textit{P. xylostella}. 
# For each directional cross, the number (No.) and proportion (p) of crosses producing eggs and adults, and the mean $\\pm$ standard error of the mean (SEM)
# numbers of progreny, are presented."

f0Cap <- "Fecundity of intra-species and reciprocal inter-species single pair crosses of \\textit{P. australiana} (\\textit{P.aus}) and
\\textit{P. xylostella} (\\textit{P.x}). Presented are the number and proportion in parentheses of replicates that produced eggs and adult offspring, 
and the mean $\\pm$ standard error of the mean number of eggs and adult offspring per replicate."

xtable(f0Latex,
       caption = f0Cap,
       align   = f0Algn,
       digits  = f0Dig,
       lab = "tab:f0crosses") %>% 
  print.xtable(floating = TRUE,
               include.rownames = FALSE,
               include.colnames = FALSE,
               add.to.row = f0Nms,
               table.placement = "p",
               caption.placement = "top",
               NA.string = "",
               scalebox = 0.9,
               booktabs = TRUE,
               sanitize.text.function = function(x){x})


# ########################
# Summarise F1 cross data
# ########################

rawF1 <- read_excel(paste0(wd, "matedCrossesv2.5Formatted.xlsx"),
                    sheet = "Mated Cross Data F1", skip = 3)

summaryF1 <- rawF1 %>%
  # for the F1 data, just use the `eggsTot` column as is (without recalculating), until checked by Kevin
  group_by(parentCrossF0, crossF1) %>%
  summarise(crossF1Type = unique(crossF1Type),
            N = length(crossF1), # num crosses set up
            Aeggs = sum(!is.na(eggsV) | !is.na(eggsD)),
            Neggs = sum(eggsV[!is.na(eggsV) > 0] | eggsD[!is.na(eggsD) > 0]), # numb crosses that produced eggs (V or D)
            Peggs = round(Neggs / Aeggs, 3),
            NeggsCount = sum(!eggsTot %in% c(NA, presVal)), # number of crosses with egg count data (for mean calc.)
            mnNrEggs = round(mean(eggsTot[!eggsTot %in% c(NA, presVal)]), 3), # mean nr eggs produced for those crosses with egg count data (not NA or presVal)
            sdNrEggs = round(sd(eggsTot[!eggsTot %in% c(NA, presVal)]), 3),
            seNrEggs = se(eggsTot[!eggsTot %in% c(NA, presVal)]),
            # just stick to presence/absence data for larvae and adults (very little count data available)
            Alarv  = sum(!is.na(larvae)),
            Nlarv  = sum(larvae[!is.na(larvae)] > 0),
            Plarv  = round(Nlarv / Alarv, 3),
            Aadul  = sum(!is.na(adults)),
            Nadul  = sum(adults[!is.na(adults)] > 0),
            Padul  = round(Nadul / Aadul, 3),
            mnNrAdul = round(mean(adults[!adults %in% c(NA, presVal)]), 5),
            sdNrAdul = round(sd(adults[!adults %in% c(NA, presVal)]), 5),
            seNrAdul = se(adults[!adults %in% c(NA, presVal)])) # mean num adults produced for crosses with count data)

#write.csv(summaryF1, "summaryF1.csv",
 #          quote = FALSE, row.names = FALSE)

# The command strings for the two hybrid F0 cross directions 
hybridF0female_Paus <- "(\\\\textit{P.aus}\\\\female{ }x \\\\textit{P.x}\\\\male)"
hybridF0female_Px <- "(\\\\textit{P.x}\\\\female{ }x \\\\textit{P.aus}\\\\male)"

# export latex code for the F1 cross table
f1Latex <- summaryF1 %>%
  ungroup() %>%
  mutate(sortOrder = c(1, 2, 4, 3, 5, 6, 7),
         # create a variable that encodes the female source of hybrids (i.e. each of the two reciprocal F0 crosses)
         crossHybrid = paste(parentCrossF0, crossF1, sep = "_"),
         crossF = gsub(
           # the paus F0 female direction
           "PafPxm_hf1hm1", paste0(hybridF0female_Paus, "\\\\female"), crossHybrid) %>%
           gsub("PafPxm_hf1Pam0", paste0(hybridF0female_Paus, "\\\\female"), .) %>%
           gsub("PafPxm_hf1Pxm0", paste0(hybridF0female_Paus, "\\\\female"), .) %>%
           gsub("PafPxm_Paf0hm1", "\\\\textit{P.aus}\\\\female", .) %>%
           gsub("PafPxm_Pxf0hmf1", "\\\\textit{P.x}\\\\female", .) %>%
           # the px F0 female direction
           gsub("PxfPam_hf1hm1", paste0(hybridF0female_Px, "\\\\female"), .) %>%
           gsub("PxfPam_hf1Pam0", paste0(hybridF0female_Px, "\\\\female"), .),
          crossM = gsub(
            # the paus F0 female direction
            "PafPxm_hf1hm1", paste0("x ",  hybridF0female_Paus, "\\\\male"), crossHybrid) %>%
           gsub("PafPxm_hf1Pam0", "x \\\\textit{P.aus}\\\\male", .) %>%
           gsub("PafPxm_hf1Pxm0", "x \\\\textit{P.x}\\\\male", .) %>%
           gsub("PafPxm_Paf0hm1", paste0("x ",  hybridF0female_Paus, "\\\\male"), .) %>%
           gsub("PafPxm_Pxf0hmf1", paste0("x ",  hybridF0female_Paus, "\\\\male"), .) %>%
           # the px F0 female direction
           gsub("PxfPam_hf1hm1", paste0("x ", hybridF0female_Px, "\\\\male"), .) %>%
           gsub("PxfPam_hf1Pam0", "x \\\\textit{P.aus}\\\\male", .),
         Peggs = wrapParenth(format(Peggs, nsmall = 2)),
         Padul    = wrapParenth(format(Padul, nsmall = 2)),
         mnNrEggs = format(round(mnNrEggs, 2), nsmall = 2),
         sdNrEggs = format(round(sdNrEggs, 2), nsmall = 2),
         seNrEggs = format(round(seNrEggs, 2), nsmall = 2),
         mnNrAdul   = format(round(mnNrAdul, 2), nsmall = 2),
         sdNrAdul   = format(round(sdNrAdul, 2), nsmall = 2),
         seNrAdul   = format(round(seNrAdul, 2), nsmall = 2),
         mnNrEggsSD = paste(mnNrEggs, sdNrEggs, sep = " \\pm "),
         mnNrEggsSE = paste(mnNrEggs, seNrEggs, sep = " \\pm "),
         mnNrAdulSD = paste(mnNrAdul, sdNrAdul, sep = " \\pm "),
         mnNrAdulSE = paste(mnNrAdul, sdNrAdul, sep = " \\pm ")) %>%
  arrange(sortOrder) %>%
  dplyr::select(
    crossF, 
    crossM,
    N, Neggs, Peggs, Nadul, Padul, mnNrEggsSE, mnNrAdulSE)

# manually replace a few values
f1Latex$mnNrEggsSE[7] <- sprintf("%.2f", 0) 
f1Latex$mnNrAdulSE[7] <- sprintf("%.2f", 0) 
f1Latex$mnNrAdulSE[1] <- "\\textendash"

# set up multicolumn headers
# set up add.to.rows for the table sub-headers
f1Nms <- list()
f1Nms$pos <- list(0, 0, 5)
f1Nms$command <- c(
  paste(c("\\multicolumn{2}{l}{\\makecell{Cross \\\\ (\\female{ } x \\male)}} & ",
          "{\\makecell{No. \\\\ replicates}} & ",
          "\\multicolumn{2}{c}{\\makecell{No. reps\\\\ eggs}} & ",
          "\\multicolumn{2}{c}{\\makecell{No. reps \\\\ adults}} & ",
          "{\\makecell{Mean $\\pm$ SEM \\\\ no. eggs}} & ",
          "{\\makecell{Mean $\\pm$ SEM \\\\ no. adults}} \\\\\n"), collapse = ""),
  "\\multicolumn{9}{l}{F0 \\textit{P.aus}\\female{} source} \\\\\n",
  "\\multicolumn{9}{l}{F0 \\textit{P.x}\\female{} source} \\\\\n")


f1Algn <- c("l",
            ">{\\quad}l", 
            "@{\\space}l",# to indent the first column (under add.to.rows subheaders)
            "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]",
            "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]",
            "@{}S[table-align-text-pre=false,table-figures-integer=1,round-mode=places,round-precision=2]",
            "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]",
            "@{}S[table-align-text-pre=false,table-figures-integer=1,round-mode=places,round-precision=2]",
            "S[separate-uncertainty,table-figures-uncertainty=1,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=2]",
            "S[separate-uncertainty,table-figures-uncertainty=1,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=2]")

f1Dig <- c(0, 
           0, 0, 0, 0, 2, 2, 0, 2, 2)
f1Cap <- "Fecundity of hybrid F1 crosses and back-crosses. 
Presented are the numbers and proportion in parentheses of replicates producing eggs and adult offspring, 
and the mean $\\pm$ standard error of the mean numbers of eggs and adults offspring per replicate.
A dash denotes an absence of count data."

xtable(f1Latex,
       caption = f1Cap,
       align   = f1Algn,
       digits  = f1Dig,
       lab = "tab:f1crosses") %>% 
  print.xtable(floating = TRUE,
               include.rownames = FALSE,
               include.colnames = FALSE,
               add.to.row = f1Nms,
               table.placement = "p",
               caption.placement = "top",
               NA.string = "",
               scalebox = 0.9,
               booktabs = TRUE,
               sanitize.text.function = function(x){x})


# #########################################################################
# 16/9/2017: 
# Make a bar plot of cross success and cross fecundity for AES presentation
# #########################################################################


f0ToPlot <- summaryF0 %>% 
  mutate(speciesCross = gsub("PafPam", "Intra-species crosses", cross) %>%
           gsub("PxfPxm", "Intra-species crosses", .) %>%
           gsub("PxfPam", "Inter-species crosses", .) %>%
           gsub("PafPxm", "Inter-species crosses", .),
         cross = gsub("PafPam", "P.aus X P.aus", cross) %>%
           gsub("PxfPxm", "P.x X P.x", .) %>%
           gsub("PafPxm", "P.aus Female X P.x Male", .) %>%
           gsub("PxfPam", "P.x Female X P.aus Male", .)) %>%
  rename(`Cross type` = speciesCross) %>%
  # reorder variable
  arrange(desc(`Cross type`), desc(cross)) %>%
  mutate(order = 1:nrow(.))
library(magrittr)
f0ToPlot$cross %<>% reorder(f0ToPlot$order)

# F0 plot
crossPal <- brewer.pal(9, "Paired")[c(1, 2)]
pdf("f0CrossSuccessPlot.pdf", height = 6, width = 4)
ggplot(f0ToPlot, 
       aes(x = cross, y = Padul, fill = `Cross type`)) +
  geom_bar(stat = 'identity') +
  ggtitle("F0 crosses producing adult offspring") +
  labs(y = "Proportion of replicates producing adult offspring") +
  theme_classic() + 
  theme(axis.title.y = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(angle = 60, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = crossPal,
                    guide = guide_legend(reverse=TRUE))
dev.off()

# F1 hybrid crosses/backcrosses from P. australiana female parent 
f1PausFemLine <- summaryF1 %>%
  filter(parentCrossF0 == "PafPxm") %>%
  mutate(speciesF1CrossType = gsub("hf1hm1", "Hybrid crosses", crossF1) %>%
           gsub("hf1Pam0", "Back-crosses", .) %>%
           gsub("hf1Pxm0", "Back-crosses", .) %>%
           gsub("Paf0hm1", "Back-crosses", .) %>%
           gsub("Pxf0hmf1", "Back-crosses", .),
         crossF1 = gsub("hf1hm1", "Hybrid X Hybrid", crossF1) %>%
           gsub("hf1Pam0", "Hybrid Female X P.aus Male", .) %>%
           gsub("hf1Pxm0", "Hybrid Female X P.x Male", .) %>%
           gsub("Paf0hm1", "P.aus Female X Hybrid Male", .) %>%
           gsub("Pxf0hmf1", "P.x Female X Hybrid Male", .)) %>%
  rename(`F1 cross type` = speciesF1CrossType) %>%
  # reorder variable
  arrange(desc(`F1 cross type`)) %>%
  mutate(order = 1:nrow(.))
f1PausFemLine$crossF1 %<>% reorder(f1PausFemLine$order)

# F1 plot - P.aus female line
crossPal2 <- brewer.pal(9, "Paired")[c(5, 6)]
pdf("f1CrossSuccessPlotPausMatline.pdf", height = 6, width = 4)
ggplot(f1PausFemLine, 
       aes(x = crossF1, y = Padul, fill = `F1 cross type`)) +
  geom_bar(stat = 'identity') +
  ggtitle("F1 crosses producing adult offspring") +
  labs(y = "Proportion of replicates producing adult offspring") +
  theme_classic() + 
  theme(axis.title.y = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(angle = 60, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0),
        legend.title = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = crossPal2,
                    guide  = guide_legend(reverse=TRUE))
dev.off()

# F1 hybrid crosses/backcrosses from P. xylostella female parent 
f1PxylFemLine <- summaryF1 %>%
  filter(parentCrossF0 == "PxfPam") %>%
  mutate(speciesF1CrossType = gsub("hf1hm1", "Hybrid crosses", crossF1) %>%
           gsub("hf1Pam0", "Back-crosses", .) %>%
           gsub("hf1Pxm0", "Back-crosses", .) %>%
           gsub("Paf0hm1", "Back-crosses", .) %>%
           gsub("Pxf0hmf1", "Back-crosses", .),
         crossF1 = gsub("hf1hm1", "Hybrid X Hybrid", crossF1) %>%
           gsub("hf1Pam0", "Hybrid Female X P.aus Male", .) %>%
           gsub("hf1Pxm0", "Hybrid Female X P.x Male", .) %>%
           gsub("Paf0hm1", "P.aus Female X Hybrid Male", .) %>%
           gsub("Pxf0hmf1", "P.x Female X Hybrid Male", .)) %>%
  rename(`F1 cross type` = speciesF1CrossType) %>%
  # reorder variable
  arrange(desc(`F1 cross type`)) %>%
  mutate(order = 1:nrow(.))
f1PxylFemLine$crossF1 %<>% reorder(f1PxylFemLine$order)

# F1 plot - P.xyl female line
pdf("f1CrossSuccessPlotPxylMatline.pdf", height = 6, width = 3)
ggplot(f1PxylFemLine, 
       aes(x = crossF1, y = Padul, fill = `F1 cross type`)) +
  geom_bar(stat = 'identity') +
  ggtitle("F1 crosses producing adult offspring") +
  labs(y = "Proportion of replicates producing adult offspring") +
  theme_classic() + 
  theme(axis.title.y = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(angle = 60, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0),
        legend.title = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = crossPal2,
                    guide  = guide_legend(reverse=TRUE))
dev.off()

