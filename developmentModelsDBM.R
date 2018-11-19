# ############################################################################
# Fit non-linear model to DBM development data
# Kym D Perry & Michael Keller 21/8/2018

# This script does:
# (1) Fits non-linear regression models to published temperature-dependent development data
# for the diamondback moth (Liu 2002) using Briere's model 2 function (Briere 1999).
# (2) Outputs the estimated parameters for Briere's model (a, Tmin, Tmax, m) for each lifestage
# (3) Plots the raw data and fitted curves. 

# Notes: 
# The main code pipeline is designed to run all models (per lifestage) in one go.
# At the bottom of the script is an alternative version for manually running models individually.
#
# ############################################################################

library(tidyverse)
library(nlreg)
library(readxl)
library(xtable)

# #########################################################
# DBM temperature-development data for each DBM lifestage #
# #########################################################

# Read in the DBM development data (Lui 2002)
devData <- file.path(
  "C:", "UserData","Kym", "PhD", "Data", "colonisation", "developmentModelsDBM", "DBM developmental models-2.xlsx"
  ) %>% read_excel(range = "AF5:AP18") %>%
  rowwise() %>%
  mutate(PupaFM = sum(c(Female_Pupa, Male_Pupa) / 2)) %>% # add a combined male-female pupa model (take mean dev rate for each)
  as_data_frame() %>%
  dplyr::select(-Total, Total) # move Total column to the end (help with order of data_frames last)
  
# ########################################
# Model fitting and prediction functions #
# ########################################

# model fitting function (df needs a column named `devRate` and a column named `TempObs`)
fitModel <- function(df, startVals = c(a = 0.05, Tmin = 7, Tmax = 36, m = 2)){
  nlreg::nlreg(devRate ~ a * TempObs * (TempObs - Tmin) * ((Tmax - TempObs) ^ (1 / m)), # Briere's model 2 (Briere 1999)
               data = df, start = startVals)
  }

# model prediction function for plotting fitted curves
predictModel <- function(mod, tempRange = c(0, 38)){
  par <- param(mod)
  data_frame(TempObs = modelr::seq_range(tempRange, n = 10000),
             devRate = par["a"] * TempObs * (TempObs - par["Tmin"]) * ((par["Tmax"] - TempObs) ^ (1 / par["m"])))
  }

# ############
# Run models #
# ############

# run models for each lifestage
fittedModels <- devData %>% # convert to `long` format (one column per variable, one row per observation)
  reshape2::melt(
    id.vars = "Temp_C",
    variable.name = "lifeStage",
    value.name = "devRate"
  ) %>%
  dplyr::rename(TempObs = Temp_C) %>% # rename to TempObs, required for the fitModel function
  filter(!lifeStage %in% c("1st_Instar", "2nd_Instar", # these models do not coverage, so exclude
                           "Female_Pupa", "Male_Pupa")) %>% # exclude individual pupal sexes, as we'll model combined dev data
  group_by(lifeStage) %>%
  nest() %>% # make a `list column` with each row containing the input data for each lifestage (single row)
  mutate(models      = map(data,   fitModel),     # iterate over the data, run a model for each lifestage, store models in a new column
         predictions = map(models, predictModel), # iterate over the models, predict data and results in a new column
         estimates   = map(models, coef)) %>%
  ungroup() %>%
  arrange(lifeStage)# interate over models, extract esimated model paramaters, store in a new column
  
# extract the estimated parameter values for each lifestage
estimatedParams <- fittedModels %>%
  unnest(estimates) %>% 
  mutate(param = rep(c("a", "Tmin", "Tmax", "m"), nrow(.) / 4)) 

# display params in a table
estimatedParams %>%
  reshape2::dcast(lifeStage ~ param, value.var = "estimates")

# To force the fitted development curves through the x-axis, plot a development rate of 0 
# for (Tmin and) Tmax temperatures, then join these data to the predicted data
zeroDev <- estimatedParams %>%
  filter(param %in% c("Tmax")) %>% # we could add Tmin here but model fits pass through x anyway, so not needed.
  mutate(devRate = 0) %>%
  dplyr::select(lifeStage, TempObs = estimates, devRate) # select and rename columns identical to predictedData (for joining)
  
# extract the predicted data for each lifestage
predictedData <- fittedModels %>% 
  unnest(predictions) %>% 
  full_join(zeroDev, by = c("lifeStage", "TempObs", "devRate"))

# plot raw data versus fitted data for each lifestage
stageNames <- c(
  "3rd_Instar" = "3rd instar",
  "4th_Instar" = "4th instar",
  "PupaFM"     = "Pupa"
  )

fittedModels %>%
  unnest(data) %>% 
  mutate(lifeStage = str_replace_all(lifeStage, stageNames)) %>%
  ggplot(aes(x = TempObs, y = devRate)) +
  geom_point() +
  geom_line(data = predictedData %>%
              mutate(lifeStage = str_replace_all(lifeStage, stageNames)),
            aes(x = TempObs, y = devRate),
            linetype = 1) +
  geom_point(data = zeroDev %>%
               mutate(lifeStage = str_replace_all(lifeStage, stageNames)), 
             aes(x = TempObs, y = devRate), colour = brewer.pal(9, "Oranges")[7]) + # red point added at Tmax to force dev through zero
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        aspect.ratio = 0.6) +
  labs(x = expression("Temperature ("*degree*C*")"), 
       y = "Development rate (1/development time) per day") +
  facet_wrap(~lifeStage, scales = "free_y", ncol = 3)
ggsave("temperatureDevelopmentCurvesDBMDraft.pdf", width = 8, height = 5.1) # write a pdf to current wd.

# #######################################
# Out a LaTex table with the parameters #
# #######################################

wrapParenth <- function(x){
  #x[!is.na(x)] <- paste0("{(}", sprintf("%.2f", x[!is.na(x)]), "{)}")
  x[!is.na(x)] <- paste0("{(}", format(x[!is.na(x)], nsmall = 2), "{)}")
  x
}

instars <- c("1st \\& 2nd Instar" = "$1^{\\\\textsc{st}}$ and $2^{\\\\textsc{nd}}$ instar",
             "3rd_Instar" = "$3^{\\\\textsc{rd}}$ instar",
             "4th_Instar" = "$4^{\\\\textsc{th}}$ instar",
             "Prepupa" = "Pre-pupa",
             "PupaFM"  = "Pupa",
             "Total"   = "Total")

# display params in a table
paramsLatex <- estimatedParams %>%
  reshape2::dcast(lifeStage ~ param, value.var = "estimates") %>%
  mutate(lifeStage = str_replace_all(lifeStage, instars)) %>% # add an escape symbol for the and
  dplyr::select(`Life stage` = lifeStage,
                `{$a$}` = a, 
                `{$m$}` = m, 
                `{$T_{\\textsc{max}}$}` = Tmax, 
                `{$T_{\\textsc{min}}$}` = Tmin)
  
# output Latex code
colAlign <- c(
  "l", "l", 
  "S[table-number-alignment=center,table-figures-decimal=5,table-figures-integer=2]\n", 
  "S[table-number-alignment=center,table-figures-decimal=5,table-figures-integer=2]\n", 
  "S[table-number-alignment=center,table-figures-decimal=5,table-figures-integer=2]\n", 
  "S[table-number-alignment=center,table-figures-decimal=5,table-figures-integer=2]\n")


colDigits <- c(0, 5, 5, 5, 5, 5)
paramsCaption <- "Briere's model parameters for \\textit{P. xylostella}"

paramsLatex %>% 
  xtable(align  = colAlign,
         digits = colDigits,
         caption = paramsCaption,
         lab = "tab:params") %>%
  print.xtable(include.rownames = FALSE,
               include.colnames = TRUE,
               table.placement = "p",
               caption.placement = "top",
               NA.string = "",
               booktabs = TRUE,
               sanitize.text.function = function(x){x})



# ====================================
# output a Latex table for publication
# ====================================

wrapParenth <- function(x){
  #x[!is.na(x)] <- paste0("{(}", sprintf("%.2f", x[!is.na(x)]), "{)}")
  #x[!is.na(x)] <- paste0("{(}", format(x[!is.na(x)], nsmall = 2), "{)}")
  paste0("{(}", x, "{)}")
}

roundNumeric <- function(x, digits = 2){
  if(is.numeric(x)){
    x <- round(x, digits)
  }
  x
  #format(x, nsmall = 2) # output a summary value of length 1. Make this character format, to force it to retain two decimal places
}

# a function to make all the minor hacks needed to the final table
cleanupLatex <- function(x){
  x[x == "NaN \\pm NA"] <- NA
  x[x == "7.00 \\pm NA"] <- "7.00"
  x[x == "1.00 \\pm NA"] <- "1.00"
  x[x == "2.00 \\pm NA"] <- "2.00"
  x[x == "NaN"] <- NA
  x
}

# for table display 
rowOrder <- c("buchan weed"  = "5", "sea rocket" = "2",
              "lincoln weed" = "1", "forage brassica" = "4",
              "nil" = "6", "vol\\. canola" = "3")

hostNames <- c("buchan weed"  = "Buchan weed, \\\\textit{Hirschfeldia incana} L., ", 
               "sea rocket" = "Sea rocket, \\\\textit{Cakile maritima} Scop.,",
               "lincoln weed" = "Lincoln weed, \\\\textit{Diplotaxis tenuifolia} L.DC", 
               "forage brassica" = "\\\\textit{Brassica} forage crops",
               "nil" = "Nil host", "vol\\. canola" = "Weedy canola, \\\\textit{Brassica napus}, L.")


dupToBlank <- function(x){
  x[duplicated(x)] <- ""
  x
}

surveySummaryLatex <- surveySummary %>%
  map_dfr(roundNumeric, 2) %>%
  mutate(
    rowOrder = str_replace_all(host, rowOrder), 
    host     = str_replace_all(host, hostNames),
    pSitesPositiveAnyStage = wrapParenth(pSitesPositiveAnyStage),
    nSitesWithLarvae = paste0(nSitesPositiveLarvae, "/", nSitesSampledLarvae),
    pSitesWithLarvae = wrapParenth(pSitesPositiveLarvae),
    nSitesWithMoths = paste0(nSitesPositiveMoths, "/", nSitesSampledMoths),
    pSitesWithMoths = wrapParenth(pSitesPositiveMoths),
    mnMothsSE       = paste(sprintf("%.2f", mnCountMoths), sprintf("%.2f", seCountMoths), sep = " \\pm "), # sprintf hack esure 2 decimals
    mnLarvaeMarchSE = paste(sprintf("%.2f", mnCountLarvaeMarch), sprintf("%.2f", seCountLarvaeMarch), sep = " \\pm "),
    mnLarvaeAprilSE = paste(sprintf("%.2f", mnCountLarvaeApril), sprintf("%.2f", seCountLarvaeApril), sep = " \\pm ")
  ) %>%
  arrange(rowOrder, host, year) %>%
  dplyr::select(-rowOrder) %>% # exclude cols and select cols must be in separate calls to dplyr::select
  dplyr::select(host, 
                year,
                nSitesSampled, nSitesPositiveAnystage, pSitesPositiveAnyStage,
                nSitesWithMoths, pSitesWithMoths, 
                nSitesWithLarvae, pSitesWithLarvae,
                mnLarvaeMarchSE, mnLarvaeAprilSE, mnMothsSE
  ) %>%
  mutate(host = dupToBlank(host)) %>%
  map_dfr(cleanupLatex) 


# output Latex code
colAlign <- c(
  "l",
  "l", # host 
  "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=4]\n", # year,
  "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2,round-mode=places,round-precision=0]\n", # sites sampled
  "S[table-number-alignment=center,table-figures-decimal=0,table-figures-integer=2]\n", # n sites positive any stage
  "@{}S[table-align-text-pre=false,table-figures-integer=1,round-mode=places,round-precision=2]\n", # p sites positive any stage,
  "c", # n sites positive larvae
  "@{}S[table-align-text-pre=false,table-figures-integer=1,round-mode=places,round-precision=2]\n", # p sites positive larvae
  "c", # n sites positive moths
  "@{}S[table-align-text-pre=false,table-figures-integer=1,round-mode=places,round-precision=2]\n", # p sites positive moths
  "S[separate-uncertainty,table-figures-uncertainty=1,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=2,round-mode=places,round-precision=2]\n", # mean larvae march
  "S[separate-uncertainty,table-figures-uncertainty=1,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=2]\n", # mean larvae april
  "S[separate-uncertainty,table-figures-uncertainty=1,table-number-alignment=center,table-figures-integer=2,table-figures-decimal=2]\n"  # mean moths
)

colDigits <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

addHeader <- list()
addHeader$pos  <- list(0, 0)
addHeader$command <- c(
  paste(c(
    # topline
    " & & & \\multicolumn{6}{c}{\\textit{Plutella} incidence} & ",
    "\\multicolumn{3}{c}{\\textit{Plutella} abundance}\\\\\n",
    "\\cmidrule(lr){4-9}\n", "\\cmidrule(lr){10-12}\n"
  ), collapse = ""),
  # second line
  paste(c(
    "Host & ",
    "{Year} & ", 
    "{$N$} & ",
    "\\multicolumn{2}{c}{$N_{\\textsc{positive}}$} & ",
    "\\multicolumn{2}{c}{$N_{\\textsc{moths}}$} & ",
    "\\multicolumn{2}{c}{$N_{\\textsc{larvae}}$} & ",
    "{\\makecell{Mean $\\pm$ SEM\\\\larvae per site\\\\(March)}} & ",
    "{\\makecell{Mean $\\pm$ SEM\\\\larvae per site\\\\(April)}} & ",
    "{\\makecell{Mean $\\pm$ SEM\\\\moths per site\\\\(March-April)}}\\\\\n"
  ), collapse = "")
)

captionSurveys <- c(
  "The incidence and abundance of \\textit{Plutella} at sites sampled during field surveys 
  conducted in South Australia during March and April in each year from 2014-2016. 
  Presented are the numbers and proportion in parentheses of: 
  Number of sites sampled ($N$), 
  sites ($N_{\\textsc{positive}}$) where any lifestage was detected in either month, 
  sites where moths ($N_{\\textsc{moths}}$) or larvae ($N_{\\textsc{larvae}}$) were detected in either month 
  and proportion in parentheses based on sites sampled for the respective lifestages,
  and the mean $\\pm$ standard error of the mean numbers of larvae collected from positive plants in each month, and moths trapped in 
  pheromone traps during the trapping period at positive sites."
)


surveySummaryLatex %>% 
  xtable(align  = colAlign,
         digits = colDigits,
         caption = captionSurveys,
         lab = "tab:autumnSurveys") %>%
  print.xtable(include.rownames = FALSE,
               include.colnames = FALSE,
               table.placement = "p",
               caption.placement = "top",
               add.to.row = addHeader,
               NA.string = "",
               #scalebox = 0.5,
               booktabs = TRUE,
               sanitize.text.function = function(x){x})

# End script
# #################