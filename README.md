# PlutellaCanola 

R scripts used for PhD thesis

These scripts were used in RAD-seq studies of Plutella xylostella and the cryptic moth, Plutella australiana (Chapters 3, 4):

1. [genotypes.R](genotypes.R)
    + makes maps, LaTex tables etc
2. [amova.R](amova.R) 
    + Runs AMOVA for P. xylostella using 1032 SNPs, outputs LaTeX table 
3. [bioassays.R](bioassays.R)
    + Performs log-logistic regression of insecticide bioassays
    + Outputs plots of dose-response curves
4. [checkGenotypes.R](checkGenotypes.R)
    + Performs PCA on technical duplicates to check the robustness of genotype calls (i.e. whether sample pairs cluster together) 
5. [collectionsPx.R](collectionsPx.R)
    + Makes LaTeX tables and plots geographic maps of P. xylostella collections
6. [crosses.R](crosses.R)
    + Makes LaTeX tables summarising experimental inter-species crosses between P. xylostella and P. australiana
    

These scripts were used in field studies of the colonisation of canola crops by P. xylostella (Chapter 5)
1. [autumnSurveys.R](autumnSurveys.R)
    + Plots geographic maps of autumn survey data
    + Outputs LaTeX table code
2. [climate.R](climate.R)
    + Extracts climate data from raster grid cells corresponding to field site for 2014 to 2016
    + Line plots of cumulative rainfall and CLIMEX Weekly Growth Index
3. [climex.R](climex.R)
    + Plots maps and spatio-temporal movies of CLIMEX "Compare Locations Years" models
4. [colonisation.R](colonisation.R)
    + Models temperature-dependent P. xylostella phenology
    + Plots spatio-temporal maps of colonisation of canola and CLIMEX Weekly Growth Index
    

Supplementary movies (Chapter 5)
1. [Plutella Climex and Canola Colonisation](https://doi.org/10.25909/5bebc11b1f1d4)
    + Spatio-temporal maps at weekly intervals of the colonisation patterns in canola crops in South Australia by P. xylostella, and the CLIMEX Weekly Growth Index for the insect. 
