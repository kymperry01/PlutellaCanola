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

These scripts were used in field studies of the colonisation of canola crops by P. xylostella (Chapter 5)
1. [autumnSurveys.R](autumnSurveys.R)
    + Plots geographic maps of autumn survey data
    + Outputs LaTeX table code

Supplementary movies (Chapter 5)
1. [Plutella Climex and Canola Colonisation](https://doi.org/10.25909/5bebc11b1f1d4)
    + Spatio-temporal maps at weekly intervals of the colonisation patterns in canola crops in South Australia by P. xylostella, and the CLIMEX Weekly Growth Index for the insect. 
