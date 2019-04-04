# RELICS
RELICS: Regulatory Element Location Identification  in CRISPR screens

`RELICS` uses a generalized linear mixed model (GLMM) to analyze CRISPR regulatory screens. Given a set of training guides (usually positive and negative controls) `RELICS` identifies regions which are more similar to the positive controls (enhancer like) and regions which are more similar to negative controls (non-regulatory like).

This work is still in progress. Please ask questions for support or [post issues](https://github.com/patfiaux/RELICS/issues).

# Installation:
RELICS uses [R](https://cran.r-project.org/bin/windows/base/). Please make sure you have R version 3.5.1 or higher

## Obtain source code
Download source code to your desired location: https://github.com/patfiaux/RELICS.git

## Install requirements
To run RELCIS you need the following packages (if you don't install them using the command after the '#'):
## R packages
dplyr # install.packages('dplyr')

ggplot2 # install.packages('ggplot2')

pROC # install.packages('pROC')

glmmTMB # install.packages('glmmTMB')

## Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")
    
IRanges # BiocManager::install("IRanges", version = "3.8")

GenomicRanges # BiocManager::install("GenomicRanges", version = "3.8")


