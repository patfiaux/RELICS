# RELICS
RELICS: Regulatory Element Location Identification  in CRISPR screens

The `RELICS` package use a generalized linear mixed model (GLMM) to analyze CRISPR regulatory screens. Given a set of training guides (usually positive and negative controls) `RELICS' identifies regions which are more similar to the positive controls (enhancer like) and regions which are more similar to negative controls (non-regulatory like).


# Install requirements
To run RELCIS you need R and have the following packages:
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


