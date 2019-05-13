# RELICS
RELICS: Regulatory Element Location Identification  in CRISPR screens

`RELICS` uses a generalized linear mixed model (GLMM) to analyze CRISPR regulatory screens. Given a set of training guides `RELICS` identifies regions which are more similar to the positive controls (enhancer like) and regions which are more similar to negative controls (non-regulatory like).

This work is continously being improved. Please ask questions for support or [post issues](https://github.com/patfiaux/RELICS/issues).

# Installation:
RELICS uses [R](https://cran.r-project.org/bin/windows/base/). Please make sure you have R version 3.5.1 or higher

## Obtain source code
Download source code to your desired location: `git clone https://github.com/patfiaux/RELICS.git`

## Install requirements
To run RELCIS you need the packages below. If you don't have them, install them using the command after the '#'):
## R packages
```
dplyr # install.packages('dplyr')

ggplot2 # install.packages('ggplot2')

pROC # install.packages('pROC')

glmmTMB # install.packages('glmmTMB')
```
## Bioconductor packages
```
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")
    
IRanges # BiocManager::install("IRanges", version = "3.8")

GenomicRanges # BiocManager::install("GenomicRanges", version = "3.8")
```

## Input data format
RELICS reqires two different files as input. One contains only the counts for each guide for each experiment. The other contains all remaining info about the guides such as targeting position and type of guide (positive control, negative control etc.)

Example count file:

| repl1_input | repl1_high | repl1_med | repl1_low | repl2_input | repl2_high | repl2_med | repl2_low |
| --- | --- | --- | --- | --- | --- |
| 11 | 9 | 12 | 11 | 152 | 119 | 189 | 102 |
68	81	39	67	360	339	280	821
96	89	109	17	3	4	5	0
104	97	116	38	190	198	194	23



## Quickstart with example data
1. source the script
`source('/path/to/script/RELICS.r')`

2. 
