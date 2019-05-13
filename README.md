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
RELICS reqires two different files as input. One contains only the counts for each guide for each experiment. Column names are necessary but the names do not matter as the user will refer to the columns by number, not by name.

Example count file: 2 replicates from a FACS experiment. Input pools was sorted into high, medium and low expression

| repl1_input | repl1_high | repl1_med | repl1_low | repl2_input | repl2_high | repl2_med | repl2_low |
|----------|----------|-------|------- |------|------|------|------|
| 11 | 9 | 12 | 11 | 152 | 119 | 189 | 102 |
| 68 | 81 | 39 | 67 | 360 | 339 | 280 | 821 |
| 96 | 89 | 109 | 17 | 3 | 4 | 5 | 0 |
| 104 | 97 | 116 | 38 | 190 | 198 | 194 | 23 |

The second file contains all remaining info about the guides such as targeting position and type of guide (positive control, negative control, exon targeting etc.). The columns specifying chromosome, guide target start, guide target end and label (chrom, start, end label) are mandatory.

| chrom | start | end | label |
|----------|----------|----------|----------|
| chr8 | 128704468 | 128704488 | chr |
| chr8 | 128704469 | 128704489 | chr |
| NA | NA | NA | neg |
| chr8 | 128704482 | 128704502 | exon |

## Quickstart with example data
### 1. source the script
`source('/path/to/script/RELICS.r')`

### 2. Setting up the analysis specification file. 
### Option 1: Modify the given template (Type_3_analysis_specs.txt)
### Option 2: Set the flags within `R` and then write them to the specs file prior to analysis
Below is an example on how to sep up the flags for the example file

Flags are set up in a list format

`analysis.specs <- list()`

Set the output name of the analysis

`analysis.specs <- 'Type_3_exampleSim'`

Give location of count and info files (easies if in same directory as the analysis is done but can also give a path to files)

```
analysis.specs$CountFileLoc <- 'Type_3_simulated_counts.csv'
analysis.specs$sgRNAInfoFileLoc <- 'Type_3_simulated_info.csv'
```

Specify the label hierarchy. This is used when labeling regions after combining overlapping guide effects. Rightmost label has highest prioirty. As an example; if a region has overlapping guides labelled as both positives ('pos') as well as targeting guides with unknown effect ('chr') then the region will be assigned the highest label from the hierarch, namely 'pos'

`analysis.specs$labelHierarchy <- c('chr', 'neg', 'pos')`

RELICS uses a GLMM and jointly analyzes all pools from each replicate. The replicates are separated by a semicolon (';') and each pool from the count file is referred to by number. 

Example: '1,2,3,4;5,6,7,8' has 2 replicates. The first four columns of the count file make up replicate 1 and column 5, 6, 7 and 8 make up replicate 2.

Note 1: Becasue of the separation by semicolon the input here is a string, not numeric!

Note 2: Analysis across multiple replicates has not been implemented yet so jointly analyzing all pools ('1,2,3,4,5,6,7,8') is not advised!

`analysis.specs$repl_groups <- '1,2,3,4;5,6,7,8'`

