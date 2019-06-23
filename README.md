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

The second file contains all remaining info about the guides such as targeting position and type of guide (positive control, negative control, exon targeting etc.). Non-targeting controls should be specified by setting `chrom`, `start`, and `end` to NA. The columns specifying chromosome, guide target start, guide target end and label (chrom, start, end label) are mandatory.

| chrom | start | end | label |
|----------|----------|----------|----------|
| chr8 | 128704468 | 128704488 | chr |
| chr8 | 128704469 | 128704489 | chr |
| NA | NA | NA | neg |
| chr8 | 128704482 | 128704502 | exon |

Row 1 in the count file should correspond to the guide in row 1 in the info file.

## Quickstart with example data
### 1. source the script
```
source('/path/to/script/RELICS.r')
```

### 2. Setting up the analysis specification file. 
There are several different parameters which have to be specified by the users before running RELICS. These parameters are set within the specification file. This has the advantage that a user can go back to an analysis at any time and see under what conditions a data set was analyzed. Below is the outline on how to set the most important flags to get the analysis going.

#### Option 1: Modify the given template in the 'Example_data' folder (Type_3_analysis_specs.txt)
There is a template specification file already set up for the example data. It contains the main flags required for successfully running RELICS. The meaning of the different flags are discussed in 'Option 2' below. 

#### Option 2: Set the flags within `R` and then write them to the specification file prior to analysis
In case of the absence of a specification file it is also possible to create it. Below is an example on how to sep up the flags for the example specification file along with their meaning. At the end of this you will be able to run RELICS on the example data

Flags are set up in a list format

```
analysis.specs <- list()
```

Set the output name of the analysis (and chose a different name from the existing file so you can compare and check that you got the same flags.

```
analysis.specs$dataName <- 'Type_3_exampleSim'
```

Give location of count and info files (easiest if in working directory but can also give a path to files)

```
analysis.specs$CountFileLoc <- 'Type_3_simulated_counts.csv'
analysis.specs$sgRNAInfoFileLoc <- 'Type_3_simulated_info.csv'
```

RELICS uses a GLMM and jointly analyzes all pools from each replicate. The replicates are separated by a semicolon (';') and each pool from the count file is referred to by number. 

Example: '1,2,3,4;5,6,7,8' has 2 replicates. The first four columns of the count file make up replicate 1 and column 5, 6, 7 and 8 make up replicate 2.

Note 1: Because of the separation by semicolon the input here is a string, not numeric!

Note 2: Analysis across multiple replicates has not been implemented yet so jointly analyzing all pools ('1,2,3,4,5,6,7,8') is not advised!

```
analysis.specs$repl_groups <- '1,2,3,4;5,6,7,8'
```

RELICS empirically estimates the GLMM parameters from the data. A set of positive and negative controls should be provided. Positive controls are usually promoter or exon targeting guides. Negative controls could be non-targeting guides. Another option is to specify everything that's not a positive control as negative control. While this increases the runtime we have observed that this reduces the noise seen in the data.

```
analysis.specs$glmm_positiveTraining <- 'exon'
analysis.specs$glmm_negativeTraining <- 'neg' 
# alternatively use: analysis.specs$glmm_negativeTraining <- c('chr', 'neg') to include everything except positives as negatives
```

Specify the method to use as RELICS
```
analysis.specs$Method <- 'RELICS-search'
```

Depending on the CRISPR system used the range of effect is different. We recommend setting the range to 20bp for `CRISPRcas9`, 1000bp for `CRISPRi` and `CRISPRa`. Note that the effect range is added to the positions specified in the info file. If the effect range is already included in the positions of the info file then it should be set to 0 here. 

In case of a `dualCRISPR` system an arbitrary `crisprEffectRange` can be specified as RELICS will automatically use the deletion range between guide 1 and guide 2 as effect range.
```
analysis.specs$crisprSystem <- 'CRISPRi' # other potions: CRISPRcas9, CRISPRa, dualCRISPR
analysis.specs$crisprEffectRange <- 1000
```

Once you have your flags set, create a specification file using the `write_specs_file()` function. The two arguments it takes are the list with flags you just set and the name of the file (.txt will be added automatically so don' include that)
```
write_specs_file(analysis.specs, 'Type_3_exampleSim_specs')
```


### 3. Run RELICS
Once you have your specification file set up simply use the `analyze_data()` function to start the RELICS analysis. For the example given it will take about 5 min, depending on your operating system.
```
analyze_data('Type_3_exampleSim_specs.txt') # or whatever you named your spec. file
```

### 3. Ouput files
RELICS will return several files. They all start with the dataName you specified above:

\* _RELICS_genomeScores.bedgraph: Contains the RELICS scores in bedgraph format and allows you to visualize your results in your preferred Genome browser.

\*_RELICS_genomeScores.csv : Contains the genome scores set up in bedgraph format. This file has 7 columns. 
> genomeScore: combined per-guide RELICS score for this region

> chrom, start, end: position of the region

> label: highest overlapping label according to the label hierarchy

> log2_FC: does not apply to RELICS but kept for backward compatibility

> nrSupportGuides: number of guide effect ranges which overlap this particular region

\*_RELICS_guideScores.csv*: Contains the per-guide RELICS scores which are combined aross regions of overlapping effects. Minimum number of columns: 8

> raw_scores, guide_scores: contain identical values. `raw_scores` was kept for backward compatibility. Use `guide_scores` when working with this file

> chrom, start, end: position of the region

> label: same as specified in the information file

> log2_FC: does not apply to RELICS but kept for backward compatibility

> replX_bf: RELICS score (Bayes Factor) for replicate `X`. Scores are reported for each replicate.

\*_RELICS_parSummary.csv: Contains RELICS parameter estimates for all replicates, for each pool, for regulatory and background model. 

\*_RELICS_replX_parEst.csv: Contains RELICS parameter estimates from the positve and negative controls for replicate `X`. For each replicate, the first pool is used as intercept, all subsequent ones defined as deviance from intercept. 

# Advanced flags for experienced RELICS users

RELICS combines information of guides which overlap with their guide effect. This can lead to scenarios where guides with different labels overlap. By default the label with fewer occurances in the data set is chosen. However, it is also possible for the user to specify the tie breaking by explicitly setting the labelHierarchy flag.
Rightmost label has highest priority. As an example; if a region has overlapping guides labeled as both exon overlapping ('exon') as well as targeting guides with unknown effect ('chr') then the region will be assigned the higher label from the hierarch, namely 'exon'.

This example data set was part of the simulations used to assess method performance and in addition to the usual 'exon' and 'neg' labels it also contains 'pos', labelling guides overlapping simulated enhancer regions. Below we will train on exon overlapping guides and negative controls to then identify the simulated positive regions ('pos').

```
analysis.specs$labelHierarchy <- c('chr', 'neg', 'exon', 'pos')
```
