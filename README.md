# RELICS :sparkles:: Regulatory Element Location Identification  in CRISPR screens

RELICS uses a generalized linear mixed model (GLMM) to analyze CRISPR regulatory screens. Given a set of guides to train on, RELICS identifies regions with high similarity to the regulatory (enhancer-like) guides and regions which are more similar to non-regulatory (negative controls) guides.

This work is continously being improved. Please ask questions or [post issues](https://github.com/patfiaux/RELICS/issues).

# Installation
RELICS runs on [R](https://cran.r-project.org/bin/windows/base/). Please make sure you have R version 3.5.1 or higher

## Obtain source code
Clone source code to your desired location with the following command: ```git clone https://github.com/patfiaux/RELICS.git```. Alternatively, download the repository.

## Install requirements
You will need the following packages to run RELICS. If you don't have them, install them using the following commands. Installations will take about 5 minutes on a standard laptop.
## R packages
dplyr
```r
install.packages('dplyr')
```
ggplot2
```r
install.packages('ggplot2')
```
pROC
```r
install.packages('pROC')
```
glmmTMB
```r
install.packages('glmmTMB')
```
## Bioconductor packages
```r
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")
``` 
IRanges
```r
BiocManager::install("IRanges", version = "3.8")
```

GenomicRanges
```r
BiocManager::install("GenomicRanges", version = "3.8")
```

## Input data format
RELICS reqires two different files as input: 
1. A **guide information file**, containing information about all the simulated guides (chromosome, start, end, label).
2. A **counts file**, containing the counts for each guide in each pool.

The counts file contains only the counts for each guide for each experiment. Column names are necessary but the names do not matter as the user will index the columns by number, not by name.

Example count file: 2 replicates from a FACS experiment. Input pools was sorted into high, medium and low expression

| repl1_input | repl1_high | repl1_med | repl1_low | repl2_input | repl2_high | repl2_med | repl2_low |
|----------|----------|-------|------- |------|------|------|------|
| 11 | 9 | 12 | 11 | 152 | 119 | 189 | 102 |
| 68 | 81 | 39 | 67 | 360 | 339 | 280 | 821 |
| 96 | 89 | 109 | 17 | 3 | 4 | 5 | 0 |
| 104 | 97 | 116 | 38 | 190 | 198 | 194 | 23 |

The guide information file contains all remaining info about the guides such as targeting position and type of guide (positive control, negative control, exon targeting etc.). Non-targeting controls should be specified by setting `chrom`, `start`, and `end` to NA. The columns specifying chromosome, guide target start, guide target end and label are mandatory and must be labelled `chrom`, `start`, `end`, and `label` respectively.

| chrom | start | end | label |
|----------|----------|----------|----------|
| chr8 | 128704468 | 128704488 | chr |
| chr8 | 128704469 | 128704489 | chr |
| NA | NA | NA | neg |
| chr8 | 128704482 | 128704502 | exon |

Row 1 in the count file should correspond to the guide in row 1 in the info file.

## Quickstart with example data
In an interactive R session:

### 1. Source the script
```r
source('/path/to/script/RELICS.r')
```

### 2. Set up the analysis specification file
Several parameters must be specified by the user before running RELICS. These parameters are set within the specification file. This allows a user to refer back to the specification file for the parameters that were used for a particular analysis. The options below  describe the most important parameters required to get the analysis going.

#### Option 1: Modify the given template in the `Example_data` folder (`Type_3_analysis_specs.txt`)
There is a template specification file already set up for the example data. It contains the main flags required to run RELICS. The meaning of the different flags are discussed in the next section. 

#### Option 2: Set the flags within `R` and save them to the specification file prior to analysis
It is also possible to create a specification file from scratch by setting parameters within R. The following steps demonstrate how to set up parameters for the example specification file. At the end of this section, you will be able to run RELICS on the example data.

1. Flags are set up in a list object

```r
analysis.specs <- list()
```

2. Set the output name of the analysis (`dataName`). Be sure to choose a different name from the existing file so that you don't overwrite the example, and you can compare and check that you got the same flags.

```r
analysis.specs$dataName <- 'Type_3_exampleSim'
```

3. Specify the path to the count and info files.

```r
analysis.specs$CountFileLoc <- './Type_3_simulated_counts.csv'
analysis.specs$sgRNAInfoFileLoc <- './Type_3_simulated_info.csv'
```

4. RELICS uses a GLMM and jointly analyzes all pools from each replicate. The replicates are separated by a semicolon (`;`) and each pool from the count file is referred to by number. 

Example: `1,2,3,4;5,6,7,8` represents two replicates. The first four columns of the count file comprise replicate 1 and the next four columns comprise replicate 2.

Note 1: The input type is a string, such as the example shown above.

Note 2: Analysis across multiple replicates has not been implemented yet, so jointly analyzing all pools (`1,2,3,4,5,6,7,8`) is not advised!

```r
analysis.specs$repl_groups <- '1,2,3,4;5,6,7,8'
```

5. The information file has a `label` column, such that each guide has a label assigned to it. Specify which guides are to be used to train the GLMM regulatory and non-regulatory parameters. Both `glmm_positiveTraining` and `glmm_negativeTraining` take either a string or a vector of strings. Regulatory guides are usually promoter- or exon-targeting guides. Non-regulatory guides can be non-targeting guides or all guides not used as regulatory guides. Although the latter option increases the runtime, we have observed that it typically reduces the noise in the results.

```r
analysis.specs$glmm_positiveTraining <- 'exon' # use all guides that overlap an exon to train the regulatory parameters
analysis.specs$glmm_negativeTraining <- 'neg' # use all negative control guides to train non-regulatory parameters
# alternatively use: analysis.specs$glmm_negativeTraining <- c('chr', 'neg') to include everything except positives as negatives
```

6. Specify the method to use as RELICS
```r
analysis.specs$Method <- 'RELICS-search'
```

7. Specify the CRISPR system used and the range of the perturbation effect. Depending on the CRISPR system used, the range of effect is different. We recommend setting the range to 20bp for `CRISPRcas9`, 1000bp for `CRISPRi` and `CRISPRa`. The `crisprEffectRange` is automatically added during the analysis so there is no need to manually extend the `start` and `end` position of the guide. However, if the range has already been accounted for in the information file, set the `crisprEffectRange` to zero.
In case of a `dualCRISPR` system, an arbitrary `crisprEffectRange` can be specified, as RELICS will automatically use the deletion range between guide 1 and guide 2 as effect range.
```r
analysis.specs$crisprSystem <- 'CRISPRi' # other options: CRISPRcas9, CRISPRa, dualCRISPR
analysis.specs$crisprEffectRange <- 1000
```

8. Once you set your flags, create a specification file using the `write_specs_file()` function. It takes two arguments: the list of flags that you have just set (`analysis.specs`) and the name of the file to write to (`.txt` will be automatically appended as the file extension):
```r
write_specs_file(analysis.specs, 'Type_3_exampleSim_specs')
```

### 3. Run RELICS
Once your specification file has been set up, simply use the `analyze_data()` function and pass it the name of the specification file (sans the `.txt` file extension) to begin the RELICS analysis. The example provided should take about 5 minutes to run on a typical desktop computer.
```r
analyze_data('Type_3_exampleSim_specs.txt') # or whatever you named your spec. file
```

### 3. Output files
RELICS will return several output files. They all start with the `dataName` specified in step 2 above:

* `{dataName}_RELICS_genomeScores.bedgraph`: This file contains the RELICS scores in bedGraph format and allows you to visualize your results in your preferred genome browser.

* `{dataName}_RELICS_genomeScores.csv`: This file contains the genome scores in bedGraph format. This file has 7 columns: 

|Column name | Column description |
|----------|----------|
| genomeScore | combined per-guide RELICS score for this region | 
| chrom | chromosome of the region |
| start | region start |
| end | region end |
| label | highest overlapping label according to the [label hierarchy](https://github.com/patfiaux/RELICS#advanced-flags) |
| log2_FC | does not apply to RELICS but kept for backward compatibility |
| nrSupportGuides | number of guide effect ranges which overlap this particular region |

* `{dataName}_RELICS_guideScores.csv`: This file contains the per-guide RELICS scores, which are combined aross regions of overlapping effects. Minimum number of columns: 8

|Column name | Column description |
|----------|----------|
| guide_score | RELICS score per guide | 
| raw_scores | identical to `guide_score`, kept for backward compatibility | 
| chrom | chromosome of the region |
| start | region start |
| end | region end |
| label | guide label |
| log2_FC | does not apply to RELICS but kept for backward compatibility |
| replX_bf | RELICS score for replicate `X`. Scores are reported for each replicate. |

* `{dataName}_RELICS_parSummary.csv`: This file contains RELICS parameter estimates for all replicates, for each pool, for regulatory and background model. 

* `{dataName}_RELICS_replX_parEst.csv`: This file contains RELICS parameter estimates from the positve and negative controls for replicate `X`. For each replicate, the first pool is used as intercept, and all subsequent ones defined as deviance from intercept. 

# Advanced flags

RELICS combines information of guides which overlap with their guide effect. This can lead to scenarios where guides with different labels overlap. By default the label with fewer occurances in the data set is chosen. However, it is also possible for the user to specify the hierarchy by explicitly setting the `labelHierarchy` flag.

The rightmost label has highest priority. Using the example below: if a region has overlapping guides labeled as both exon overlapping (`exon`) as well as targeting guides with unknown effect (`chr`), then the region will be assigned the label with higher priority in the hierarchy - in this case being `exon`.

This example dataset was part of the simulations used to assess method performance. In addition to the usual `exon` and `neg` labels, it also contains a `pos` label for guides overlapping simulated enhancer regions. If specifying the `labelHierarchy`, all labels should be provided. Guides with labels that were not included will not be properly used the the analysis.

```r
analysis.specs$labelHierarchy <- c('chr', 'neg', 'exon', 'pos')
```
