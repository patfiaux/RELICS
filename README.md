# RELICS :sparkles:: Regulatory Element Location Identification  in CRISPR screens

RELICS is a method of analyzing tiling CRISPR screens for detecting functional sequences. The current version (v.2.0) of RELICS uses a Bayesian hierarchical model and considers the overlapping effects of multiple guides, can jointly analyze multiple pools per replicate and estimates the number of functional sequences supported by the data.

This work is continuously being improved. Please ask questions or [post issues](https://github.com/patfiaux/RELICS/issues).

# Installation
RELICS runs on [R](https://cran.r-project.org/bin/windows/base/). Please make sure you have R version 3.5.1 or higher

## Obtain source code
Clone source code to your desired location with the following command: ```git clone https://github.com/patfiaux/RELICS.git```. Alternatively, download the repository.

## Install requirements
You will need the following packages to run RELICS. If you don't have them, install them using the following commands. Installations will take about 5 minutes on a standard laptop.
## R packages
ggplot2 (for plotting)
```r
install.packages('ggplot2')
```
gridExtra (for combining multiple plots in one figure)
```r
install.packages('gridExtra')
```
poibin (for calculating the Poisson-Binomial)
```r
install.packages('poibin')
```
extraDistr (for calculating the Dirichlet-Multinomial)
```r
install.packages('extraDistr')
```
gtools (for calculating combinations)
```r
install.packages('gtools')
```
## Bioconductor packages
```r
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")
``` 
IRanges (for handling genomic coordinates)
```r
BiocManager::install("IRanges", version = "3.8")
```

GenomicRanges (for handling genomic coordinates)
```r
BiocManager::install("GenomicRanges", version = "3.8")
```

## Input data format
RELICS requires an input file conatining the sgRNA targets and the corresponding counts in the different pools. 
The required columns must have the following information: chromosome, sgRNA target start, sgRNA target end, sgRNA label, ...(sgRNA counts in different pools)...

The columns specifying chromosome, guide target start, guide target end and label are mandatory and must be labelled `chrom`, `start`, `end`, and `label` respectively.

For the columns containing sgRNA counts, names are necessary but the names do not matter as the user will index the columns by number, not by name.

Example count file: 2 replicates from a FACS experiment. Input pools was sorted into high, medium and low expression

| chrom | start | end | label | repl1_input | repl1_high | repl1_med | repl1_low | repl2_input | repl2_high | repl2_med | repl2_low |
|----------|----------|----------|----------|----------|----------|-------|------- |------|------|------|------|
| chr8 | 128704468 | 128704488 | chr | 11 | 9 | 12 | 11 | 152 | 119 | 189 | 102 |
| chr8 | 128704469 | 128704489 | chr | 68 | 81 | 39 | 67 | 360 | 339 | 280 | 821 |
| NA | NA | NA | neg | 96 | 89 | 109 | 17 | 3 | 4 | 5 | 0 |
| chr8 | 128704482 | 128704502 | exon | 104 | 97 | 116 | 38 | 190 | 198 | 194 | 23 |

## Quickstart with example data
In an interactive R session:

### 1. Source the script
```r
source('/path/to/script/RELICS.v2.r')
```

### 2. Set up the analysis specification file
Several parameters must be specified by the user before running RELICS. Tutorial walkthrough below describes the most important parameters required to get the analysis going. In the process you will analyze a subset of the CD69 CRISPRa screen from (Simeonov et al.)[https://www.nature.com/articles/nature23875].

#### Option 1: Modify the given template in the `RELICS_tutorial` folder (`Example_analysis_specifications.txt`)
There is a template specification file already set up for the example data. It contains the main flags required to run RELICS. The meaning of the different flags are discussed in the next section. 

#### Option 2: Set the flags within `R` and save them to the specification file prior to analysis
It is also possible to create a specification file from scratch by setting parameters within R. The following steps demonstrate how to set up parameters for the example specification file. At the end of this section, you will be able to run RELICS on the example data.

1. Flags are set up in a list object

```r
relics.parameters <- list()
```

2. Set the output name of the analysis (`dataName`). Be sure to choose a different name from the existing file so that you don't overwrite the example, and you can compare and check that you got the same flags.

```r
relics.parameters$dataName <- 'CD69_example_analysis'
```

3. Specify the data file.

```r
relics.parameters$DataInputFileLoc <- './Example_data/CD69_data_example.csv'
```

4. RELICS uses a Dirichlet-Multinomial and jointly analyzes all pools from each replicate. In the specification file the replicates are separated by a semicolon (`;`) and referred to by column position in the data file. Within R, the replicates are each an element within a `list()`

Example: `5,7,9,11,13;6,8,10,12,14` represents two replicates. In this example, the pools of each replicate are altering.

Note 1: The input type is a string, such as the example shown above.

Note 2: Analysis across multiple replicates has not been implemented yet, so jointly analyzing all pools (`5,7,9,11,13,6,8,10,12,14`) is not advised!

```r
repl.pools <- list()
repl.pools[[1]] <- c(5, 7, 9, 11, 13)
repl.pools[[2]] <- c(6, 8, 10, 12, 14)

relics.parameters$repl_groups <- repl.pools
```

5. The `label` column assigns each sgRNA to a category. These categories are used in the beginning to specify the training sets. By default a set of sgRNAs overlapping known functional sequences (FSs) have to be provided to the `FS0_label` flag. In this case these are guides overlapping the CD69 promoter. All other sgRNAs are used to train the background parameters. The counts of the sgRNAs overlapping each FS detected are iteratively added to the set of counts used to determine the Dirichlet parameters of the FS.

As an option, it is also possible to specify the sgRNA labels to be used as the background. This is done with the `background_label` flag.

In both cases the flags are given either as string or as vector of strings.

```r
relics.parameters$FS0_label <- 'CD69_promoter' # use all sgRNAs that overlap the CD69 promoter are used to initially train the FS parameters

# option: specify the background parameters
# relics.parameters$glmm_negativeTraining <- c('chr', 'exon') # specify what sgRNAs to use to initially train the background parameters
```

6. Specify the minimum number of functional sequences to look for. It is possible that the data will have more or less FS than you specify. In this example, RELICS will find 4. Because is reaches convergence after 4 FS it will stop. In case RELICS does not converge we recommend you increase the allowed number of FS.
```r
relics.parameters$min_FS_nr <- 8
```

7. Specify the CRISPR system used. Depending on the CRISPR system used, the area of effect is different. By default RELICS assumes that it's 20bp for `CRISPRcas9` and 400bp for both `CRISPRi` and `CRISPRa`. In case of a `dualCRISPR` system RELICS will automatically use the deletion range between guide 1 and guide 2 as effect range. It is also possible to manually set the area of effect using the `crisprEffectRange` flag.
```r
relics.parameters$crisprSystem <- 'CRISPRa' # other options: CRISPRcas9, CRISPRi, dualCRISPR

# optional: specify the area of effect for your CRISPR system
# relics.parameters$crisprEffectRange <- 200
```

8. Give the loacation of the output directory by setting the `out_dir` flag. Either reference to full path or the path from the current working directory. In this example we will do the latter and assume you are in the `RELICS_tutorial` folder. We recommend you create a new file in which the results are saved. Nore, RELICS will NOT create non-existent files for you. In this example, first create the `CD69_tutorial_output` folder, then set the flag:
```r
relics.parameters$out_dir <- 'CD69_tutorial_output/'
```

### 3. Run RELICS
Once your specification file has been set up, simply use the `analyze_data()` function and pass it the name of the specification file (sans the `.txt` file extension) to begin the RELICS analysis. The example provided should take about 5 minutes to run on a typical desktop computer.
```r
RELICS(input.parameter.list = relics.parameters)

# option: use the parameter file instead
# RELICS('Example_analysis_specifications.txt')
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

* `{dataName}_RELICS_guideScores.csv`: This file contains the per-guide RELICS scores, which are combined across regions of overlapping effects. Minimum number of columns: 8

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

* `{dataName}_RELICS_replX_parEst.csv`: This file contains RELICS parameter estimates from the positive and negative controls for replicate `X`. For each replicate, the first pool is used as intercept, and all subsequent ones defined as deviance from intercept. 

# Advanced flags

RELICS combines information of guides which overlap with their guide effect. This can lead to scenarios where guides with different labels overlap. By default the label with fewer occurrences in the data set is chosen. However, it is also possible for the user to specify the hierarchy by explicitly setting the `labelHierarchy` flag.

The rightmost label has highest priority. Using the example below: if a region has overlapping guides labeled as both exon overlapping (`exon`) as well as targeting guides with unknown effect (`chr`), then the region will be assigned the label with higher priority in the hierarchy - in this case being `exon`.

This example dataset was part of the simulations used to assess method performance. In addition to the usual `exon` and `neg` labels, it also contains a `pos` label for guides overlapping simulated enhancer regions. If specifying the `labelHierarchy`, all labels should be provided. Guides with labels that were not included will not be properly used the analysis.

```r
analysis.specs$labelHierarchy <- c('chr', 'neg', 'exon', 'pos')
```

## Input data format contd. (for backward compatibility)
Instead of providing one joint file containing both coordinates and counts it is also possible to supply them separately. In this case the format is the following:
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
