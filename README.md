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
Once you have set up your parameters you can run RELICS by directly giving it the list we set up above, or by first saving it to a `.txt` file. In the latter case, the flags and their values should be separated by a colon (`:`, see `Example_analysis_specifications.txt`).
The CD69 example provided should take about 10 minutes to run on a typical desktop computer.
```r
RELICS(input.parameter.list = relics.parameters)

# option: use the parameter file instead
# RELICS('Example_analysis_specifications.txt')
```

### 4. Output files
RELICS will return several output files. They all start with the `dataName` specified in step 2 above. By default, RELICS will give you the genome segments that were used, as well as the files associated with finding the last functional sequence before convergence:

* `{dataName}_segmentInfo.csv`: This file contains the segments used by RELICS. It contians the information of chromosome, start and end location of the segment, as well as the label of the segment.
|Column name | Column description |
|----------|----------|
| chrom | chromosome of the region |
| start | region start |
| end | region end |
| label | highest overlapping label according to the [label hierarchy](https://github.com/patfiaux/RELICS#advanced-flags) |

In all subsequent file names, the pattern `_kX_` refers to `X` functional sequences detected.

* `{dataName}_final_kX_total_pp.bedgraph`: This file contains the sum of posteriors across all functional sequences detected.

* `{dataName}_final_kX_FS_locations.bed`: This file contains all genome segments part of the functional sequences detected.

* `{dataName}_final_kX.csv`: This file contains the functional sequence probabilities of all functional sequences detected. Each column corresponds to a genome segment, ordered as in `{dataName}_segmentInfo.csv`. Each row correponds to the functional sequence probabilities of a particular functional sequence. The first row corresponds to FS0, the second to FS1 etc.

* `{dataName}_final_kX_ll_progression.csv`: This file keeps track of the -log-likelihood model improvement with each additional functioanl sequence detected. Correctly detecting an additional functional sequence should improve the model fit if this is supported by the data. The initial controbutions are usually quite large and then start plateauing as all functional sequences are detected.
|Column name | Column description |
|----------|----------|
| FS | the functional sequence which is included in the overall model |
| FS_ll | the -log-likelihood of the model by including all FS functional sequences |

* `{dataName}_final_kX_perFS_LLcontributions.csv`: This file contains the per-functional sequence contribution to the model improvement.
|Column name | Column description |
|----------|----------|
| FS | the functional sequence |
| ll | the number of functional sequences included in the model |
| nr_fs | the number of genome segments considered to be part of the functional sequence as defined by the functional sequence threshold |

* `{dataName}_final_k4_alphas.csv`: This file contains the dirichlet sorting parameters for the background (`alpha0`) and for functional sequences (`alpha1`) for all replicates.


* `{dataName}_final_allCorrs_k(X+1).csv`: This file contains all pairwise correlation between the functional sequence probabilities up to and including the functional sequence which resulted in a correlation above the convergence threshold and leading to RELICS to stop.


* `{dataName}__summaryStatPlots.pdf`: Plot showing the overall model -log-likelihood progression and the correlation of the functional sequence probabilities.


* `{dataName}_final_kX.tiff`: This file contains the plots of the individual functional sequence probabilities and notes how many segments are contained within each. Segemnts which are above the functional sequence threshold are labelled in purple. The 'Sum of Poseteriors' shows the sum of all functional sequence probabilities.
* `{dataName}_final_kX_FS_locations.bed`: This file contains all genome segments part of the functional sequences detected.


# Advanced flags

RELICS combines information of guides which overlap with their guide effect. This can lead to scenarios where guides with different labels overlap. By default the label with fewer occurrences in the data set is chosen. However, it is also possible for the user to specify the hierarchy by explicitly setting the `labelHierarchy` flag.

The rightmost label has highest priority. Using the example below: if a region has overlapping guides labeled as both promoter overlapping (`CD69_promoter`) as well as targeting guides with unknown effect (`chr`), then the region will be assigned the label with higher priority in the hierarchy - in this case being `CD69_promoter`.

If specifying the `labelHierarchy`, all labels should be provided. Guides with labels that were not included will not be properly used the analysis. In this case this means also specifying the `exon` label fo rguides overlapping CD69 exons.

```r
relics.parameters$labelHierarchy <- c('chr', 'exon', 'CD69_promoter')
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
