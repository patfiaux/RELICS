#source('RELICS_2.R')
source('/iblm/netapp/home/pfiaux/RELICS_2/dev_Code_30/RELICS.v2.R')

# initiate the relics flgas list
analysis.specs <- list()

# name of the analysis. All output files will start with this
analysis.specs$dataName <- 'CD69_example_analysis'

# location of the data file
analysis.specs$DataInputFileLoc <- 'Example_data/CD69_data_guideScored_filtered.csv'

# within the data file, specify the location of the guide counts. Each replicate is a new list.
repl.pools <- list()
repl.pools[[1]] <- c(5, 7, 9, 11, 13)
repl.pools[[2]] <- c(6, 8, 10, 12, 14)
analysis.specs$repl_groups <- repl.pools

# specify what label to use as known functional sequence (FS0)
analysis.specs$FS0_label <- 'CD69_promoter'

# specify the expected number of functional sequences and how many to look for in total
analysis.specs$max_fs_nr <- 15
analysis.specs$expected_fs_nr <- 5 # expected based on previous findings by Simeonov et al., 2017 and Fiaux et al., 2020

# specify the CRISPR system
relics.parameters$crisprSystem <- 'CRISPRa' # other options: CRISPRcas9, CRISPRi, dualCRISPR

# ouput file
relics.parameters$out_dir <- 'CD69_tutorial_output/'

# specify the number of bins to group the guide counts into and the degrees of freedom of the spline function for each replicate
analysis.specs$nr_disp_bins <- 15
analysis.specs$repl_spline_df <- list(repl_1 = 3, repl_2 = 3)


#####
# additional flags

# specify the high and low pools to compute the log2 fold change
l2fc.repl.groups <- list()
l2fc.repl.groups[[1]] <- c(7, 13)
l2fc.repl.groups[[2]] <- c(8, 14)
analysis.specs$l2fc_groups <- l2fc.repl.groups

# specify the names of the different pools. Helps with naming hyperparameters in the output files
all.repl.names <- list()
all.repl.names[[1]] <- c("back", "baseline", "low", "medium", "high")
all.repl.names[[2]] <- c("back", "baseline", "low", "medium", "high")
analysis.specs$pool_names <- all.repl.names

# the order of importance of labels can be specified. It is also possible to remove a set of guides with a label by leaving said class out from the label hierarchy
analysis.specs$labelHierarchy <- c('chr', 'exon', 'CD69_promoter')

# guide efficiency scores can be used to further improve results. The guide scores can either be provided in a separate file or
# be part of the data file. In this example only one set of scores has been provided (at column 21) but it's possible to provide multiple scores
analysis.specs$guide_efficiency_loc <- analysis.specs$DataInputFileLoc
analysis.specs$guide_efficiency_cols <- c(21)

analysis.specs$record_pp <- TRUE
analysis.specs$record_cs_pp <- TRUE
analysis.specs$record_deltas <- TRUE

analysis.specs$fixed_ge_coeff <- FALSE

# additional flags end
#####

RELICS(input.parameter.list = analysis.specs, record.all.fs = TRUE, data.file.split = F)

