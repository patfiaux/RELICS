source('../Code//RELICS.v2.R')

# initiate the relics flgas list
relics.parameters <- list()

# name of the replics. All output files will start with this
relics.parameters$dataName <- 'CD69_example_analysis'

# location of the data file
relics.parameters$DataInputFileLoc <- 'Example_data/CD69_data_guideScored_filtered.csv'

# within the data file, specify the location of the guide counts. Each replicate is a new list.
repl.pools <- list()
repl.pools[[1]] <- c(5, 7, 9, 11, 13)
repl.pools[[2]] <- c(6, 8, 10, 12, 14)
relics.parameters$repl_groups <- repl.pools

# specify what label to use as known functional sequence (FS0)
relics.parameters$FS0_label <- 'CD69_promoter'

# specify the expected number of functional sequences and how many to look for in total
relics.parameters$max_fs_nr <- 15
relics.parameters$expected_fs_nr <- 5 # expected based on previous findings by Simeonov et al., 2017 and Fiaux et al., 2020

# specify the CRISPR system
relics.parameters$crisprSystem <- 'CRISPRa' # other options: CRISPRcas9, CRISPRi, dualCRISPR

# ouput file
relics.parameters$out_dir <- 'CD69_example_analysis/'

# specify the number of bins to group the guide counts into and the degrees of freedom of the spline function for each replicate
relics.parameters$nr_disp_bins <- 15
relics.parameters$repl_spline_df <- list(repl_1 = 3, repl_2 = 3)

#########################
# additional flags

# specify the names of the different pools. Helps with naming hyperparameters in the output files
all.repl.names <- list()
all.repl.names[[1]] <- c("back", "baseline", "low", "medium", "high")
all.repl.names[[2]] <- c("back", "baseline", "low", "medium", "high")
relics.parameters$pool_names <- all.repl.names

# guide efficiency scores can be used to further improve results. The guide scores can either be provided in a separate file or
# be part of the data file. In this example only one set of scores has been provided (at column 21) but it's possible to provide multiple scores
relics.parameters$guide_efficiency_loc <- relics.parameters$DataInputFileLoc
relics.parameters$guide_efficiency_cols <- c(21)

relics.parameters$record_pp <- TRUE
relics.parameters$record_cs_pp <- TRUE

# additional flags end
#########################

# this should take ~6-7 min to run on a laptop
RELICS(input.parameter.list = relics.parameters)
