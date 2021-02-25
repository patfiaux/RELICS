
# This script contains two examples for how to run RELICS on the example data set
# Both ways return the same result

# The first option is to provide a .txt file with the analysis parameters to use
RELICS('Example_analysis_specifications.txt')

# The second option is to provide the parameters directly in R, as a list
repl.pools <- list()
repl.pools[[1]] <- c(5, 7, 9, 11, 13)
repl.pools[[2]] <- c(6, 8, 10, 12, 14)

relics.parameters <- list(dataName = 'CD69_example_analysis',
                          repl_groups = repl.pools, 
                          DataInputFileLoc = 'Example_data/CD69_data_example.csv',
                          min_FS_nr = 8,
                          FS0_label = 'CD69_promoter',
                          crisprSystem = 'CRISPRa',
                          out_dir = 'CD69_example_analysis/')

RELICS(input.parameter.list = relics.parameters)
