
print('Loading packages...')
suppressPackageStartupMessages(require(IRanges))
suppressPackageStartupMessages(require(GenomicRanges))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(glmmTMB))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(pROC))

# master analysis function
# input: spec.file
#     specification file contains instructions about analysis types to use, plotting and evaluations
#   optional argument: scoreFile, labelFile
#     score.file: contains the location of the file containing already computed scores (avoids re-calculation)
#     label.file: contains the genomic coordinates for specific labels
#     data.dir: contains the location of the data directory. Used in case data only referenced by name, not loaction)
# output: none
#   There is no direct output to this function
#   Analyzing, file saving and plotting are all done during the function call
# per analysis run, only one type of data re-arrangement can be made
# per analysis run, only one type of filtering can be made
analyze_data <- function(spec.file = NULL, label.file = NULL, data.dir = NULL){

  # read in specification file
  if(is.null(spec.file)) stop('Error in "analyze_data()": specification file has to be provided')
  analysis.specs <- read_analysis_specs(spec.file, data.dir)

  if('noPostScoreAnalysis' %in% names(analysis.specs)){
    analysis.specs$postScoringAnalysis <- NULL
  }
  analysis.specs.names <- names(analysis.specs)

  final.score.list <- list()  #instantiate the final score list (specifiy format?)

  final.score.list <- calculate_scores(analysis.specs, label.file = label.file)

  print('Genome scores obtained ...')

  if('postScoringAnalysis' %in% analysis.specs.names){
    print('Saving genome scores ...')
    save_all_scores(final.score.list, analysis.specs)
  }

  if('createBedGraph' %in% analysis.specs.names){
    print('Creating Bedgraphs ...')
    create_bedgraphs(final.score.list, analysis.specs)
  }

  if('scorePlotting' %in% analysis.specs.names){
    print('Score plotting ...')
    score_plotting(final.score.list, analysis.specs)
  }

}

#output: list. Each analysis specification is stored within for donwstream instructions
# contains many unused flags, but will keep here to keep compatability with
# RELICS-perform repository
read_analysis_specs <- function(in_specs_loc, data.dir = NULL){
  raw_specs <- scan(in_specs_loc,what='character')
  out_specs_list <- list()

  data.loc <- ''  #directory where the data is located (used in case data file is only referenced by name and not entrie path)
  if(! is.null(data.dir)){
    data.loc <- data.dir
  }
  #added during analysis:
  # analysis.specs$labelFile

  out_specs_list$InitialDispersion <- 0.5
  out_specs_list$RRAmaxPercent <- 0.1
  out_specs_list$CRISPReffectRange <- 20
  # default initial dispersion and zero-mass proportion estimate:
  out_specs_list$InitialDispersion <- 0.5
  out_specs_list$InitialEta <- 0.5
  out_specs_list$ZINBminThresh <- 0.1
  # default RRA properties:
  out_specs_list$rraPermutNr <- 5
  out_specs_list$rraNaBin <- 10
  out_specs_list$scoreThresh <- 0.3  # threshold used for calclating alpha-RRA, scores below this won't be considered

  # DESeq2 default is for paired end analysis
  out_specs_list$DESeq2_type <- 'paired-end'

  out_specs_list$dataFilter <- 'none'
  out_specs_list$saveGuideScores <- 'yes'
  out_specs_list$createBedGraph <- 'yes'
  out_specs_list$postScoringAnalysis <- 'yes'
  out_specs_list$RELICS_genomeScoring <- 'yes'
  out_specs_list$dataArranging <- 'none'

  for(spec in raw_specs){
    spec_id <- strsplit(spec,':')[[1]][1]
    if(spec_id == 'dataArranging'){
      out_specs_list$pcrDuplicates <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'CountFileLoc'){
      out_specs_list$CountFileLoc <- paste(data.loc, paste(strsplit(spec,':')[[1]][2:length(strsplit(spec,':')[[1]])],collapse = ':'), sep = '')
      print(paste('countFile: ', out_specs_list$CountFileLoc, sep = ''))
    }
    if(spec_id == 'sgRNAInfoFileLoc'){
      out_specs_list$sgRNAInfoFileLoc <- paste(data.loc, paste(strsplit(spec,':')[[1]][2:length(strsplit(spec,':')[[1]])],collapse = ':'), sep = '')
    }
    if(spec_id == 'sgTargetPosFileLoc'){
      out_specs_list$sgTargetPosFileLoc <- paste(data.loc, paste(strsplit(spec,':')[[1]][2:length(strsplit(spec,':')[[1]])],collapse = ':'), sep = '')
    }
    if(spec_id == 'Method'){
      out_specs_list$Method <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'foldChangePaired'){
      out_specs_list$foldChangePaired <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'Group1'){
      out_specs_list$Group1 <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'Group2'){
      out_specs_list$Group2 <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'Group1Comb'){
      out_specs_list$Group1Comb <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'Group2Comb'){
      out_specs_list$Group2Comb <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'unsortedGroups'){
      out_specs_list$unsortedGroups <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'dataArranging'){
      out_specs_list$dataArranging <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'dataFilter'){
      out_specs_list$dataFilter <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'CPMs'){
      out_specs_list$dataFilter <- 'CPMs'
      out_specs_list$cpm <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'dataName'){
      out_specs_list$dataName <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'MAGeCKrra'){
      out_specs_list$MAGeCKrra <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'CRESTrra'){
      out_specs_list$CRESTrra <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'repl_groups'){
      out_specs_list$repl_groups <- lapply(strsplit(strsplit(strsplit(spec,':')[[1]][2],';')[[1]],','), as.numeric)
      names(out_specs_list$repl_groups) <- paste0('repl_',c(1:length(out_specs_list$repl_groups)))
    }
    if(spec_id == 'multivar_repl_groups'){
      out_specs_list$multivar_repl_groups <- lapply(strsplit(strsplit(strsplit(spec,':')[[1]][2],';')[[1]],','), as.numeric)
      names(out_specs_list$multivar_repl_groups) <- paste0('repl_',c(1:length(out_specs_list$multivar_repl_groups)))
    }
    if(spec_id == 'crisprSystem'){
      out_specs_list$crisprSystem <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'crisprEffectRange'){
      out_specs_list$crisprEffectRange <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'deletionSize'){
      out_specs_list$deletionSize <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'RELICS_genomeScoring'){
      out_specs_list$RELICS_genomeScoring <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'noPostScoreAnalysis'){
      out_specs_list$noPostScoreAnalysis <- strsplit(spec,':')[[1]][2]
    }

    #parameters with default settings:
    if(spec_id == 'saveAnalysisCounts'){
      out_specs_list$saveAnalysisCounts <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'saveAnalysisInfo'){
      out_specs_list$saveAnalysisInfo <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'saveAnalysisRowsRemoved'){
      out_specs_list$saveAnalysisRowsRemoved <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'InitialDispersion'){
      out_specs_list$InitialDispersion <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'InitialEta'){
      out_specs_list$InitialEta <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'ZINBminThresh'){
      out_specs_list$ZINBminThresh <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'binSize'){
      out_specs_list$binSize <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'preBinning'){
      out_specs_list$preBinning <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'preOverlapBinning'){
      out_specs_list$preOverlapBinning <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'binBy'){
      out_specs_list$binBy <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'labelFile'){
      out_specs_list$labelFile <- paste(data.loc, strsplit(spec,':')[[1]][2], sep = '')
    }
    if(spec_id == 'labelHierarchy'){
      out_specs_list$labelHierarchy <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'postScoringAnalysis'){
      out_specs_list$postScoringAnalysis <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'postScoringCREST'){
      out_specs_list$postScoringCREST <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'postScoreRRA'){
      out_specs_list$postScoreRRA <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'postScoreAlphaRRA'){
      out_specs_list$postScoreAlphaRRA <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'postScoreSlidingWindow'){
      out_specs_list$postScoreSlidingWindow <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'guidePerSlidingWindow'){
      out_specs_list$guidePerSlidingWindow <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'maxWindowSize'){
      out_specs_list$maxWindowSize <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'postScoringMAGeCK'){
      out_specs_list$postScoringMAGeCK <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'evaluate_perGuide_Performance'){
      out_specs_list$evaluate_perGuide_Performance <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'evaluate_perRegion_Performance'){
      out_specs_list$evaluate_perRegion_Performance <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'evaluate_perElement_Performance'){
      out_specs_list$evaluate_perElement_Performance <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'pos_regions'){
      out_specs_list$pos_regions <- paste(data.loc, paste(strsplit(spec,':')[[1]][2:length(strsplit(spec,':')[[1]])],collapse = ':'), sep = '')
    }
    if(spec_id == 'pos_region_size'){
      out_specs_list$pos_region_size <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'saveGuideScores'){
      out_specs_list$saveGuideScores <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'save_dirichlet_info'){
      out_specs_list$save_dirichlet_info <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'save_glmmUntrained_info'){
      out_specs_list$save_glmmUntrained_info <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'positiveLabels'){
      out_specs_list$positiveLabels <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'negativeLabels'){
      out_specs_list$negativeLabels <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'genePosFileLoc'){
      out_specs_list$plotGenes <- 'yes'
      plotGenesDataLoc <- paste(data.loc, paste(strsplit(spec,':')[[1]][2:length(strsplit(spec,':')[[1]])],collapse = ':'), sep = '')
      out_specs_list$plotGenesData <- read.csv(plotGenesDataLoc,as.is = TRUE,header=TRUE)
    }
    if(spec_id == 'zoomRange'){
      out_specs_list$zoomRange <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
      out_specs_list$zoomChromosomePlot <- 'yes'
    }
    if(spec_id == 'tileZoomRange'){
      out_specs_list$tileZoomRange <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
      out_specs_list$tileZoomChromosomePlot <- 'yes'
    }
    if(spec_id == 'plotAllChrom'){
      out_specs_list$plotAllChrom <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'plotChroms'){
      out_specs_list$plotChroms <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'plotSeparateChroms'){
      out_specs_list$plotSeparateChroms <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'RRAmaxPercent'){
      out_specs_list$RRAmaxPercent <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'CRISPReffectRange'){
      out_specs_list$CRISPReffectRange <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'rraPermutNr'){
      out_specs_list$rraPermutNr <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'rraNaBin'){
      out_specs_list$rraNaBin <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'scoreThresh'){
      out_specs_list$scoreThresh <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'createBedGraph'){
      out_specs_list$createBedGraph <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'saveNB_LLR'){
      out_specs_list$saveNB_LLR <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'postScoreGuideWindowRRA'){
      out_specs_list$postScoreGuideWindowRRA <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'guidePerFixedWindow'){
      out_specs_list$guidePerFixedWindow <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'crestEdgeRfilter'){
      out_specs_list$crestEdgeRfilter <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'simulated_data'){
      out_specs_list$simulated_data <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'glmm_positiveTraining'){
      out_specs_list$glmm_positiveTraining <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'glmm_negativeTraining'){
      out_specs_list$glmm_negativeTraining <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'glmm_untrainedNull'){
      out_specs_list$glmm_untrainedNull <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'dirichlet_null'){
      out_specs_list$dirichlet_null <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'glmm_allData_positiveTraining'){
      out_specs_list$glmm_allData_positiveTraining <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'glmm_allData_negativeTraining'){
      out_specs_list$glmm_allData_negativeTraining <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'glmm_reducedData_positiveTraining'){
      out_specs_list$glmm_reducedData_positiveTraining <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'glmm_reducedData_negativeTraining'){
      out_specs_list$glmm_reducedData_negativeTraining <- strsplit(strsplit(spec,':')[[1]][2],',')[[1]]
    }
    if(spec_id == 'glmm_neg_trainingFraction'){
      out_specs_list$glmm_neg_trainingFraction <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'searchSilencers'){
      out_specs_list$searchSilencers <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'hmmStates'){
      out_specs_list$hmmStates <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'transition_prob'){
      out_specs_list$transition_prob <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'hmmStartStateProbs'){
      out_specs_list$hmmStartStateProbs <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'transition_type'){
      out_specs_list$transition_type <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'dataDir'){
      out_specs_list$dataDir <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'simName'){
      out_specs_list$simName <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'scorePlotting'){
      out_specs_list$scorePlotting <- strsplit(spec,':')[[1]][2]
    }
    if(spec_id == 'nrSims'){
      out_specs_list$nrSims <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'categoryProbs'){
      out_specs_list$categoryProbs <- as.numeric(strsplit(strsplit(spec,':')[[1]][2],',')[[1]])
    }
    if(spec_id == 'nrFreeCores'){
      out_specs_list$nrFreeCores <- as.numeric(strsplit(spec,':')[[1]][2])
    }
    if(spec_id == 'parallelOutfile'){
      out_specs_list$parallelOutfile <- strsplit(spec,':')[[1]][2]
    }
  }
  return(out_specs_list)
}

# given the specifications, calculate the scores for the specified methods
# input: analysis specification file, label file (optional)
#   for further details on input files, see analyze_data() documentation
# output: list containing dataframes for each analysis method, each dataframe having the following columns:
#   $method:
#     rawScores: calculated probabilities of p-values
#     formatScores: adjusted scores such that they are signed by ratio and can be ranked from most to least significant
#     log2_rate_ratio: log2 ratio used for sign-adjusting scores, either based on count ratio or calculated rates
#     chromosome, labels, sgRNA targets (can be 1 or 2 columns) (initial values from the info file)
calculate_scores <- function(analysis.specs, label.file = NULL){

  analysis.specs.names <- names(analysis.specs)
  count.data <- read.csv(analysis.specs$CountFileLoc,as.is = TRUE)
  data.info <- read.csv(analysis.specs$sgRNAInfoFileLoc,as.is = TRUE, stringsAsFactors = FALSE)

  #convert NA in chromosome column to string 'NA'
  if(length(which(is.na(data.info[,1]))) > 0){
    data.info[which(is.na(data.info[,1])),1] <- 'NA'
  }
  #if a label file is supplied, add corresponding note to the specification file
  if(! is.null(label.file)){
    analysis.specs$labelFile <- label.file
  }

  #extract samples of interest
  # output: list
  #   $counts, contains counts of relevant samples
  #   $info, contains chromsome and positional information
  #   $specs, contains specs and added $subGroup1 and $subGroup2
  sample.list <- extract_samples(count.data,data.info,analysis.specs)
  sample.counts <- sample.list$counts
  sample.info <- sample.list$info
  sample.specs <- sample.list$specs

  # re-arrange data if specified
  # output: list
  #   $counts, contains re-arranged counts
  #   $info, contains chromsome and positional information about re-arranged data
  #   $specs, contains specs about re-arranged data
  arranged.data.list <- reArrange_data(sample.counts, sample.info, sample.specs)
  arranged.data.counts <- arranged.data.list$counts
  arranged.data.info <- arranged.data.list$info
  arranged.data.specs <- arranged.data.list$specs

  #filter data if specified
  filtered.data.list <- filtering_counts(arranged.data.counts, arranged.data.info, arranged.data.specs)
  filtered.data.counts <- filtered.data.list$counts
  filtered.data.info <- filtered.data.list$info
  filtered.data.specs <- filtered.data.list$specs

  #added this line
  filtered.data.specs$new_grp1 <- filtered.data.specs$subGroup1
  filtered.data.specs$new_grp2 <- filtered.data.specs$subGroup2

  # if specified; counts, info and rows removed will be saved for reproducability
  save_analysisData(filtered.data.counts, filtered.data.info, filtered.data.specs)

  #output scores from specified methods
  method.score.list <- score_calculation(filtered.data.counts, filtered.data.info, filtered.data.specs)

  print('Guide scores obtained ...')
  if('saveGuideScores' %in% names(filtered.data.specs)){
    print('Saving guide scores ...')
    save_all_guide_scores(method.score.list, filtered.data.specs)
  }
  #apply additional methods on the output scores for high-level scores
  post.method.score.list <- post_score_calculation(method.score.list, filtered.data.specs)

  return(post.method.score.list)
}

# extract the samples of interest
# input: count df, info df, specification list
# output: list
#   $counts, contains counts of relevant samples
#   $info, contains chromsome and positional information
#   $specs, contains specs and added $subGroup1, $subGroup2 and subUnsorted (latter if Model_X is used for analysis)
extract_samples <- function(input.counts, input.info, input.specs){

  columns.to.extract <- c()
  # columns.to.extract <- c(input.specs$Group1, input.specs$Group2)
  if('Group1Comb' %in% names(input.specs)){
    input.specs$Group1 <- vector(mode = 'numeric', length = length(input.specs$Group1Comb))
    input.specs$Group2 <- vector(mode = 'numeric', length = length(input.specs$Group2Comb))
    for(i in 1:length(input.specs$Group1Comb)){
      if(grepl('+', input.specs$Group1Comb[i], fixed = T)){
        temp.split <- strsplit(input.specs$Group1Comb[i],'+')[[1]]
        extract.cols <- seq(1,length(temp.split),2)
        temp.cols <- as.numeric(temp.split[extract.cols])
        col.sum <- rowSums(input.counts[, temp.cols])
        added.col.index <- ncol(input.counts) + 1
        input.specs$Group1[i] <- added.col.index
        input.counts[[paste0(temp.cols, collapse = '_')]] <- col.sum
      } else {
        input.specs$Group1[i] <- as.numeric(input.specs$Group1Comb[i])
      }
    }
    for(i in 1:length(input.specs$Group2Comb)){
      if(grepl('+', input.specs$Group2Comb[i], fixed = T)){
        temp.split <- strsplit(input.specs$Group2Comb[i],'+')[[1]]
        extract.cols <- seq(1,length(temp.split),2)
        temp.cols <- as.numeric(temp.split[extract.cols])
        col.sum <- rowSums(input.counts[, temp.cols])
        added.col.index <- ncol(input.counts) + 1
        input.specs$Group2[i] <- added.col.index
        input.counts[[paste0(temp.cols, collapse = '_')]] <- col.sum
      } else {
        input.specs$Group2[i] <- as.numeric(input.specs$Group2Comb[i])
      }
    }
  }

  columns.to.extract <- c(input.specs$Group1, input.specs$Group2)

  # if 'repl_groups' are specified, add them to the columns which have to be extracted
  if('repl_groups' %in% names(input.specs)){
    columns.to.extract <- unique(c(columns.to.extract, as.numeric(unlist(input.specs$repl_groups))))
  }

  if('multivar_repl_groups' %in% names(input.specs)){
    columns.to.extract <- unique(c(columns.to.extract, as.numeric(unlist(input.specs$multivar_repl_groups))))
  }

  # find the mapping to the above...

  if(max(columns.to.extract) > ncol(input.counts)) stop('Error in "extract_samples": gropus are referencing non-existent columns')
  output.list <- list()
  output.list$counts <- input.counts[,columns.to.extract]
  final.info <- input.info
  if(! all(c('chrom','label','start') %in% colnames(final.info) == TRUE)){
    print("Missing column labels; either 'chrom', 'label' or 'start' ")
  }
  # colnames(final.info)[1:4] <- c('chrom','label','start','end')
  output.list$info <- final.info
  updated.specs <- input.specs

  group1.pos <- which(columns.to.extract %in% input.specs$Group1)
  group2.pos <- which(columns.to.extract %in% input.specs$Group2)

  # adjust the referenced columns from the replicate gropus
  if('repl_groups' %in% names(input.specs)) {
    # need to preserve the ordering to have meaningful intercepts in GLMM
    repl.grp.pos <- lapply(input.specs$repl_groups, function(x){
      out.order <- c()
      for(i in 1:length(x)){
        out.order <- c(out.order, which(columns.to.extract == x[i]))
      }
      out.order
      })
    updated.specs$repl_groups <- repl.grp.pos
  }

  if('multivar_repl_groups' %in% names(input.specs)) {
    # need to preserve the ordering to have meaningful intercepts in GLMM
    mv.repl.grp.pos <- lapply(input.specs$multivar_repl_groups, function(x){
      mv.out.order <- c()
      for(i in 1:length(x)){
        mv.out.order <- c(mv.out.order, which(columns.to.extract == x[i]))
      }
      mv.out.order
      })
    updated.specs$multivar_repl_groups <- mv.repl.grp.pos
  }


  updated.specs$subGroup1 <- group1.pos # c(1:length(input.specs$Group1))
  updated.specs$subGroup2 <- group2.pos # c(1:length(input.specs$Group2)) + length(input.specs$Group1)

  output.list$specs <- updated.specs
  return(output.list)
}

# re-arrange data
# input: count df, info df, specification list
# output: list
#   $counts, contains re-arranged counts
#   $info, contains chromsome and positional information about re-arranged data
#   $specs, contains specs about re-arranged data
reArrange_data <- function(input.counts, input.info, input.specs){
  output.list <- list()
  if(input.specs$dataArranging == 'none'){
    print('No data rearrangements were made.')
    output.list$counts <- input.counts
    output.list$info <- input.info
    output.list$specs <- input.specs
    return(output.list)
  }
  if(input.specs$dataArranging == 'yes'){
    if('preBinning' %in% names(input.specs)){
      binning.list <- per_sample_binning(input.counts, input.info, input.specs, 'preBinning')
      output.list$counts <- binning.list$counts
      output.list$info <- binning.list$info
      output.list$specs <- input.specs
      return(output.list)
    }
    if('preOverlapBinning' %in% names(input.specs)){
      binning.list <- per_sample_binning(input.counts, input.info, input.specs, 'preOverlapBinning')
      output.list$counts <- binning.list$counts
      output.list$info <- binning.list$info
      output.list$specs <- input.specs
      return(output.list)
    }
  }
}

# filter counts
# input: count df, info df, specification list
# output: list
#   $counts, contains re-arranged counts
#   $info, contains chromsome and positional information about re-arranged data
#   $specs, contains specs about re-arranged data
#   $rowsRemoved, rows which were removed
filtering_counts <- function(input.counts, input.info, input.specs){
  out.list <- list()

  #if no filtering is applied
  if(input.specs$dataFilter == 'none'){
    out.list$counts <- input.counts
    out.list$info <- input.info
    out.list$specs <- input.specs
    return(out.list)
  }

  #if zero rows should be removed
  if(input.specs$dataFilter == 'zeroRows'){
    rows.to.remove <- which(rowSums(input.counts) == 0)
    if(length(rows.to.remove) > 0){
      out.list$counts <- input.counts[-rows.to.remove,]
      out.list$info <- input.info[-rows.to.remove,]
      out.list$specs <- input.specs
    } else {
      out.list$counts <- input.counts
      out.list$info <- input.info
      out.list$specs <- input.specs
    }
    return(out.list)
  }

  #if CPMs below indicated shoul dbe removed
  if(input.specs$dataFilter == 'CPMs'){
    cpm.cutoff <- input.specs$cpm
    rows.to.remove <- which(rowSums(cpm(input.counts)) < cpm.cutoff)
    if(length(rows.to.remove) > 0){
      print(paste('Removed ',length(rows.to.remove),' with CPM less than ', cpm.cutoff, sep = ''))
      out.list$counts <- input.counts[-rows.to.remove,]
      out.list$info <- input.info[-rows.to.remove,]
      out.list$specs <- input.specs
    } else {
      out.list$counts <- input.counts
      out.list$info <- input.info
      out.list$specs <- input.specs
    }
    return(out.list)
  }

}

# saving data to be analyzed
# input: count df, info df, specification list
# output: csv tables is specified
#   counts to be analyzed
#   info to be analyzed
#   rows removed in the filtering process
save_analysisData <- function(input.counts, input.info, input.specs){
  if('saveAnalysisCounts' %in% names(input.specs)){
    write.csv(input.counts, file = 'filtered_analysis_counts.csv',row.names = FALSE)
  }
  if('saveAnalysisInfo' %in% names(input.specs)){
    write.csv(input.info, file = 'filtered_analysis_info.csv',row.names = FALSE)
  }
  if('saveAnalysisRowsRemoved' %in% names(input.specs)){
    write.csv(input.specs, file = 'filtered_analysis_rowsRemoved.csv',row.names = FALSE)
  }
}

# obtain the scores for all specified analysis methods
# input: count df, info df, specification list
# output: list, each method adds the following to the
# reduced version: does not calculate anything else besides RELCIS scores
#   $method, where each $mthod contains the following rows:
#     rawScores: calculated probabilities of p-values
#     formatScores: adjusted scores such that they are signed by ratio and can be ranked from most to least significant
#     log2_rate_ratio: log2 ratio used for sign-adjusting scores, either based on count ratio or calculated rates
#     chromosome, labels, sgRNA targets (can be 1 or 2 columns) (all values from the info file, can be altered in methods like CREST analysis)
score_calculation <- function(input.counts, input.info, input.specs){
  out.list <- list()
  if('RELICS-search' %in% input.specs$Method){
    print('#####################################')
    print('Running RELICS-search')
    out.list$RELICS <- RELICS_search(input.counts, input.info, input.specs)
    print('Finished RELICS-search')
  }
  return(out.list)
}

# RELICS-search: using a GLMM a Bayes factor is calculated for each guide
# input: count df, info df, specification list
# output: df, each method adds the following to the
#   rawScores, formatScores, log2_rate_ratio, chromosome, labels, sgRNA targets (can be 1 or 2 columns)
RELICS_search <- function(input.counts, input.info, input.specs){

  out.score.list <- list()

  glmm.info <- input.info
  glmm.counts <- input.counts

  pos.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_positiveTraining), ]
  neg.counts <- glmm.counts[which(glmm.info$label %in% input.specs$glmm_negativeTraining), ]

  nr.pos.crossv <- round(nrow(pos.counts) / 3)
  nr.neg.crossv <- round(nrow(neg.counts) / 3)

  par.list <- list()

  print('Establishing RELICS parameters...')
  set.seed(13)  # for reproducibility of the sampling
  for(i in 1:length(input.specs$repl_groups)){
    repl.pars <- c() # all fold estimates of parameters
    # estimate parameters, AUC and prAUC from all the data
    temp.label.counts <- rbind(pos.counts[, input.specs$repl_groups[[i]]],
      neg.counts[, input.specs$repl_groups[[i]]])
    guide.type <- c(rep('pos', nrow(pos.counts)), rep('neg', nrow(neg.counts)))
    temp.par.list <- RELICS_parameterEstimation(temp.label.counts, guide.type, input.specs)
    temp.sample.bayes.factor <- RELICS_likelihood_calculation(temp.label.counts, list(temp.par.list$out_list),
      list(repl_groups = list(1:length(input.specs$repl_groups[[i]]))) )
    temp.bayes.factor <- rowSums(temp.sample.bayes.factor)
    temp.auc.out <- AUC_x_y(temp.bayes.factor, guide.type, 'pos', 'neg')
    repl.pars <- rbind(repl.pars, c(temp.par.list$out_df[1,], temp.auc.out$num_AUC))

    for(j in 1:10){
      temp.label.counts <- rbind(pos.counts[sample(1:nrow(pos.counts), nr.pos.crossv),
        input.specs$repl_groups[[i]]], neg.counts[sample(1:nrow(neg.counts), nr.neg.crossv),
        input.specs$repl_groups[[i]]])
      guide.type <- c(rep('pos', nr.pos.crossv), rep('neg', nr.neg.crossv))
      temp.par.list <- RELICS_parameterEstimation(temp.label.counts, guide.type, input.specs)
      temp.sample.bayes.factor <- RELICS_likelihood_calculation(temp.label.counts, list(temp.par.list$out_list),
        list(repl_groups = list(1:length(input.specs$repl_groups[[i]]))) )
      temp.bayes.factor <- rowSums(temp.sample.bayes.factor)
      temp.auc.out <- AUC_x_y(temp.bayes.factor, guide.type, 'pos', 'neg')
      repl.pars <- rbind(repl.pars, c(temp.par.list$out_df[1,], temp.auc.out$num_AUC))
    }

    colnames(repl.pars) <- c(names(temp.par.list$out_df[1,]), 'AUC')
    repl.pars <- as.data.frame(repl.pars, stringsAsFactors = F)
    # check if parameters can distinguish between positive and negatie controls
    temp.mean.auc <- mean(repl.pars$AUC[2:11])
    temp.sd.auc <- sd(repl.pars$AUC[2:11])
    if((temp.mean.auc - 2*temp.sd.auc) < 0.5){
      print('Warning: parameter estimates can not easily distinguish between positive and negative controls')
    }

    # nr betas is equivalen to number of pools - 1 due to intercept
    par.list[[i]] <- RELICS_extract_median_pars(repl.pars[c(2:11),], (length(input.specs$repl_groups[[i]]) - 1))
    write.csv(repl.pars, file = paste0(input.specs$dataName, 'RELICS_repl', i, '_parEst.csv'), row.names = F)
  }

  guide.sample.bayes.factor <- RELICS_likelihood_calculation(glmm.counts, par.list, input.specs)  #total.guide.pos.ll - total.guide.neg.ll
  guide.bayes.factor <- rowSums(guide.sample.bayes.factor)

  par.summary <- RELICS_summary(par.list)
  write.csv(par.summary, file = paste0(input.specs$dataName,'_RELICS_parSummary.csv') , row.names = F)

  out.df <- cbind.data.frame(rawScores = guide.bayes.factor, formatScores = guide.bayes.factor,
    log2_rate_ratio = rep(1, length(guide.bayes.factor)), stringsAsFactors = FALSE)

  for(i in 1:ncol(guide.sample.bayes.factor)){
    out.df[[paste0('repl', i, '_bf')]] <- guide.sample.bayes.factor[,i]
  }
  #add the info part
  out.df <- cbind(out.df, glmm.info, stringsAsFactors = FALSE)

  return(out.df)

}

# this estimates the parameters for the positive and negative controls
# input.type describes if the guide is from a positive pool or not
# input.counts contain first all positive, then all negative guides!
RELICS_parameterEstimation <- function(input.counts, input.type, input.specs){

  pool.IDs <- colnames(input.counts)
  guide.nrs <- c(1:nrow(input.counts))
  guide.names <- paste0('guide', guide.nrs)
  pool.nrs <- c(1:ncol(input.counts))
  pool.names <- paste0('pool', pool.nrs)
  positive.types <- c(1:(ncol(input.counts) - 1))

  #change to : intercept, pos1, pos
  pos.type.label <- c('interc', paste0('pos', positive.types))
  nr.positives <- length(which(input.type == 'pos'))

  neg.type.label <- c('interc', paste0('neg', positive.types))
  nr.negatives <- length(which(input.type == 'neg'))

  guide.labels <- rep(guide.names, each = length(pool.names))
  pool.labels <- rep(pool.names, length(guide.names))
  type.labels <- c(rep(pos.type.label, nr.positives), rep(neg.type.label, nr.negatives))

  long.df <- cbind.data.frame(observed = c(t(input.counts)), guide_labels = guide.labels,
    pool_labels = pool.labels, type_labels = type.labels, stringsAsFactors = T)

  tmb.nb2.model <- glmmTMB(observed ~ type_labels + (1|guide_labels),
    data = long.df, family=nbinom2)

  model.dispersion <- sigma(tmb.nb2.model) # residual variance, theta corresponding to the NB variance
  model.tau <- exp(getME(tmb.nb2.model, name = 'theta'))  # SD of the random effect
  #  no idea why this works, no idea how else to get random effect SD

  nr.beta.neg <- ncol(input.counts) - 1
  nr.beta.pos <- ncol(input.counts) - 1
  nr.total.beta <- nr.beta.neg + nr.beta.pos
  tmb.betas.orig <- tmb.nb2.model$fit$par[2:(nr.total.beta + 1)]
  beta.neg <- tmb.betas.orig[1:nr.beta.neg] #+ c(0, rep(tmb.betas.orig[1], nr.beta.neg - 1))
  beta.pos <- tmb.betas.orig[(nr.beta.neg + 1):nr.total.beta] #beta.neg + c(0, tmb.betas.orig[(nr.beta.neg + 1):nr.total.beta])

  if('searchSilencers' %in% names(input.specs)){
    beta.pos <- beta.neg - (beta.pos - beta.neg)  # inverse the effect of the positives
  }

  tmb.neg_logLik <- logLik(tmb.nb2.model)[1]

  out.parameters <- c(tmb.nb2.model$fit$par[1], beta.neg, beta.pos, model.tau, model.dispersion)
  out.df <- t(as.data.frame(out.parameters))
  colnames(out.df) <- c('beta_intercept', paste('beta_neg_pred', c(1:length(beta.neg)), sep = '_'),
    paste('beta_pos_pred', c(1:length(beta.neg)), sep = '_'), 'tau', 'disp')

  out.list <- list(betas = tmb.betas.orig[2:(nr.beta.neg + 1)], beta_pos_pred = beta.pos,
    beta_neg_pred = beta.neg, model_tau = model.tau, model_dispersion = model.dispersion,
    neg_logLik = tmb.neg_logLik, pool_IDs = pool.IDs, beta_intercept = tmb.nb2.model$fit$par[1])

  final.list <- list(out_df = out.df, out_list = out.list)

  return(final.list)

}

# calculate hte likelihood of observing the data counts for each guide
# given the positive and negative control parameters
RELICS_likelihood_calculation <- function(input.counts, input.par.list, input.specs){
  guide.pos.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(input.counts)))
  guide.neg.ll <- data.frame(matrix(ncol = length(input.specs$repl_groups), nrow = nrow(input.counts)))
  for(i in 1:nrow(input.counts)){
    for(j in 1:length(input.specs$repl_groups)){
      temp.guide.counts <- unlist(input.counts[i,input.specs$repl_groups[[j]]])
      temp_pos_integration_function <- RELICS_function_generator(temp.guide.counts,
        c(0, input.par.list[[j]]$beta_pos_pred), input.par.list[[j]]$model_tau,
        input.par.list[[j]]$model_dispersion, input.par.list[[j]]$beta_intercept)
      temp_neg_integration_function <- RELICS_function_generator(temp.guide.counts,
        c(0, input.par.list[[j]]$beta_neg_pred), input.par.list[[j]]$model_tau,
        input.par.list[[j]]$model_dispersion, input.par.list[[j]]$beta_intercept)
      guide.pos.ll[i,j] <- log(integrate(temp_pos_integration_function,
        lower = (-5*input.par.list[[j]]$model_tau + input.par.list[[j]]$beta_intercept),
        upper = (5*input.par.list[[j]]$model_tau + input.par.list[[j]]$beta_intercept) )$value[1])
      guide.neg.ll[i,j] <- log(integrate(temp_neg_integration_function,
        lower = (-5*input.par.list[[j]]$model_tau + input.par.list[[j]]$beta_intercept),
        upper = (5*input.par.list[[j]]$model_tau + input.par.list[[j]]$beta_intercept) )$value[1])
    }
  }

  adjusted.guide.pos.ll <- guide.pos.ll
  adjusted.guide.neg.ll <- guide.neg.ll

  # check below when debugging: how many fall in this category?
  for(i in 1:ncol(adjusted.guide.pos.ll)){
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == -Inf),i] <- min(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] > -Inf), i])
    adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] == Inf),i] <- max(adjusted.guide.pos.ll[which(adjusted.guide.pos.ll[,i] < Inf), i])

    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == -Inf),i] <- min(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] > -Inf), i])
    adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] == Inf),i] <- max(adjusted.guide.neg.ll[which(adjusted.guide.neg.ll[,i] < Inf), i])
  }

  sample.bayes.factor <- adjusted.guide.pos.ll - adjusted.guide.neg.ll

  return(sample.bayes.factor)

}

# generic function generator; given input parameters a function for liikelihood calculation is created
RELICS_function_generator <- function(input.obs, input.beta, input.tau,
  input.dispersion, input.mean){
    # x is ranging around the mean
  function(x){
    likelihood <- dnorm(x, mean = input.mean, sd = input.tau)
    for(i in 1:length(input.obs)){
      likelihood <- likelihood * dnbinom(input.obs[i], mu = exp(input.beta[i] + x),
        size = input.dispersion)
    }
    likelihood
  }
}

# calculating area under the curve using regular AUC, looking at %negative control vs %positive control included
# input:
#   formatted.pvals: formatted p-values as input. Largest one is treated as most significant one
#   in_label: vector of equal length as formatted.pvals, assigns leable with each p-value
#   pos_labels: label identifier for positive labels (ex.: 'pos', 'exon', ..)
#   neg_labels: label identifier for negative labels (ex.: 'neg', 'nonTargeting', ..)
#output: list containing: xAxis, yAxis, AUC (already numeric)
AUC_x_y <- function(formatted.pvals,in_label,pos_labels,neg_labels){

  comb_df <- cbind.data.frame(-formatted.pvals,in_label,stringsAsFactors = FALSE)
  sorted_comb_df <- comb_df[order(comb_df[,1]),]

  #x-value represents the number of false positives incorporated
  x_axis <- 0
  #y-value represents the % of total true positives found
  y_axis <- 0

  total_nr_pos <- length(which(in_label %in% pos_labels))  #count_elemnt_occurance(pos_labels,in_label)
  total_nr_neg <- length(which(in_label %in% neg_labels))  #$count_elemnt_occurance(neg_labels,in_label)

  #not all positives have been found
  all_pos_NOT_found <- TRUE
  posit_tracker <- 1
  nr_pos_found <- 0
  nr_neg_found <- 0

  #category is either true positive (1) or false positive (0)
  auc_category <- c()

  while(all_pos_NOT_found){

    temp_row <- sorted_comb_df[posit_tracker,]

    if(temp_row[2] %in% pos_labels){
      nr_pos_found <- nr_pos_found + 1
      x_axis <- c(x_axis,nr_neg_found)
      y_axis <- c(y_axis,nr_pos_found)
      auc_category <- c(auc_category, 1)
      if(nr_pos_found == total_nr_pos){
        all_pos_NOT_found <- FALSE
        break
      }
    }
    else if(temp_row[2] %in% neg_labels){
      nr_neg_found <- nr_neg_found + 1
      x_axis <- c(x_axis,nr_neg_found)
      y_axis <- c(y_axis,nr_pos_found)
      auc_category <- c(auc_category, 0)
    }
    posit_tracker <- posit_tracker + 1
  }

  if(nr_neg_found < total_nr_neg){
    x_axis <- c(x_axis,c((nr_neg_found+1):total_nr_neg))
    y_axis <- c(y_axis,rep(total_nr_pos,length(c((nr_neg_found+1):total_nr_neg))))
  }

  #following example on r bloggers
  add_negatives <- rep(0,total_nr_neg-nr_neg_found)
  auc_category <- c(auc_category,add_negatives)
  auc_prediction <- rev(seq_along(auc_category))
  roc_obj <- roc(auc_category, auc_prediction)

  return(list(xAxis = (x_axis/total_nr_neg), yAxis = (y_axis/total_nr_pos), num_AUC = round(pROC::auc(roc_obj),3)))
}

# extract the median parameter from the 10-fold corssvalidation
RELICS_extract_median_pars <- function(input.df, nr.betas){
  out.list <- list()
  out.list$beta_intercept <- median(input.df$beta_intercept)

  temp.beta.neg <- c()
  for(i in 2:(nr.betas + 1)){
    temp.beta.neg <- c(temp.beta.neg, median(input.df[,i]))
  }
  out.list$beta_neg_pred <- temp.beta.neg

  temp.beta.pos <- c()
  for(i in (nr.betas + 2):(nr.betas*2 + 1)){
    temp.beta.pos <- c(temp.beta.pos, median(input.df[,i]))
  }
  out.list$beta_pos_pred <- temp.beta.pos

  out.list$model_tau <- median(input.df[,(nr.betas*2 + 2)])
  out.list$model_dispersion <- median(input.df[,(nr.betas*2 + 3)])

  return(out.list)
}

# summarize RELICS parameters and write to output file
# 'repl_ID', 'pool_name', 'neg_beta', 'pos_beta', 'beta_diff', 'tau', 'theta'
RELICS_summary <- function(input.list){

  nr.betas <- 0
  pool.names <- c()
  for(i in 1:length(input.list)){
    nr.betas <- nr.betas + length(input.list[[i]]$beta_neg_pred) + 1
    pool.names <- c(pool.names, input.list[[i]]$pool_IDs)
  }

  out.df <- data.frame(matrix(ncol = 7, nrow = nr.betas ))
  colnames(out.df) <- c('repl_ID', 'pool_name', 'neg_beta', 'pos_beta', 'beta_diff', 'tau', 'theta')

  temp.pos.betas <- c()
  temp.neg.betas <- c()
  temp.tau <- c()
  temp.theta <- c()
  temp.rownames <- c()
  for(i in 1:length(input.list)){
    temp.name <- paste0('repl1_',i)
    temp.rownames <- c(temp.rownames, 'intercept', paste0(temp.name,'_beta_', c(1:length(input.list[[i]]$beta_neg_pred) )))
    temp.neg.betas <- c(temp.neg.betas, input.list[[i]]$beta_intercept, input.list[[i]]$beta_neg_pred)
    temp.pos.betas <- c(temp.pos.betas, input.list[[i]]$beta_intercept, input.list[[i]]$beta_pos_pred)
    temp.tau <- c(temp.tau, rep(input.list[[i]]$model_tau, length(input.list[[i]]$beta_neg_pred) + 1))
    temp.theta <- c(temp.theta, rep(input.list[[i]]$model_dispersion, length(input.list[[i]]$beta_neg_pred) + 1))
  }
  beta.diff <- temp.pos.betas - temp.neg.betas

  out.df$repl_ID <- temp.rownames
  out.df$pool_name <- pool.names
  out.df$neg_beta <- round(temp.neg.betas,3)
  out.df$pos_beta <- round(temp.pos.betas,3)
  out.df$beta_diff <- round(beta.diff,3)
  out.df$tau <- round(temp.tau,3)
  out.df$theta <- round(temp.theta,3)
  return(out.df)
}

#apply sliding window to all score datasets
#input:
#   input.list: all method-specifi results to be analyzed using a specific RRA
#   input.specs: specification list
#output: list
#   rawScores, formatScores, log2_rate_ratio, chromosome, label, sgRNA targets (can be 1 or 2 columns)
RELICS_genomeScoring_wrapper <- function(input.list, input.specs){
  out.list <- list()
  list.names <- names(input.list)

  for(i in 1:length(input.list)){
    if(! list.names[i] %in% c('viterbi', 'frwdBkwd', 'dirichletMultinom', 'nbGlmmUntr')){
      temp.list.name <- paste(list.names[i], 'genomeScores', sep = '_')
      out.list[[temp.list.name]] <- RELICS_genomeScoring(input.list[[i]], input.specs, list.names[i])
    } else {
      temp.list.name <- list.names[i]
      out.list[[temp.list.name]] <- input.list[[i]]
    }
  }

  return(out.list)
}

# slide window along scores and claculate sum or average, depending on type of analysis
# input.window, input.ll.df, input.label.hierarchy, input.gene, input.pdf.name
RELICS_genomeScoring <- function(input.df, input.specs, analysis.name){

  # generate a window of effect for each guide
  if(input.specs$crisprSystem %in% c('CRISPRcas9', 'CRISPRi', 'CRISPRa')){
    input.df$start <- input.df$start - round(input.specs$crisprEffectRange / 2)
    input.df$end <- input.df$end + round(input.specs$crisprEffectRange / 2)
  }

  # need to sort the data according to start of guide
  targeting.window.scores <- c()

  #separate out chromsome from non-chromosome guides
  targeting.df <- c()  # data frame containing only valid chromosomes
  avg.overlaps <- c()
  non.targeting.df <- c()
  if('NA' %in% input.df$chrom){
    targeting.df <- filter(input.df, grepl('chr',chrom))
    non.targeting.df <- filter(input.df, grepl('NA',chrom))
  } else {
    targeting.df <- input.df
  }

  scoring.type <- c()

  if(analysis.name %in% c('RELICS','nbGlmm','nbGlmmDev', 'allDataNbGlmmCategoryParallel',
    'allDataNbGlmm', 'numIntPoisGlmm', 'uniPoisGlmm', 'multiPoisGlmm', 'allDataMultiPoisGlmmPerRepl',
    'allDataMultiPoisGlmmAcrossRepl', 'reducedDataMultiPoisGlmmPerRepl', 'reducedDataMultiPoisGlmmAcrossRepl')){
      scoring.type <- 1
    } else if(analysis.name %in% 'FoldChange'){
      scoring.type <- 2
    } else { # DESeq2, edgeR, ... p-value based
      scoring.type <- 3
    }

  if(length(unique(targeting.df$chrom)) == 1){
    sorted.targeting.chrom <- targeting.df[order(targeting.df$start),]
    targeting.window.scores <- score_regions(sorted.targeting.chrom,
      input.specs$labelHierarchy, scoring.type)
  } else {
    all.chroms <- unique(targeting.df$chrom)
    for(i in 1:length(all.chroms)){
      targeting.chrom <- targeting.df[(which(targeting.df$chrom == all.chroms[i])),]
      sorted.targeting.chrom <- targeting.chrom[order(targeting.chrom$start),]
      targeting.window.scores <- score_regions(sorted.targeting.chrom,
        input.specs$labelHierarchy, scoring.type)
      # avg.overlaps <- c(avg.overlaps, temp.score.list$avg_overlaps)
      targeting.window.scores <- rbind(targeting.window.scores, temp.targeting.window.scores)
    }
  }
  if(nrow(targeting.df) != nrow(input.df)){
    avg.overlap <- round(median(targeting.window.scores$nrSupportGuides))
    # avg.overlap <- round(mean(avg.overlaps))
    non.targeting.window.scores <- score_nonTarget_regions(scoring.type,
      non.targeting.df, avg.overlap, non.targeting.df$label[1])
    targeting.window.scores <- rbind(targeting.window.scores, non.targeting.window.scores)
  }

  return(targeting.window.scores)
}

# function to score regions:
# breaks genome into windows of non-changing information and assigns scores
# based on input scores (summation for bayes factors, average for fold change,
#  fisher's method for p-values)
score_regions <- function(input.df, input.label.hierarchy, scoringType){

  all.breaks <- unique(c(input.df$start, input.df$end))
  sorted.breaks <- sort(all.breaks)
  all.region.breaks <- sorted.breaks
  start.regions <- all.region.breaks[c(1:(length(all.region.breaks) - 1))]
  end.regions <- all.region.breaks[c(2:length(all.region.breaks))] - 1

  region.ranges <- GRanges(seqnames = rep(input.df$chrom[1], length(start.regions)),
    ranges = IRanges(start.regions, end.regions))

  # adjust score ranges such that they fall within, not onto, boarders of genomic ranges
  score.ranges <- GRanges(seqnames = input.df$chrom,
    ranges = IRanges(input.df$start+1, input.df$end-1))

  score.overlaps <- as.data.frame(findOverlaps(region.ranges, score.ranges, type = 'any'))

  region.overlap.list <- split(score.overlaps, score.overlaps$queryHits)

  # all.region.scores <- vector('numeric', length = length(start.regions))
  all.region.rawScores <- vector('numeric', length = length(start.regions))
  all.region.formatScores <- vector('numeric', length = length(start.regions))
  all.region.foldChange <- vector('numeric', length = length(start.regions))

  all.region.guide.supports <- vector('numeric', length = length(start.regions))
  unique.regions <- unique(score.overlaps$queryHits)

  region.labels <- unlist(lapply(region.overlap.list, function(x){
    temp.label <- input.df$label[x$subjectHits]
    temp.labels.present <- which(input.label.hierarchy %in% temp.label)
    out.label <- input.label.hierarchy[max(temp.labels.present)]
    return(out.label)
    }))

  avg.region.overlaps <- unlist(lapply(region.overlap.list, function(x){
    temp.overlaps <- length(input.df$formatScores[x$subjectHits])
    return(temp.overlaps)
    }))

  # if the input is based on bayes factor, add the scores, else take the mean
  region.scores <- c()

  if(scoringType == 1) {
    region.scores <- do.call(rbind, lapply(region.overlap.list, function(x){
      temp.scores <- input.df$formatScores[x$subjectHits]
      out.score <- sum(temp.scores)
      out.values <- c(out.score, out.score, rep(1, length(out.score)) )
      return(out.values)
      }))
  } else if(scoringType == 2){ # fold change
    region.scores <- do.call(rbind, lapply(region.overlap.list, function(x){
      temp.df <- input.df[x$subjectHits,]
      out.score <- mean(temp.df$formatScores)
      out.values <- c(out.score, out.score, out.score)
      return(out.values)
      }))
  } else { # p-value based result, use Fisher's method to combine
    region.scores <- do.call(rbind, lapply(region.overlap.list, function(x){
      temp.df <- input.df[x$subjectHits,]
      avg.sign <- mean(temp.df$log2_rate_ratio)
      temp.pval <- pchisq(-2 * sum(log(temp.df$rawScores)), df = 2 * nrow(temp.df), lower.tail=F)
      out.pval <- -log10(temp.pval) * sign(avg.sign)
      out.values <- c(temp.pval, out.pval, avg.sign)
      return(out.values)
      }))
  }

  all.region.rawScores[unique.regions] <- region.scores[,1]
  all.region.formatScores[unique.regions] <- region.scores[,2]
  all.region.foldChange[unique.regions] <- region.scores[,3]

  all.region.labels <- rep('chr', length(all.region.rawScores))
  all.region.labels[unique.regions] <- region.labels
  all.region.guide.supports[unique.regions] <- avg.region.overlaps

  out.df <- data.frame(rawScores = all.region.rawScores, formatScores = all.region.formatScores,
    log2_rate_ratio = all.region.foldChange,
    chrom = rep(input.df$chrom[1], length(all.region.rawScores)),
    label = all.region.labels, start = start.regions, end = end.regions,
    nrSupportGuides = all.region.guide.supports, stringsAsFactors = F)

  return(out.df)
}

# score when conmbining information of non-targeting guides
score_nonTarget_regions <- function(scoringType, input.df, avg.overlap, input.label){

  row.seq <- seq(1, nrow(input.df), by = avg.overlap)
  out.raw.scores <- vector('numeric', length = length(row.seq))
  out.format.scores <- vector('numeric', length = length(row.seq))
  out.log2.fc <- vector('numeric', length = length(row.seq))

  overlapping <- c(0:(avg.overlap - 1))
  for(i in 1:(length(row.seq) - 1)){
    temp.rows <- row.seq[i] + overlapping
    if(scoringType == 1){ # llr based
      out.raw.scores[i] <- sum(input.df$formatScores[temp.rows])
      out.format.scores[i] <- sum(input.df$formatScores[temp.rows])
      out.log2.fc[i] <- 1
    } else if(scoringType == 2){
      out.raw.scores[i] <- mean(input.df$formatScores[temp.rows])
      out.format.scores[i] <- mean(input.df$formatScores[temp.rows])
      out.log2.fc[i] <- mean(input.df$formatScores[temp.rows])
    } else {
      temp.df <- input.df[temp.rows,]
      avg.sign <- mean(temp.df$log2_rate_ratio)
      temp.pval <- pchisq(-2 * sum(log(temp.df$rawScores)), df = 2 * nrow(temp.df), lower.tail=F)

      out.raw.scores[i] <- temp.pval
      out.format.scores[i] <- -log10(temp.pval) * sign(avg.sign)
      out.log2.fc[i] <- avg.sign
    }
  }
  last.index <- length(row.seq)
  temp.rows <- c(row.seq[last.index]:nrow(input.df))
  if(scoringType == 1){ # llr based
    out.raw.scores[last.index] <- sum(input.df$formatScores[temp.rows])
    out.format.scores[last.index] <- sum(input.df$formatScores[temp.rows])
    out.log2.fc[last.index] <- 1
  } else if(scoringType == 2){
    out.raw.scores[last.index] <- mean(input.df$formatScores[temp.rows])
    out.format.scores[last.index] <- mean(input.df$formatScores[temp.rows])
    out.log2.fc[last.index] <- mean(input.df$formatScores[temp.rows])
  } else {
    temp.df <- input.df[temp.rows,]
    avg.sign <- mean(temp.df$log2_rate_ratio)
    temp.pval <- pchisq(-2 * sum(log(temp.df$rawScores)), df = 2 * nrow(temp.df), lower.tail=F)

    out.raw.scores[last.index] <- temp.pval
    out.format.scores[last.index] <- -log10(temp.pval) * sign(avg.sign)
    out.log2.fc[last.index] <- avg.sign
  }

  out.df <- data.frame(rawScores = out.raw.scores, formatScores = out.format.scores,
    log2_rate_ratio = rep(1, length(out.raw.scores)), chrom = rep('NA', length(out.raw.scores)),
    label = rep(input.label, length(out.raw.scores)), start = rep(NA, length(out.raw.scores)),
    end = rep(NA, length(out.raw.scores)), nrSupportGuides = c(rep(avg.overlap,
      (length(out.raw.scores)-1)), length(temp.rows)), stringsAsFactors = F)

  return(out.df)
}

# save te final scores used for plotting, AUC etc.
# input: list with scores per method, input specifications
# output: .csv for each methods scores
save_all_guide_scores <- function(input.score.list, input.specs){
  score.names <- names(input.score.list)
  for(i in 1:length(score.names)){
    temp.score.df <- input.score.list[[i]]
    temp.score.name <- score.names[i]
    method.identifier <- strsplit(temp.score.name, '_')[[1]][1]
    out.csv <- temp.score.df
    out.names <- names(out.csv)
    out.names[which(out.names == 'formatScores')] <- 'guideScore'
    names(out.csv) <- out.names

    write.csv(out.csv, file = paste(input.specs$dataName, method.identifier, 'guideScores.csv', sep = '_'), row.names = FALSE)
  }
}

#for all scores, apply additional post-scoring methods
# obtain the scores for all specified analysis methods
# input: list, containing the following:
#   rawScores, formatScores, log2_rate_ratio, chromosome, label, sgRNA targets (can be 1 or 2 columns)
# output: list, each method adds the following to the
#   rawScores, formatScores, log2_rate_ratio, chromosome, label, sgRNA targets (can be 1 or 2 columns)
#   for more description see: score_calculation
post_score_calculation <- function(input.list, input.specs){
  if('postScoringAnalysis' %in% names(input.specs)){
    out.list <- list()
    if('RELICS_genomeScoring' %in% names(input.specs)){
      out.list <- c(out.list, RELICS_genomeScoring_wrapper(input.list, input.specs))
    }
    return(out.list)
  } else {
    return(input.list)
  }
}

# save the final scores used for plotting, AUC etc.
# input: list with scores per method, input specifications
# output: .csv for each methods scores
save_all_scores <- function(input.score.list, input.specs){
  score.names <- names(input.score.list)
  for(i in 1:length(score.names)){
    temp.score.df <- input.score.list[[i]]
    temp.score.name <- score.names[i]
    method.type <- strsplit(temp.score.name, '_')[[1]]
    method.identifier <- method.type[1]
    if(length(method.type) > 1){
      if(method.type[2] == 'genomeScores'){
        out.csv <- data.frame(genomeScore = temp.score.df$formatScores, chrom = temp.score.df$chrom,
          start = temp.score.df$start, end = temp.score.df$end, label = temp.score.df$label,
          log2_FC = temp.score.df$log2_rate_ratio, nrSupportGuides = temp.score.df$nrSupportGuides,
          stringsAsFactors = F)
      } else {
        out.csv <- data.frame(genomeScore = temp.score.df$formatScores, chrom = temp.score.df$chrom,
          start = temp.score.df$start, end = temp.score.df$end, label = temp.score.df$label,
          log2_FC = temp.score.df$log2_rate_ratio, stringsAsFactors = F)
      }
      write.csv(out.csv, file = paste0(input.specs$dataName,  '_', method.identifier,  '_', method.type[2], '.csv'), row.names = FALSE)
    }
  }
}

# given a list of scores, create bedgraphs for each one of them
create_bedgraphs <- function(input.score.list, input.specs){
  score.names <- names(input.score.list)
  nr.bg <- length(score.names)  # create unique colors for each bedgraph to plot
  bedgraph.colors <- col2rgb(rainbow(nr.bg,s = 1, v = 1, start = 0, end =  max(1, nr.bg - 1)/nr.bg, alpha = 1))
  for(i in 1:length(score.names)){
    intit.temp.score.df <- input.score.list[[i]]
    temp.score.df <- intit.temp.score.df[which(! is.na(intit.temp.score.df$start)),]
    temp.score.name <- score.names[i]
    print(paste0('Creating Bedgraph for: ', temp.score.name))
    temp.bg.color <- paste(bedgraph.colors[1,i], bedgraph.colors[2,i], bedgraph.colors[3,i], sep = ',')
    first.chrom <- temp.score.df[which(temp.score.df$chrom == temp.score.df$chrom[1]),]
    temp.header1 <- paste0('browser position ', first.chrom$chrom[1], ':', min(first.chrom$start, na.rm = T),
      '-', max(first.chrom$end, na.rm = T), '\n')
    temp.header2 <- paste0("track type=bedGraph name='", input.specs$dataName, temp.score.name, "_bg' description='",
      input.specs$dataName, temp.score.name, "_bg' visibility=full color=", temp.bg.color, '\n')

    temp.score.gr <- GRanges(seqnames = temp.score.df$chrom, ranges = IRanges(temp.score.df$start, temp.score.df$end),
      score = temp.score.df$formatScores)
    temp.score.gr <- sortSeqlevels(temp.score.gr)
    temp.score.gr <- sort(temp.score.gr)
    df.temp.score <- as.data.frame(temp.score.gr)
    df.temp.score.final <- df.temp.score[,c(1,2,3,6)]

    tmp.file <- paste0(input.specs$dataName, '_',temp.score.name, ".bedgraph")
    cat(temp.header1, file = tmp.file)
    cat(temp.header2, file = tmp.file, append = TRUE)

    write.table(df.temp.score.final, file = tmp.file, append = TRUE, row.names = F, col.names = F, quote = F, sep = '\t')

    temp.method.name <- strsplit(score.names[i],'_')[[1]]

  }
}

# wrapper function for various plotting options
score_plotting <- function(input.scores, input.specs){

  scores.to.plot <- input.scores

  # below is depricated: FDR scores are currently not calculated (conflict with HMM mehtods)
  # # for alphaRRA we want both p-values and FDR andjusted p-values
  # # need to duplicate results from alphaRRA and convert formatScores to -log10 of rawScores
  # if('postScoringAnalysis' %in% names(input.specs)){
  #   scores.to.duplicate <- which(unlist(lapply(strsplit(names(scores.to.plot), '_'), function(x){x[2] == 'alphaRRA'})) == TRUE)
  #   if(length(scores.to.duplicate) >= 1){
  #     score.names <- names(scores.to.plot)
  #     for(i in 1:length(scores.to.duplicate)){
  #       temp.name <- paste0(score.names[scores.to.duplicate[i]], '_FDR')
  #       temp.scores <- scores.to.plot[[scores.to.duplicate[i]]]
  #       temp.scores$formatScores <- -log10(temp.scores$rawScores)
  #       scores.to.plot[[temp.name]]<- temp.scores
  #     }
  #   }
  # }
  score.names <- names(scores.to.plot)

  # need to account for viterbi labelling
  orig.label.hierarchy <- input.specs$labelHierarchy

  for(i in 1:length(score.names)){
    raw.methd.scores <- scores.to.plot[[i]]
    na.rows <- which(is.na(raw.methd.scores$rawScores))

    # if not enough reads present, DESeq2 will return NA
    # remove these rows before plotting
    methd.scores <- c()
    if(length(na.rows) > 0){
      methd.scores <- raw.methd.scores[-na.rows,]
    } else {
      methd.scores <- raw.methd.scores
    }

    methd.name <- score.names[i]
    if('viterbiLabels' %in% colnames(methd.scores)){
      input.specs$labelHierarchy <- c(paste(orig.label.hierarchy, 'nonEnh', sep = '_'), paste(orig.label.hierarchy, 'enh', sep = '_'))
    } else {
      input.specs$labelHierarchy <- orig.label.hierarchy
    }
    if('plotAllChrom' %in% names(input.specs)){
      all.chroms <- unique(methd.scores$chrom)
      plotMultiChrom(methd.scores, input.specs, methd.name, all.chroms, 'plotAllChrom')
    }
    if('plotChroms' %in% names(input.specs)){
      chroms.to.plot <- input.specs$plotChroms
      plotMultiChrom(methd.scores, input.specs, methd.name, chroms.to.plot, 'plotChroms')
    }
    if('plotSeparateChroms' %in% names(input.specs)){
      chroms.to.plot <- unique(methd.scores$chrom)

      # if there is a 'NA' chromosome, remove it since there is no point in
      # plotting non-targeting guides by themselves
      if(NA %in% chroms.to.plot){
        chroms.to.plot <- chroms.to.plot[-which(is.na(chroms.to.plot))]
      }
      if('NA' %in% chroms.to.plot){
        chroms.to.plot <- chroms.to.plot[-which(chroms.to.plot == 'NA')]
      }
      plotSingleChrom(methd.scores, input.specs, methd.name, chroms.to.plot, 'plotSeparateChroms')
    }
    if('zoomChromosomePlot' %in% names(input.specs)){
      zoom.file <- read.csv(input.specs$zoomRange, as.is = TRUE, stringsAsFactors = FALSE)
      plotZoomChrom(methd.scores, input.specs, methd.name, chroms.to.plot, 'plotZoomChroms', zoom.file)
    }
    if('tileZoomChromosomePlot' %in% names(input.specs)){
      zoom.file <- read.csv(input.specs$tileZoomRange, as.is = TRUE, stringsAsFactors = FALSE)
      plotTileZoomChrom(methd.scores, input.specs, methd.name, chroms.to.plot, 'plotTileZoomChroms', zoom.file)
    }
  }
}

# plot multiple chromosome scores next to each other
# chromosome coordinates will be adjusted
# input:
#   input.scores: data frame containing scores (rawScores, formatScores, log2_rate_ratio, chromosome, label, start, end)
#   input.specs: specification list
#   input.name: name of scores
#   input.chroms: chromosomes to be plotted next to each other
#   plot.type: specifying whether it's plotting all chromosomes or a subset of them
# output:
#   pdf containing scores per chromosome
plotMultiChrom <- function(input.scores, input.specs, input.name, input.chroms, plot.type){

  chrom.df <- filter(input.scores, grepl('chr',chrom))
  all.chroms <- input.chroms[grep('chr', input.chroms, perl = TRUE)]
  na.df <- c()

  if(nrow(chrom.df) < nrow(input.scores)){
    na.df <- filter(input.scores, grepl('NA',chrom))
  }

  x.pos <- c()
  y.pos <- c()
  final.labels <- c()
  chrom.range <- c()
  inter.chr.space <- c()
  gene.pos.df <- c()  # data frame containing the coordinates for genes which are to be plotted

  for(chrom in all.chroms){
    temp.chrom.df <- chrom.df[chrom.df$chrom == chrom,]
    sorted.temp.chrom.df <- temp.chrom.df[order(temp.chrom.df$start),]
    temp.x.pos <- (sorted.temp.chrom.df$start + sorted.temp.chrom.df$end) / 2
    temp.y.pos <- sorted.temp.chrom.df$formatScores
    temp.labels <- sorted.temp.chrom.df$label

    if(length(x.pos) == 0){
      #have x-position be 1, and keep spacing between guides
      #define spacing between chromosomes
      shift.temp.x.pos <- temp.x.pos - temp.x.pos[1] + 1
      # inter.chr.space <- round(0.25*mean(shift.temp.x.pos))
      inter.chr.space <- round(1/8* (max(temp.x.pos) - min(temp.x.pos)))
      #if specific gene are to be plotted, check if for current chromosome a
      # gene is to be plotted, adjust gene position relative to the shift made
      gene.pos.df <- get_gene_Position(input.specs, chrom, shift.temp.x.pos, gene.pos.df, temp.x.pos[1], 0)

    } else {
      #have x-position be 1, and keep spacing between guides
      # start where the last chromosome left off with plus the interchromosome space to visually keep them apart
      coord.shift <- inter.chr.space + max(x.pos)
      shift.temp.x.pos <- temp.x.pos - temp.x.pos[1] + 1 + coord.shift
      gene.pos.df <- get_gene_Position(input.specs, chrom, shift.temp.x.pos, gene.pos.df, temp.x.pos[1], coord.shift)
    }
    x.pos <- c(x.pos, shift.temp.x.pos)
    y.pos <- c(y.pos, temp.y.pos)
    final.labels <- c(final.labels, temp.labels)

    chrom.range <- get_chrom_range(chrom, shift.temp.x.pos, gene.pos.df, chrom.range)
  }

  if(nrow(chrom.df) < nrow(input.scores)){
    all.na.labels <- unique(na.df$label)
    for(lab in all.na.labels){
      na.lab.df <- na.df[na.df$label == lab,]
      na.lab.start <- max(c(chrom.range[,3],x.pos)) + inter.chr.space
      na.lab.stepSize <- median(sort(x.pos)[2:length(x.pos)] - sort(x.pos)[1:(length(x.pos) - 1)])  # adjust step size to median step size of targeting guides
      na.lab.x.pos <- na.lab.start + na.lab.stepSize*c(1:nrow(na.lab.df))

      x.pos <- c(x.pos, na.lab.x.pos)
      y.pos <- c(y.pos, na.lab.df$formatScores)
      final.labels <- c(final.labels, na.lab.df$label)
    }
  }

  if('showFDR' %in% names(input.specs)){
    fdrs <- c()
    for(i in 1:length(input.specs$showFDR)){
      # amongs the rows which are below a certain FDR in the rawScore, return the minimum raw
      temp.fdr <- min(input.scores$formatScores[which(input.scores$rawScores <= input.specs$showFDR[i])])
      fdrs <- c(fdrs, temp.fdr)
    }
    input.specs$FDRcutoffs <- fdrs
  }

  pdf(paste(input.specs$dataName, input.name, plot.type, '.pdf',sep = '_'),width = 14, height = 7)
  create_chromosome_scorePlot(x.pos, y.pos, final.labels, gene.pos.df, chrom.range, input.specs)
  dev.off()
}

# plot chromosomes individually, option to zoom in into specific region
# input:
#   input.scores: data frame containing scores (rawScores, formatScores, log2_rate_ratio, chromosome, label, start, end)
#   input.specs: specification list
#   input.name: name of scores
#   input.chroms: chromosomes to be plotted next to each other
#   plot.type: specifying whether it's plotting all chromosomes or a subset of them
# output:
#   pdf containing scores per chromosome
plotSingleChrom <- function(input.scores, input.specs, input.name, input.chroms, plot.type){

  chrom.df <- filter(input.scores, grepl('chr',chrom))
  chrom.midpoints <- (chrom.df$start + chrom.df$end) / 2
  all.chroms <- input.chroms[grep('chr', input.chroms, perl = TRUE)]
  na.df <- c()

  if(nrow(chrom.df) < nrow(input.scores)){
    na.df <- filter(input.scores, grepl('NA',chrom))
  }

  for(chrom in all.chroms){
    temp.chrom.df <- chrom.df[chrom.df$chrom == chrom,]
    temp.x.pos <- chrom.midpoints[chrom.df$chrom == chrom]
    temp.y.pos <- temp.chrom.df$formatScores
    temp.labels <- temp.chrom.df$label
    temp.inter.chr.space <- round(1/8* (max(temp.x.pos) - min(temp.x.pos))) + 1

    gene.pos.df <- c()
    gene.pos.df <- get_gene_Position(input.specs, chrom, temp.x.pos, gene.pos.df, 0, 0)

    if(! is.null(na.df) ){
      all.na.labels <- unique(na.df$label)
      for(lab in all.na.labels){
        na.lab.df <- na.df[na.df$label == lab,]
        na.lab.start <- max(temp.x.pos) + temp.inter.chr.space
        na.lab.stepSize <- median(sort(temp.x.pos)[2:length(temp.x.pos)] - sort(temp.x.pos)[1:(length(temp.x.pos) - 1)])  # adjust step size to median step size of targeting guides
        if(is.na(na.lab.stepSize)){
          na.lab.stepSize <- 1
        }
        na.lab.x.pos <- na.lab.start + na.lab.stepSize*c(1:nrow(na.lab.df))

        temp.x.pos <- c(temp.x.pos, na.lab.x.pos)
        temp.y.pos <- c(temp.y.pos, na.lab.df$formatScores)
        temp.labels <- c(temp.labels, na.lab.df$label)
      }
    }
    chrom.range <- NULL

    if('showFDR' %in% names(input.specs)){
      fdrs <- c()
      for(i in 1:length(input.specs$showFDR)){
        # amongs the rows which are below a certain FDR in the rawScore, return the minimum raw
        temp.fdr <- min(input.scores$formatScores[which(input.scores$rawScores <= input.specs$showFDR[i])])
        fdrs <- c(fdrs, temp.fdr)
      }
      input.specs$FDRcutoffs <- fdrs
    }

    pdf(paste(input.specs$dataName, input.name, plot.type, chrom, '.pdf',sep = '_'),width = 14, height = 7)
    create_chromosome_scorePlot(temp.x.pos, temp.y.pos, temp.labels, gene.pos.df, chrom.range, input.specs)
    dev.off()
  }
}

# given a chromosome to plot, obtaing the gene coordinates, if available, and adjust accordingly
# input:
#   input.specs: specification list
#   input.chrom: chromsome for which gene coordinates are to be obtained (ex: 'chr6')
#   input.coords: chromosome coordinates of all scores (used to adjust gene coordinates)
#   input.gene.pos: exisiting gene coordinates for other chromsomes, already adjusted (cols: chrom, start, end)
#   input.initial.shift: if multiple chromosomes are plotted next to each other , this shoul dbe the minimum x-value, to start the coordinates at 1
#   input.additional.coord.shift: number by which genes should be shifted after setting coordinates using initial.shift, used if multiple chromosomes are plotted at the same time
# output: data dataframe
#   cols. chrom, start, end
get_gene_Position <- function(input.specs, input.chrom, input.coords, input.gene.pos, input.initial.shift, input.additional.coord.shift){
  sorted.input.coords <- sort(input.coords)
  if('plotGenes' %in% names(input.specs)){
    chrom.gene.rows <- which(input.specs$plotGenesData[,1] == input.chrom)
    if(length(chrom.gene.rows) > 0){
      temp.gene.pos <- input.specs$plotGenesData[chrom.gene.rows,]
      new.gene.pos <- cbind.data.frame(chrom = temp.gene.pos[,1],
        start = (temp.gene.pos[,2] - input.initial.shift + 1 + input.additional.coord.shift),
        end = (temp.gene.pos[,3] - input.initial.shift + 1 + input.additional.coord.shift), stringsAsFactors = FALSE)
      gene.pos.df <- rbind(input.gene.pos, new.gene.pos)
    } else {  # if no gene is specified for said chromosome simply give a minimal range
      new.gene.pos <- cbind.data.frame(chrom = input.chrom,
        start = sorted.input.coords[1], end = sorted.input.coords[1] + 1, stringsAsFactors = FALSE)
      gene.pos.df <- rbind(input.gene.pos, new.gene.pos)
    }
  } else {
    new.gene.pos <- cbind.data.frame(chrom = input.chrom,
      start = sorted.input.coords[1], end = sorted.input.coords[1] + 1, stringsAsFactors = FALSE)
    gene.pos.df <- rbind(input.gene.pos, new.gene.pos)
  }
  return(gene.pos.df)
}

# for a given chromosoem, create the sector delineating start and end (for plotting)
# input:
#   input.chrom: chromsome for which sector to be obtained (ex: 'chr6')
#   input.coords: chromosome coordinates of all scores (used to adjust gene coordinates)
#   input.gene.pos: exisiting gene coordinates for other chromsomes, already adjusted (cols: chrom, start, end)
#   input.range: existing sectors from other chromosomes
# output: data frame
#   cols: chromRange, start, end
get_chrom_range <- function(input.chrom, input.coords, input.gene.pos, input.range){
  chrom.gene.pos <- input.gene.pos[input.gene.pos$chrom == input.chrom,]
  chrom.range <- cbind.data.frame(chromRange = input.chrom,
    start = min(min(input.coords), min(chrom.gene.pos[,2])),
    end =  max(max(input.coords), max(chrom.gene.pos[,3])), stringsAsFactors = FALSE)
  out.chr.range <- rbind(input.range, chrom.range)
  return(out.chr.range)
}

#plot the p-values along the chromosome of interest
create_chromosome_scorePlot <- function(input.x.values, input.y.values, input.labels, input.gene.pos, input.chrom.range, input.specs){
  chrom.info <- cbind.data.frame(coords = input.x.values, scores = input.y.values, labels = input.labels)

  out.chrom.plot <- ggplot()+
    geom_point(data = chrom.info[chrom.info$label == input.specs$labelHierarchy[1], ], aes(x = coords, y = scores,
      colour = labels), size = 1)

  for(i in 2:length(input.specs$labelHierarchy)){
    if(input.specs$labelHierarchy[i] %in% unique(input.labels)){
      out.chrom.plot <- out.chrom.plot + geom_point(data = chrom.info[chrom.info$label == input.specs$labelHierarchy[i], ],
        aes(x = coords, y = scores, colour = labels), size = 1.5)
    }
  }

  # keep chromosome range in a list to avoid overrriding of level specific variables
  if(! is.null(input.chrom.range)){
    chrom.range.list <- list()
    for(i in 1:nrow(input.chrom.range)){
      chrom.range <- cbind.data.frame(pos = c(input.chrom.range[i,2],input.chrom.range[i,3]),
      val = c((min(input.y.values)-2),(min(input.y.values)-2)))
      chrom.range.list[[i]] <- as.data.frame(chrom.range)
      out.chrom.plot <- out.chrom.plot + geom_line(data = chrom.range.list[[i]], aes(x = pos,y = val, size = 3),colour = 'red')
    }
  }

  # keep gene positions in a list to avoid overrriding of level specific variables
  if(! is.null(input.gene.pos)){
  gene.pos.list <- list()
    for(i in 1:nrow(input.gene.pos)){
      gene.pos <- cbind.data.frame(pos = c(input.gene.pos[i,2],input.gene.pos[i,3]), val = c(0,0))
      gene.pos.list[[i]] <- as.data.frame(gene.pos)
      out.chrom.plot <- out.chrom.plot + geom_line(data = gene.pos.list[[i]], aes(x = pos,y = val, size = 4),colour = 'orange')
    }
  }

  if('showFDR' %in% names(input.specs)){
    fdr.cutoff.list <- list()
    for(i in 1:length(input.specs$showFDR)){
      fdr.line <- cbind.data.frame(pos = c(min(chrom.info$coords),max(chrom.info$coords)),
        FDRcutoff = c(input.specs$FDRcutoffs[i], input.specs$FDRcutoffs[i]))
      fdr.cutoff.list[[i]] <- fdr.line
      out.chrom.plot <- out.chrom.plot + geom_line(data = fdr.cutoff.list[[i]],
        aes(x = pos, y = FDRcutoff, size = 0.8), colour = 'blue')
    }
  }

  out.chrom.plot <- out.chrom.plot + theme_bw() + scale_colour_discrete(name = 'Scores')
  print(out.chrom.plot)
}

# plot each zoom regions individually
# input:
#   input.scores: data frame containing scores (rawScores, formatScores, log2_rate_ratio, chromosome, label, start, end)
#   input.specs: specification list
#   input.name: name of scores
#   input.chroms: chromosomes to be plotted next to each other
#   plot.type: specifying whether it's plotting all chromosomes or a subset of them
#   zoom.regions: regions which shuld be zoomed into (cols; chrom, start, end)
# output:
#   pdf containing scores per chromosome
plotZoomChrom <- function(input.scores, input.specs, input.name, input.chroms, plot.type, zoom.regions){

  chrom.df <- filter(input.scores, grepl('chr',chrom))
  na.df <- c()

  if(nrow(chrom.df) < nrow(input.scores)){
    na.df <- filter(input.scores, grepl('NA',chrom))
  }

  for(i in 1:nrow(zoom.regions)){
    zoom.df <- input.scores[input.scores$chrom == zoom.regions$chrom[i] & input.scores$start >= zoom.regions$start[i] & input.scores$end <= zoom.regions$end[i], ]
    zoom.x.pos <- (zoom.df$start + zoom.df$end) / 2
    zoom.y.pos <- zoom.df$formatScores
    zoom.labels <- zoom.df$label
    zoom.inter.chr.space <- round(1/8* (max(zoom.x.pos) - min(zoom.x.pos)))

    gene.pos.df <- c()
    if(! is.null(input.specs$plotGenesData)){
      all.gene.pos <- input.specs$plotGenesData
      zoom.gene.df <- all.gene.pos[all.gene.pos$chrom == zoom.regions$chrom[i] & all.gene.pos$start >= zoom.regions$start[i] & all.gene.pos$end <= zoom.regions$end[i], ]
      if(nrow(zoom.gene.df) > 0){
        gene.pos.df <- zoom.gene.df
      }
    }

    if(nrow(chrom.df) < nrow(input.scores)){
      all.na.labels <- unique(na.df$label)
      for(lab in all.na.labels){
        na.lab.df <- na.df[na.df$label == lab,]
        na.lab.start <- max(zoom.x.pos) + zoom.inter.chr.space
        na.lab.stepSize <- median(sort(zoom.x.pos)[2:length(zoom.x.pos)] - sort(zoom.x.pos)[1:(length(zoom.x.pos) - 1)])  # adjust step size to median step size of targeting guides
        na.lab.x.pos <- na.lab.start + na.lab.stepSize*c(1:nrow(na.lab.df))

        zoom.x.pos <- c(zoom.x.pos, na.lab.x.pos)
        zoom.y.pos <- c(zoom.y.pos, na.lab.df$formatScores)
        zoom.labels <- c(zoom.labels, na.lab.df$label)
      }
    }
    chrom.range <- NULL
    pdf(paste(input.specs$dataName, input.name, plot.type, 'zoom', i, '.pdf',sep = '_'),width = 14, height = 7)
    create_chromosome_scorePlot(zoom.x.pos, zoom.y.pos, zoom.labels, gene.pos.df, chrom.range, input.specs)
    dev.off()
  }
}

# based on list input, create the analysis specification file:
write_specs_file <- function(input.specs.list, input.filename){

  specs.names <- names(input.specs.list)
  for(name in specs.names){
    temp_line <- paste(name,':',paste(input.specs.list[[name]],collapse = ','),sep='')
    write(temp_line,file = paste(input.filename, '.txt', sep = ''),append = TRUE)
  }

}
