suppressMessages(library(GenomicRanges))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))

suppressMessages(library(poibin)) # to calculate Poisson-Binomial

suppressMessages(library(extraDistr)) # Dirichlet-Multinomial distributions

suppressMessages(library(gtools)) # for getting all combinations when computing correlations between PPs


#' @title RELICS 2.0 analysis function. Uses IBSS to return a set of functional sequences FS for CRISPR regulatory screens
#' @param input.parameter.file: location of file containing all analysis parameters
#' @param input.parameter.list: RELICS analysis parameters in list format. Default = NULL
#' @param record.all.fs: logical, if information of all intermediate FS should be recorded, in addition to the final setof FS
#' @return list of final per-layer posteriors
#' @export RELICS()

RELICS <- function(input.parameter.file, input.parameter.list = NULL, data.file.split = FALSE, 
                   record.all.fs = FALSE, return.init.hypers = FALSE){

  analysis.parameters <- list()

  if(!is.null(input.parameter.list)){
    print('Using provided parameter list')
    analysis.parameters <- check_parameter_list(input.parameter.list, data.file.split)
  } else {
    print('Reading in provided parameter file')
    analysis.parameters <- read_analysis_parameters(input.parameter.file)
    analysis.parameters <- check_parameter_list(analysis.parameters, data.file.split)
  }

  # if parameters are missing, then do not run the analysis
  if(is.logical(analysis.parameters)){
    print('Please fill in the missing parameters before running RELICS')
    break()
  }

  data.setup <- set_up_RELICS_data(analysis.parameters, data.file.split,
                                   guide.offset = analysis.parameters$crisprEffectRange,
                                   repl_pools = analysis.parameters$repl_groups,
                                   labelHierarchy = analysis.parameters$labelHierarchy,
                                   fs0.label = analysis.parameters$FS0_label,
                                   file.save.dir = paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName),
                                   save.files = analysis.parameters$save_input_files,
                                   min.seg.dist = analysis.parameters$seg_dist)
  
  # the whole section below should be redundant
  # # if guide efficiency scores are provided, calculate guide efficiency and include in the model
  # if(! is.null(data.setup$guide_efficiency_scores)){
  #   
  #   # fixed_ge_coeff
  #   data.setup$guide_efficiency <- 1 / (1 + exp(-(data.setup$guide_efficiency_scores %*% rep(1, ncol(data.setup$guide_efficiency_scores)))))
  #         # 
  #   # data.setup$fix_guideEfficiency <- analysis.parameters$fix_guideEfficiency
  #   # data.setup$estimate_ge_alphaDiffScale <- analysis.parameters$estimate_ge_alphaDiffScale
  #   # data.setup$ge_beta_estimation <- analysis.parameters$ge_beta_estimation
  #   # 
  #   # # unde the conditions below, it is assumed that the data.setup$guide_efficiency_scores matrix is only one column
  #   # # therefore convert to: 
  #   # if(analysis.parameters$fix_guideEfficiency & ! analysis.parameters$estimate_ge_alphaDiffScale){
  #   #   data.setup$guide_efficiency <- data.setup$guide_efficiency_scores[,1]
  #   # } else {
  #   #   print('non-fixed guide efficiency or reestimation of alpha diff scale is not yet implemented')
  #   #   # immediate prining of the betas
  #   #   # also, might want to see if 'ge_beta_estimation' is required here
  #   # }
  # }

  # check if hyper parameters are provided, otherwise estimate from data
  if(! 'hyper_pars' %in% names(analysis.parameters)){

    #labels for background
    background.labels <- c()
    if(analysis.parameters$background_label_specified){
      background.labels <- analysis.parameters$background_label
    } else {
      background.labels <- analysis.parameters$labelHierarchy[-which(analysis.parameters$labelHierarchy %in% analysis.parameters$FS0_label)]
    }
    
    fs0.alphas <- estimate_hyper_parameters(data.setup,
                                            analysis.parameters,
                                            #input.guide.offset = analysis.parameters$crisprEffectRange,
                                            input.repl.pools = analysis.parameters$repl_groups,
                                            #input.labelHierarchy = analysis.parameters$FS0_label,
                                            fs0.label = analysis.parameters$FS0_label,
                                            analysis.parameters$one_dispersion
                                            #min.seg.dist = analysis.parameters$seg_dist
                                            )

    # background.alpha0 <- estimate_hyper_parameters(analysis.parameters,
    #                                                data.file.split,
    #                                                 input.guide.offset = analysis.parameters$crisprEffectRange,
    #                                                 input.repl.pools = analysis.parameters$repl_groups,
    #                                                 input.labelHierarchy = background.labels,
    #                                                fs0.label = background.labels,
    #                                                 min.seg.dist = analysis.parameters$seg_dist)
    # 
    # fs0.alpha1 <- estimate_hyper_parameters(analysis.parameters,
    #                                         data.file.split,
    #                                         input.guide.offset = analysis.parameters$crisprEffectRange,
    #                                         input.repl.pools = analysis.parameters$repl_groups,
    #                                         input.labelHierarchy = analysis.parameters$FS0_label,
    #                                         fs0.label = analysis.parameters$FS0_label,
    #                                         min.seg.dist = analysis.parameters$seg_dist)

    analysis.parameters$hyper_pars <- fs0.alphas #list(alpha0 = background.alpha0, alpha1 = fs0.alpha1, L = 1)
  }
  
  if(return.init.hypers){
    # record alphas used for posterior calculation
    alpha.out.dir <- paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName)
    record_alphas(background.alpha0, fs0.alpha1, alpha.out.dir, 'init')
    break()
  }
  
  data.setup$fixed_ge_coeff <- analysis.parameters$fixed_ge_coeff
       
  # if guide efficiency scores are provided, calculate guide efficiency and include in the model
  if(! is.null(data.setup$guide_efficiency_scores)){
    
    if(! analysis.parameters$fixed_ge_coeff){
      
      temp.relics.params <- init_relics_param(analysis.parameters$hyper_pars, data.setup)

      ge.list <- recompute_ge_coefficients(temp.relics.params,
                                           analysis.parameters$hyper_pars,
                                           data.setup$data,
                                           data.setup$guide_to_seg_lst,
                                           data.setup$guide_efficiency_scores,
                                           c(0, rep(1, ncol(data.setup$guide_efficiency_scores))))
      
      data.setup$guide_efficiency <- ge.list$guide_efficiency
      data.setup$ge_coeff <- ge.list$ge_coeff
    }
    
    # data.setup$fix_guideEfficiency <- analysis.parameters$fix_guideEfficiency
    # data.setup$estimate_ge_alphaDiffScale <- analysis.parameters$estimate_ge_alphaDiffScale
    # data.setup$ge_beta_estimation <- analysis.parameters$ge_beta_estimation
    # 
    # # unde the conditions below, it is assumed that the data.setup$guide_efficiency_scores matrix is only one column
    # # therefore convert to: 
    # if(analysis.parameters$fix_guideEfficiency & ! analysis.parameters$estimate_ge_alphaDiffScale){
    #   data.setup$guide_efficiency <- data.setup$guide_efficiency_scores[,1]
    # } else {
    #   print('non-fixed guide efficiency or reestimation of alpha diff scale is not yet implemented')
    #   # immediate prining of the betas
    #   # also, might want to see if 'ge_beta_estimation' is required here
    # }
  }

  run_RELICS_2(input.data = data.setup,
               final.layer.nr = analysis.parameters$min_FS_nr,
               out.dir = paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName),
               input.hypers = analysis.parameters$hyper_pars,
               fix.hypers = analysis.parameters$fix_hypers,
               nr.segs = analysis.parameters$nr_segs,
               geom.p = analysis.parameters$geom_p,
               fs.correlation.cutoff = analysis.parameters$fs_correlation_cutoff,
               input.min.rs.pp = analysis.parameters$min_fs_pp,
               auto.stop = analysis.parameters$auto_stop,
               record.all.fs = record.all.fs,
               input.convergence.tol = analysis.parameters$convergence_tol,
               adjust.tol = analysis.parameters$adjust_tol,
               one.dispersion = analysis.parameters$one_dispersion,
               recompute.fs0 = analysis.parameters$recompute_fs0,
               local.max = analysis.parameters$local_max, 
               local.max.range = analysis.parameters$local_max_range)


}


#' @title Ensure that all required parameters for RELICS are correct and add default parameters where not specified
#' @param input.parameter.list: list: required: 'dataName','repl_groups', 'CountFileLoc', 'sgRNAInfoFileLoc', 'min_FS_nr', 'crisprSystem', 'FS0_label'
#' @return list with parameters for running RELICS or logical FALSE if required parameter is missing
#' @export check_parameter_list()

check_parameter_list <- function(input.parameter.list, data.file.split){

  out.parameter.list <- input.parameter.list

  par.given <- names(input.parameter.list)

  # set default to TRUE. But if not background is provided, then switch to false
  out.parameter.list$background_label_specified <- TRUE

  if(! 'out_dir' %in% par.given){
    out.parameter.list$out_dir <- getwd()
  }
  if(! 'adjust_tol' %in% par.given){
    out.parameter.list$adjust_tol <- FALSE #TRUE
  }
  if(! 'convergence_tol' %in% par.given){
    out.parameter.list$convergence_tol <- 0.1
  }
  if(! 'fix_hypers' %in% par.given){
    out.parameter.list$fix_hypers <- FALSE
  }
  if(! 'recompute_fs0' %in% par.given){
    out.parameter.list$recompute_fs0 <- FALSE
  }
  if(! 'local_max' %in% par.given){
    out.parameter.list$local_max <- FALSE
  }
  if(! 'local_max_range' %in% par.given){
    out.parameter.list$local_max_range <- 5
  }
  if(! 'iterative_hyper_est' %in% par.given){
    out.parameter.list$iterative_hyper_est <- FALSE
  }
  if(! 'nr_segs' %in% par.given){
    out.parameter.list$nr_segs <- 10
  }
  if(! 'geom_p' %in% par.given){
    #probability of success for geometric distribution
    out.parameter.list$geom_p <- 0.1
    #geom.norm.factr <- pgeom(nr.segs, geom.p)
  }
  # if(! 'seg_length_dir' %in% par.given){
  #   out.parameter.list$seg_length_dir <- 'continous-norm'
  # }
  if(! 'fs_correlation_cutoff' %in% par.given){
    out.parameter.list$fs_correlation_cutoff <- 0.1
  }
  if(! 'min_fs_pp' %in% par.given){
    out.parameter.list$min_fs_pp <- 0.1
    # input.min.rs.pp = 0.1,
  }
  if(! 'auto_stop' %in% par.given){
    out.parameter.list$auto_stop <- TRUE
  }
  if(! 'seg_dist' %in% par.given){
    out.parameter.list$seg_dist <- 100
  }
  if(! 'save_input_files' %in% par.given){
    out.parameter.list$save_input_files <- FALSE
  }
  if(! 'background_label' %in% par.given){
    out.parameter.list$background_label_specified <- FALSE
  }
  if(! 'one_dispersion' %in% par.given){
    out.parameter.list$one_dispersion <- TRUE
  }
  if(! 'dualToSingle' %in% par.given){
    out.parameter.list$dualToSingle <- FALSE
  }
  
  # guide efficiency related parameters
  if('guide_efficiency_loc' %in% par.given){
    temp.ge.file <- read.csv(input.parameter.list$guide_efficiency_loc, stringsAsFactors = F)
    if(! 'guide_efficiency_cols' %in% par.given){
      print(paste0("Error: Guide efficiency file provided but no column index given. Please provide a numeric input for 'guide_efficiency_cols'."))
      missing.parameters <- TRUE
    } else {
      out.parameter.list$guide_efficiency_scores <- temp.ge.file[,input.parameter.list$guide_efficiency_cols, drop = FALSE]
    }
    
    # check if fixed
    if(! 'fixed_ge_coeff' %in% par.given){
      out.parameter.list$fixed_ge_coeff <- FALSE
      # check if betas are provided, else give default
    } else {
      out.parameter.list$fixed_ge_coeff <- input.parameter.list$fixed_ge_coeff
    }
    
    # if(! 'fix_guideEfficiency' %in% par.given){
    #   out.parameter.list$fix_guideEfficiency <- TRUE
    #   # check if betas are provided, else give default
    # } else {
    #   out.parameter.list$fix_guideEfficiency <- input.parameter.list$fix_guideEfficiency
    # }
    
    if(! 'ge_betas' %in% par.given){
      print('Default guide efficiency betas are not yet implemented')
    } else {
      out.parameter.list$ge_betas <- input.parameter.list$ge_betas
    }
    
    if(! 'ge_beta_estimation' %in% par.given){
      out.parameter.list$ge_beta_estimation <- FALSE
    } else {
      out.parameter.list$ge_beta_estimation <- input.parameter.list$ge_beta_estimation
    }
    
    # #if the GE is fixed but there are multiple metics, then the alpha diff scaling has to be reestimated
    # # this flag is used for first time estimation only
    # if(out.parameter.list$fix_guideEfficiency & ncol(out.parameter.list$guide_efficiency_scores) < 2){
    #   out.parameter.list$estimate_ge_alphaDiffScale <- FALSE
    # } else { # if there are multiple scores the guide efficiency has to be roughly estimated from the 2+ scoring schemes
    #   out.parameter.list$estimate_ge_alphaDiffScale <- TRUE
    # }
    
    
    # ge_beta_estimation (wether the scales have to be re-estimated)
    # ToDo
    
  } else {
    out.parameter.list$guide_efficiency_scores <- NULL
    out.parameter.list$fixed_ge_coeff <- NULL
  }

  minimum.parameters <- c()
  if(data.file.split){
    minimum.parameters <- c('dataName','repl_groups', 'CountFileLoc', 'sgRNAInfoFileLoc', 'min_FS_nr', 'crisprSystem', 'FS0_label')
  } else {
    minimum.parameters <- c('dataName','repl_groups', 'DataInputFileLoc', 'min_FS_nr', 'crisprSystem', 'FS0_label')
  }
  # crisprSystem options: CRISPRcas9, CRISPRa, dualCRISPR

  par.given <- names(input.parameter.list)
  missing.parameters <- FALSE

  for(i in 1:length(minimum.parameters)){

    if(! minimum.parameters[i] %in% par.given){

      print(paste0('Missing parameter: ', minimum.parameters[i]))
      missing.parameters <- TRUE
    }

    # for the info file, make sure that the guides are ordered and contain no NAs
    if('sgRNAInfoFileLoc' == minimum.parameters[i] | 'DataInputFileLoc' == minimum.parameters[i]){
      test.info <- c()
      if('sgRNAInfoFileLoc' == minimum.parameters[i]){
        test.info <- read.csv(input.parameter.list$sgRNAInfoFileLoc, stringsAsFactors = F)
      } else {
        test.info <- read.csv(input.parameter.list$DataInputFileLoc, stringsAsFactors = F)
      }

      # for the info file, check if the label hierarch is given, else establish it from the info file
      if(!'labelHierarchy' %in% par.given){
        temp.info <- test.info #read.csv(input.parameter.list$sgRNAInfoFileLoc, stringsAsFactors = F)
        out.parameter.list$labelHierarchy <- c()
        counted.lables <- table(temp.info$label)

        while(length(counted.lables) > 1){
          max.lab.pos <- which(counted.lables == max(counted.lables)[1])[1]
          max.lab.name <- names(counted.lables)[max.lab.pos]

          out.parameter.list$labelHierarchy <- c(out.parameter.list$labelHierarchy, max.lab.name)

          counted.lables <- counted.lables[-max.lab.pos]
        }

        # then last part
        out.parameter.list$labelHierarchy <- c(out.parameter.list$labelHierarchy, names(counted.lables))
      }
    }

    # for the crisprSystem, check if effect range is given, else use default
    if(minimum.parameters[i] == 'crisprSystem'){
      if(! 'crisprEffectRange' %in% par.given){
        if(input.parameter.list$crisprSystem == 'CRISPRi'){
          out.parameter.list$crisprEffectRange <- 200
        } else if(input.parameter.list$crisprSystem == 'CRISPRa'){
          out.parameter.list$crisprEffectRange <- 200
        } else if(input.parameter.list$crisprSystem == 'CRISPRcas9'){
          out.parameter.list$crisprEffectRange <- 20
        } else if(input.parameter.list$crisprSystem == 'dualCRISPR'){
          out.parameter.list$crisprEffectRange <- 0
        } else{
          print(paste0('Please provide a valid CRISPR system: CRISPRi, CRISPRa, CRISPRcas9 or dualCRISPR'))
          missing.parameters <- TRUE
        }
      } else if(! out.parameter.list$crisprSystem %in% c('CRISPRi', 'CRISPRa', 'CRISPRcas9', 'dualCRISPR')){
        print(paste0('Please provide a valid CRISPR system: CRISPRi, CRISPRa, CRISPRcas9 or dualCRISPR'))
        missing.parameters <- TRUE
      }
    }
  }

  if(missing.parameters){
    return(FALSE)
  } else {
    return(out.parameter.list)
  }
}


#' @title Read in the RELICS parameters from an input file
#' @param parameter.file.loc: string, location of the parameter file
#' @return list with parameters for running RELICS
#' @export read_analysis_parameters()

read_analysis_parameters <- function(parameter.file.loc){
  raw.parameters <- scan(parameter.file.loc,what='character')
  out.parameter.list <- list()

  for(parameter in raw.parameters){

    parameter.id <- strsplit(parameter,':')[[1]][1]

    if('out_dir' == parameter.id){
      out.parameter.list$out_dir <- strsplit(parameter,':')[[1]][2]
    }
    if('adjust_tol' == parameter.id){
      out.parameter.list$adjust_tol <- as.logical(strsplit(parameter,':')[[1]][2])
    }
    if('convergence_tol' == parameter.id){
      out.parameter.list$convergence_tol <- as.numeric(strsplit(parameter,':')[[1]][2])
    }
    if('fix_hypers' == parameter.id){
      out.parameter.list$fix_hypers <- as.logical(strsplit(parameter,':')[[1]][2])
    }
    if('recompute_fs0' == parameter.id){
      out.parameter.list$recompute_fs0 <- as.logical(strsplit(parameter,':')[[1]][2])
    }
    if('local_max' == parameter.id){
      out.parameter.list$local_max <- as.logical(strsplit(parameter,':')[[1]][2])
    }
    if('local_max_range' == parameter.id){
      out.parameter.list$local_max_range <- as.numeric(strsplit(parameter,':')[[1]][2])
    }
    if('iterative_hyper_est' == parameter.id){
      out.parameter.list$iterative_hyper_est <- as.logical(strsplit(parameter,':')[[1]][2])
    }
    if('nr_segs' == parameter.id){
      out.parameter.list$nr_segs <- as.numeric(strsplit(parameter,':')[[1]][2])
    }
    if('geom_p' == parameter.id){
      #probability of success for geometric distribution
      out.parameter.list$geom_p <- as.numeric(strsplit(parameter,':')[[1]][2])
      #geom.norm.factr <- pgeom(nr.segs, geom.p)
    }
    # if('seg_length_dir' == parameter.id){
    #   out.parameter.list$seg_length_dir <- 'continous-norm'
    # }
    if('fs_correlation_cutoff' == parameter.id){
      out.parameter.list$fs_correlation_cutoff <- as.numeric(strsplit(parameter,':')[[1]][2])
    }
    if('min_fs_pp' == parameter.id){
      out.parameter.list$min_fs_pp <- as.numeric(strsplit(parameter,':')[[1]][2])
      # input.min.rs.pp = 0.1,
    }
    if('auto_stop' == parameter.id){
      out.parameter.list$auto_stop <- as.logical(strsplit(parameter,':')[[1]][2])
    }
    if(parameter.id == 'labelHierarchy'){
      out.parameter.list$labelHierarchy <- strsplit(strsplit(parameter.id,':')[[1]][2],',')[[1]]
    }
    if('dataName' == parameter.id){
      out.parameter.list$dataName <- strsplit(parameter,':')[[1]][2]
    }
    if('repl_groups' == parameter.id){
      out.parameter.list$repl_groups <- lapply(strsplit(strsplit(strsplit(parameter,':')[[1]][2],';')[[1]],','), as.numeric)
      names(out.parameter.list$repl_groups) <- paste0('repl_',c(1:length(out.parameter.list$repl_groups)))
    }
    if('DataInputFileLoc' == parameter.id){
      out.parameter.list$DataInputFileLoc <- strsplit(parameter,':')[[1]][2]
    }
    if('CountFileLoc' == parameter.id){
      out.parameter.list$CountFileLoc <- strsplit(parameter,':')[[1]][2]
    }
    if('sgRNAInfoFileLoc' == parameter.id){
      out.parameter.list$sgRNAInfoFileLoc <- strsplit(parameter,':')[[1]][2]
    }
    if('min_FS_nr' == parameter.id){
      out.parameter.list$min_FS_nr <- as.numeric(strsplit(parameter,':')[[1]][2])
    }
    if('crisprSystem' == parameter.id){
      out.parameter.list$crisprSystem <- strsplit(parameter,':')[[1]][2]
    }
    if('dualToSingle' == parameter.id){
      out.parameter.list$dualToSingle <- as.logical(strsplit(parameter,':')[[1]][2])
    }
    if('FS0_label' == parameter.id){
      out.parameter.list$FS0_label <- strsplit(parameter,':')[[1]][2]
    }
    if('background_label' == parameter.id){
      out.parameter.list$background_label <- strsplit(parameter,':')[[1]][2]
    }
    if('alpha0' == parameter.id){
      out.parameter.list$alpha0 <- lapply(strsplit(strsplit(strsplit(parameter,':')[[1]][2],';')[[1]],','), as.numeric)
      names(out.parameter.list$alpha0) <- paste0('repl_',c(1:length(out.parameter.list$alpha0)))
    }
    if('alpha1' == parameter.id){
      out.parameter.list$alpha1 <- lapply(strsplit(strsplit(strsplit(parameter,':')[[1]][2],';')[[1]],','), as.numeric)
      names(out.parameter.list$alpha1) <- paste0('repl_',c(1:length(out.parameter.list$alpha1)))
    }
    if('seg_dist' == parameter.id){
      out.parameter.list$seg_dist <- as.numeric(strsplit(parameter,':')[[1]][2])
    }

  }
  # aggregate alpha0 and alpha1 if provided
  if('alpha0' %in% names(out.parameter.list) & 'alpha1' %in% names(out.parameter.list)){
    out.parameter.list$hyper_pars <- list(alpha0 = out.parameter.list$alpha0, alpha1 = out.parameter.list$alpha1, L = 1)
  }
  if('one_dispersion' == parameter.id){
    out.parameter.list$one_dispersion <- as.logical(strsplit(parameter,':')[[1]][2])
  }

  ##############################################
  # guide efficiency related parameters
  # location of file (can be same as info, or total combines, or separate)
  if('guide_efficiency_loc' == parameter.id){
    out.parameter.list$guide_efficiency_loc <- strsplit(parameter,':')[[1]][2]
  }
  # cloumns to use for GE
  if('guide_efficiency_cols' == parameter.id){
    out.parameter.list$guide_efficiency_cols <- as.numeric(strsplit(parameter,':')[[1]][2])
  }
  # boolean whether to simply use given guide coefficients or trying to estimate them
  if('fixed_ge_coeff' == parameter.id){
    out.parameter.list$fixed_ge_coeff <- as.logical(strsplit(parameter,':')[[1]][2])
  }
  
  
  # # boolean whether to simply use given guide efficiencies or trying to estimate a scaling
  # if('fix_guideEfficiency' == parameter.id){
  #   out.parameter.list$fix_guideEfficiency <- as.logical(strsplit(parameter,':')[[1]][2])
  # }
  
  
  # to do: ge_betas, ge_beta_estimation (this serves for estimation of the scaling factors, logical)

  return(out.parameter.list)
}


#' @title set up data structures for RELICS analysis
#' @param input.parameter.list: parameter list either containing $CountFileLoc and $sgRNAInfoFileLoc or $DataInputFileLoc
#' @param data.file.split: whether the data is split into two files
#' @param guide.offset: range by which each guide has an effect beyon the target site (default = 500)
#' @param repl_pools: list, each element is a set of columns which correspond to a replicate
#' @param labelHierarchy: ordering lof guide labels, left to right signifies ordering hor hierarchy to override in case of multiple labels overlapping
#' @param fs0.label: label used for generating FS0
#' @param file.save.dir: directory in which to save the used counts, info and segment info files
#' @param save.files: logical, by default, do not save the files used for analysis.
#' @param min.seg.dist: minimum distance of a segment. Default = 100
#' @return list: data, fs0_idx, segLabels, seg_info, guide_to_seg_lst, seg_to_guide_lst, next_guide_lst
#' @export set_up_RELICS_data()

set_up_RELICS_data <- function(input.parameter.list, data.file.split, guide.offset = 500,
                               repl_pools, labelHierarchy = c("chr", "neg", "pos", "exon"),
                               fs0.label = 'exon', file.save.dir, save.files = FALSE, 
                               min.seg.dist = 100){

  sim.counts <- c()
  sim.info <- c()
  if(data.file.split){
    sim.counts <- read.csv(input.parameter.list$CountFileLoc, stringsAsFactors = F)
    sim.info <- read.csv(input.parameter.list$sgRNAInfoFileLoc, stringsAsFactors = F)
  } else {
    all.data <- read.csv(input.parameter.list$DataInputFileLoc, stringsAsFactors = F)

    # separate info from counts
    sim.info <- all.data[, c('chrom', 'start', 'end', 'label')]
    count.cols <- unlist(repl_pools)
    sim.counts <- all.data[, count.cols]

    # adjust the replicate pool idx
    repl_pools.original <- repl_pools

    pool.idx <- 1
    for(i in 1:length(repl_pools.original)){
      repl_pools[[i]] <- c(pool.idx:(pool.idx + length(repl_pools.original[[i]]) - 1) )
      pool.idx <- pool.idx + length(repl_pools.original[[i]])
    }

  }

  #if guide efficiency is used, add the scores to info file
  ges.cols <- c() # guide efficiency scores columns
  if(! is.null(input.parameter.list$guide_efficiency_scores)){
    cols.in.info <- ncol(sim.info)
    ges.cols <- c((cols.in.info + 1):((cols.in.info + 1) + (ncol(input.parameter.list$guide_efficiency_scores) - 1)))
    
    # adjust 0 and 1 values to avoid -Inf and Inf during logit/expit
    guide.eff.filt <- input.parameter.list$guide_efficiency_scores
    guide.eff.filt[guide.eff.filt == 0] <- 0.0001
    guide.eff.filt[guide.eff.filt == 1] <- 0.9999
    
    sim.info <- cbind(sim.info, guide.eff.filt)
  }

  if(length(which(! sim.info$label %in% labelHierarchy)) > 0){
    #print(paste0("Removing labels not present in the 'labelHierarchy': ", unique(sim.info$label[which(! sim.info$label %in% labelHierarchy)])))
    sim.counts <- sim.counts[-which(! sim.info$label %in% labelHierarchy),]
    sim.info <- sim.info[-which(! sim.info$label %in% labelHierarchy),]
    #save.files <- TRUE
  }

  if(length(which(is.na(sim.info$start))) > 0){
    print('Removing rows due to NAs in the guide target sites')
    sim.counts <- sim.counts[-which(is.na(sim.info$start)),]
    sim.info <- sim.info[-which(is.na(sim.info$start)),]
    #save.files <- TRUE
  }
  
  if(length(unique(sim.info$chrom)) > 1){
    print(paste0("Currently RELICS only processes one chromosome per analysis. Multi-chromosome analysis is in active development and will hopefully be deployed soon."))
    missing.parameters <- TRUE
  }

  if(is.unsorted(sim.info$start)){
    print('Guide targets are unsorted. Adjusting them!')
    sim.counts <- sim.counts[order(sim.info$start),]
    sim.info <- sim.info[order(sim.info$start),]
    #save.files <- TRUE
  }
  
  # if GE scores exist, extract them from the info file
  filtered.ges <- NULL
  if(! is.null(input.parameter.list$guide_efficiency_scores)){
    processed.ge.scores <- as.matrix(sim.info[, ges.cols, drop = F])
    filtered.ges <- processed.ge.scores
  }
  
  # initialize varaible
  sim.guide.seg.list <- c()
  
  if(input.parameter.list$crisprSystem == 'dualCRISPR' & input.parameter.list$dualToSingle){
    
    # set up two data frames, one for each guide
    sim.info.g1 <- sim.info
    sim.info.g2 <- sim.info
    
    sim.info.g1$start <- sim.info$start - guide.offset
    sim.info.g1$end <- sim.info$start + guide.offset
    
    sim.info.g2$start <- sim.info$end - guide.offset
    sim.info.g2$end <- sim.info$end + guide.offset
    
    sim.guide.seg.list <- adapt_data_to_regionFormat_forDualCRISPR(sim.counts, repl_pools, sim.info.g1, 
                                                                   sim.info.g2, labelHierarchy, min.seg.dist)
    
  } else {
    # adjust the target positions according to CRISPRi
    sim.info$start <- sim.info$start - guide.offset
    sim.info$end <- sim.info$end + guide.offset
    
    # generate the guide-segment matrix and the per-segment labels
    # seg_info, guide_to_seg_lst, seg_to_guide_lst, counts
    sim.guide.seg.list <- adapt_data_to_regionFormat(sim.counts, repl_pools, sim.info, labelHierarchy, min.seg.dist)
    
  }
  
  sim.seg.info <- sim.guide.seg.list$seg_info

  # set up the delta vector and the vector containing positions of positive controls
  segment.labels <- as.numeric(as.factor(sim.seg.info$label))

  # next guie index list
  next.guide.list <- generate_next_guide_list(sim.guide.seg.list$seg_to_guide_lst)

  format.data.beta <- list(data = sim.guide.seg.list$counts,
                           fs0_idx = which(sim.seg.info$label %in% fs0.label),
                           segLabels = segment.labels,
                           seg_info = sim.seg.info,
                           guide_to_seg_lst = sim.guide.seg.list$guide_to_seg_lst,
                           seg_to_guide_lst = sim.guide.seg.list$seg_to_guide_lst,
                           next_guide_lst = next.guide.list,
                           guide_efficiency_scores = filtered.ges)

  write.csv(sim.seg.info, file = paste0(file.save.dir, '_segmentInfo.csv'), row.names = F)

  if(save.files){
    write.csv(sim.seg.info, file = paste0(file.save.dir, '_segmentInfo.csv'), row.names = F)
    write.csv(sim.info, file = paste0(file.save.dir, '_guideInfo.csv'), row.names = F)
    write.csv(sim.counts, file = paste0(file.save.dir, '_guideCounts.csv'), row.names = F)
  }
  
  return(format.data.beta)

}

#' @title Adapt an info file to a overlap matrix and corresponding data structures and filter the counts, return by replicate. Adapted for dual guide screens
#' @param input.counts: data.frame, rows are guides, columns are pools
#' @param replicate.list: each element is a replicate
#' @param input.info.g1: data frame, for first guide: each row contains chrom, start, end, a label
#' @param input.info.g2: data frame, for second guide: each row contains chrom, start, end, a label
#' @param input.label.hierarchy: ordering lof guide labels, left to right signifies ordering hor hierarchy to override in case of multiple labels overlapping
#' @param min.seg.dist: minimum distance of a segment. Default = 100
#' @return list: seg_info, guide_to_seg_lst, seg_to_guide_lst, counts (list, per replicate)
#' @export adapt_data_to_regionFormat_forDualCRISPR()

adapt_data_to_regionFormat_forDualCRISPR <- function(input.counts, replicate.list, input.info.g1, 
                                                     input.info.g2, input.label.hierarchy, min.seg.dist = 100){
  
  info.targeting.df.g1 <- c()
  info.targeting.df.g2 <- c()
  counts.targeting.df <- c()
  if(length(which(is.na(input.info.g1$start))) > 0 | length(which(is.na(input.info.g2$start))) > 0){
    print('Error: data was not properly filtered!')
    break()
  }
  
  
  info.targeting.df.g1 <- input.info.g1
  info.targeting.df.g2 <- input.info.g2
  counts.targeting.df <- input.counts
  
  # create guide-segment matrix for targeting guides in dual guide design: seg_info, guide_to_seg_lst, seg_to_guide_lst
  targeting.gs.list <- create_targeting_guide_segment_matrix_forDualCRISPR(info.targeting.df.g1, info.targeting.df.g2, 
                                                                           input.label.hierarchy, min.seg.dist)
  
  targeting.gs.list$counts <- list()
  for(i in 1:length(replicate.list)){
    temp.counts <- counts.targeting.df[, replicate.list[[i]]]
    colnames(temp.counts) <- paste('pool', c(1:ncol(temp.counts)), sep = '_')
    temp.counts$n <- rowSums(temp.counts)
    targeting.gs.list$counts[[i]] <- temp.counts
  }
  
  return(targeting.gs.list)
  
}


#' @title Generate the guide to segment and segment to guide mappings along with the segment info file, for gual guide design
#' @note Modified by ont trimming guides at start and end but only end. Avoid having information-less segments
#' @param input.targeting.info.g1: data.frame for first guide, start, end, chrom, label
#'  @param input.targeting.info.g2: data.frame for second guide, start, end, chrom, label
#' @param input.label.hierarchy: order in which labels are assigned to segments
#' @param min.seg.dist: minimum distance of a segment. Default = 100
#' @return list: seg_info, guide_to_seg_lst, seg_to_guide_lst
#' @export create_targeting_guide_segment_matrix_forDualCRISPR()

create_targeting_guide_segment_matrix_forDualCRISPR <- function(input.targeting.info.g1, input.targeting.info.g2, 
                                                                input.label.hierarchy, min.seg.dist = 100){
  
  all.breaks <- unique(c(input.targeting.info.g1$start, input.targeting.info.g1$end, 
                         input.targeting.info.g2$start, input.targeting.info.g2$end))
  sorted.breaks <- sort(all.breaks)
  
  bin.start <- c() #sorted.breaks[1]
  bin.ends <- c()
  
  break.idx <- 1
  next.counter <- 0
  proxy.start <- sorted.breaks[1]
  while(break.idx < length(sorted.breaks)){
    
    next.counter <- 1
    while(sorted.breaks[break.idx + next.counter] - proxy.start < min.seg.dist){
      next.counter <- next.counter + 1
      if(break.idx + next.counter > length(sorted.breaks)){
        bin.start <- c(bin.start, proxy.start)
        bin.ends <- c(bin.ends, sorted.breaks[break.idx + next.counter - 1])
        break.idx <- break.idx + next.counter
        break()
      }
    }
    if(break.idx + next.counter > length(sorted.breaks)){
      break()
    }
    
    bin.start <- c(bin.start, proxy.start)
    bin.ends <- c(bin.ends, sorted.breaks[break.idx + next.counter])
    
    proxy.start <- sorted.breaks[break.idx + next.counter] + 1
    break.idx <- break.idx + next.counter
    
  }
  
  region.df <- data.frame(chrom = rep(input.targeting.info.g1$chrom[1], length(bin.start)),
                          start = bin.start, end = bin.ends, stringsAsFactors = F)
  
  region.ranges <- GRanges(seqnames = region.df$chrom,
                           ranges = IRanges(region.df$start, region.df$end))
  
  # adjust score ranges such that they fall within, not onto, boarders of genomic ranges
  #score.ranges
  g1.ranges <- GRanges(seqnames = input.targeting.info.g1$chrom,
                          ranges = IRanges(input.targeting.info.g1$start+ 1, input.targeting.info.g1$end-1))
  g2.ranges <- GRanges(seqnames = input.targeting.info.g2$chrom,
                          ranges = IRanges(input.targeting.info.g2$start+ 1, input.targeting.info.g2$end-1))
  
  g1.overlaps <- as.data.frame(findOverlaps(region.ranges, g1.ranges, type = 'any'))
  g2.overlaps <- as.data.frame(findOverlaps(region.ranges, g2.ranges, type = 'any'))
  
  rows.to.keep <- unique(sort(c(g1.overlaps$queryHits, g2.overlaps$queryHits)))
  
  region.df.filtered <- region.df[rows.to.keep,]
  region.ranges.filtered <- region.ranges[rows.to.keep]
  
  g1.overlaps.filtered <- as.data.frame(findOverlaps(region.ranges.filtered, g1.ranges, type = 'any'))
  g2.overlaps.filtered <- as.data.frame(findOverlaps(region.ranges.filtered, g2.ranges, type = 'any'))
  
  if(length(unique(c(g1.overlaps.filtered$queryHits, g2.overlaps.filtered$queryHits))) != length(rows.to.keep)){
    print('removing segments not overlapping guides failed!')
    break
  }
  if(max(c(g1.overlaps.filtered$queryHits, g2.overlaps.filtered$queryHits)) != nrow(region.df.filtered)){
    print('removing segments not overlapping guides failed! Unequal max nr.')
    break
  }
  
  # combine the per-bin overlaps of both guides
  guide.overlaps.filtered <- rbind(g1.overlaps.filtered, g2.overlaps.filtered)
  region.overlap.list <- split(guide.overlaps.filtered, guide.overlaps.filtered$queryHits)
  seg.to.guide.list <- generate_seg_to_guide_list(region.overlap.list, c(1:nrow(input.targeting.info.g1)))
  
  g1.to.seg.overlaps <- as.data.frame(findOverlaps(g1.ranges, region.ranges.filtered, type = 'any'))
  g2.to.seg.overlaps <- as.data.frame(findOverlaps(g2.ranges, region.ranges.filtered, type = 'any'))
  
  guide.to.seg.overlaps <- rbind(g1.to.seg.overlaps, g2.to.seg.overlaps)
  guide.overlap.list <- split(guide.to.seg.overlaps, guide.to.seg.overlaps$queryHits)
  guide.to.seg.list <- generate_guide_to_seg_list(guide.overlap.list)
  
  # set the labels for each segment
  region.labels <- unlist(lapply(region.overlap.list, function(x){
    temp.label <- input.targeting.info.g1$label[x$subjectHits]
    temp.labels.present <- which(input.label.hierarchy %in% temp.label)
    out.label <- input.label.hierarchy[max(temp.labels.present)]
    return(out.label)
  }))
  
  region.df.filtered$label <- region.labels
  
  
  out.list <- list(seg_info = region.df.filtered,
                   guide_to_seg_lst = guide.to.seg.list,
                   seg_to_guide_lst = seg.to.guide.list)
  
  return(out.list)
}



#' @title Adapt an info file to a overlap matrix and corresponding data structures and filter the counts, return by replicate
#' @param input.counts: data.frame, rows are guides, columns are pools
#' @param replicate.list: each element is a replicate
#' @param input.info: data frame: each row contains chrom, start, end, a label
#' @param input.label.hierarchy: ordering lof guide labels, left to right signifies ordering hor hierarchy to override in case of multiple labels overlapping
#' @param min.seg.dist: minimum distance of a segment. Default = 100
#' @return list: seg_info, guide_to_seg_lst, seg_to_guide_lst, counts (list, per replicate)
#' @export adapt_data_to_regionFormat()

adapt_data_to_regionFormat <- function(input.counts, replicate.list, input.info, input.label.hierarchy, min.seg.dist = 100){

  info.targeting.df <- c()
  counts.targeting.df <- c()
  if(length(which(is.na(input.info$start))) > 0){
    print('Error: data was not properly filtered!')
    break()
  }

  info.targeting.df <- input.info
  counts.targeting.df <- input.counts

  # first create guide-segment matrix for targeting guides: seg_info, guide_to_seg_lst, seg_to_guide_lst
  targeting.gs.list <- create_targeting_guide_segment_matrix(info.targeting.df, input.label.hierarchy, min.seg.dist)

  targeting.gs.list$counts <- list()
  for(i in 1:length(replicate.list)){
    temp.counts <- counts.targeting.df[, replicate.list[[i]]]
    colnames(temp.counts) <- paste('pool', c(1:ncol(temp.counts)), sep = '_')
    temp.counts$n <- rowSums(temp.counts)
    targeting.gs.list$counts[[i]] <- temp.counts
  }

  return(targeting.gs.list)

}


#' @title Generate the guide to segment and segment to guide mappings along with the segment info file
#' @note Modified by ont trimming guides at start and end but only end. Avoid having information-less segments
#' @param input.targeting.info: data.frame, start, end, chrom, label
#' @param input.label.hierarchy: order in which labels are assigned to segments
#' @param min.seg.dist: minimum distance of a segment. Default = 100
#' @return list: seg_info, guide_to_seg_lst, seg_to_guide_lst
#' @export create_targeting_guide_segment_matrix()

create_targeting_guide_segment_matrix <- function(input.targeting.info, input.label.hierarchy, min.seg.dist = 100){

  all.breaks <- unique(c(input.targeting.info$start, input.targeting.info$end))
  sorted.breaks <- sort(all.breaks)

  bin.start <- c() #sorted.breaks[1]
  bin.ends <- c()

  break.idx <- 1
  next.counter <- 0
  proxy.start <- sorted.breaks[1]
  while(break.idx < length(sorted.breaks)){

    next.counter <- 1
    while(sorted.breaks[break.idx + next.counter] - proxy.start < min.seg.dist){
      next.counter <- next.counter + 1
      if(break.idx + next.counter > length(sorted.breaks)){
        bin.start <- c(bin.start, proxy.start)
        bin.ends <- c(bin.ends, sorted.breaks[break.idx + next.counter - 1])
        break.idx <- break.idx + next.counter
        break()
      }
    }
    if(break.idx + next.counter > length(sorted.breaks)){
      break()
    }

    bin.start <- c(bin.start, proxy.start)
    bin.ends <- c(bin.ends, sorted.breaks[break.idx + next.counter])

    proxy.start <- sorted.breaks[break.idx + next.counter] + 1
    break.idx <- break.idx + next.counter

  }

  region.df <- data.frame(chrom = rep(input.targeting.info$chrom[1], length(bin.start)),
                          start = bin.start, end = bin.ends, stringsAsFactors = F)

  region.ranges <- GRanges(seqnames = region.df$chrom,
                           ranges = IRanges(region.df$start, region.df$end))

  # adjust score ranges such that they fall within, not onto, boarders of genomic ranges
  score.ranges <- GRanges(seqnames = input.targeting.info$chrom,
                          ranges = IRanges(input.targeting.info$start+ 1, input.targeting.info$end-1))

  score.overlaps <- as.data.frame(findOverlaps(region.ranges, score.ranges, type = 'any'))

  rows.to.keep <- unique(score.overlaps$queryHits)

  region.df.filtered <- region.df[rows.to.keep,]
  region.ranges.filtered <- region.ranges[rows.to.keep]

  score.overlaps.filtered <- as.data.frame(findOverlaps(region.ranges.filtered, score.ranges, type = 'any'))

  if(length(unique(score.overlaps.filtered$queryHits)) != length(rows.to.keep)){
    print('removing segments not overlapping guides failed!')
    break
  }
  if(max(score.overlaps.filtered$queryHits) != nrow(region.df.filtered)){
    print('removing segments not overlapping guides failed! Unequal max nr.')
    break
  }

  region.overlap.list <- split(score.overlaps.filtered, score.overlaps.filtered$queryHits)
  seg.to.guide.list <- generate_seg_to_guide_list(region.overlap.list, c(1:nrow(input.targeting.info)))

  guide.to.seg.overlaps <- as.data.frame(findOverlaps(score.ranges, region.ranges.filtered, type = 'any'))
  guide.overlap.list <- split(guide.to.seg.overlaps, guide.to.seg.overlaps$queryHits)
  guide.to.seg.list <- generate_guide_to_seg_list(guide.overlap.list)

  # set the labels for each segment
  region.labels <- unlist(lapply(region.overlap.list, function(x){
    temp.label <- input.targeting.info$label[x$subjectHits]
    temp.labels.present <- which(input.label.hierarchy %in% temp.label)
    out.label <- input.label.hierarchy[max(temp.labels.present)]
    return(out.label)
  }))

  region.df.filtered$label <- region.labels


  out.list <- list(seg_info = region.df.filtered,
                   guide_to_seg_lst = guide.to.seg.list,
                   seg_to_guide_lst = seg.to.guide.list)

  return(out.list)
}


#' @title Establish the mapping from segments to guides
#' @param input.seg.overlaps: list: each element is the segment, containing the indexes of the guides with which it overlaps
#' @param guide.idx: vector equivalent to the total number of guides
#' @return list. each element corresponds to a segment, containg a list: guide_idx, nonGuide_idx
#' @export generate_seg_to_guide_list()

generate_seg_to_guide_list <- function(input.seg.overlaps, guide.idx){

  out.seg.to.guide.lst <- lapply(input.seg.overlaps, function(x){

    seg.lst <- list(guide_idx = x$subjectHits,
                    nonGuide_idx = guide.idx[-x$subjectHits])
  })

  return(out.seg.to.guide.lst)
}


#' @title Establish the mapping from guides to segments
#' @param input.score.overlaps: list: each element is the guide, containing the indexes of the segments with which it overlaps
#' @return list. each element is a guide, contaiing a vecotr of indices, indexing the segemnts overlapped by it
#' @export generate_guide_to_seg_list()

generate_guide_to_seg_list <- function(input.guide.overlaps){

  out.guide.to.seg.lst <- lapply(input.guide.overlaps, function(x){
    x$subjectHits
  })

  return(out.guide.to.seg.lst)
}


#' @title given a segment, what guide indeces are added compared to the previous one
#' @param input.seg.to.guide.list: list: each element is the segment, containing the indexes of the guides with which it overlaps and doesn't overlap
#' @return list. each element corresponds to a segment, containg a list: next_guide_idx
#' @export generate_next_guide_list()

generate_next_guide_list <- function(input.seg.to.guide.list){

  out.next.guide.list<- list()
  out.next.guide.list[[1]] <- list()
  out.next.guide.list[[1]]$next_guide_idx <- c()

  for(i in 2:length(input.seg.to.guide.list)){
    prev.idx <- input.seg.to.guide.list[[i-1]]$guide_idx
    cur.idx <- input.seg.to.guide.list[[i]]$guide_idx
    out.next.guide.list[[i]] <- list()
    out.next.guide.list[[i]]$next_guide_idx<- cur.idx[which(! cur.idx %in% prev.idx)]
  }

  return(out.next.guide.list)

}


#' @title extract hyper parameters given a label for multiple replicates. Deprecated. Uses the old setup of individually estimating alpha0 and alpha1
#' @param input.counts.loc: string location of the count file
#' @param input.info.loc: string location of the info file
#' @param input.guide.offset: range by which each guide has an effect beyon the target site (default = 500)
#' @param input.repl.pools: list, each element is a set of columns which correspond to a replicate
#' @param input.labelHierarchy: ordering lof guide labels, left to right signifies ordering hor hierarchy to override in case of multiple labels overlapping
#' @param fs0.label: label used for generating FS0
#' @param file.save.dir: directory in which to save the used counts, info and segment info files
#' @param save.files: logical, by default, do not save the files used for analysis. However, if a filtering step is introduced, file saving is switched on
#' @param min.seg.dist: minimum distance of a segment. Default = 100
#' @return list: hyper parameter estimates for each replicate per list element, $alpha0, alpha1, and $L
#' @export estimate_hyper_parameters_deprecated()

estimate_hyper_parameters_deprecated <- function(analysis.parameters, data.file.split, input.guide.offset,
                                     input.repl.pools, input.labelHierarchy, fs0.label, min.seg.dist = 100){

  data.par.list <- set_up_RELICS_data(analysis.parameters,
                                      data.file.split,
                                      input.guide.offset,
                                      input.repl.pools,
                                      input.labelHierarchy,
                                      fs0.label,
                                      min.seg.dist)
  
  # most basic implementation of GE
  if(! is.null(analysis.parameters$guide_efficiency_scores)){
    
    data.par.list$guide_efficiency <- 1 / (1 + exp(-(data.par.list$guide_efficiency_scores %*% rep(1, ncol(data.par.list$guide_efficiency_scores)))))

  }

  # # if guide efficiency scores are provided, calculate guide efficiency and include in the model
  # if(! is.null(analysis.parameters$guide_efficiency_scores)){
  #   
  #   # estimate_ge_alphaDiffScale should work for one-time calculation
  #   # same as in 'RELICS', assumption is that data.par.list$guide_efficiency_scores is n x 1
  #   if(analysis.parameters$fix_guideEfficiency & ! analysis.parameters$estimate_ge_alphaDiffScale){ 
  #     data.par.list$guide_efficiency <- data.par.list$guide_efficiency_scores[,1]
  #   } else {
  #     # ge_beta_estimation
  #     print('non-fixed guide efficiency or reestimation of alpha diff scale is not yet implemented')
  #     
  #     # if scaling is calculated, then immediately print it!
  #   }
  # }

  # make sure all guides are considered to be the same category
  dirichlet.guide.ll <- compute_perGuide_fs_ll(rep(1, nrow(data.par.list$seg_info)), data.par.list$guide_to_seg_lst, hyper.setup = TRUE)

  final.alpha <- list()

  for(i in 1:length(data.par.list$data)){
    temp.drch.hypers <- list(alpha0 = rep(0.5, length(input.repl.pools[[i]])),
                             alpha1 = rep(0.5, length(input.repl.pools[[i]])))

    temp.hyper.params <- c(sqrt(temp.drch.hypers$alpha0), sqrt(temp.drch.hypers$alpha1))
    temp.alpha0.idx <- c(1:length(temp.drch.hypers$alpha0))
    temp.alpha1.idx <- c(1:length(temp.drch.hypers$alpha0)) + max(temp.alpha0.idx)
    
    temp.res.drch <- c()
    
    if(analysis.parameters$one_dispersion){
      temp.res.drch <- optim(temp.hyper.params, prior_dirichlet_ll_singleDisp, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
                             data = data.par.list$data[[i]], region.ll.list = dirichlet.guide.ll,
                             alpha0.idx = temp.alpha0.idx, alpha1.idx = temp.alpha1.idx, 
                             guide.efficiency = data.par.list$guide_efficiency)
      
      temp.alpha0s <- temp.res.drch$par[temp.alpha0.idx]**2
      temp.alpha1s <- temp.res.drch$par[temp.alpha1.idx]**2
      
      temp.alpha1s.norm <- temp.alpha1s / sum(temp.alpha1s)
      temp.alpha1s.adj <- temp.alpha1s.norm * sum(temp.alpha0s)

      final.alpha[[i]] <- round(temp.alpha1s.adj, 3)
    } else {
      temp.res.drch <- optim(temp.hyper.params, prior_dirichlet_ll, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
                             data = data.par.list$data[[i]], region.ll.list = dirichlet.guide.ll,
                             alpha0.idx = temp.alpha0.idx, alpha1.idx = temp.alpha1.idx, 
                             guide.efficiency = data.par.list$guide_efficiency)
      
      final.alpha[[i]] <- round(temp.res.drch$par[temp.alpha1.idx]**2, 3)
    }


  }

  return(final.alpha)

}


#' @title extract hyper parameters given a label for multiple replicates
#' @param data.par.list: list, contains the elements of the processed and formatted data
#' @param analysis.par.list: list, contains all analysis flags
#' @param input.repl.pools: list, each element is a set of columns which correspond to a replicate
#' @param fs0.label: label used for generating FS0
#' @param one.dispersion: logical, whether analysis is happening with one ro two dispersions
#' @return list: hyper parameter estimates for each replicate per list element, $alpha0, alpha1, and $L
#' @export estimate_hyper_parameters()

estimate_hyper_parameters <- function(data.par.list, analysis.par.list, input.repl.pools,  fs0.label, one.dispersion){

  # most basic implementation of GE
  if(! is.null(data.par.list$guide_efficiency_scores)){
    if(! analysis.par.list$fixed_ge_coeff){
      #data.par.list$guide_efficiency <- 1 / (1 + exp(-(data.par.list$guide_efficiency_scores %*% rep(1, ncol(data.par.list$guide_efficiency_scores)))))
      
      logit.ge <- apply(data.par.list$guide_efficiency_scores, 2, function(x){
        log(x / (1 - x))
      })
      data.par.list$guide_efficiency <- 1 / (1 + exp(-(logit.ge %*% rep(1, ncol(data.par.list$guide_efficiency_scores)))))
    } else {
      data.par.list$guide_efficiency <- data.par.list$guide_efficiency_scores
    }
  }
  
  # make sure all guides are considered to be the same category
  fs.assignment <- rep(0, nrow(data.par.list$seg_info))
  fs.assignment[which(data.par.list$seg_info$label %in% fs0.label)] <- 1
  
  dirichlet.guide.ll <- compute_perGuide_fs_ll(fs.assignment, data.par.list$guide_to_seg_lst, hyper.setup = TRUE)
  
  final.alpha <- list()
  final.alpha$alpha0 <- list()
  final.alpha$alpha1 <- list()
  
  for(i in 1:length(data.par.list$data)){
    temp.drch.hypers <- list(alpha0 = rep(1 / length(input.repl.pools[[i]]), length(input.repl.pools[[i]])),
                             alpha1 = rep(1 / length(input.repl.pools[[i]]), length(input.repl.pools[[i]])))
    
    temp.hyper.params <- c(sqrt(temp.drch.hypers$alpha0), sqrt(temp.drch.hypers$alpha1))
    temp.alpha0.idx <- c(1:length(temp.drch.hypers$alpha0))
    temp.alpha1.idx <- c(1:length(temp.drch.hypers$alpha0)) + max(temp.alpha0.idx)
    
    temp.res.drch <- c()
    
    if(one.dispersion){
      temp.res.drch <- optim(temp.hyper.params, prior_dirichlet_ll_singleDisp, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
                             data = data.par.list$data[[i]], region.ll.list = dirichlet.guide.ll,
                             alpha0.idx = temp.alpha0.idx, alpha1.idx = temp.alpha1.idx, 
                             guide.efficiency = data.par.list$guide_efficiency)
      
      temp.alpha0s <- temp.res.drch$par[temp.alpha0.idx]**2
      temp.alpha1s <- temp.res.drch$par[temp.alpha1.idx]**2
      
      temp.alpha1s.norm <- temp.alpha1s / sum(temp.alpha1s)
      temp.alpha1s.adj <- temp.alpha1s.norm * sum(temp.alpha0s)
      
      final.alpha$alpha0[[i]] <- temp.alpha0s
      final.alpha$alpha1[[i]] <- temp.alpha1s.adj
      
      # final.alpha[[i]] <- round(temp.alpha1s.adj, 3)
    } else {
      temp.res.drch <- optim(temp.hyper.params, prior_dirichlet_ll, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
                             data = data.par.list$data[[i]], region.ll.list = dirichlet.guide.ll,
                             alpha0.idx = temp.alpha0.idx, alpha1.idx = temp.alpha1.idx, 
                             guide.efficiency = data.par.list$guide_efficiency)
      
      temp.alpha0s <- temp.res.drch$par[temp.alpha0.idx]**2
      temp.alpha1s <- temp.res.drch$par[temp.alpha1.idx]**2
      
      final.alpha$alpha0[[i]] <- temp.alpha0s
      final.alpha$alpha1[[i]] <- temp.alpha1s
      
    }
    
    
  }
  
  final.alpha$L <- 1
  
  return(final.alpha)
  
}


#' @title optimize the hyper parameters, only one dispersion across the two distributions
#' @param hyper.param: hyperparameters
#' @param data: data, consists of: $pool1, pool2... and $n
#' @param region.ll.list: list containing the ll and the indexes of the null alternative and both
#' @param alpha0.idx: positions in the par vector of the null alphas, using this dispersion
#' @param alpha1.idx: positions in the par vector of the alternative alphas
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @return sum of the -log likelihood across all guides
#' @export prior_dirichlet_ll_singleDisp()

prior_dirichlet_ll_singleDisp <- function(hyper.param, data, region.ll.list, alpha0.idx, alpha1.idx, guide.efficiency) {
  
  alpha0s <- hyper.param[alpha0.idx]**2
  alpha1s <- hyper.param[alpha1.idx]**2
  
  alpha1s.norm <- alpha1s / sum(alpha1s)
  alpha1s.adj <- alpha1s.norm * sum(alpha0s)
  
  # cat(u0, u1, v0, v1, "\n")
  hyper <- list(alpha0 = alpha0s,
                alpha1 = alpha1s.adj)
  
  out.sg.ll <- estimate_relics_sgrna_log_like(hyper, data, region.ll.list, guide.efficiency)
  
  -sum(out.sg.ll$total_guide_ll)
}


#' @title optimize the hyper parameters, each distribution with it's own variance
#' @param hyper.param: hyperparameters
#' @param data: data, consists of: $pool1, pool2... and $n
#' @param region.ll.list: list containing the ll and the indexes of the null alternative and both
#' @param alpha0.idx: positions in the par vector of the null alphas
#' @param alpha1.idx: positions in the par vector of the alternative alphas
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @return sum of the -log likelihood across all guides
#' @export prior_dirichlet_ll()

prior_dirichlet_ll <- function(hyper.param, data, region.ll.list, alpha0.idx, alpha1.idx, guide.efficiency) {

  alpha0s <- hyper.param[alpha0.idx]**2
  alpha1s <- hyper.param[alpha1.idx]**2

  # cat(u0, u1, v0, v1, "\n")
  hyper <- list(alpha0 = alpha0s,
                alpha1 = alpha1s)

  out.sg.ll <- estimate_relics_sgrna_log_like(hyper, data, region.ll.list, guide.efficiency)
  
  -sum(out.sg.ll$total_guide_ll)
  
}


#' @title Compute log likelihoods of observed counts for each sgRNA, given. Formerly 'compute_dirichlet_sgrna_log_like'
#' @param hyper: hyperparameters, $alpha0, alpha1
#' @param data: data, consists of: pool1, pool2... and $n
#' @param region.ll.list: list containing the ll and the indexes of the null alternative and both
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @return data frame: total_guide_ll, alt_only_ll
#' @export estimate_relics_sgrna_log_like()

estimate_relics_sgrna_log_like <- function(hyper, data, region.ll.list, guide.efficiency, return.model.ll = FALSE){

  pool.cols <- c(1:(ncol(data) - 1))

  sgrna.null.log.like <- ddirmnom(data[,pool.cols], size = data[,ncol(data)], alpha = hyper$alpha0, log = T)
  
  sgrna.alt.log.like <- c()
  if(is.null(guide.efficiency)){
    sgrna.alt.log.like <- ddirmnom(data[,pool.cols], size = data[,ncol(data)], alpha = hyper$alpha1, log = T)
  } else {
    alpha0.matrix <- t(apply(matrix(0, ncol = length(hyper$alpha0), nrow = length(guide.efficiency)), 1, function(x){x + hyper$alpha0}))
    alpha.diffs <- hyper$alpha0 - hyper$alpha1
    
    alpha1.matrix <- matrix(0, ncol = length(hyper$alpha0), nrow = length(guide.efficiency))
    
    for(i in 1:length(guide.efficiency)){
      alpha1.matrix[i,] <- alpha0.matrix[i,] - alpha.diffs * guide.efficiency[i]
    }
    
    sgrna.alt.log.like <- ddirmnom(data[,pool.cols], size = data[,ncol(data)], alpha = alpha1.matrix, log = T)
  }
  
  if(return.model.ll){
    return(data.frame(null_only_ll = sgrna.null.log.like, alt_only_ll = sgrna.alt.log.like))
  }

  out.comb <- vector('numeric', length = length(sgrna.null.log.like))

  out.comb[region.ll.list$null_only_idx] <- sgrna.null.log.like[region.ll.list$null_only_idx] + region.ll.list$null_only_lls
  out.comb[region.ll.list$alt_only_idx] <- sgrna.alt.log.like[region.ll.list$alt_only_idx] + region.ll.list$alt_only_lls

  comb.df <- data.frame(null_guide = sgrna.null.log.like[region.ll.list$both_idx],
                        null_ll = region.ll.list$both_lls[,1],
                        alt_guide = sgrna.alt.log.like[region.ll.list$both_idx],
                        alt_ll = region.ll.list$both_lls[,2])
  comb.lls <- apply(comb.df, 1, function(x){addlogs(x[1] + x[2], x[3] + x[4])})

  out.comb[region.ll.list$both_idx] <- comb.lls
  
  # output is now a data frame, where the first column is the total guide ll
  # the second one is the alternate ll (FS ll)
  out.df <- data.frame(total_guide_ll = out.comb, alt_only_ll = sgrna.alt.log.like)

  return(out.df)

}


#' @title Compute log likelihoods of each sgRNA to overlap a FS. Formerly 'compute_region_dpoibin_fast1'
#' @param cumulative.pp: sum of all previous posteriors
#' @param guide.seg.idx.lst: list, each element is a guide and contains the indexes of the segments it overlaps
#' @param hyper.setup: logical, if this function is called during the hyper parameter setup or not.
#' @return list: null_only_idx, null_only_lls, alt_only_idx, alt_only_lls, both_idx, both_lls (last one is a df)
#' @export compute_perGuide_fs_ll()

compute_perGuide_fs_ll <- function(cumulative.pp, guide.seg.idx.lst, hyper.setup = FALSE){

  n.sgrna <- length(guide.seg.idx.lst)

  null.ll <- vector('numeric', length = n.sgrna)
  alt.ll <- vector('numeric', length = n.sgrna)
  null.indicator <- vector('numeric', length = n.sgrna)
  alt.indicator <- vector('numeric', length = n.sgrna)

  null.only.idx <- vector('numeric', length = n.sgrna)
  null.only.counter <- 1
  alt.only.idx <- vector('numeric', length = n.sgrna)
  alt.only.counter <- 1
  both.idx <- vector('numeric', length = n.sgrna)
  both.counter <- 1


  for(j in 1:n.sgrna) {
    region.pp <- cumulative.pp[guide.seg.idx.lst[[j]] ]

    # compute probability that number of regulatory regions overlapped by this sgRNA is 0 or >0, given posterior
    # probabilities, using poisson binomial probability mass function
    p.k.eq.0 <- dpoibin(0, region.pp)
    p.k.gt.0 <- 1-p.k.eq.0

    null.ll[j] <- log(p.k.eq.0)
    alt.ll[j] <- log(p.k.gt.0)

    if(p.k.eq.0 >= 1.0) {
      null.only.idx[null.only.counter] <- j
      null.only.counter <- null.only.counter + 1
    }
    else if(p.k.gt.0 >= 1.0) {
      alt.only.idx[alt.only.counter] <- j
      alt.only.counter <- alt.only.counter + 1
    } else {
      both.idx[both.counter] <- j
      both.counter <- both.counter + 1
    }
  }

  null.idxs <- null.only.idx[1:(null.only.counter - 1)]
  alt.idxs <-alt.only.idx[1:(alt.only.counter - 1)]
  both.idxs <- both.idx[1:(both.counter - 1)]

  # if there are overlaps between the indeces that's a problem!
  if(sum(c(length(null.idxs), length(alt.idxs), length(both.idxs))) != length(unique(c(null.idxs, alt.idxs, both.idxs))) & !hyper.setup){
    print('Potential issue with guide likelihood assignments!')
  }

  out.list <- list(null_only_idx = null.idxs,
                   null_only_lls = null.ll[null.idxs],
                   alt_only_idx = alt.idxs,
                   alt_only_lls = alt.ll[alt.idxs],
                   both_idx = both.idxs,
                   both_lls = data.frame(null_lls = null.ll[both.idxs], alt_lls = alt.ll[both.idxs]))

  return(out.list)

}

#' @title Run RELICS. Iteratively detect FS until convergence is reached
#' @param input.data: list: $guide_to_seg_lst, $data, $fs0_idx, $overlapMat, $guide_to_seg_lst, $seg_to_guide_lst
#' @param final.layer.nr: number of layers to go up to
#' @param out.dir: directory to which all files are to be written. Default = NULL
#' @param adjust.tol: whether convergence tolerance should be adjusted dynamically as layers are computed (increasingly more rigorous)
#' @param fix.hypers: whether the input parameters remain unchanged or are altered after calculation of posteriors. Default = FALSE
#' @param iterative.hyper.est: whether hyper parameters are to be estimated iteratively with the posteriors (TRUE) or after convergence of posteriors (FALSE)
#' @param input.hypers: list, $alpha0, $alpha1, both elements contain lists of length wequal to nr. replicates (input.hypers$alpha0[[1]])
#' @param nr.segs, number of segments to consider for the length of a regulatory element
#' @param geom.p: proababilty of the genometric distribution to penalize for enhancers of increasing length
#' @param fs.correlation.cutoff: what correlation between layers determines stopping of layer calculation (default: 0.1)
#' @param input.min.rs.pp: minimum posterior required to be part of a regulatory set
#' @param auto.stop: whether or not computations should be stopped after recomended stopping point
#' @param record.all.fs: logical, if information of all intermediate FS should be recorded, in addition to the final set of FS
#' @param one.dispersion: whether there should be one or 2 dispersions estimated for the 2 hyper parameters
#' @param recompute.fs0: logical, whether FS0 is to be set to 0 or keep the initial assignment
#' @param local.max: logical, whether a local maximum should be computed
#' @param local.max.range: number of segments to include in addition to the the ones with highest PP (one way, so multiply by 2 for total nr of segs)
#' @return list of final per-layer posteriors
#' @export run_RELICS_2()

run_RELICS_2 <- function(input.data, final.layer.nr, out.dir = NULL,
                     fix.hypers = FALSE,
                     iterative.hyper.est = FALSE,
                     input.hypers = NULL,
                     nr.segs = 10,
                     geom.p = 0.1,
                     fs.correlation.cutoff = 0.1,
                     input.min.rs.pp = 0.1,
                     auto.stop = TRUE,
                     record.all.fs,
                     input.convergence.tol = 0.01,
                     adjust.tol = TRUE,
                     one.dispersion = TRUE,
                     recompute.fs0,
                     local.max, local.max.range){

  final.layer.posterior <- list()
  final.layer.alpha0 <- list()
  final.layer.alpha1 <- list()
  final.layer.ll <- c()
  final.per.layer.ll <- list()
  layer.corr.df <- c()
  layer.corr.df.list <- list()
  final.all.layers.ll <- list() # record the overal model ll progression
  ge.coeff.list <- list()
  ge.coeff.list[[1]] <- input.data$ge_coeff
  final.fs.ll.rt <- list()

  no.convergence.layer <- c() # during what layers whas there no convergence
  no.convergence.lls <- c() #model ll when convergence was not reached
  no.convergence.bool <- c() # boolean to keep track if convergence was reached

  layer.times <- vector('numeric', length = final.layer.nr)

  relics.hyper <- c()

  if(is.null(input.hypers)){

    print('Hyper parameters must be provided')
    break()

  } else {
    relics.hyper <- input.hypers
  }

  relics.param <- init_relics_param(relics.hyper, input.data, local.max)
  if(recompute.fs0){
    relics.param$delta.pp[1,] <- 0
  }

  coverg.tol <- rep(input.convergence.tol, final.layer.nr)

  if(adjust.tol){

    first.third.stop <- round(0.3 * final.layer.nr)
    second.third.stop <- round(0.7 * final.layer.nr)

    if(length(unique(c(1, first.third.stop, second.third.stop, final.layer.nr))) < 4){
      print('problem with adjutable tolerance!')
    }
    coverg.tol[1:first.third.stop] <- 0.05
    coverg.tol[(first.third.stop + 1):second.third.stop] <- 0.01
    coverg.tol[(second.third.stop + 1): final.layer.nr] <- 0.001
  }

  for(i in 1:final.layer.nr){
    print(paste0('Computing FS: ', i))
    layer.time <- proc.time()
    layer.data <- relics_compute_FS_k(input.param = relics.param,
                                      input.hyper = relics.hyper,
                                      input.data.list = input.data,
                                      input.tol = coverg.tol[i],
                                      fix.hypers, iterative.hyper.est, nr.segs, geom.p,
                                      min.pp = input.min.rs.pp, input.data$guide_efficiency,
                                      one.dispersion,
                                      local.max, local.max.range)

    layer.time.final <- proc.time() - layer.time
    layer.times[i] <- layer.time.final[3] / 60
    layer.size <- length(layer.data$posterior_trace_list)

    layer.posteriors <- layer.data$posterior_trace_list[[layer.size]]
    layer.posteriors[layer.posteriors > 1] <- 1
    
    fs.ll.rt <- list()
    if(local.max){
      fs.ll.rt <- layer.data$fs_ll_rt_trace[[layer.size]]
    }
    
    # order the layers according to their model contributions
    order.pps.lst <- order_pps(layer.posteriors, layer.data$ll_tract[[layer.size]],
                               input.data, layer.data$alpha0[[layer.size]], layer.data$alpha1[[layer.size]], 
                               input.data$guide_efficiency, local.max, fs.ll.rt)

    final.layer.posterior[[i]] <- order.pps.lst$pp_ordered
    final.layer.alpha0[[i]] <- layer.data$alpha0[[layer.size]]
    final.layer.alpha1[[i]] <- layer.data$alpha1[[layer.size]]
    final.layer.ll <- c(final.layer.ll, layer.data$ll_tract[[layer.size]])
    final.fs.ll.rt[[i]] <- order.pps.lst$fs_ll_rt_ordered

    #$nr_rs, $rs_prob, $training_overlap, $prct_overlap
    pp.stat.lst <- pps_stats(order.pps.lst$pp_ordered, input.min.rs.pp)

    # record the per-layer contribution to the model
    per.layer.ll.contribution <- order.pps.lst$layer_ll_ordered
    per.layer.ll.df <- data.frame(FS = c('total_model_ll', paste0('k_', c(1:(length(per.layer.ll.contribution)))) ),
                                  ll = c(layer.data$ll_tract[[layer.size]], per.layer.ll.contribution),
                                  nr_fs = c(0, pp.stat.lst$nr_rs),
                                  stringsAsFactors = F)
    final.per.layer.ll[[i]] <- per.layer.ll.df
    #per.layer.ll.df$layer_Iteration <- rep(as.character(i), nrow(per.layer.ll.df))


    # prep the parameters for the next iteration (both the posteriors and hypers)
    relics.hyper$L <- relics.hyper$L + 1
    relics.param$delta.pp <- rbind(order.pps.lst$pp_ordered, rep(0, ncol(relics.param$delta.pp)))
    relics.hyper$alpha0 <- layer.data$alpha0[[layer.size]]
    relics.hyper$alpha1 <- layer.data$alpha1[[layer.size]]
    
    if(local.max){
      relics.param$ll_rt <- rbind(order.pps.lst$fs_ll_rt_ordered, rep(0, ncol(relics.param$delta.pp)))
    }

    # if there was an issue with convergence, then record it
    no.convergence.bool <- c(no.convergence.bool, layer.data$max_iter_reached)
    if(layer.data$max_iter_reached){
      no.convergence.layer <- c(no.convergence.layer, i) # during what layers whas there no convergence
      no.convergence.lls <- c(no.convergence.lls, layer.data$ll_tract[[layer.size]])
    }

    # # record the PP correlations and overlaps
    layer.corrs <- pps_corr(order.pps.lst$pp_ordered, i)
    layer.corr.df <- rbind(layer.corr.df, layer.corrs)
    layer.corr.df.list[[i]] <- layer.corr.df

    # record sum of posteriors as bedgraph
    total.per.seg.posterior <- colSums(order.pps.lst$pp_ordered)
    total.per.seg.posterior[which(total.per.seg.posterior > 1)] <- 1
    out.bedgraph <- input.data$seg_info
    out.bedgraph$score <- total.per.seg.posterior
    to.bg.list <- list(total_pp = out.bedgraph)
    
    if(local.max){
      to.bg.list$total_ll <- input.data$seg_info
      to.bg.list$total_ll$score <- colSums(order.pps.lst$fs_ll_rt_ordered)
    }

    # record the segments with FS probabilities above the threshold.
    all.seg.fs.df <- extract_fs_locations(order.pps.lst$pp_ordered, input.data$seg_info, input.min.rs.pp)

    # record the overal model ll progression
    all.layers.ll <- data.frame(FS = c(1:length(final.layer.ll)),
                                FS_ll = final.layer.ll,
                                stringsAsFactors = F)
    final.all.layers.ll[[i]] <- all.layers.ll

    if(record.all.fs){
      # record posteriors
      write.csv(order.pps.lst$pp_ordered, file = paste0(out.dir, '_k',i,'.csv'), row.names = F)
      
      if(local.max){
        write.csv(order.pps.lst$fs_ll_rt_ordered, file = paste0(out.dir, '_k',i,'_llRt.csv'), row.names = F)
      }

      # record bedgraph
      create_bedgraphs(to.bg.list, paste0(out.dir, '_k',i) )

      # save FS locations
      write.table(all.seg.fs.df, file = paste0(out.dir, '_k',i,'_FS_locations.bed'),
                  sep = '\t', quote = F, row.names = F, col.names = F)

      # record alphas used for posterior calculation
      record_alphas(final.layer.alpha0[[i]], final.layer.alpha1[[i]], out.dir, i)

      # record the per-layer contribution to the model
      write.csv(per.layer.ll.df, file = paste0(out.dir, '_k',i,'_perFS_LLcontributions.csv'), row.names = F)

      # model ll progression
      write.csv(all.layers.ll, file = paste0(out.dir, '_k',i,'_ll_progression.csv'), row.names = F)

      # record the per-layer correlations
      write.layer.corrs <- data.frame(FSpair_1 = layer.corrs$layer_1, FSpair_2 = layer.corrs$layer_2,
                                      corr = layer.corrs$corr, nr_FS_compared = layer.corrs$layer_Iteration, stringsAsFactors = F)
      write.csv(write.layer.corrs, file = paste0(out.dir, '_k',i,'_correlations.csv'), row.names = F)
      
      # guide efficiency vars:
      if(! is.null(input.data$fixed_ge_coeff) && ! input.data$fixed_ge_coeff){
        ge.coff.df <- data.frame(ge_coeff = c('beta0', paste0('beta', 1:ncol(input.data$guide_efficiency_scores))), 
                                 ge_coeff_scores = round(input.data$ge_coeff, 3))
        write.csv(ge.coff.df, file = paste0(out.dir, '_k',i,'_ge_coeff.csv'), row.names = F)
      }
      
      

    }

    # plot the outputs
    if(i > 1){

      if(record.all.fs){
        display_relics_fs_as_tiff(final.layer.posterior[[i]],
                               input.data$segLabels,
                               paste0(out.dir, '_k_', i),
                               input.min.rs.pp)

        # plot the ll progression and mark layers where convergence failed
        # plot the correlation between PP
        plot_fs_stats(all.layers.ll, layer.corr.df, out.dir, i, fs.correlation.cutoff)
      }

      fs.correlation.cutoff.list <- determine_FS_nr_cutoff(layer.corrs, fs.correlation.cutoff)

      if(fs.correlation.cutoff.list$need_to_stop){
        print(fs.correlation.cutoff.list$why_to_stop)

        final.pp.out <- final.layer.posterior[[i - 1]]
        all.layers.ll.out <- final.all.layers.ll[[i - 1]]
        final.per.layer.ll.out <- final.per.layer.ll[[i - 1]]
        final.ll.rt.out <- final.fs.ll.rt[[i - 1]]

        # record sum of posteriors as bedgraph
        total.per.seg.posterior <- colSums(final.pp.out)
        total.per.seg.posterior[which(total.per.seg.posterior > 1)] <- 1
        out.bedgraph <- input.data$seg_info
        out.bedgraph$score <- total.per.seg.posterior
        to.bg.list <- list(total_pp = out.bedgraph)
        
        if(local.max){
          to.bg.list$total_ll <- input.data$seg_info
          to.bg.list$total_ll$score <- colSums(final.ll.rt.out)
        }

        # record the segments with FS probabilities above the threshold.
        all.seg.fs.df.final <- extract_fs_locations(final.pp.out, input.data$seg_info, input.min.rs.pp)

        if(auto.stop){
          # record bedgraph
          create_bedgraphs(to.bg.list, paste0(out.dir, '_final_k', i - 1) )

          # save FS locations
          write.table(all.seg.fs.df.final, file = paste0(out.dir, '_final_k', i - 1,'_FS_locations.bed'),
                      sep = '\t', quote = F, row.names = F, col.names = F)

          # record alphas used for posterior calculation
          record_alphas(final.layer.alpha0[[i - 1]], final.layer.alpha1[[i - 1]], paste0(out.dir, '_final'), i - 1)

          # record the per-layer contribution to the model
          write.csv(final.per.layer.ll.out, file = paste0(out.dir, '_final_k',i - 1,'_perFS_LLcontributions.csv'), row.names = F)

          # model ll progression
          write.csv(all.layers.ll.out, file = paste0(out.dir, '_final_k',i - 1,'_ll_progression.csv'), row.names = F)

          # record posteriors
          write.csv(final.pp.out, file = paste0(out.dir, '_final_k', i - 1,'.csv'), row.names = F)
          
          if(local.max){
            write.csv(final.ll.rt.out, file = paste0(out.dir, 'final_k', i-1, '_llRt.csv'), row.names = F)
          }

          # record all correlations
          write.layer.corrs.df <- data.frame(FSpair_1 = layer.corr.df$layer_1, FSpair_2 = layer.corr.df$layer_2,
                                          corr = layer.corr.df$corr, nr_FS_compared = layer.corr.df$layer_Iteration, stringsAsFactors = F)
          write.csv(write.layer.corrs.df, file = paste0(out.dir, '_final_allCorrs_k', i,'.csv'), row.names = F)

          display_relics_fs_as_tiff(final.layer.posterior[[i]],
                                    input.data$segLabels,
                                    paste0(out.dir, '_atFinal_k', i),
                                    input.min.rs.pp)

          display_relics_fs_as_tiff(final.pp.out,
                                    input.data$segLabels,
                                    paste0(out.dir, '_final_k', i-1),
                                    input.min.rs.pp)

          # # plot the ll progression and mark layers where convergence failed
          # # plot the correlation between PP
          plot_fs_stats(all.layers.ll.out, layer.corr.df.list[[i-1]], paste0(out.dir, '_final'), i-1, fs.correlation.cutoff)
          #
          # also plot summary at the stopping point so see how the correlation changed
          plot_fs_stats(all.layers.ll, layer.corr.df, paste0(out.dir, '_atCutoff'), i, fs.correlation.cutoff)
          
          # guide efficiency vars:
          if(! is.null(input.data$fixed_ge_coeff) && ! input.data$fixed_ge_coeff){
            ge.coff.df <- data.frame(ge_coeff = c('beta0', paste0('beta', 1:ncol(input.data$guide_efficiency_scores))), 
                                     ge_coeff_scores = round(ge.coeff.list[[i - 1]], 3))
            write.csv(ge.coff.df, file = paste0(out.dir, '_final_k', i-1, '_ge_coeff.csv'), row.names = F)
          }

          break()
        } else {
          print(paste0('Recommend stopping RELICS 2'))
          print(fs.correlation.cutoff.list$why_to_stop)
          print('Continuing to run because auto.stop was set to FALSE')

          # record bedgraph
          create_bedgraphs(to.bg.list, paste0(out.dir, '_recommendedFinal_k', i - 1) )

          # save FS locations
          write.table(all.seg.fs.df.final, file = paste0(out.dir, '_recommendedFinal_k', i - 1,'_FS_locations.bed'),
                      sep = '\t', quote = F, row.names = F, col.names = F)

          # record alphas used for posterior calculation
          record_alphas(final.layer.alpha0[[i - 1]], final.layer.alpha1[[i - 1]], paste0(out.dir, '_recommendedFinal'), i - 1)

          # record the per-layer contribution to the model
          write.csv(final.per.layer.ll.out, file = paste0(out.dir, '_recommendedFinal_k',i - 1,'_perFS_LLcontributions.csv'), row.names = F)

          # model ll progression
          write.csv(all.layers.ll.out, file = paste0(out.dir, '_recommendedFinal_k',i - 1,'_ll_progression.csv'), row.names = F)

          # record posteriors
          write.csv(final.pp.out, file = paste0(out.dir, '_recommendedFinal_k', i - 1,'.csv'), row.names = F)
          
          if(local.max){
            write.csv(final.ll.rt.out, file = paste0(out.dir, '_recommendedFinal_k', i-1, '_llRt.csv'), row.names = F)
          }

          # record posteriors
          write.csv(final.pp.out, file = paste0(out.dir, '_recommendedFinal_k', i - 1,'.csv'), row.names = F)

          # record all correlations
          write.layer.corrs.df <- data.frame(FSpair_1 = layer.corr.df$layer_1, FSpair_2 = layer.corr.df$layer_2,
                                             corr = layer.corr.df$corr, nr_FS_compared = layer.corr.df$layer_Iteration, stringsAsFactors = F)
          write.csv(write.layer.corrs.df, file = paste0(out.dir, '_recommendedFinal_allCorrs_k', i,'.csv'), row.names = F)

          display_relics_fs_as_tiff(final.layer.posterior[[i]],
                                    input.data$segLabels,
                                    paste0(out.dir, '_atRecommendedFinal_k', i),
                                    input.min.rs.pp)

          display_relics_fs_as_tiff(final.pp.out,
                                    input.data$segLabels,
                                    paste0(out.dir, '_recommendedFinal_k', i-1),
                                    input.min.rs.pp)

          # plot the ll progression and mark layers where convergence failed
          # plot the correlation between PP
          plot_fs_stats(all.layers.ll.out, layer.corr.df, paste0(out.dir, '_recommendedFinal'), i-1, fs.correlation.cutoff)

          # also plot summary at the stopping point so see how the correlation changed
          plot_fs_stats(all.layers.ll, layer.corr.df, paste0(out.dir, '_atRecommendedCutoff'), i, fs.correlation.cutoff)
          
          # guide efficiency vars:
          if(! is.null(input.data$fixed_ge_coeff) && ! input.data$fixed_ge_coeff){
            ge.coff.df <- data.frame(ge_coeff = c('beta0', paste0('beta', 1:ncol(input.data$guide_efficiency_scores))), 
                                     ge_coeff_scores = round(ge.coeff.list[[i - 1]], 3))
            write.csv(ge.coff.df, file = paste0(out.dir, '_recommendedFinal_k',i,'_ge_coeff.csv'), row.names = F)
          }

        }

      }

    }

    # if the hyper parameters are reestimated after convergence of the posteriors
    if(! fix.hypers & !iterative.hyper.est){
      relics.hyper <- recompute_hyper_parameters(relics.param,
                                                 relics.hyper,
                                                 input.data$data,
                                                 input.data$guide_to_seg_lst,
                                                 input.data$guide_efficiency,
                                                 one.dispersion)
    }
    
    if(! is.null(input.data$guide_efficiency_scores)){
      if(! input.data$fixed_ge_coeff){
        ge.list <- recompute_ge_coefficients(relics.param,
                                             relics.hyper,
                                             input.data$data,
                                             input.data$guide_to_seg_lst,
                                             input.data$guide_efficiency_scores,
                                             input.data$ge_coeff)
        
        input.data$guide_efficiency <- ge.list$guide_efficiency
        input.data$ge_coeff <- ge.list$ge_coeff
        ge.coeff.list[[i + 1]] <- ge.list$ge_coeff
      }

    }

  }

}


#' @title wrapper which returns the updated guide efficiency and the corresponding coefficients used to calculate it
#' @param hyper: hyperparameters
#' @param param: matrix of all posterior probs
#' @param data: data, consists of: $y1, $y2 and $n
#' @param guide.seg.idx.lst: list: each element is a guide, contaiing the indexes of the segments it overlaps
#' @param guide.efficiency.scores: matrix, each column containing a different set of scores per guide
#' @param ge.coeff vector, containing the guide efficiency coefficients
#' @param one.dispersion. logical, if there should be one or two dispersions for the hyper parameters
#' @return list: $guide_efficiency, $ge_coeff
#' @export recompute_ge_coefficients()

recompute_ge_coefficients <- function(param, hyper, data, guide.seg.idx.lst, guide.efficiency.scores, ge.coeff) {
  cumulative.pp <- colSums(param$delta.pp) #apply(param$delta.pp, 2, sum)
  cumulative.pp[cumulative.pp > 1] <- 1
  
  guide.lls.list <- compute_perGuide_fs_ll(cumulative.pp, guide.seg.idx.lst)
  
  ge.coeff.param <- ge.coeff
  
  res <- optim(ge.coeff.param, guide_coeff_ll, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
               data=data, region.ll.list = guide.lls.list,
               alpha0.input = hyper$alpha0, alpha1.input = hyper$alpha1, 
               guide.efficiency.scores = guide.efficiency.scores)
  
  if(res$convergence==0) {
    
    out.list <- list()
    # out.list$guide_efficiency <- 1 / (1 + exp(-(res$par[1] + guide.efficiency.scores %*% res$par[2:length(res$par)])))

    guide.efficiency.scores.logit <- apply(guide.efficiency.scores, 2, function(x){
      log(x / (1-x))
    })
    
    out.list$guide_efficiency <- 1 / (1 + exp(-(res$par[1] + guide.efficiency.scores.logit %*% res$par[2:length(res$par)])))
    out.list$ge_coeff <- res$par
    
  } else {
    warning("estimation of hyperparameters failed to converge")
  }

  return(out.list)
}


#' @title optimize the guide efficiency coefficients
#' @param ge.coeff.param: guide efficiency coefficients (beta0, beta1, ...)
#' @param data: data, consists of: $pool1, pool2... and $n
#' @param region.ll.list: list containing the ll and the indexes of the null alternative and both
#' @param alpha0.input: hyper parameters of the background
#' @param alpha1.input: hyper parameters of the functional sequences
#' @param guide.efficiency.scores: matrix, each column containing a different set of scores per guide
#' @return sum of the -log likelihood across all guides
#' @export guide_coeff_ll()

guide_coeff_ll <- function(ge.coeff.param, data, region.ll.list, alpha0.input, alpha1.input, guide.efficiency.scores) {
  
  #guide.efficiency <- 1 / (1 + exp(-(ge.coeff.param[1] + guide.efficiency.scores %*% ge.coeff.param[2:length(ge.coeff.param)])))
  
  guide.efficiency.scores.logit <- apply(guide.efficiency.scores, 2, function(x){
    log(x / (1-x))
  })
  guide.efficiency <- 1 / (1 + exp(-(ge.coeff.param[1] + guide.efficiency.scores.logit %*% ge.coeff.param[2:length(ge.coeff.param)])))
  
  total.neg.ll <- 0
  
  for(i in 1:length(data)){
    hyper <- list(alpha0 = alpha0.input[[i]],
                  alpha1 = alpha1.input[[i]])
    
    temp.neg.ll <- estimate_relics_sgrna_log_like(hyper, data[[i]], region.ll.list, guide.efficiency)
    
    total.neg.ll <- total.neg.ll + -sum(temp.neg.ll$total_guide_ll)
    
  }
  
  total.neg.ll
}


#' @title Based on correlations, determine of another FS can be calculated; all segments above a cutoff. Formerly 'determine_fs_correlation_cutoff'
#' @param input.fs.pp: data.frame, contains the per FS posteriors
#' @param input.seg.info: location of each segment
#' @param fs.threshold: cutoff for each FS
#' @return data.frame: 'chrom', 'start', 'end', 'label', 'score'
#' @export extract_fs_locations()

extract_fs_locations <- function(input.fs.pp, input.seg.info, fs.threshold){

  if(nrow(input.seg.info) != ncol(input.fs.pp)){
    print('Error! Number of segments and posterior values do not match up!')
    break()
  }

  fs.label <- ''
  fs.df <- c()
  for(i in 1:nrow(input.fs.pp)){
    fs.label <- paste0('FS', i-1)

    fs.idx <- which(input.fs.pp[i,] > fs.threshold)
    
    if(length(fs.idx) > 0){
      fs.scores <- input.fs.pp[i, which(input.fs.pp[i,] > fs.threshold)]
      
      temp.fs.position <- input.seg.info[fs.idx,]
      temp.fs.position$label <- fs.label
      temp.fs.position$score <- fs.scores
      
      fs.df <- rbind(fs.df, temp.fs.position)
    }
  
  }

  out.fs.df <- fs.df[,c('chrom', 'start', 'end', 'label', 'score')]
  
  # zero-index the output
  out.fs.df$start <- out.fs.df$start - 1
  out.fs.df$end <- out.fs.df$end - 1

  return(out.fs.df)
}


#' @title Based on correlations, determine of another FS can be calculated; all segments above a cutoff. Formerly 'determine_fs_correlation_cutoff'
#' @param input.layer.corrs: data frame, max. correlations across all layers: $layer_1, $layer_2, $corr, $layer_comb, $layer_Iteration
#' @param fs.correlation.cutoff: if correlation between layers exceeds fs.correlation.cutoff, recommend stopping
#' @param out.dir: directory to which to write
#' @return list: $need_to_stop, $flags_boolean, $flags_layer, (if stopping: $why_to_stop)
#' @export determine_FS_nr_cutoff()

determine_FS_nr_cutoff <- function(input.layer.corrs, fs.correlation.cutoff){

  need.to.stop <- FALSE
  why.to.stop <- c()

  out.list <- list(need_to_stop = FALSE)

  if(length(which(input.layer.corrs$corr > fs.correlation.cutoff)) > 0){

    all.above.corrs <- which(input.layer.corrs$corr > fs.correlation.cutoff)
    out.list$need_to_stop <- TRUE

    for(i in 1:length(all.above.corrs)){
      why.to.stop <- c(why.to.stop, paste0('Correlation between FS',
                                           input.layer.corrs$layer_1[all.above.corrs[i]],
                                           ' and FS', input.layer.corrs$layer_2[all.above.corrs[i]],
                                           ' exceeds cutoff of ', fs.correlation.cutoff))
    }

    out.list$why_to_stop <- why.to.stop

  }

  return(out.list)

}


#' @title Display the credible sets along with cumulative posteriors. Formerly 'display_mrvr_norm_tiff'
#' @param input.L matrix of posteriors
#' @param input.labels: labelling of all segments
#' @param tiff.name: name for plot
#' @param fs.threshold: threshold for defining functional sequences
#' @return One plot for each selection step, credible sets colored in purple
#' @export display_relics_fs_as_tiff()

display_relics_fs_as_tiff <- function(input.L, input.labels, tiff.name, fs.threshold){

  #cs.mat <- compute_PP_RS(input.L, fs.threshold)

  # attempt to put 8 per page
  if(nrow(input.L) <= 9){
    tiff(paste0(tiff.name, '.tiff'), width = 11520 / 2, height = 5760 / 2, units = "px", res = 400)
    par(mfrow = c(3,3))

    for(i in 1:nrow(input.L)){
      cs.col <- rep("black", ncol(input.L))
      cs.col[input.L[i,] > fs.threshold] <- "purple"

      plot(input.L[i,], pch=21,
           main=paste("FS", i-1, ", Nr. segments > ", fs.threshold, " = ", length(which(input.L[i,] > fs.threshold)), sep=""), col=cs.col,
           ylab = 'PP', xlab = 'Genome Segment')
    }
    p.sum <- c()
    if(nrow(input.L) > 2){
      p.sum <- colSums(input.L[1:nrow(input.L),])
    } else {
      p.sum <- input.L[1:2,]
    }
    p.sum[p.sum > 1] <- 1
    plot(p.sum, pch=21,
         main="Sum of Posteriors",
         col="black", bg = input.labels, ylab = 'PP', xlab = 'Genome Segment')

    dev.off()
  } else {
    nr.pages <- ceiling(nrow(input.L) / 8)

    for(pg in 1:nr.pages){
      tiff(paste0(tiff.name, '_', pg, '.tiff'), width = 11520 / 2, height = 5760 / 2, units = "px", res = 400)
      par(mfrow = c(3,3))

      for(i in 1:8){

        row.index <- (pg - 1)*8 + i

        if(row.index > nrow(input.L)){
          p.sum <- c()
          if(nrow(input.L) > 2){
            p.sum <- colSums(input.L[1:nrow(input.L),])
          } else {
            p.sum <- input.L[1:2,]
          }
          p.sum[p.sum > 1] <- 1
          plot(p.sum, pch=21,
               main="Sum of Posteriors",
               col="black", bg = input.labels, ylab = 'PP', xlab = 'Genome Segment')
          #dev.off()
          break()
        } else{

          cs.col <- rep("black", ncol(input.L))
          cs.col[input.L[row.index,] > fs.threshold] <- "purple"

          # ', RS-pp: ', temp.cs.pp,
          plot(input.L[row.index,], pch=21,
               main=paste("FS", row.index-1, ", Nr. segments > ", fs.threshold, " = ", length(which(input.L[row.index,] > fs.threshold)),sep=""), col=cs.col,
               ylab = 'PP', xlab = 'Genome Segment')
          if(i == 8){
            p.sum <- c()
            if(nrow(input.L) > 2){
              p.sum <- colSums(input.L[1:nrow(input.L),])
            } else {
              p.sum <- input.L[1:2,]
            }
            p.sum[p.sum > 1] <- 1
            plot(p.sum, pch=21,
                 main="Sum of Posteriors",
                 col="black", bg = input.labels, ylab = 'PP', xlab = 'Genome Segment')
          }
        }
      }
      dev.off()
    }

  }

}

#' @title Compute regulatory sets for all segments above a cutoff. Formerly 'compute_mrvr_norm_CS'
#' @param input.L: matrix, nrow = number of selections performed, ncol = number of elements
#' @param threshold: double, max probability of the sum of all elements in this CS
#' @return matrix: for each row, elements with 1 are part of credible set, otherwise not part (0)
#' @export compute_PP_RS()

compute_PP_RS <- function(input.L, threshold){

  L <- nrow(input.L)

  # make credible set matrix with 1s indicating the regions that are part of the
  # credible set for each row
  cs.matrix <- matrix(0, nrow=L, ncol=ncol(input.L))

  # skip first row, which was the predefined
  for(i in 2:L) {
    cs.matrix[i,which(input.L[i,] >= threshold)] <- 1
  }

  return(cs.matrix)

}


#' @title Compute for the credible sets for the segment-length normalized posteriors; all segments above a cutoff. Formerly 'plot_layer_stats'
#' @param input.layer.ll.df: data frame, ll progression through layers
#' @param layer.corr.df: btween layer correlations
#' @param out.dir: directory to which to write
#' @return list: $pp_ordered, $layer_ll_ordered
#' @export plot_fs_stats()

plot_fs_stats <- function(input.layer.ll.df, layer.corr.df, out.dir, layer.nr, fs.correlation.cutoff){

  pdf(paste0(out.dir, '_k', layer.nr, '_summaryStatPlots.pdf'))

  p.ll.prog <- ggplot() + geom_line(data = input.layer.ll.df, aes(x = FS, y = FS_ll)) +
    ggtitle('Log-Likelihood progression') + theme_bw() + labs(x='Number of FSs (K)', y='Log-Likelihood')

  layer.corr.df.mod <- layer.corr.df
  layer.corr.df.mod$cutoff <- 'below'
  layer.corr.df.mod$cutoff[which(layer.corr.df.mod$corr > fs.correlation.cutoff)] <- 'above'
  layer.corr.df.mod$layer_comb[which(layer.corr.df.mod$corr < fs.correlation.cutoff)] <- layer.corr.df.mod$layer_1[which(layer.corr.df.mod$corr < fs.correlation.cutoff)]

  # attempt to arrange the factors
  layer.corr.df.mod$layer_Iteration <- factor(layer.corr.df.mod$layer_Iteration, levels = c(1:max(as.numeric(layer.corr.df.mod$layer_Iteration))))

  layer.corr.df.mod.iterSplit <- split(layer.corr.df.mod, layer.corr.df.mod$layer_Iteration)
  layer.corr.df.mod.max.corrs <- do.call(rbind, lapply(layer.corr.df.mod.iterSplit, function(x){

    max.row <- which(x$corr == max(x$corr))[1]
    x[max.row,]
  }))

  p.corr <- ggplot() +
    geom_boxplot(data = layer.corr.df.mod, aes(y = corr, x = (layer_Iteration) ),position = position_dodge(0.7)) +
    #geom_point(data = fSim21.all.ap.corrs, aes(y = corr, x = (layer_Iteration) )) +
    geom_line(data = layer.corr.df.mod.max.corrs, aes(x = as.numeric(layer_Iteration), y = corr, color = 'Max. corr'), linetype = "dashed") +
    ggtitle('FS correlations (Boxplots)') + theme_classic() + geom_hline(yintercept=0.1, linetype="dotted", color = "black") +
    labs(x='Number of FSs (K)', y='Correlation of FS prob.', col="", shape="Cutoff = 0.1")

  grid.arrange(p.ll.prog, p.corr, nrow = 2)

  dev.off()
}


#' @title Record the hyper parameters for the background and the functional sorting probabilities
#' @param input.alpha0.list: list of the null-alphas resulting the max lolg-lik. for that layer
#' @param input.alpha1.list: list of the alternative-alphas resulting the max lolg-lik. for that layer
#' @param input.alpha.outDir: directory to which the alphas are to be written
#' @param layer.nr: number of functional sequences recorded
#' @return .csv file
#' @export record_alphas()

record_alphas <- function(input.alpha0.list, input.alpha1.list, input.alpha.outDir, layer.nr){

  total.rows <- length(input.alpha0.list) + length(input.alpha1.list)
  total.cols <- max(c( unlist(lapply(input.alpha0.list, function(x){length(x)})),
                       unlist(lapply(input.alpha1.list, function(x){length(x)}))))

  alpha.matrix <- matrix(0, nrow = total.rows, ncol = total.cols + 1)

  for(i in 1:length(input.alpha0.list)){
    alpha.matrix[i, c(1:length(input.alpha0.list[[i]]))] <- round(input.alpha0.list[[i]] / sum(input.alpha0.list[[i]]), 3)
    alpha.matrix[i, length(input.alpha0.list[[i]]) + 1] <- round(sum(input.alpha0.list[[i]]), 3)
  }

  for(j in (length(input.alpha0.list) + 1):total.rows){
    alpha.matrix[j, c(1:length(input.alpha1.list[[j - length(input.alpha0.list)]]))] <- round(input.alpha1.list[[j - length(input.alpha0.list)]] / sum(input.alpha1.list[[j - length(input.alpha0.list)]]), 3)
    alpha.matrix[j, length(input.alpha1.list[[j - length(input.alpha0.list)]]) + 1] <- round(sum(input.alpha1.list[[j - length(input.alpha0.list)]]), 3)
  }

  alpha.names <- c(paste(rep('alpha0', length(input.alpha0.list)), 'r', c(1:length(input.alpha0.list)), sep = '_'),
                   paste(rep('alpha1', length(input.alpha1.list)), 'r', c(1:length(input.alpha1.list)), sep = '_'))


  out.alpha.df <- cbind(alpha.names, alpha.matrix)

  colnames(out.alpha.df) <- c('alpha_type', paste('pool', c(1:total.cols)), 'dispersion')

  write.csv(out.alpha.df, file = paste0(input.alpha.outDir, '_k', layer.nr, '_alphas.csv'), row.names = F, quote = F)

}


#' @title layer exons such that they don't overlap
#' @param input.score.list: list, each list element contains the name of the bedgraph to be plotted and a data frame with: $chrom, $start, $end, $formatScores
#' @param bg.name: name to add to bedgraph
#' @return bedgraph file
#' @export create_bedgraphs()

create_bedgraphs <- function(input.score.list, bg.name){
  score.names <- names(input.score.list)
  nr.bg <- length(score.names)  # create unique colors for each bedgraph to plot
  bedgraph.colors <- col2rgb(rainbow(nr.bg,s = 1, v = 1, start = 0, end =  max(1, nr.bg - 1)/nr.bg, alpha = 1))
  for(i in 1:length(score.names)){
    intit.temp.score.df <- input.score.list[[i]]
    temp.score.df <- intit.temp.score.df[which(! is.na(intit.temp.score.df$start)),]

    if(nrow(temp.score.df) == 0){
      print('bedraph plotting failed. No non-targeting guides to remove. Fix code to account for that')
      break()
    }

    if(! 'formatScores' %in% names(temp.score.df)){
      if('genomeScore' %in% names(temp.score.df)){
        temp.score.df$formatScores <- temp.score.df$genomeScore
      } else if('guideScore' %in% names(temp.score.df)){
        temp.score.df$formatScores <- temp.score.df$guideScore
      } else if('score' %in% names(temp.score.df)){
        temp.score.df$formatScores <- temp.score.df$score
      } else{
        print('failed to identify score column')
        break()
      }
    }

    temp.score.name <- score.names[i]
    #print(paste0('Creating Bedgraph')) # for: ', temp.score.name))
    temp.bg.color <- paste(bedgraph.colors[1,i], bedgraph.colors[2,i], bedgraph.colors[3,i], sep = ',')
    first.chrom <- temp.score.df[which(temp.score.df$chrom == temp.score.df$chrom[1]),]
    temp.header1 <- paste0('browser position ', first.chrom$chrom[1], ':', min(first.chrom$start, na.rm = T),
                           '-', max(first.chrom$end, na.rm = T), '\n')
    temp.header2 <- paste0("track type=bedGraph name='", bg.name, temp.score.name, "_bg' description='",
                           bg.name, temp.score.name, "_bg' visibility=full color=", temp.bg.color, '\n')

    temp.score.gr <- GRanges(seqnames = temp.score.df$chrom, ranges = IRanges(temp.score.df$start, temp.score.df$end),
                             score = temp.score.df$formatScores)
    temp.score.gr <- sortSeqlevels(temp.score.gr)
    temp.score.gr <- sort(temp.score.gr)
    df.temp.score <- as.data.frame(temp.score.gr)
    df.temp.score.final <- df.temp.score[,c(1,2,3,6)]
    
    # zero-index the files
    df.temp.score.final[,2] <- df.temp.score.final[,2] - 1
    df.temp.score.final[,3] <- df.temp.score.final[,3] - 1

    tmp.file <- paste0(bg.name, '_',temp.score.name, ".bedgraph")
    cat(temp.header1, file = tmp.file)
    cat(temp.header2, file = tmp.file, append = TRUE)

    write.table(df.temp.score.final, file = tmp.file, append = TRUE, row.names = F, col.names = F, quote = F, sep = '\t')

    temp.method.name <- strsplit(score.names[i],'_')[[1]]

  }
}


#' @title Compute the correlation and the overlaps between fs rows
#' @param input.pp: posteriors from convergence
#' @param layer.nr: what layer iteration is being processed
#' @return data frame: $layer_1, $layer_2, $corr, $layer_comb, $layer_Iteration
#' @export pps_corr()

pps_corr <- function(input.pp, layer.nr){

  if(length(unique(input.pp[1,])) > 2){
    print('Posterior matrix is not standard format')
    return()
  }

  layer.cor.tests <- combinations(n=nrow(input.pp), r=2,v=c(1:nrow(input.pp)), repeats.allowed=F)
  layer.cors <- matrix(0, nrow = nrow(layer.cor.tests), ncol = ncol(layer.cor.tests) + 1)
  layer.cors[,1] <- layer.cor.tests[,1] - 1
  layer.cors[,2] <- layer.cor.tests[,2] - 1

  #layers.to.test <- unique(layer.cor.tests[,1])
  for(i in 1:nrow(layer.cor.tests)){

    layer.cors[i,3] <- cor(input.pp[layer.cor.tests[i,1],], input.pp[layer.cor.tests[i,2],])

  }

  # for(i in layers.to.test){
  #
  #   layer.tests <- layer.cor.tests[which(layer.cor.tests[,1] == i), , drop = F]
  #
  #   max.layer.corr <- -3
  #   max.layer.corr.nr <- -1
  #
  #   for(j in 1:nrow(layer.tests)){
  #
  #     temp.corr <- cor(input.pp[layer.tests[j,1],], input.pp[layer.tests[j,2],])
  #
  #     if(temp.corr > max.layer.corr){
  #       max.layer.corr <- temp.corr
  #       max.layer.corr.nr <- layer.tests[j,2]
  #     }
  #
  #   }
  #
  #   layer.cors[i,1] <- i - 1
  #   layer.cors[i,2] <- max.layer.corr.nr - 1
  #   layer.cors[i,3] <- max.layer.corr
  #
  # }

  layer.cors <- as.data.frame(layer.cors)
  colnames(layer.cors) <- c('layer_1', 'layer_2', 'corr')
  layer.cors$layer_comb <- paste(layer.cors$layer_1, layer.cors$layer_2, sep = '_')
  layer.cors$layer_Iteration <- rep(as.character(layer.nr), nrow(layer.cors))


  #layer.cors <- layer.cors[1:layers.to.test,]

  return(layer.cors)

}


#' @title Compute for the regulatory sets for the segment-length normalized posteriors; all segments above a cutoff
#' @param input.pp: posteriors from convergence
#' @param rs.cutoff: ctuoff for defining a regulatory set
#' @return list: $nr_rs, $rs_prob, $training_overlap, $prct_overlap
#' @export pps_stats()

pps_stats <- function(input.pp, rs.cutoff){

  nr.rs <- c()
  rs.prob <- c()
  training.overlap <- c()
  prct.training.overl <- c()

  training.idx <- which(input.pp[1,] == 1)

  for(i in 2:nrow(input.pp)){
    temp.rs.idx <-which(input.pp[i,] > rs.cutoff)

    if(length(temp.rs.idx) == 0){
      nr.rs <- c(nr.rs, 0)
      rs.prob <- c(rs.prob, 0)
      training.overlap <- c(training.overlap, 2) # use 2 as 'don't know'
      prct.training.overl <- c(prct.training.overl, 2)
    } else {
      nr.rs <- c(nr.rs, length(temp.rs.idx))
      rs.prob <- c(rs.prob, round(1 - dpoibin(0, input.pp[i, temp.rs.idx]), 3))

      if(length(which(temp.rs.idx %in% training.idx)) > 0){
        training.overlap <- c(training.overlap, 1)
        prct.training.overl <- c(prct.training.overl, round(length(which(temp.rs.idx %in% training.idx)) / length(temp.rs.idx), 3))
      } else {
        training.overlap <- c(training.overlap, 0)
        prct.training.overl <- c(prct.training.overl, 0)
      }
    }

  }

  out.list <- list()
  out.list$nr_rs <- nr.rs
  out.list$rs_prob <- rs.prob
  out.list$training_overlap <- training.overlap
  out.list$prct_overlap <- prct.training.overl

  return(out.list)
}


#' @title Arrange the fs rows in order of their likelihood model contribution
#' @param input.pp: posteriors from convergence
#' @param input.total.ll: overall log-likelihood of the model
#' @param input.data: list: $guide_to_seg_lst,
#' @param input.alpha0: list, alpha 0s for the replicates
#' @param input.alpha1: list, alpha 1s for the replicates
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @param local.max: logical, whether a local max was computed
#' @param fs.ll.rt: per functional sequence, per segment log-likelihood ratio. Is empty list if local.max == FALSE
#' @return list: $pp_ordered, $layer_ll_ordered, $fs_ll_rt_ordered
#' @export order_pps()

order_pps <- function(input.pp, input.total.ll, input.data, input.alpha0, input.alpha1, 
                      guide.efficiency, local.max, fs.ll.rt){

  final.ll <- input.total.ll
  final.pp <- input.pp
  layer.ll.diff <- c()

  for(f in 2:nrow(final.pp)){
    temp.final.pp <- final.pp
    temp.final.pp[f,] <- 0

    temp.final.pp.sum <- colSums(temp.final.pp)
    temp.final.pp.sum[temp.final.pp.sum > 1] <- 1

    temp.guide.ll <- compute_perGuide_fs_ll(temp.final.pp.sum, input.data$guide_to_seg_lst)

    temp.ll <- 0
    for(r in 1:length(input.data$data)){
      temp.hypers <- list(alpha0 = input.alpha0[[r]], alpha1 = input.alpha1[[r]])
      temp.sgRNa.ll <- estimate_relics_sgrna_log_like(temp.hypers,
                                                                input.data$data[[r]],
                                                                temp.guide.ll,
                                                              guide.efficiency)

      temp.dirichlet.ll <- sum(temp.sgRNa.ll$total_guide_ll)
      temp.ll <- temp.ll + temp.dirichlet.ll
    }

    layer.ll.diff <- c(layer.ll.diff, temp.ll - final.ll)
  }

  out.list <- list()

  if(nrow(input.pp) == 2){
    out.list$pp_ordered <- input.pp
    out.list$layer_ll_ordered <- layer.ll.diff
    out.list$fs_ll_rt_ordered <- fs.ll.rt

    return(out.list)
  } else {
    pp.to.order <- input.pp[2:nrow(input.pp),]
    pp.ordered <- pp.to.order[order(layer.ll.diff),]

    out.pp <- rbind(input.pp[1,], pp.ordered)

    layer.ll.diff.ordered <- layer.ll.diff[order(layer.ll.diff)]
    
    fs.ll.rt.ordered  <- fs.ll.rt
    if(local.max){
      fs.ll.rt.to.order <- fs.ll.rt[2:nrow(input.pp),]
      fs.ll.rt.ordered <- fs.ll.rt.to.order[order(layer.ll.diff),]
      fs.ll.rt.ordered <- rbind(fs.ll.rt[1, ], fs.ll.rt.ordered)
    }

    out.list$pp_ordered <- out.pp
    out.list$layer_ll_ordered <- layer.ll.diff.ordered
    out.list$fs_ll_rt_ordered <- fs.ll.rt.ordered

  }

  return(out.list)
}


#' @title Set initial hyper parameters, generate the fs posterior matrix
#' @param hyper: hyperparameters
#' @param in.data.list: list, $seg_info, $fs0_idx
#' @param known.reg: position of known regions
#' @return list
#' @export init_relics_param()

init_relics_param <- function(hyper, in.data.list, local.max) {
  param <- list()

  # delta posterior probabilities
  param$delta.pp <- matrix(0, nrow=hyper$L+1, ncol= nrow(in.data.list$seg_info)) #nrow(in.data.list$overlapMat))

  # first row is special and contains posterior probs of known regulatory elements
  param$delta.pp[1, in.data.list$fs0_idx] <- 1
  
  param$ll_rt <- list()
  
  if(local.max){
    param$ll_rt <- matrix(0, nrow=hyper$L+1, ncol= nrow(in.data.list$seg_info))
  }

  return(param)
}

#' @title Display the credible sets along with cumulative posteriors. Formerly 'dirichlet_addLayers_mrvr_norm'
#' @param input.param, list containing the matrix for the posteriors
#' @param input.hyper: hyper parameters
#' @param input.data.list, list: $guide_to_seg_lst, $data, $true_pos_seg, $overlapMat, $guide_to_seg_lst, $seg_to_guide_lst
#' @param input.tol: nr of layers to create
#' @param fix.hypers, whetehr the hyper parameters are fixed
#' @param iterative.hyper.est: whether the hyper parameters should be re-calculated after each round of estimating the posteriors. Default = TRUE
#' @param nr.segs: max. length a functional segment is considered to have
#' @param geom.p: proababilty of the genometric distribution to penalize for enhancers of increasing length
#' @param min.pp: minimum posterior required to be part of a regulatory set
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @param one.dispersion: whether the hyper parameters should be estimated with one or two dispersions
#' @param local.max: logical, whether a local maximum should be computed
#' @param local.max.range: number of segments to include in addition to the the ones with highest PP (one way, so multiply by 2 for total nr of segs)
#' @return ll_tract, posterior_trace_list, alpha0_est, alpha1_est, max_iter_reached, fs_ll_rt_trace
#' @export relics_compute_FS_k()

relics_compute_FS_k <- function(input.param,
                          input.hyper,
                          input.data.list,
                          input.tol = 0.001,
                          fix.hypers = FALSE,
                          iterative.hyper.est = TRUE,
                          nr.segs = 10, geom.p = 0.1,
                          min.pp = 0.1,
                          guide.efficiency,
                          one.dispersion,
                          local.max, local.max.range){


  dirichlet.hyper <- input.hyper
  dirichlet.param <- input.param

  dirichlet.pp <- colSums(dirichlet.param$delta.pp)
  dirichlet.pp[dirichlet.pp > 1] <- 1

  tol <- input.tol
  max.iter <- 20
  max.iter.reached <- FALSE
  max.diff <- 1e6
  iter <- 0

  dirichlet.ll.trace <- c()
  dirichlet.posterior.trace <- list()
  max.diff.trace <- c()
  alpha0.est <- list()
  alpha1.est <- list()
  fs.ll.rt.trace <- list() # per functional sequence, per segment log-likelihood ratio of local max

  too.low.pps <- 3 # 3 strikes to get posteriors higher than lowest possibility

  while((max.diff > tol) & (iter < max.iter) & too.low.pps > 0) {

    iter <- iter + 1

    # print(paste('alpha0:', lapply(dirichlet.hyper$alpha0, function(x){round(x, 3)})))
    # print(paste('alpha1:', lapply(dirichlet.hyper$alpha1, function(x){round(x, 3)})))

    alpha0.est[[iter]] <- dirichlet.hyper$alpha0
    alpha1.est[[iter]] <- dirichlet.hyper$alpha1

    cur.delta <- dirichlet.param$delta.pp

    dirichlet.param <- relics_estimate_pp(dirichlet.param,
                                           dirichlet.hyper,
                                           input.data.list$data,
                                           input.data.list$true_pos_seg,
                                           input.data.list$guide_to_seg_lst,
                                           input.data.list$seg_to_guide_lst,
                                           input.data.list$next_guide_lst,
                                           nr.segs, geom.p, guide.efficiency,
                                          local.max, local.max.range)


    # keep track of changes in posterior estimates for delta
    delta.diff <- cur.delta - dirichlet.param$delta.pp
    max.diff <- max(abs(delta.diff))

    # print(paste0('diff in posterior: ', max.diff))

    max.diff.trace <- c(max.diff.trace, max.diff)
    dirichlet.posterior.trace[[iter]] <- dirichlet.param$delta.pp
    
    if(local.max){
      fs.ll.rt.trace[[iter]] <- dirichlet.param$ll_rt
    }

    if(!fix.hypers & iterative.hyper.est){
      dirichlet.hyper <- recompute_hyper_parameters(dirichlet.param,
                                                     dirichlet.hyper,
                                                     input.data.list$data,
                                                     input.data.list$guide_to_seg_lst,
                                                    guide.efficiency,
                                                    one.dispersion)
    }

    dirichlet.pp <- colSums(dirichlet.param$delta.pp)
    dirichlet.pp[dirichlet.pp > 1] <- 1

    dirichlet.guide.ll <- compute_perGuide_fs_ll(dirichlet.pp, input.data.list$guide_to_seg_lst)
    dirichlet.ll <- 0
    for(i in 1:length(input.data.list$data)){
      temp.hypers <- list(alpha0 = dirichlet.hyper$alpha0[[i]], alpha1 = dirichlet.hyper$alpha1[[i]])
      temp.sgRNa.ll <- estimate_relics_sgrna_log_like(temp.hypers,
                                                                input.data.list$data[[i]],
                                                                dirichlet.guide.ll,
                                                              guide.efficiency)
      
      temp.dirichlet.ll <- sum(temp.sgRNa.ll$total_guide_ll)
      dirichlet.ll <- dirichlet.ll + temp.dirichlet.ll
    }

    #print(paste0('Updated ll: ', sum(dirichlet.ll)))
    dirichlet.ll.trace <- rbind(dirichlet.ll.trace, dirichlet.ll)

    # check if the last layer has a posterior below a certain threshold
    if(nrow(dirichlet.param$delta.pp) > 2){
      iter.pps <- dirichlet.param$delta.pp[2:nrow(dirichlet.param$delta.pp),]
      iter.pps.max <- apply(iter.pps, 1, max)

      # three attempts to not have a low layer
      if(length(which(iter.pps.max < min.pp)) > 0){
        too.low.pps <- too.low.pps - 1
      }
    }

  }

  # if max.inter has been reached, look back 5 iterations, select the one with best ll to continue with
  if(iter == max.iter){

    max.iter.reached = TRUE

    # if there still was a constant decrease in the posterior difference then keep the last iteration, otherwise pick the one with min-ll
    last.iter.diffs <- round(max.diff.trace[iter:(iter-5)],3)

    # if either there was no constand decrease in posterior difference or that the derease is too small
    if(is.unsorted(last.iter.diffs) | sum(last.iter.diffs - last.iter.diffs[1]) < 0.01){
      max.ll.idx <- which(dirichlet.ll.trace[iter:(iter-5)] == max(dirichlet.ll.trace[iter:(iter-5)]))[1]
      dirichlet.posterior.trace[[iter]] <- dirichlet.posterior.trace[[iter - max.ll.idx + 1]]
      dirichlet.ll.trace[iter] <- dirichlet.ll.trace[iter - max.ll.idx + 1]
      alpha0.est[[iter]] <- alpha0.est[[iter - max.ll.idx + 1]]
      alpha1.est[[iter]] <- alpha1.est[[iter - max.ll.idx + 1]]
      
      if(local.max){
        fs.ll.rt.trace[[iter]] <- fs.ll.rt.trace[[iter - max.ll.idx + 1]]
      }
      
    }


  }

  out.list <- list(ll_tract = dirichlet.ll.trace,
                   posterior_trace_list = dirichlet.posterior.trace,
                   alpha0_est = alpha0.est,
                   alpha1_est = alpha1.est,
                   max_iter_reached = max.iter.reached,
                   fs_ll_rt_trace = fs.ll.rt.trace)


  return(out.list)

}


#' @title wrapper which returns -log likelihood given provided prior distribution hyperparameters, and current estimates of the cumulative posteriors. Formerly 'fit_prior_dirichlet_distr_mvr'
#' @param hyper: hyperparameters
#' @param param: matrix of all posterior probs
#' @param data: data, consists of: $y1, $y2 and $n
#' @param guide.seg.idx.lst: list: each element is a guide, contaiing the indexes of the segments it overlaps
#' @param guide.efficiency: either a vector of guide efficiency or NULL
#' @param one.dispersion. logical, if there should be one or two dispersions for the hyper parameters
#' @return list of the re-estimated hyperparameters
#' @export recompute_hyper_parameters()

recompute_hyper_parameters <- function(param, hyper, data, guide.seg.idx.lst, guide.efficiency, one.dispersion) {
  cumulative.pp <- colSums(param$delta.pp) #apply(param$delta.pp, 2, sum)
  cumulative.pp[cumulative.pp > 1] <- 1

  guide.lls.list <- compute_perGuide_fs_ll(cumulative.pp, guide.seg.idx.lst)

  for(i in 1:length(data)){
    # need alphas to be greater than 0, taking square root now to exponentiate after
    hyper.param <- c(sqrt(hyper$alpha0[[i]]), sqrt(hyper$alpha1[[i]]))
    alpha0.idx <- c(1:length(hyper$alpha0[[i]]))
    alpha1.idx <- c(1:length(hyper$alpha0[[i]])) + max(alpha0.idx)
    
    res <- c()
    
    if(one.dispersion){
      res <- optim(hyper.param, prior_dirichlet_ll_singleDisp, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
                   data=data[[i]], region.ll.list = guide.lls.list,
                   alpha0.idx = alpha0.idx, alpha1.idx = alpha1.idx, 
                   guide.efficiency = guide.efficiency)
    } else {
      res <- optim(hyper.param, prior_dirichlet_ll, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
                   data=data[[i]], region.ll.list = guide.lls.list,
                   alpha0.idx = alpha0.idx, alpha1.idx = alpha1.idx, 
                   guide.efficiency = guide.efficiency)
    }


    if(res$convergence==0) {
      # return new estimates of hyperparamers
      
      if(one.dispersion){
        alpha0s <- res$par[alpha0.idx]**2
        alpha1s <- res$par[alpha1.idx]**2
        
        alpha1s.norm <- alpha1s / sum(alpha1s)
        alpha1s.adj <- alpha1s.norm * sum(alpha0s)
        
        hyper$alpha0[[i]] <- alpha0s
        hyper$alpha1[[i]] <- alpha1s.adj
      } else {
        hyper$alpha0[[i]] <- res$par[alpha0.idx]**2
        hyper$alpha1[[i]] <- res$par[alpha1.idx]**2
      }


    } else {
      warning("estimation of hyperparameters failed to converge")
    }
  }


  return(hyper)
}


#' @title estimate the posterior probabilities for placing a specified number of functional sequences over a specified number of segments. Formerly 'fit_dirichlet_model_fast1_mrvr_norm'
#' @param hyper: hyperparameters
#' @param param: matrix of posterior probabilities
#' @param data: list, each element is a replicate
#' @param known.reg: region positions of the known regions
#' @param guide.to.seg.lst: mapping of guides to segments
#' @param seg.to.guide.lst: mapping of segments to guides, both $guide_idx, $nonGuide_idx
#' @param next.guide.lst: list that contains the index of $next_guide_idx and next_nonGuide_idx
#' @param nr.segs: max number of segemnts that can be combined
#' @param geom.p: probability of the geometic distribution used to calculate the probability of there being x segments in the regulatory element
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @param local.max: logical, whether a local maximum should be computed
#' @param local.max.range: number of segments to include in addition to the the ones with highest PP (one way, so multiply by 2 for total nr of segs)
#' @return log likelihood for each region
#' @export relics_estimate_pp()

relics_estimate_pp <- function(param, hyper, data, known.reg,
                                guide.to.seg.lst, seg.to.guide.lst,
                                next.guide.lst, nr.segs = 10, 
                               geom.p = 0.1, guide.efficiency,
                               local.max, local.max.range) {
  n.sgrna <- length(guide.to.seg.lst)
  n.region <- length(seg.to.guide.lst)

  for(l in 2:(hyper$L+1)) {
    # set posteriors for this row to 0
    param$delta.pp[l,] <- rep(0, n.region)

    # sum posteriors from other rows
    pp <- apply(param$delta.pp, 2, sum)

    # cumulative posterior probability can become greater than 1,
    # So we restrict it to 1.  doesn't seem like the ideal solution...
    pp[pp > 1.0] <- 1.0

    # Compute the ll for each RE length for both replicates, combine the lls and then add to the ll total
    # set up the data structues as lists to access the different replicates
    layer.guide.ll <- compute_perGuide_fs_ll(pp, guide.to.seg.lst) # processed.pp
    sgrna.log.like.list <- list()
    data.mat.list <- list()
    data.total.list <- list()

    for(repl in 1:length(data)){
      temp.hypers <- list(alpha0 = hyper$alpha0[[repl]], alpha1 = hyper$alpha1[[repl]])
      temp.data <- data[[repl]]
      temp.data.counts <- temp.data[, c(1:(ncol(temp.data) - 1))]
      temp.data.totals <- temp.data[, ncol(temp.data)]

      data.mat.list[[repl]] <- temp.data.counts
      data.total.list[[repl]] <- temp.data.totals

      temp.sgrna.log.like <- estimate_relics_sgrna_log_like(temp.hypers, temp.data, layer.guide.ll, guide.efficiency)
      sgrna.log.like.list[[repl]] <- temp.sgrna.log.like

    }

    delta.pps <- estimate_fs_pp(#hyper, data.mat.list, data.total.list, 
                                seg.to.guide.lst, next.guide.lst, nr.segs, geom.p, sgrna.log.like.list)
    
    # if the signal is to be located at a local region
    if(local.max){
      
      # identify segment with highest posterior, pull our segments +/- local.max.region
      local.max.idx <- extract_local_max_idx(delta.pps, local.max.range) # toDo, add 'local.max.range's
      
      # temp.local.seg.to.guide.lst <- seg.to.guide.lst[local.max.idx]
      # # adjust $nonGuide_idx to match local region
      # local.seg.to.guide.lst <- extract_local_seg_to_guide(temp.local.seg.to.guide.lst, local.max.idx)
      local.seg.to.guide.lst <- seg.to.guide.lst[local.max.idx]
      
      # local.seg.to.guide.lst <- seg.to.guide.lst[local.max.idx] # toDo (pull out the ones needed)
      local.next.guide.lst <- next.guide.lst[local.max.idx] # toDo
      # sg.ll.list # can be left, as the index for the guide remains unchanged
      
      local.delta.pps <- estimate_fs_pp(local.seg.to.guide.lst, local.next.guide.lst, 
                                        nr.segs, geom.p, sgrna.log.like.list)
      
      delta.pps <- rep(0, length(delta.pps))
      delta.pps[local.max.idx] <- local.delta.pps
      
      # compute log-lik. ratio for each segment
      ll.rt <- rep(0, length(delta.pps))
      
      for(repl in 1:length(data)){
        temp.hypers <- list(alpha0 = hyper$alpha0[[repl]], alpha1 = hyper$alpha1[[repl]])
        temp.data <- data[[repl]]
        temp.data.counts <- temp.data[, c(1:(ncol(temp.data) - 1))]
        temp.data.totals <- temp.data[, ncol(temp.data)]
        
        temp.guide.model.ll <- estimate_relics_sgrna_log_like(temp.hypers, temp.data, layer.guide.ll, guide.efficiency, return.model.ll = TRUE)
        
        local.ll.rt <- compute_local_ll_ratio(temp.guide.model.ll, local.seg.to.guide.lst)
        ll.rt[local.max.idx] <- ll.rt[local.max.idx] + local.ll.rt
        
      }

      param$ll_rt[l, ] <- ll.rt
      
    }

    param$delta.pp[l, ] <- delta.pps
  }

  return(param)
}


#' @title Identifying the segment with the highest PP and surrounding segments (focuses on the first segment if multiple max PPs)
#' @param input.pps: vector, contains the PPs
#' @param local.max.range: number of segments to include in addition to the the ones with highest PP (one way, so multiply by 2 for total nr of segs)
#' @return vector with idx of region to look for signal
#' @export extract_local_max_idx()

extract_local_max_idx <- function(input.pps, local.max.range){
  
  center.idx <- which(input.pps == max(input.pps))[1]
  min.idx <- center.idx - local.max.range
  max.idx <- center.idx + local.max.range
  
  if(min.idx < 1){
    min.idx <- 1
  }
  if(max.idx > length(input.pps)){
    max.idx <- length(input.pps)
  }
  
  return(c(min.idx:max.idx))
  
}


#' @title adjust the non_guide index according to the window
#' @param seg.to.guide.lst: list, $guide_idx, $nonGuide_idx
#' @param local.max.range: number of segments to include in addition to the the ones with highest PP (one way, so multiply by 2 for total nr of segs)
#' @return vector with idx of region to look for signal
#' @export extract_local_seg_to_guide()

extract_local_seg_to_guide <- function(seg.to.guide.lst, local.max.idx){
  
  all.local.guides <- unique(unlist(lapply(seg.to.guide.lst, function(x){
    x$guide_idx
  })))
  
  local.seg.to.guide.lst <- lapply(seg.to.guide.lst, function(x){
    x$nonGuide_idx <- all.local.guides[which(! all.local.guides %in% x$guide_idx)]
    
    x
  })

  return(local.seg.to.guide.lst)
  
}


#' @title comput log-lik. ratio for the selected idx
#' @param local.seg.to.guide.lst: list, $guide_idx, $nonGuide_idx
#' @param guide.ll.df: data.frame, null_only_ll, alt_only_ll
#' @return vector with idx of region to look for signal
#' @export compute_local_ll_ratio()

compute_local_ll_ratio <- function(guide.ll.df, local.seg.to.guide.lst){
  
  fs.ll <- unlist(lapply(local.seg.to.guide.lst, function(x){
    sum(guide.ll.df$alt_only_ll[x$guide_idx])
  }))
  
  background.ll <- unlist(lapply(local.seg.to.guide.lst, function(x){
    sum(guide.ll.df$null_only_ll[x$guide_idx])
  }))
  
  out.ll.rt <- -2 * (background.ll - fs.ll)
  
  return(out.ll.rt)
  
}


#' @title Compute the log likelihood of each possible configuration of the placement of a fs amonst a set f segments. Formerly 'compute_dirichlet_delta_ll_mrvr_nromSegLLMat'
#' @param hyper: list, hyperparameters
#' @param in.data.list: list, each element contains a replicate, the total column has already been removed
#' @param in.data.totals: list, each element contains the per-guide totals for a replicate
#' @param seg.to.guide.lst: list, each element is a segment, contains a list with guide_idx, nonGuide_idx
#' @param next.guide.lst: list that contains the index of $next_guide_idx and next_nonGuide_idx
#' @param nr.segs: max number of segemnts that can be combined
#' @param geom.p: probability of the geometic distribution used to calculate the probability of there being x segments in the regulatory element
#' @param sg.ll.list: list, each element is a data frame per.guide log likelihood given the posteriors for a replicate and alt. ll (total_guide_ll, alt_only_ll)
#' @return log likelihood for each segment
#' @export estimate_fs_pp()

estimate_fs_pp <- function(#hyper, in.data.list, in.data.totals,
                           seg.to.guide.lst, next.guide.lst, nr.segs = 10, geom.p = 0.1,
                           sg.ll.list) {

  # set up a list, the length of the number of segments
  # each element contains a vector, which contains the ll of that segment being overlapped by a regulatory element
  segment.ll.list <- list()

  total.segs <- length(seg.to.guide.lst)
  region.priors <- rep(1/total.segs, total.segs)
  geom.norm.factr <- pgeom(nr.segs, geom.p)

  # initialize the total likelihood with the first, one-segment likelihood
  initial.segment.ll <- 0
  initial.guide.idx<- seg.to.guide.lst[[1]]$guide_idx
  initial.nonGuide.idx<- seg.to.guide.lst[[1]]$nonGuide_idx

  for(repl in 1:length(sg.ll.list)){ # in.data.list
    # initial.guide.ll <- sum(ddirmnom(in.data.list[[repl]][initial.guide.idx,],
    #                                  size = in.data.totals[[repl]][initial.guide.idx],
    #                                  alpha = hyper$alpha1[[repl]], log = T))
    
    initial.guide.ll <- sum(sg.ll.list[[repl]][initial.guide.idx,2])
    initial.nonGuide.ll <- sum(sg.ll.list[[repl]][initial.nonGuide.idx,1])
    initial.segment.ll <- initial.segment.ll + initial.guide.ll + initial.nonGuide.ll + log(dgeom(1, geom.p) / geom.norm.factr)
  }

  segment.ll.list[[1]] <- initial.segment.ll
  total.ll <- initial.segment.ll

  # make an initial pass with one-segment length RE placements
  for(seg in 2:total.segs){
    temp.segment.ll <- 0
    temp.guide.idx<- seg.to.guide.lst[[seg]]$guide_idx
    temp.nonGuide.idx<- seg.to.guide.lst[[seg]]$nonGuide_idx

    for(repl in 1:length(sg.ll.list)){ # in.data.list
      # temp.guide.ll <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx,],
      #                               size = in.data.totals[[repl]][temp.guide.idx],
      #                               alpha = hyper$alpha1[[repl]], log = T))
      
      temp.guide.ll <- sum(sg.ll.list[[repl]][temp.guide.idx,2])
      temp.nonGuide.ll <- sum(sg.ll.list[[repl]][temp.nonGuide.idx,1])
      temp.segment.ll <- temp.segment.ll + temp.guide.ll + temp.nonGuide.ll + log(dgeom(1, geom.p) / geom.norm.factr)
    }
    segment.ll.list[[seg]] <- temp.segment.ll
    total.ll <- addlogs(total.ll, temp.segment.ll)
  }

  all.guide.idx <- c(1:length(sg.ll.list[[1]][,1])) # helps to differentitate between guides for null vs alternative

  # now make a pass through all 2 to nr.segs length segments
  for(ns in 1:(nr.segs-1) ){

    for(seg in 1:(total.segs - nr.segs)){

      temp.stretch.ll <- 0
      temp.stretch.segs <- c(seg:(seg + ns)) # index of segments to be included as part of the stretch and be updated in this iteration
      temp.guide.idx <- c(seg.to.guide.lst[[seg]]$guide_idx, unlist(next.guide.lst[(seg + 1):(seg + ns)]))
      temp.nonGuide.idx <- all.guide.idx[-temp.guide.idx]

      for(repl in 1:length(sg.ll.list)){ # in.data.list
        # temp.guide.ll <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx,],
        #                               size = in.data.totals[[repl]][temp.guide.idx],
        #                               alpha = hyper$alpha1[[repl]], log = T))
        
        temp.guide.ll <- sum(sg.ll.list[[repl]][temp.guide.idx,2])
        temp.nonGuide.ll <- sum(sg.ll.list[[repl]][temp.nonGuide.idx,1])

        # likelihood of this continous stretch of bins to contain a regulatory element
        temp.stretch.ll <- temp.stretch.ll + temp.guide.ll + temp.nonGuide.ll + log(dgeom(1 + ns, geom.p) / geom.norm.factr)

      }

      total.ll <- addlogs(total.ll, temp.stretch.ll)

      for(span in 1:length(temp.stretch.segs)){
        segment.ll.list[[temp.stretch.segs[span]]] <- c(segment.ll.list[[temp.stretch.segs[span]]], temp.stretch.ll)
      }
    }
  }

  # normalize the list by the total, and add the logs across the segments
  segment.lls <- unlist(lapply(segment.ll.list, function(x){
    x.norm <- x - total.ll

    x.norm.sum <- combine_segment_ll(x.norm)
    x.norm.sum
  }))

  segment.posteriors <- exp(segment.lls)

  return(segment.posteriors)

}


#' @title given the already computed segment log-likelihoods, update the bins specified with the addition of the log-lik.
#' @param input.segment.lls: vector, contains log-likelihoods of the different regulalotry elements overlapping it
#' @return list of processed posteriors and vector of length equivalent to the number of segments handled
#' @export combine_segment_ll()

combine_segment_ll <- function(input.segment.lls){

  out.ll <- input.segment.lls[1]

  if(length(input.segment.lls) > 1){
    for(i in 2:length(input.segment.lls)){
      out.ll <- addlogs(out.ll, input.segment.lls[i])
    }
  }

  return(out.ll)

}

#' @title Numerically stable addition of logs
#' @param loga: first log value
#' @param logb: second log value
#' @return log addition
#' @export addlogs()

addlogs <- function(loga, logb) {
  max(loga, logb) + log(1 + exp(-abs(loga - logb)))
}

