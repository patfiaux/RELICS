suppressMessages(library(GenomicRanges))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))

suppressMessages(library(poibin)) # to calculate Poisson-Binomial

suppressMessages(library(extraDistr)) # Dirichlet-Multinomial distributions

suppressMessages(library(gtools)) # for getting all combinations when computing correlations between PPs

suppressMessages(library(splines)) # enables computing of splines



#' @title RELICS 2.0 analysis function. Uses IBSS to return a set of functional sequences FS for CRISPR regulatory screens
#' @param input.parameter.file: location of file containing all analysis parameters
#' @param input.parameter.list: RELICS analysis parameters in list format. Default = NULL
#' @param data.file.split: logical, whether counts are separate from info, default = FALSE
#' @param record.all.fs: logical, if information of all intermediate FS should be recorded, in addition to the final setof FS
#' @param return.init.hypers: logical, whether the hyperparamters should be returned after initial computing
#' @return list of final per-layer posteriors
#' @export RELICS()

RELICS <- function(input.parameter.file, input.parameter.list = NULL, data.file.split = FALSE,
                   record.all.fs = FALSE, return.init.hypers = FALSE, input.hyper.prop = NULL){

  analysis.parameters <- read_parameters(input.parameter.list, input.parameter.file, data.file.split)

  data.setup <- set_up_RELICS_data(analysis.parameters, data.file.split,
                                   guide.offset = analysis.parameters$crisprEffectRange,
                                   repl_pools = analysis.parameters$repl_groups,
                                   labelHierarchy = analysis.parameters$labelHierarchy,
                                   fs0.label = analysis.parameters$FS0_label,
                                   file.save.dir = paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName),
                                   save.files = analysis.parameters$save_input_files,
                                   min.seg.dist = analysis.parameters$seg_dist)

  
  
  if(is.null(input.hyper.prop)){
    analysis.parameters <- compute_hyperparameters(analysis.parameters, data.setup)
    if(return.init.hypers){
      out.hyper.list <- record_hyperparameters(input.bkg.alpha = analysis.parameters$hyper_par_components$bkg_alpha, 
                             input.fs.alpha = analysis.parameters$hyper_par_components$FS_alpha, 
                             input.bkg.disp = analysis.parameters$hyper_par_components$dispersion, 
                             paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName), 
                             0, analysis.parameters$pool_names, print.out = TRUE)
      return(out.hyper.list)
    }
  } else {
    analysis.parameters <- integrate_hyperparameters(analysis.parameters, data.setup, input.hyper.prop)
  }

  data.setup$fixed_ge_coeff <- analysis.parameters$fixed_ge_coeff
  # # if guide efficiency scores are provided, calculate guide efficiency and include in the model
  if(! is.null(data.setup$guide_efficiency_scores)){

    if(! analysis.parameters$fixed_ge_coeff){

      temp.relics.params <- init_relics_param(analysis.parameters$hyper_pars, data.setup, analysis.parameters$local_max, FALSE)

      ge.list <- recompute_ge_coefficients(temp.relics.params,
                                           analysis.parameters$hyper_par_components, #analysis.parameters$hyper_pars,
                                           data.setup$data,
                                           data.setup$guide_to_seg_lst,
                                           data.setup$guide_efficiency_scores,
                                           c(0, rep(1, ncol(data.setup$guide_efficiency_scores))))

      data.setup$guide_efficiency <- ge.list$guide_efficiency
      data.setup$ge_coeff <- ge.list$ge_coeff
      
      analysis.parameters$hyper_pars <- incorporate_ge(analysis.parameters$hyper_pars,
                                                       analysis.parameters$hyper_par_components,
                                                       data.setup$data,
                                                       data.setup$guide_efficiency,
                                                       analysis.parameters)
    } else {
      logit.ge <- apply(data.setup$guide_efficiency_scores, 2, function(x){
        log(x / (1 - x))
      })
      guide.efficiency <- 1 / (1 + exp(-(logit.ge %*% rep(1, ncol(data.setup$guide_efficiency_scores)))))
      data.setup$guide_efficiency <- guide.efficiency
      data.setup$ge_coeff <- rep(1, ncol(data.setup$guide_efficiency_scores))
    }

  }

  # plot the per-segment ll-ratio
  out.pars <- list(out_dir = paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName),
                   iter = 0)
  # issue is that this does not take into account the AoE
  # an option is to to run the full pp ll calculation and use that instead
  # as safety I should check if using the normal AoE will result in the same as the current implemetation
  # while the new one with AoE should mirror the acutal results more....
  record_ll_ratio(analysis.parameters$hyper_pars, data.setup, out.pars, analysis.parameters, '')
  record_guide_ll_ratio(analysis.parameters$hyper_pars, data.setup, out.pars, analysis.parameters, '')
  
  analysis.parameters <- set_up_fs_priors(analysis.parameters, data.setup)
  
  run_RELICS_7(input.data = data.setup,
               final.layer.nr = analysis.parameters$max_fs_nr,
               out.dir = paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName),
               input.hypers = analysis.parameters$hyper_pars,
               input.hyper.components = analysis.parameters$hyper_par_components,
               fix.hypers = analysis.parameters$fix_hypers,
               nr.segs = analysis.parameters$nr_segs,
               geom.p = analysis.parameters$geom_p,
               input.min.rs.pp = analysis.parameters$min_fs_pp,
               auto.stop = analysis.parameters$auto_stop,
               record.all.fs = record.all.fs,
               input.convergence.tol = analysis.parameters$convergence_tol,
               one.dispersion = analysis.parameters$one_dispersion,
               recompute.fs0 = analysis.parameters$recompute_fs0,
               local.max = analysis.parameters$local_max,
               local.max.range = analysis.parameters$local_max_range,
               pool.names = analysis.parameters$pool_names,
               mean.var.type = analysis.parameters$mean_var_type,
               analysis.parameters = analysis.parameters,
               guide.dist.to.seg = data.setup$guide_dist_to_seg)
  
}


#' @title helper function, manages the reading in of the analysis parameters
#' @param analysis.parameters: list, contains necessary flags and elements
#' @param data.setup: data fromatted for RELICS analysis
#' @param input.hyper.prop: list of lists, containing new hyperparameter sorting proportions
#' @return list: all parameters embedded in a list
#' @export integrate_hyperparameters()

integrate_hyperparameters <- function(analysis.parameters, data.setup, input.hyper.prop){
  
  if(analysis.parameters$model_dispersion){
    
    if(analysis.parameters$estimateSpline){
      print('not implemented yet')
      # repl.splines <- identify_splines()
      break()
    } else {
      repl.spline.df <- analysis.parameters$repl_spline_df
      repl.disp <- disp_from_spline(repl.spline.df, data.setup$data,
                                    paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName),
                                    analysis.parameters$nr_disp_bins,
                                    analysis.parameters$FS0_label, 
                                    data.setup$guide_info)
      analysis.parameters$repl_disp <- repl.disp
    }
    
    fs0.alphas <- reverse_hyperparameters(data.setup, analysis.parameters, input.hyper.prop, analysis.parameters$FS0_label)
    
  } else {
    print('not implemented yet')
  }
  
  analysis.parameters$hyper_pars <- fs0.alphas$hyper_pars
  analysis.parameters$hyper_par_components <- fs0.alphas$hyper_par_components
  analysis.parameters$init_model_ll <- fs0.alphas$init_model_ll
  
  # if 'intermediate' signal selected:
  if(analysis.parameters$hyper_adj < 1){
    record_orig_fs0_alphas(fs0.alphas$orig_hyper_par_components, analysis.parameters, analysis.parameters$pool_names) # record just the fs.alphas
  }
  
  return(analysis.parameters)
  
}


#' @title helper function, manages the reading in of the analysis parameters
#' @param data.par.list: list, contains the elements of the processed and formatted data
#' @param analysis.par.list: list, contains all analysis flags
#' @param input.hyper.prop: list of lists, containing new hyperparameter sorting proportions
#' @param fs0.label: label of FS0
#' @return list: all parameters embedded in a list
#' @export reverse_hyperparameters()

reverse_hyperparameters <- function(data.par.list, analysis.par.list, input.hyper.prop, fs0.label){
  
  model.ll <- 0
  
  orig.fs.par <- list()
  
  final.dirichlet.pars <- list()
  final.dirichlet.pars$bkg_alpha <- list()
  final.dirichlet.pars$FS_alpha <- list()
  final.dirichlet.pars$dispersion <- list()
  
  final.alpha <- list()
  final.alpha$alpha0 <- list()
  final.alpha$alpha1 <- list()
  
  # make sure all guides are considered to be the same category
  fs.assignment <- rep(0, nrow(data.par.list$seg_info))
  fs.assignment[which(data.par.list$seg_info$label %in% fs0.label)] <- 1
  
  dirichlet.guide.ll <- compute_perGuide_fs_ll(fs.assignment, data.par.list$guide_to_seg_lst, hyper.setup = TRUE)
  
  for(i in 1:length(data.par.list$data)){
    
    bkg.prop <- input.hyper.prop$bkg_prop[[i]]
    fs.prop <- input.hyper.prop$fs_prop[[i]]
    bkg.multiplier <- 1/bkg.prop[1]
    fs.multiplier <- 1/fs.prop[1]
    bkg.par <- bkg.prop[2:length(bkg.prop)] * bkg.multiplier
    fs.par <- fs.prop[2:length(fs.prop)] * fs.multiplier
    
    temp.hyper.params <- c(bkg.par, fs.par)
    temp.bkg.idx <- c(1:length(bkg.par))
    temp.fs.idx <- c(1:length(fs.par)) + max(temp.bkg.idx)
    
    repl.data <- data.par.list$data[[i]]
    temp.repl.disp <- analysis.par.list$repl_disp[[i]]$repl_disp
    
    model.ll <- model.ll + prior_dirichlet_proportions(temp.hyper.params, data = repl.data, 
                                                       region.ll.list = dirichlet.guide.ll,
                                                       bkg.idx = temp.bkg.idx, 
                                                       fs.idx = temp.fs.idx, 
                                                       guide.efficiency = data.par.list$guide_efficiency, 
                                                       repl.disp = temp.repl.disp,
                                                       model.disp = analysis.par.list$model_dispersion)
    
    orig.fs.par[[i]] <- fs.par
    
    if(analysis.par.list$hyper_adj < 1){
      fs.par <- adjust_hypers(bkg.par, fs.par, analysis.par.list$hyper_adj)
    }
    
    temp.bkg.alpha <- c(1, bkg.par) / sum(c(1, bkg.par))
    temp.fs.alpha <- c(1, fs.par) / sum(c(1, fs.par))
    
    final.dirichlet.pars$bkg_alpha[[i]] <- bkg.par
    final.dirichlet.pars$FS_alpha[[i]] <- fs.par
    final.dirichlet.pars$dispersion[[i]] <- temp.repl.disp
    
    temp.hyper <- reparameterize_hypers(temp.bkg.alpha, temp.fs.alpha, temp.repl.disp, data.par.list$guide_efficiency, analysis.par.list$model_dispersion)
    final.alpha$alpha0[[i]] <- temp.hyper$alpha0
    final.alpha$alpha1[[i]] <- temp.hyper$alpha1
  }
  
  final.alpha$L <- 1
  
  out.list <- list(hyper_pars = final.alpha,
                   hyper_par_components = final.dirichlet.pars,
                   init_model_ll = model.ll,
                   orig_hyper_par_components = orig.fs.par)
  return(out.list)
}


#' @title helper function, manages the reading in of the analysis parameters
#' @param input.parameter.file: location of file containing all analysis parameters
#' @param input.parameter.list: RELICS analysis parameters in list format. Default = NULL
#' @param data.file.split: logical, whether counts are separate from info, default = FALSE
#' @return list: all parameters embedded in a list
#' @export read_parameters()

read_parameters <- function(input.parameter.list, input.parameter.file, data.file.split){
  
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
  
  return(analysis.parameters)
}


#' @title given hyper parameters and guide efficiency (optional), report the per-segment ll ratio
#' @param input.hypers: list: hyper parameters
#' @param input.data: list: $seg_to_guide_lst, $guide_efficiency_scores (and $guide_efficiency if former not NULL)
#' @param input.ge: guide efficiency
#' @param input.seg.to.guide.lst: list: hyper parameters
#' @return list with parameters for running RELICS or logical FALSE if required parameter is missing
#' @export record_ll_ratio()

record_ll_ratio <- function(input.hypers, input.data, out.pars, input.parameters, file.extension){
  
  not.used <- c()
  guide.efficiency <- NULL
  if(! is.null(input.data$guide_efficiency_scores)){
    guide.efficiency <- input.data$guide_efficiency
  }
  
  ll.rt <- vector('numeric', length = length(input.data$seg_to_guide_lst))
  for(i in 1:length(input.data$data)){
    temp.data <- input.data$data[[i]]
    temp.hyper <- list(alpha0 = input.hypers$alpha0[[i]],
                       alpha1 = input.hypers$alpha1[[i]])
    
    guide.model.ll <- estimate_relics_sgrna_log_like(temp.hyper, temp.data, not.used, guide.efficiency, return.model.ll = TRUE)
    temp.ll.rt <- compute_seg_ll_ratio(guide.model.ll, input.data$seg_to_guide_lst, 
                                       input.data$guide_dist_to_seg, input.parameters$areaOfEffect_type)
    ll.rt <- ll.rt + temp.ll.rt
  }
  
  # combine with segment info
  segment.info <- input.data$seg_info
  segment.info$score <- round(ll.rt, 3)
  
  # create bedgraph
  to.bg.list <- list(seg_llRt = segment.info)
  
  # write bedgraph to output
  # out.dir <- paste0(out.pars$out_dir, '_FS', out.pars$iter)
  out.dir <- paste0(out.pars$out_dir, file.extension, '_FS', out.pars$iter)
  create_bedgraphs(to.bg.list, out.dir)
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
  
  
  if(! 'fs_signif' %in% par.given){
    out.parameter.list$fs_signif <- 0.05
  }
  if(! 'plot_raw_pp' %in% par.given){
    out.parameter.list$plot_raw_pp <- FALSE
  }
  if(! 'plot_raw_cs_pp' %in% par.given){
    out.parameter.list$plot_raw_cs_pp <- FALSE
  }
  if(! 'plot_cs_pp' %in% par.given){
    out.parameter.list$plot_cs_pp <- FALSE
  }
  if(! 'plot_segment_lls' %in% par.given){
    out.parameter.list$plot_segment_lls <- FALSE
  }
  if(! 'record_pp' %in% par.given){
    out.parameter.list$record_pp <- FALSE
  }
  if(! 'record_cs_pp' %in% par.given){
    out.parameter.list$record_cs_pp <- FALSE
  }
  if(! 'record_deltas' %in% par.given){
    out.parameter.list$record_deltas <- FALSE
  }
  
  if(! 'hyper_adj' %in% par.given){
    out.parameter.list$hyper_adj <- 1
  }
  out.parameter.list$cs_params <- list()
  if(! 'cs_sw_continous' %in% par.given){
    out.parameter.list$cs_params$cs_sw_continous <- TRUE
  } else {
    out.parameter.list$cs_params$cs_sw_continous <- input.parameter.list$cs_sw_continous
  }
  if(! 'cs_sw_size' %in% par.given){
    out.parameter.list$cs_params$cs_sw_size <- 10
  } else {
    out.parameter.list$cs_params$cs_sw_size <- input.parameter.list$cs_sw_size
  }
  if(! 'min_cs_sum' %in% par.given){
    out.parameter.list$cs_params$min_cs_sum <- 0.01
  } else {
    out.parameter.list$cs_params$min_cs_sum <- input.parameter.list$min_cs_sum
  }
  if(! 'cs_threshold' %in% par.given){
    out.parameter.list$cs_params$cs_threshold <- 0.9
  } else {
    out.parameter.list$cs_params$cs_threshold <- input.parameter.list$cs_threshold
  }
  # if(! 'max_cs_seg' %in% par.given){
  #   out.parameter.list$cs_params$max_cs_seg <- 10
  # } else {
  #   out.parameter.list$cs_params$max_cs_seg <- input.parameter.list$max_cs_seg
  # }
  # if(! 'max_mu' %in% par.given){
  #   out.parameter.list$max_mu <- 50
  # }
  # if(! 'max_size' %in% par.given){
  #   out.parameter.list$max_size <- 100
  # }
  # if(! 'record_posteriors_w_priors' %in% par.given){
  #   out.parameter.list$record_posteriors_w_priors <- FALSE
  # }
  # if(! 'record_posteriors' %in% par.given){
  #   out.parameter.list$record_posteriors <- FALSE
  # }
  # if(! 'record_segment_lls' %in% par.given){
  #   out.parameter.list$record_segment_lls <- FALSE
  # }
  if(! 'model_dispersion' %in% par.given){
    if(! 'mean_var_type' %in% par.given){
      out.parameter.list$mean_var_type <- 'spline'
      out.parameter.list$model_dispersion <- TRUE
    } else {
      out.parameter.list$model_dispersion <- FALSE
    }
  }
  if('mean_var_type' %in% par.given){
    if(out.parameter.list$model_dispersion & !out.parameter.list$mean_var_type == 'spline'){
      out.parameter.list$mean_var_type <- 'spline'
    }
  }
  if(!'mean_var_type' %in% par.given){
    out.parameter.list$mean_var_type <- 'spline'
  }
  if(! 'estimateSpline' %in% par.given){
    out.parameter.list$estimateSpline <- FALSE
  }
  if(! 'nr_disp_bins' %in% par.given){
    out.parameter.list$nr_disp_bins <- 20
  }
  if(! 'pp_calculation' %in% par.given){
    out.parameter.list$pp_calculation <- 'v4'
  }
  if(! 'out_dir' %in% par.given){
    out.parameter.list$out_dir <- getwd()
  }
  if(! 'convergence_tol' %in% par.given){
    out.parameter.list$convergence_tol <- 0.1
  }
  if(! 'fix_hypers' %in% par.given){
    out.parameter.list$fix_hypers <- TRUE
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
  # if(! 'iterative_hyper_est' %in% par.given){
  #   out.parameter.list$iterative_hyper_est <- FALSE
  # }
  if(! 'nr_segs' %in% par.given){
    out.parameter.list$nr_segs <- 10
  }
  if(! 'geom_p' %in% par.given){
    #probability of success for geometric distribution
    out.parameter.list$geom_p <- 0.1
  }
  if(! 'fs_correlation_cutoff' %in% par.given){
    out.parameter.list$fs_correlation_cutoff <- 0.1
  }
  if(! 'min_fs_pp' %in% par.given){
    out.parameter.list$min_fs_pp <- 0.1
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
  if(! 'pool_names' %in% par.given){
    out.parameter.list$pool_names <- NULL
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
  } else {
    out.parameter.list$guide_efficiency_scores <- NULL
    out.parameter.list$fixed_ge_coeff <- TRUE
  }
  
  # l2fc related parameters
  if('l2fc_groups' %in% par.given){
    out.parameter.list$l2fc_groups <- input.parameter.list$l2fc_groups
    if('l2fc_slidWind' %in% par.given){
      out.parameter.list$l2fc_slidWind <- input.parameter.list$l2fc_slidWind
    } else {
      out.parameter.list$l2fc_slidWind <- TRUE
    }
    
    if(out.parameter.list$l2fc_slidWind){
      if('l2fc_slidWind_guideNr' %in% par.given){
        out.parameter.list$l2fc_slidWind_guideNr <- input.parameter.list$l2fc_slidWind_guideNr
      } else {
        out.parameter.list$l2fc_slidWind_guideNr <- 10
      }
      
      if('l2fc_slidWind_maxGap' %in% par.given){
        out.parameter.list$l2fc_slidWind_maxGap <- input.parameter.list$l2fc_slidWind_maxGap
      } else {
        out.parameter.list$l2fc_slidWind_maxGap <- 3000
      }
    }
  }
  
  # type of area of effect. options are 'normal' (default), 'uniform', 'slab_and_spike'
  if(! 'areaOfEffect_type' %in% par.given){
    out.parameter.list$areaOfEffect_type <- 'normal'
  }
  
  minimum.parameters <- c()
  if(data.file.split){
    minimum.parameters <- c('dataName','repl_groups', 'CountFileLoc', 'sgRNAInfoFileLoc', 'crisprSystem', 'FS0_label')
  } else {
    minimum.parameters <- c('dataName','repl_groups', 'DataInputFileLoc', 'crisprSystem', 'FS0_label')
  }
  
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
    
    # # for the crisprSystem, check if effect range is given, else use default
    if(minimum.parameters[i] == 'crisprSystem'){
      if(! out.parameter.list$crisprSystem %in% c('CRISPRi', 'CRISPRa', 'Cas9', 'CRISPRcas9', 'dualCRISPR', 'dualCRISPRi', 'dualCRISPRa')){
        print(paste0('Please provide a valid CRISPR system: CRISPRi, CRISPRa, Cas9 (or CRISPRcas9), dualCRISPR, dualCRISPRi or dualCRISPRa'))
        missing.parameters <- TRUE
      }
    }
  }
  
  if(out.parameter.list$areaOfEffect_type == 'uniform'){
    if(! 'crisprEffectRange' %in% par.given){
      if(out.parameter.list$crisprSystem == 'CRISPRi' | out.parameter.list$crisprSystem == 'CRISPRa'){
        out.parameter.list$crisprEffectRange <- 500
      } else if(out.parameter.list$crisprSystem %in% c('Cas9','CRISPRcas9') ){
        out.parameter.list$crisprEffectRange <- 10
      } else if(out.parameter.list$crisprSystem == 'dualCRISPR'){
        out.parameter.list$crisprEffectRange <- 0
      } else {
        print("Error: please specify a valid CRISPR system for the 'uniform' area of effect (Cas9 (or CRISPRcas9), CRISPRi, CRISPRa, dualCRISPR)")
        break()
      }
    }
  }
  
  # set default area of effects sd and crisprEffectRange if not provided for normal AoE
  if (out.parameter.list$areaOfEffect_type == 'normal') {
      
      # set defaults for area of effect sd based on CRISPR system
      if (! 'normal_areaOfEffect_sd' %in% par.given) {
          if (out.parameter.list$crisprSystem %in% c('CRISPRi', 'CRISPRa', 'dualCRISPRi', 'dualCRISPRa')) {
              out.parameter.list$normal_areaOfEffect_sd <- 170
          } else if (out.parameter.list$crisprSystem %in% c('Cas9', 'CRISPRcas9', 'dualCRISPR')) {
              out.parameter.list$normal_areaOfEffect_sd <- 8.5
          } else {
              print("ERROR: please specify a valid CRISPR system for 'normal' are of effect")
              print("Cas9, CRISPRcas9, CRISPRi, CRISPRa, dualCRISPR, dualCRISPRa")
              break()
          }
      }
      
      # set defaults for crisprEffectRange based on CRISPR system
      if (! 'crisprEffectRange' %in% par.given) {
          if (out.parameter.list$crisprSystem %in% c('CRISPRi', 'CRISPRa', 'dualCRISPRi', 'dualCRISPRa')) {
              out.parameter.list$crisprEffectRange <- 415
          } else if (out.parameter.list$crisprSystem %in% c('Cas9', 'CRISPRcas9', 'dualCRISPR')) {
              out.parameter.list$crisprEffectRange <- 21
          } else {
              print("ERROR: please specify a valid CRISPR system for 'normal' are of effect")
              print("Cas9, CRISPRcas9, CRISPRi, CRISPRa, dualCRISPR, dualCRISPRa")
              break()
          }
      }
  }
  
  # set default parameters for slab and spike AoE
  if (out.parameter.list$areaOfEffect_type == "slab_and_spike") {
    
    # set default flanking distance to 500bp
    if (! 'flanking.distance' %in% par.given) {
      out.parameter.list$flanking.distance <- 500
    }
    
    # set default slab AoE to 0.1
    if (! 'slab.aoe' %in% par.given) {
      out.parameter.list$slab.aoe <- 0.1
    }
  }
  
          
#   if(out.parameter.list$areaOfEffect_type == 'normal'){
#     if(! 'normal_areaOfEffect_sd' %in% par.given){
#       if(out.parameter.list$crisprSystem %in% c('CRISPRi', 'CRISPRa', 'dualCRISPRi', 'dualCRISPRa') ){
#         out.parameter.list$normal_areaOfEffect_sd <- 170
#         if(! 'crisprEffectRange' %in% par.given){
#           out.parameter.list$crisprEffectRange <- 415
#         }
#       } else if(out.parameter.list$crisprSystem %in% c('Cas9','CRISPRcas9', 'dualCRISPR') ){
#         out.parameter.list$normal_areaOfEffect_sd <- 8.5
#         if(! 'crisprEffectRange' %in% par.given){
#           out.parameter.list$crisprEffectRange <- 21
#         }
#       } else {
#         print("Error: please specify a valid CRISPR system for the 'uniform' area of effect (Cas9 (or CRISPRcas9), CRISPRi, CRISPRa, dualCRISPR)")
#         break()
#       }
#     }
#   }
  
  if(out.parameter.list$crisprSystem == 'dualCRISPR' & (! 'deletionProb' %in% par.given)){
    out.parameter.list$deletionProb <- 0.1
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
    
    if('fs_signif' == parameter.id){
      out.parameter.list$fs_signif <- as.numeric(strsplit(parameter,':')[[1]][2])
    }
    if('plot_raw_pp' == parameter.id){
      out.parameter.list$plot_raw_pp <- as.logical(strsplit(parameter,':')[[1]][2])
    }
    if('plot_raw_cs_pp' == parameter.id){
      out.parameter.list$plot_raw_cs_pp <- as.logical(strsplit(parameter,':')[[1]][2])
    }
    if('plot_cs_pp' == parameter.id){
      out.parameter.list$plot_cs_pp <- as.logical(strsplit(parameter,':')[[1]][2])
    }
    if('plot_segment_lls' == parameter.id){
      out.parameter.list$plot_segment_lls <- as.logical(strsplit(parameter,':')[[1]][2])
    }
    if('record_pp' == parameter.id){
      out.parameter.list$record_pp <- as.logical(strsplit(parameter,':')[[1]][2])
    }
    if('record_cs_pp' == parameter.id){
      out.parameter.list$record_cs_pp <- as.logical(strsplit(parameter,':')[[1]][2])
    }
    if('record_deltas' == parameter.id){
      out.parameter.list$record_deltas <- as.logical(strsplit(parameter,':')[[1]][2])
    }
    if('expected_fs_nr' == parameter.id){
      out.parameter.list$expected_fs_nr <- as.numeric(strsplit(parameter,':')[[1]][2])
    }
    if('max_fs_nr' == parameter.id){
      out.parameter.list$max_fs_nr <- as.numeric(strsplit(parameter,':')[[1]][2])
    }
    
    if('hyper_adj' == parameter.id){
      out.parameter.list$hyper_adj <- as.numeric(strsplit(parameter,':')[[1]][2])
    }
    if('cs_sw_continous' == parameter.id){
      out.parameter.list$cs_sw_continous <- as.logical(strsplit(parameter,':')[[1]][2])
    }
    if('cs_sw_size' == parameter.id){
      out.parameter.list$cs_sw_size <- as.numeric(strsplit(parameter,':')[[1]][2])
    }
    if('min_cs_sum' == parameter.id){
      out.parameter.list$min_cs_sum <- as.numeric(strsplit(parameter,':')[[1]][2])
    }
    if('cs_threshold' == parameter.id){
      out.parameter.list$cs_threshold <- as.numeric(strsplit(parameter,':')[[1]][2])
    }
    # if('max_cs_seg' == parameter.id){
    #   out.parameter.list$max_cs_seg <- as.numeric(strsplit(parameter,':')[[1]][2])
    # }
    # if('max_mu' == parameter.id){
    #   out.parameter.list$max_mu <- as.numeric(strsplit(parameter,':')[[1]][2])
    # }
    # if('max_size' == parameter.id){
    #   out.parameter.list$max_size <- as.numeric(strsplit(parameter,':')[[1]][2])
    # }
    # if('record_posteriors' == parameter.id){
    #   out.parameter.list$record_posteriors <- as.logical(strsplit(parameter,':')[[1]][2])
    # }
    # if('record_segment_lls' == parameter.id){
    #   out.parameter.list$record_segment_lls <- as.logical(strsplit(parameter,':')[[1]][2])
    # }
    if('fs_prior' == parameter.id){
      out.parameter.list$fs_prior <- lapply(strsplit(strsplit(strsplit(parameter,':')[[1]][2],';')[[1]],','), as.numeric)
      names(out.parameter.list$repl_groups) <- c('mean', 'sd')
    }
    if('model_dispersion' == parameter.id){
      out.parameter.list$model_dispersion <- as.logical(strsplit(parameter,':')[[1]][2])
    }
    if('repl_spline_df' == parameter.id){
      out.parameter.list$repl_spline_df <- lapply(strsplit(strsplit(strsplit(parameter,':')[[1]][2],';')[[1]],','), as.numeric)
      names(out.parameter.list$repl_spline_df) <- paste0('df_repl_',c(1:length(out.parameter.list$repl_groups)))
    }
    if('estimateSpline' == parameter.id){
      out.parameter.list$estimateSpline <- as.logical(strsplit(parameter,':')[[1]][2])
    }
    if('nr_disp_bins' == parameter.id){
      out.parameter.list$nr_disp_bins <- as.numeric(strsplit(parameter,':')[[1]][2])
    }
    if('pp_calculation' == parameter.id){
      out.parameter.list$pp_calculation <- strsplit(parameter,':')[[1]][2]
    }
    if('out_dir' == parameter.id){
      out.parameter.list$out_dir <- strsplit(parameter,':')[[1]][2]
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
    # if('iterative_hyper_est' == parameter.id){
    #   out.parameter.list$iterative_hyper_est <- as.logical(strsplit(parameter,':')[[1]][2])
    # }
    if('nr_segs' == parameter.id){
      out.parameter.list$nr_segs <- as.numeric(strsplit(parameter,':')[[1]][2])
    }
    if('geom_p' == parameter.id){
      #probability of success for geometric distribution
      out.parameter.list$geom_p <- as.numeric(strsplit(parameter,':')[[1]][2])
      #geom.norm.factr <- pgeom(nr.segs, geom.p)
    }
    if('fs_correlation_cutoff' == parameter.id){
      out.parameter.list$fs_correlation_cutoff <- as.numeric(strsplit(parameter,':')[[1]][2])
    }
    if('min_fs_pp' == parameter.id){
      out.parameter.list$min_fs_pp <- as.numeric(strsplit(parameter,':')[[1]][2])
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
    if('pool_names' == parameter.id){
      out.parameter.list$pool_names <- lapply(strsplit(strsplit(strsplit(parameter,':')[[1]][2],';')[[1]],','), as.numeric)
      names(out.parameter.list$pool_names) <- paste0('repl_',c(1:length(out.parameter.list$pool_names)))
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
  
  
  # l2fc parameters
  if('l2fc_groups' == parameter.id){
    out.parameter.list$l2fc_groups <- lapply(strsplit(strsplit(strsplit(parameter,':')[[1]][2],';')[[1]],','), as.numeric)
    names(out.parameter.list$l2fc_groups) <- paste0('repl_',c(1:length(out.parameter.list$l2fc_groups)))
  }
  if('l2fc_slidWind' == parameter.id){
    out.parameter.list$l2fc_slidWind <- as.logical(strsplit(parameter,':')[[1]][2])
  }
  if('l2fc_slidWind_guideNr' == parameter.id){
    out.parameter.list$l2fc_slidWind_guideNr <- as.numeric(strsplit(parameter,':')[[1]][2])
  }
  if('l2fc_slidWind_maxGap' == parameter.id){
    out.parameter.list$l2fc_slidWind_maxGap <- as.numeric(strsplit(parameter,':')[[1]][2])
  }
  
  # area of effect parameters
  if('areaOfEffect_type' == parameter.id){
    out.parameter.list$areaOfEffect_type <- as.numeric(strsplit(parameter,':')[[1]][2])
  }
  if('normal_areaOfEffect_sd' == parameter.id){
    out.parameter.list$normal_areaOfEffect_sd <- as.numeric(strsplit(parameter,':')[[1]][2])
  }
  
  if('mean_var_type' == parameter.id){
    out.parameter.list$mean_var_type <- strsplit(parameter,':')[[1]][2]
  }
  
  if('deletionProb' == parameter.id){
    out.parameter.list$deletionProb <- as.numeric(strsplit(parameter,':')[[1]][2])
  }
  
  return(out.parameter.list)
}


#' @title helper function, manages the computing of the hyperparameters
#' @param analysis.parameters: list, contains necessary flags and elements
#' @param data.setup: data fromatted for RELICS analysis
#' @return list: all parameters embedded in a list
#' @export compute_hyperparameters()

compute_hyperparameters <- function(analysis.parameters, data.setup){
  
  print('Estimating initial hyperparameters ...')
  #labels for background
  background.labels <- c()
  if(analysis.parameters$background_label_specified){
    background.labels <- analysis.parameters$background_label
  } else {
    background.labels <- analysis.parameters$labelHierarchy[-which(analysis.parameters$labelHierarchy %in% analysis.parameters$FS0_label)]
  }
  
  fs0.alphas <- c()
  # either dispersion is modeled via splines or independent of the total counts
  if(analysis.parameters$model_dispersion){
    # # if no spline has yet been specified, identify the optimal one
    # # To Do 5-fold crossvalidation
    repl.splines <- c() # spline df for each replicate
    if(analysis.parameters$estimateSpline){
      print('not implemented yet')
      # repl.splines <- identify_splines()
      break()
    } else {
      repl.spline.df <- analysis.parameters$repl_spline_df
      repl.disp <- disp_from_spline(repl.spline.df, data.setup$data,
                                    paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName),
                                    analysis.parameters$nr_disp_bins,
                                    analysis.parameters$FS0_label, 
                                    data.setup$guide_info)
      analysis.parameters$repl_disp <- repl.disp
    }
    
    fs0.alphas <- estimate_dirichlet_proportions(data.setup,
                                                 analysis.parameters,
                                                 input.repl.pools = analysis.parameters$repl_groups,
                                                 fs0.label = analysis.parameters$FS0_label)
  } else {
    fs0.alphas <- estimate_hyper_parameters(data.setup,
                                            analysis.parameters,
                                            input.repl.pools = analysis.parameters$repl_groups,
                                            fs0.label = analysis.parameters$FS0_label)
  }
  
  analysis.parameters$hyper_pars <- fs0.alphas$hyper_pars
  analysis.parameters$hyper_par_components <- fs0.alphas$hyper_par_components
  analysis.parameters$init_model_ll <- fs0.alphas$init_model_ll
  
  # if 'intermediate' signal selected:
  if(analysis.parameters$hyper_adj < 1){
    record_orig_fs0_alphas(fs0.alphas$orig_hyper_par_components, analysis.parameters, analysis.parameters$pool_names) # record just the fs.alphas
  }
  
  return(analysis.parameters)
}


#' @title Record the original hyperparameters 
#' @param input.fs0.list: list of FS0 hyperparameters
#' @param analysis.parameters: list of all analysis paramters
#' @param pool.names: names of the different pools
#' @return list of final per-layer posteriors
#' @export record_orig_fs0_alphas()

record_orig_fs0_alphas <- function(input.fs0.list, analysis.parameters, pool.names){
  
  out.dir <- paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName, '_orig_FS0_hyperPars.csv')
  
  out.alpha.df <- c()
  
  total.rows <- length(input.fs0.list)
  total.cols <- max(c( unlist(lapply(input.fs0.list, function(x){length(x) + 1}))))
  
  alpha.names <- c(paste(rep('FS', length(input.fs0.list)), 'r', c(1:length(input.fs0.list)), sep = '_'))
  
  alpha.matrix <- matrix(0, nrow = total.rows, ncol = total.cols)
  
  if(! is.null(pool.names)){
    temp.pool.names <- unique(unlist(pool.names))
    
    for(i in 1:length(input.fs0.list)){
      alpha1.scores <- rep(0, length(length(temp.pool.names)))
      
      temp.fs.alpha <- c(1,input.fs0.list[[i]]) / sum(c(1, input.fs0.list[[i]]))
      
      alpha1.scores[match(pool.names[[i]], temp.pool.names)] <- round(temp.fs.alpha, 3)
      
      alpha.matrix[i,] <- alpha1.scores
      
      out.alpha.df <- cbind(alpha.names, alpha.matrix)
      
      colnames(out.alpha.df) <- c('hyperPar_type', temp.pool.names)
      
    }
    
  } else {
    
    for(i in 1:length(input.fs0.list)){
      temp.fs0.alpha <- c(1,input.fs0.list[[i]]) / sum(c(1, input.fs0.list[[i]]))
      alpha.matrix[i, c(1:length(temp.fs0.alpha))] <- round(temp.fs0.alpha, 3)
    }
    
    out.alpha.df <- cbind(alpha.names, alpha.matrix)
    
    colnames(out.alpha.df) <- c('hyperPar_type', paste0('pool', c(1:(nrow(out.alpha.df) - 2))))
  }
  
  write.csv(out.alpha.df, file = out.dir, row.names = F, quote = F)
  
}


#' @title plot the ll results
#' @param input.fs.results: list, with all the results
#' @param fs.nr: how many FS have been calculated
#' @param out.dir: directory to which to write
#' @return list: $pp_ordered, $layer_ll_ordered
#' @export set_up_fs_priors()

set_up_fs_priors <- function(input.params, input.data.setup){
  
  param.names <- names(input.params)
  out.params <- input.params
  if(! 'expected_fs_nr' %in% param.names){
    data.segs <- input.data.setup$seg_info
    avg.seg.size <- input.params$seg_dist
    
    total.range <- nrow(data.segs) * avg.seg.size
    expected.fs.nr <- round(total.range / 50000)
    out.params$expected_fs_nr <- expected.fs.nr
    print(paste0("Setting 'expected_fs_nr' to ",out.params$expected_fs_nr) )
  }
  if(! 'max_fs_nr' %in% param.names){
    out.params$max_fs_nr <- out.params$expected_fs_nr + max(3, round(0.3 * out.params$expected_fs_nr))
    print(paste0("Setting 'max_fs_nr' to ",out.params$max_fs_nr) )
  }
  if(out.params$expected_fs_nr > out.params$max_fs_nr ){
    out.params$max_fs_nr <- out.params$expected_fs_nr + 3
    print(paste0("Adjusting 'max_fs_nr' to ",out.params$max_fs_nr) )
  }
  if(! 'fs_prior' %in% param.names){
    out.params$fs_prior <- dbinom(c(0:(out.params$max_fs_nr)), 
                                  size = out.params$max_fs_nr, 
                                  prob = out.params$expected_fs_nr / out.params$max_fs_nr)
  }
  if(! 'fs_ll_signif' %in% param.names){
    out.params$fs_ll_signif <- qchisq(out.params$fs_signif, df = 1, lower.tail = F)
  }
  
  return(out.params)
}


#' @title Run RELICS v7. Iteratively place FS. Calculate the CS by selecting best FS length. Terminate based on max FS nr or min PP
#' @param input.data: list: $guide_to_seg_lst, $data, $fs0_idx, $overlapMat, $guide_to_seg_lst, $seg_to_guide_lst
#' @param final.layer.nr: number of layers to go up to
#' @param out.dir: directory to which all files are to be written. Default = NULL
#' @param adjust.tol: whether convergence tolerance should be adjusted dynamically as layers are computed (increasingly more rigorous)
#' @param fix.hypers: whether the input parameters remain unchanged or are altered after calculation of posteriors. Default = FALSE
#' @param input.hypers: list, $alpha0, $alpha1, both elements contain lists of length wequal to nr. replicates (input.hypers$alpha0[[1]])
#' @param nr.segs, number of segments to consider for the length of a regulatory element
#' @param geom.p: proababilty of the genometric distribution to penalize for enhancers of increasing length
#' @param input.min.rs.pp: minimum posterior required to be part of a regulatory set
#' @param auto.stop: whether or not computations should be stopped after recomended stopping point
#' @param record.all.fs: logical, if information of all intermediate FS should be recorded, in addition to the final set of FS
#' @param one.dispersion: whether there should be one or 2 dispersions estimated for the 2 hyper parameters
#' @param recompute.fs0: logical, whether FS0 is to be set to 0 or keep the initial assignment
#' @param local.max: logical, whether a local maximum should be computed
#' @param local.max.range: number of segments to include in addition to the the ones with highest PP (one way, so multiply by 2 for total nr of segs)
#' @param pool.names, if not NULL, the names of the pools to be used when recording output
#' @param mean.var.type: type of mean-variance relationship
#' @param analysis.parameters: list, contains various elements used for analysis
#' @param guide.dist.to.seg, list of lists, pre-computed distances of guides to FS
#' @return list of final per-layer posteriors
#' @export run_RELICS_7()

run_RELICS_7 <- function(input.data, final.layer.nr, out.dir = NULL,
                         fix.hypers = FALSE,
                         input.hypers = NULL,
                         input.hyper.components,
                         nr.segs = 10,
                         geom.p = 0.1,
                         input.min.rs.pp = 0.1,
                         auto.stop = TRUE,
                         record.all.fs,
                         input.convergence.tol = 0.01,
                         adjust.tol = TRUE,
                         one.dispersion = TRUE,
                         recompute.fs0,
                         local.max, local.max.range,
                         pool.names,
                         mean.var.type,
                         analysis.parameters,
                         guide.dist.to.seg){
  
  relics.hyper <- input.hypers #c()
  hyper.components <- input.hyper.components
  
  relics.param <- initialize_relics_mtxs(relics.hyper, input.data, local.max, recompute.fs0)
  
  coverg.tol <- rep(input.convergence.tol, final.layer.nr)
  
  fs.result.lists <- initialize_fs_result_list(relics.param$deltas, 
                                               relics.param$pps, 
                                               input.data, 
                                               input.data$guide_efficiency,
                                               relics.hyper, 
                                               fix.hypers)
  
  for(i in 1:analysis.parameters$max_fs_nr){
    print(paste0('Computing FS: ', i))
    
    fs.data <- compute_FS_placement(input.param = relics.param,
                                    input.hyper = relics.hyper,
                                    input.data.list = input.data,
                                    input.tol = coverg.tol[i],
                                    nr.segs, geom.p,
                                    min.pp = input.min.rs.pp, input.data$guide_efficiency,
                                    one.dispersion,
                                    local.max, local.max.range, analysis.parameters$areaOfEffect_type,
                                    guide.dist.to.seg,
                                    fix.hypers,
                                    fs.result.lists,
                                    analysis.parameters$fs_prior,
                                    analysis.parameters$cs_params,
                                    analysis.parameters$fs_ll_signif)
    
    # fs.data already contains the most recent additions
    if(record.all.fs){
      
      record_relics_results(fs.data, analysis.parameters, input.data, i, relics.hyper, hyper.components, '')
      
    }
    
    if(sum(fs.data$deltas[i+1,]) < 1){
      print('No more FS to detect. Max. FS nr. not reached.')
      record_relics_results(fs.data, analysis.parameters, input.data, i, relics.hyper, hyper.components, '_finalFS')
      break()
    }
    if(i == final.layer.nr){
      print('Max. FS nr. reached.')
      record_relics_results(fs.data, analysis.parameters, input.data, i, relics.hyper, hyper.components, '_finalFS')

      break()
    }
    
    # update the results
    fs.result.lists <- fs.data
    
    # prep the parameters for the next iteration (both the posteriors and hypers)
    relics.hyper$L <- relics.hyper$L + 1
    relics.param$pps <- rbind(fs.result.lists$posteriors, rep(0, ncol(relics.param$pps)))
    relics.param$deltas <- rbind(fs.result.lists$deltas, rep(0, ncol(relics.param$deltas)))
    relics.param$cs_pps <- rbind(fs.result.lists$cs_posteriors, rep(0, ncol(relics.param$pps)))
    
    # if the hyper parameters are re-estimated after convergence of the posteriors
    # or if the guide efficiency scores are re-estimated
    if(! fix.hypers){
      
      relics.hyper.list <- c()
      if(analysis.parameters$model_dispersion){
        relics.hyper.list <- recompute_dirichlet_hypers_w_deltas(relics.param, 
                                                                 input.data$guide_to_seg_lst, 
                                                                 hyper.components, 
                                                                 input.data$data, 
                                                                 input.data$guide_efficiency, 
                                                                 relics.hyper,
                                                                 analysis.parameters)
      } else {
        relics.hyper.list <- recompute_hyper_parameters_w_deltas(relics.param,
                                                                 relics.hyper,
                                                                 hyper.components,
                                                                 input.data$data,
                                                                 input.data$guide_to_seg_lst,
                                                                 input.data$guide_efficiency,
                                                                 analysis.parameters)
      }
      
      relics.hyper <- relics.hyper.list$hyper_pars
      hyper.components <- relics.hyper.list$hyper_par_components
    }
    
    if(! is.null(input.data$guide_efficiency_scores)){
      if(! input.data$fixed_ge_coeff){
        ge.list <- recompute_ge_coefficients_w_deltas(relics.param,
                                                      hyper.components, #relics.hyper,
                                                      input.data$data,
                                                      input.data$guide_to_seg_lst,
                                                      input.data$guide_efficiency_scores,
                                                      input.data$ge_coeff)
        
        input.data$guide_efficiency <- ge.list$guide_efficiency
        input.data$ge_coeff <- ge.list$ge_coeff
        
        relics.hyper <- incorporate_ge(relics.hyper,
                                           hyper.components,
                                           input.data$data,
                                           input.data$guide_efficiency,
                                           analysis.parameters)
      }
    }
  }
  
}


#' @title use re-computed guide-efficiency scores to adjust the alphas
#' @param hyper.components: individual components of the hyper parameters
#' @param data: list, contains the elements of the processed and formatted data
#' @param guide.efficiency: either a vector of guide efficiency or NULL
#' @param hyper: hyperparameters, already divided into alpha0 and alpha1
#' @param analysis.parameters: list of analysis paramters
#' @return alphas
#' @export incorporate_ge()

incorporate_ge <- function(hyper, hyper.components, data, guide.efficiency, analysis.parameters){
  
  for(i in 1:length(data)){
    
    temp.bkg.alpha <- c(1, hyper.components$bkg_alpha[[i]]) / sum(c(1, hyper.components$bkg_alpha[[i]]))
    temp.fs.alpha <- c(1, hyper.components$FS_alpha[[i]]) / sum(c(1, hyper.components$FS_alpha[[i]]))
    
    if(analysis.parameters$model_dispersion){
      temp.repl.disp <- hyper.components$dispersion[[i]]
      temp.hyper <- reparameterize_hypers(temp.bkg.alpha, temp.fs.alpha, temp.repl.disp, 
                                          guide.efficiency, analysis.parameters$model_dispersion)
      hyper$alpha0[[i]] <- temp.hyper$alpha0
      hyper$alpha1[[i]] <- temp.hyper$alpha1
    } else {
      temp.disp <- hyper.components$bkg_dispersion[[i]]
      temp.hyper <- reparameterize_hypers(temp.bkg.alpha, temp.fs.alpha, temp.disp, guide.efficiency, TRUE)
      hyper$alpha0[[i]] <- temp.hyper$alpha0
      hyper$alpha1[[i]] <- temp.hyper$alpha1
    }
  }
  
  return(hyper)
}


#' @title Set initial hyper parameters, generate the fs posterior matrix and delta matrix
#' @param hyper: hyperparameters
#' @param in.data.list: list, $seg_info, $fs0_idx
#' @param local.max: whether local.max is computed
#' @param recompute.fs0: logical, whether FS0 should be recomputed
#' @return list
#' @export initialize_relics_mtxs()

initialize_relics_mtxs <- function(hyper, in.data.list, local.max, recompute.fs0) {
  param <- list()
  
  #  posterior probabilities
  param$pps <- matrix(0, nrow=hyper$L+1, ncol= nrow(in.data.list$seg_info)) #nrow(in.data.list$overlapMat))
  
  # first row is special and contains posterior probs of known regulatory elements
  param$pps[1, in.data.list$fs0_idx] <- 1
  
  param$deltas <- param$pps
  param$cs_pps <- param$pps
  
  param$ll_rt <- list()
  
  if(local.max){
    param$ll_rt <- matrix(0, nrow=hyper$L+1, ncol= nrow(in.data.list$seg_info))
  }
  
  if(recompute.fs0){
    param$pps[1,] <- 0
    param$deltas[1,] <- 0
    param$cs_pps[1,] <- 0
  }
  
  return(param)
}


#' @title Display the credible sets along with cumulative posteriors. Formerly 'display_mrvr_norm_tiff'
#' @param input.L matrix of posteriors
#' @param input.labels: labelling of all segments
#' @param tiff.name: name for plot
#' @param seg.info: contains info for each segment, specicially chromosome information
#' @param delta.mtx: matrix containin the location of the FS
#' @return One plot for each selection step, credible sets colored in purple
#' @export plot_relics_pp_as_tiff()

plot_relics_pp_as_tiff <- function(input.L, input.labels, tiff.name, seg.info, delta.mtx){
  
  cs.col.orig <- rep("darkgrey", ncol(input.L))
  
  if(length(unique(seg.info$chrom)) > 1){
    even.chroms <- unique(seg.info$chrom)[c(FALSE,TRUE)]
    even.chroms.idx <- which(seg.info$chrom %in% even.chroms)
    cs.col.orig[even.chroms.idx] <- 'lightgrey'
  }
  
  # attempt to put 8 per page
  if(nrow(input.L) <= 8){
    tiff(paste0(tiff.name, '.tiff'), width = 11520 / 2, height = 5760 / 2, units = "px", res = 400)
    par(mfrow = c(3,3))
    
    for(i in 1:nrow(input.L)){
      cs.col <- cs.col.orig
      cs.col[which(delta.mtx[i,] == 1)] <- "purple"
      
      plot(input.L[i,], pch=21,
           main=paste("FS", i-1, ", Nr. segments in FS: ", length(which(delta.mtx[i,] == 1)), sep=""), col=cs.col,
           ylab = 'PP', xlab = 'Genome Segment')
    }
    p.sum <- c()
    if(nrow(input.L) > 2){
      p.sum <- colSums(input.L[1:nrow(input.L),])
    } else {
      p.sum <- colSums(input.L[1:2,])
    }
    p.sum[p.sum > 1] <- 1
    plot(p.sum, pch=21,
         main="Sum of Posteriors",
         col=cs.col.orig, bg = input.labels, ylab = 'PP', xlab = 'Genome Segment')
    
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
            p.sum <- colSums(input.L[1:2,])
          }
          p.sum[p.sum > 1] <- 1
          plot(p.sum, pch=21,
               main="Sum of Posteriors",
               col=cs.col.orig, bg = input.labels, ylab = 'PP', xlab = 'Genome Segment')
          #dev.off()
          break()
        } else{
          
          cs.col <- cs.col.orig #rep("black", ncol(input.L))
          cs.col[which(delta.mtx[row.index,] == 1)] <- "purple"
          
          # ', RS-pp: ', temp.cs.pp,
          plot(input.L[row.index,], pch=21,
               main=paste("FS", row.index-1, ", Nr. segments in FS: ", length(which(delta.mtx[row.index,] == 1)),sep=""), col=cs.col,
               ylab = 'PP', xlab = 'Genome Segment')
          if(i == 8){
            p.sum <- c()
            if(nrow(input.L) > 2){
              p.sum <- colSums(input.L[1:nrow(input.L),])
            } else {
              p.sum <- colSums(input.L[1:2,])
            }
            p.sum[p.sum > 1] <- 1
            plot(p.sum, pch=21,
                 main="Sum of Posteriors",
                 col=cs.col.orig, bg = input.labels, ylab = 'PP', xlab = 'Genome Segment') #"black"
          }
        }
      }
      dev.off()
    }
    
  }
  
}


#' @title plot the ll results
#' @param input.fs.results: list, with all the results
#' @param fs.nr: how many FS have been calculated
#' @param out.dir: directory to which to write
#' @param file.extension: additional file labels
#' @return list: $pp_ordered, $layer_ll_ordered
#' @export plot_fs_model_ll_stats()

plot_fs_model_ll_stats <- function(out.dir, file.extension, fs.nr, input.fs.results){
  
  pdf(paste0(out.dir, file.extension, '_k', fs.nr, '_summaryStatPlots.pdf'), useDingbats = FALSE)
  par(mfrow = c(2,2))
  
  # plot the progression of the total model ll with each FS
  plot(x = c(0:(length(input.fs.results$total_model_ll_w_prior) - 1) ), y = round(input.fs.results$total_model_ll_w_prior, 2), 
       main = 'Model ll progression', xlab = 'Nr. FS', ylab = 'll')
  
  # plot the conditional ll improvement with each FS added
  plot(x = c(1:length(input.fs.results$conditional_fs_ll_w_prior)), y = round(input.fs.results$conditional_fs_ll_w_prior, 2), 
       main = 'll contribution of each FS', xlab = 'Nr. FS', ylab = 'll')
  
  # plot the progression of the total model ll with each FS
  plot(x = c(0:(length(input.fs.results$total_model_ll) - 1) ), y = round(input.fs.results$total_model_ll, 2), 
       main = 'Raw Model ll progression', xlab = 'Nr. FS', ylab = 'll')
  
  # plot the conditional ll improvement with each FS added
  plot(x = c(1:length(input.fs.results$conditional_fs_ll)), y = round(input.fs.results$conditional_fs_ll, 2), 
       main = 'Raw ll contribution of each FS', xlab = 'Nr. FS', ylab = 'll')

  dev.off()
}


#' @title Wrapper function for recording the output from the most recent run
#' @param input.results.list: lest, each element keeps track of results from a number of FS identified
#' @param analysis.parameters: list containing all the analysis parameters
#' @param input.data: list containing the counts
#' @param relics.hyper: list, RELICS hyperparameters, alpha0 and alpha1, each list of df, one for each repl
#' @param hyper.components: hyperpparameter components (alpha0 and alpha1 proportions and dispersions)
#' @param file.extension: extension for files, either '', '_final', or '_recommendedFinal'
#' @param fs.iter: current FS interation
#' @return list with updated input.results.list
#' @export record_relics_results()

record_relics_results <- function(input.results.list, analysis.parameters, input.data, fs.iter, relics.hyper, hyper.components, file.extension){
  
  out.dir <- paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName)
  
  segment.lls <- input.results.list$segment_lls
  
  # plot the outputs
  if(analysis.parameters$plot_raw_pp){
    plot_relics_pp_as_tiff(input.results.list$posteriors,
                           input.data$segLabels,
                           paste0(out.dir, file.extension, '_k_', fs.iter, '_raw_PP'),
                           input.data$seg_info,
                           input.results.list$deltas)
  }
  if(analysis.parameters$plot_raw_cs_pp){
    plot_relics_pp_as_tiff(input.results.list$cs_posteriors,
                           input.data$segLabels,
                           paste0(out.dir, file.extension, '_k_', fs.iter, '_raw_cs_PP'),
                           input.data$seg_info,
                           input.results.list$deltas)
  }
  if(analysis.parameters$plot_cs_pp){
    plot_relics_pp_as_tiff(input.results.list$posteriors_w_prior,
                           input.data$segLabels,
                           paste0(out.dir, file.extension, '_k_', fs.iter, '_cs_PP'),
                           input.data$seg_info,
                           input.results.list$deltas)
  }
  if(analysis.parameters$plot_segment_lls){
    display_relics_segLLs_as_tiff(segment.lls,
                                  input.data$segLabels,
                                  paste0(out.dir, file.extension, '_k_', fs.iter, '_segLLs'),
                                  analysis.parameters$min_fs_pp,
                                  input.data$seg_info)
  }
  
  plot_fs_model_ll_stats(out.dir, file.extension, fs.iter, input.results.list)
  
  ll_ratio_recording(input.data, analysis.parameters, relics.hyper, fs.iter, file.extension)
  
  if(fs.iter == 1 |file.extension == '_finalFS'){
    hyperparameter_recording(list(bkg_hyper = hyper.components$bkg_alpha, fs_hyper = hyper.components$FS_alpha, bkg_disp = hyper.components$bkg_dispersion), 
                             fs.iter, hyper.components, analysis.parameters, file.extension) 
    
    # record dispersion
    record_dispersion(fs.iter, hyper.components, analysis.parameters, file.extension)
  } else if(! analysis.parameters$fix_hypers){
    hyperparameter_recording(list(bkg_hyper = hyper.components$bkg_alpha, fs_hyper = hyper.components$FS_alpha, bkg_disp = hyper.components$bkg_dispersion), 
                             fs.iter, hyper.components, analysis.parameters, file.extension) 
    # record dispersion
    record_dispersion(fs.iter, hyper.components, analysis.parameters, file.extension)
  }
  
  # guide efficiency vars:
  if(! is.null(input.data$fixed_ge_coeff) && ! input.data$fixed_ge_coeff){
    record.ge.coeff <- input.data$ge_coeff#[[fs.iter]]
    ge.coff.df <- data.frame(ge_coeff = c('beta0', paste0('beta', 1:ncol(input.data$guide_efficiency_scores))),
                             ge_coeff_scores = round(record.ge.coeff, 3))
    write.csv(ge.coff.df, file = paste0(out.dir, file.extension, '_k', fs.iter, '_ge_coeff.csv'), row.names = F)
  }
  
  write.table(input.results.list$posteriors_w_prior, file = paste0(out.dir, file.extension, '_k', fs.iter,'_cs_pp.csv'), row.names = F, col.names = F, sep = ',')
  record.pp.w.priors <- colSums(input.results.list$posteriors_w_prior)
  record.pp.w.priors[which(record.pp.w.priors > 1)] <- 1
  out.bedgraph <- input.data$seg_info
  out.bedgraph$score <- record.pp.w.priors
  to.bg.list <- list(cs_pp = out.bedgraph)
  
  if(analysis.parameters$record_pp){
    write.table(input.results.list$posteriors, file = paste0(out.dir, file.extension, '_k', fs.iter,'_raw_pp.csv'), row.names = F, col.names = F, sep = ',')
    record.pp <- colSums(input.results.list$posteriors)
    record.pp[which(record.pp > 1)] <- 1
    record.pp.bedgraph <- input.data$seg_info
    record.pp.bedgraph$score <- record.pp
    to.bg.list$raw_pp <- record.pp.bedgraph
  }
  
  if(analysis.parameters$record_cs_pp){
    write.table(input.results.list$cs_posteriors, file = paste0(out.dir, file.extension, '_k', fs.iter,'_raw_cs_pp.csv'), row.names = F, col.names = F, sep = ',')
    record.cs.pp <- colSums(input.results.list$cs_posteriors)
    record.cs.pp[which(record.cs.pp > 1)] <- 1
    record.cs.pp.bedgraph <- input.data$seg_info
    record.cs.pp.bedgraph$score <- record.cs.pp
    to.bg.list$raw_cs_pp <- record.cs.pp.bedgraph
  }
  
  if(analysis.parameters$record_deltas){
    write.table(input.results.list$deltas, file = paste0(out.dir, file.extension, '_k', fs.iter,'_raw_deltas.csv'), row.names = F, col.names = F, sep = ',')
    all.seg.fs.df.final <- extract_fs_locations(input.results.list$deltas, input.data$seg_info, analysis.parameters$min_fs_pp)
    write.table(all.seg.fs.df.final, file = paste0(out.dir, file.extension, '_k', fs.iter,'_raw_FS_locations.bed'),
                sep = '\t', quote = F, row.names = F, col.names = F)
  }
  
  # record the segments with FS probabilities above the threshold.
  write.table(input.results.list$deltas_w_prior, file = paste0(out.dir, file.extension, '_k', fs.iter,'_deltas.csv'), row.names = F, col.names = F, sep = ',')
  all.seg.fs.df.final <- extract_fs_locations(input.results.list$deltas_w_prior, input.data$seg_info, analysis.parameters$min_fs_pp)
  write.table(all.seg.fs.df.final, file = paste0(out.dir, file.extension, '_k', fs.iter,'_FS_locations.bed'),
              sep = '\t', quote = F, row.names = F, col.names = F)
  
  ll.df <- data.frame(FS = c(0:(length(input.results.list$total_model_ll)-1)),
                      fs_ll_contrib = c(0, input.results.list$conditional_fs_ll_w_prior),
                      raw_fs_ll_contrib = c(0, input.results.list$conditional_fs_ll),
                      model_ll = input.results.list$total_model_ll_w_prior,
                      raw_model_ll = input.results.list$total_model_ll,
                      cs_raw_pp = input.results.list$cs_total_pp,
                      stringsAsFactors = F)
  
  # record all the model lls
  write.csv(ll.df, file = paste0(out.dir, file.extension, '_k', fs.iter, '_model_lls.csv'), row.names = F)
  
  create_bedgraphs(to.bg.list, paste0(out.dir, file.extension, '_k', fs.iter) )

}


#' @title initialize the results list
#' @param input.delta: matrix, contains placement of FS
#' @param input.pp, matrix of posterior probabilities
#' @param dirichlet.hyper: hyper parameters
#' @param input.data.list, list: $guide_to_seg_lst, $data, $true_pos_seg, $overlapMat, $guide_to_seg_lst, $seg_to_guide_lst
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @param fix.hypers: logical, whether the hyperparameters are fixed
#' @return list: posteriors, total_model_ll, conditional_fs_ll, model_fs_ll, fs_placement_ll, fs_total_ll
#' @export record_dispersion()

record_dispersion <- function(fs.iter, hyper.components, analysis.parameters, file.extension){
  
  out.dir <- paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName)
  
  hyper.names <- names(hyper.components)
  if("dispersion" %in% hyper.names){
    disp.df <- data.frame(disp_repl1 = 1/hyper.components$dispersion[[1]], stringsAsFactors = F)
    if(length(hyper.components$dispersion) > 1){
      for(i in 2:length(hyper.components$dispersion)){
        disp.df[paste0('disp_repl',i)] <- 1/hyper.components$dispersion[[i]]
      }
    }
    
    write.csv(disp.df, file = paste0(out.dir, file.extension, '_k', fs.iter, '_disp.csv'), row.names = F)
  }
  
}


#' @title initialize the results list
#' @param input.delta: matrix, contains placement of FS
#' @param input.pp, matrix of posterior probabilities
#' @param dirichlet.hyper: hyper parameters
#' @param input.data.list, list: $guide_to_seg_lst, $data, $true_pos_seg, $overlapMat, $guide_to_seg_lst, $seg_to_guide_lst
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @param fix.hypers: logical, whether the hyperparameters are fixed
#' @return list: posteriors, total_model_ll, conditional_fs_ll, model_fs_ll, fs_placement_ll, fs_total_ll
#' @export initialize_fs_result_list()

initialize_fs_result_list <- function(input.delta, input.pp, input.data.list, guide.efficiency, 
                                           dirichlet.hyper, fix.hypers){
  
  fs.result.lists <- list()
  fs.result.lists$posteriors <- input.pp
  fs.result.lists$deltas <- input.delta
  fs.result.lists$cs_posteriors <- input.pp
  
  dirichlet.guide.ll <- compute_perGuide_fs_ll(colSums(input.delta), input.data.list$guide_to_seg_lst)
  total.model.ll <- 0
  for(i in 1:length(input.data.list$data)){
    temp.hypers <- list(alpha0 = dirichlet.hyper$alpha0[[i]], alpha1 = dirichlet.hyper$alpha1[[i]])
    temp.sgRNa.ll <- estimate_relics_sgrna_log_like(temp.hypers,
                                                    input.data.list$data[[i]],
                                                    dirichlet.guide.ll,
                                                    guide.efficiency)
    
    temp.dirichlet.ll <- sum(temp.sgRNa.ll$total_guide_ll)
    total.model.ll <- total.model.ll + temp.dirichlet.ll
  }
  
  # technically delta, but doesn't matter in this situation
  temp.pp <- input.pp
  temp.pp[1,] <- 0
  temp.pp.sum <- colSums(temp.pp)
  temp.pp.sum[temp.pp.sum > 1] <- 1
  condit.guide.ll <- compute_perGuide_fs_ll(temp.pp.sum, input.data.list$guide_to_seg_lst)
  
  temp.cond.fs.ll <- 0
  
  for(r in 1:length(input.data.list$data)){
    temp.hypers <- list(alpha0 = dirichlet.hyper$alpha0[[r]], alpha1 = dirichlet.hyper$alpha1[[r]])
    temp.cond.sgRNa.ll <- estimate_relics_sgrna_log_like(temp.hypers,
                                                         input.data.list$data[[r]],
                                                         condit.guide.ll,
                                                         guide.efficiency)
    temp.cond.fs.ll <- temp.cond.fs.ll +  sum(temp.cond.sgRNa.ll$total_guide_ll)
    
  }
  
  fs.result.lists$total_model_ll <- total.model.ll
  fs.result.lists$conditional_fs_ll <- temp.cond.fs.ll - total.model.ll # (the ll contribution of a single FS conditional on all other FS present) - (the total model ll)
  fs.result.lists$total_model_ll_w_prior <- total.model.ll
  fs.result.lists$conditional_fs_ll_w_prior <- temp.cond.fs.ll - total.model.ll
  fs.result.lists$cs_total_pp <- 0 # the posterior probability of placing an FS at said FS
  fs.result.lists$segment_lls <- c()
  
  return(fs.result.lists)
}


#' @title Compute the placement of functional sequences and credible sets while accounting for the prior of adding fs
#' @param input.param, list containing the matrix for the posteriors
#' @param input.hyper: hyper parameters
#' @param input.data.list, list: $guide_to_seg_lst, $data, $true_pos_seg, $overlapMat, $guide_to_seg_lst, $seg_to_guide_lst
#' @param input.tol: nr of layers to create
#' @param nr.segs: max. length a functional segment is considered to have
#' @param geom.p: proababilty of the genometric distribution to penalize for enhancers of increasing length
#' @param min.pp: minimum posterior required to be part of a regulatory set
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @param one.dispersion: whether the hyper parameters should be estimated with one or two dispersions
#' @param local.max: logical, whether a local maximum should be computed
#' @param local.max.range: number of segments to include in addition to the the ones with highest PP (one way, so multiply by 2 for total nr of segs)
#' @param area.of.effect: type of area of effect. Currently either 'normal' or 'uniform'
#' @param guide.dist.to.seg, list of lists, pre-computed distances of guides to FS
#' @param fix.hypers: logical, whether hyperparameters are to be fixed
#' @param fs.ll.signif: chi-square max. value for significance
#' @return ll_tract, posterior_trace_list, alpha0_est, alpha1_est, max_iter_reached, fs_ll_rt_trace
#' @export compute_FS_placement()

compute_FS_placement <- function(input.param,
                                 input.hyper,
                                 input.data.list,
                                 input.tol = 0.001,
                                 nr.segs = 10, geom.p = 0.1,
                                 min.pp = 0.1,
                                 guide.efficiency,
                                 one.dispersion,
                                 local.max, local.max.range,  
                                 area.of.effect,
                                 guide.dist.to.seg, 
                                 fix.hypers,
                                 fs.result.lists,
                                 fs.prior,
                                 cs.params,
                                 fs.ll.signif){
  
  
  dirichlet.hyper <- input.hyper
  posteriors <- input.param$pps
  deltas <- input.param$deltas
  cs.posteriors <- input.param$cs_pps
  
  fs.model.list <- estimate_FS_CS(posteriors,
                                  dirichlet.hyper,
                                  input.data.list$data,
                                  input.data.list$guide_to_seg_lst,
                                  input.data.list$seg_to_guide_lst,
                                  input.data.list$next_guide_lst,
                                  nr.segs, geom.p, guide.efficiency, 
                                  area.of.effect,
                                  guide.dist.to.seg, fix.hypers,
                                  fs.result.lists,
                                  deltas, cs.params, cs.posteriors)

  model.lls <- compute_model_fit(fs.model.list$cs_posteriors, input.data.list, 
                                 guide.efficiency, dirichlet.hyper, 
                                 fs.result.lists$total_model_ll)

  ordered.fs.list <- order_FS(fs.model.list, model.lls)
  
  model.lls.w.prior <- compute_model_ll_fit_w_priors(ordered.fs.list$cs_posteriors, input.data.list, 
                                                     guide.efficiency, dirichlet.hyper, fs.prior,
                                                     fs.result.lists$total_model_ll_w_prior,
                                                     ordered.fs.list$deltas,
                                                     fs.ll.signif)
  
  out.list <- list(posteriors = ordered.fs.list$posteriors,
                   deltas = ordered.fs.list$deltas,
                   cs_posteriors = ordered.fs.list$cs_posteriors,
                   segment_lls = ordered.fs.list$segment_lls,
                   cs_total_pp = ordered.fs.list$cs_total_pp,
                   total_model_ll = ordered.fs.list$total_model_ll,
                   conditional_fs_ll = ordered.fs.list$conditional_fs_ll,
                   total_model_ll_w_prior = model.lls.w.prior$total_model_ll_w_prior,
                   conditional_fs_ll_w_prior = model.lls.w.prior$conditional_fs_ll_w_prior,
                   posteriors_w_prior = model.lls.w.prior$posteriors_w_prior,
                   deltas_w_prior = model.lls.w.prior$deltas_w_prior)
  
  return(out.list)
}


#' @title compute the model lls while accounting for the number of priors
#' @param input.pp, matrix of posterior probabilities
#' @param dirichlet.hyper: hyper parameters
#' @param input.data.list, list: $guide_to_seg_lst, $data, $true_pos_seg, $overlapMat, $guide_to_seg_lst, $seg_to_guide_lst
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @param fs.prior: prior distribution
#' @param input.model.ll: vector of ll for a given model
#' @param input.deltas: matrix of FS placement
#' @param fs.ll.signif: chi-square max. value for significance
#' @return list: total_model_ll, per_fs_ll (the model.ll increase for each FS, taking into account all other FS), model_ll_increase (the step-by-step increase in model ll with each FS)
#' @export compute_model_ll_fit_w_priors()

compute_model_ll_fit_w_priors <- function(input.pp, input.data.list, guide.efficiency, 
                                          dirichlet.hyper, 
                                          fs.prior, input.model.ll, input.deltas,
                                          fs.ll.signif){
  
  total.model.ll <- input.model.ll
  
  # first compute the weighted pp matrix
  combined.pp.mtx <- matrix(0, ncol = ncol(input.pp), nrow = nrow(input.pp))
  norm.weigth <- fs.prior[1:nrow(input.pp)] / sum(fs.prior[1:nrow(input.pp)])
  for(fs in 1:nrow(input.pp)){
    combined.pp.mtx[fs,] <- input.pp[fs,] * sum(norm.weigth[fs:length(norm.weigth)])
  }
  
  # take the sum across all weighted FS placements
  # compute the model ll
  combined.pp <- colSums(combined.pp.mtx)
  combined.pp[combined.pp > 1] <- 1
  dirichlet.guide.ll <- compute_perGuide_fs_ll(combined.pp, input.data.list$guide_to_seg_lst)
  
  temp.total.model.ll <- 0
  
  for(i in 1:length(input.data.list$data)){
    temp.hypers <- list(alpha0 = dirichlet.hyper$alpha0[[i]], alpha1 = dirichlet.hyper$alpha1[[i]])
    temp.sgRNa.ll <- estimate_relics_sgrna_log_like(temp.hypers,
                                                    input.data.list$data[[i]],
                                                    dirichlet.guide.ll,
                                                    guide.efficiency)
    temp.dirichlet.ll <- sum(temp.sgRNa.ll$total_guide_ll)
    temp.total.model.ll <- temp.total.model.ll + temp.dirichlet.ll
  }
  
  total.model.ll <- c(total.model.ll, temp.total.model.ll)
  
  # compute the conditional contribution of each FS to the model (but not FS0)
  conditional.fs.ll <- c()
  cond.pp <- matrix(0, ncol = ncol(combined.pp.mtx), nrow = nrow(combined.pp.mtx))
  for(fs in 1:nrow(combined.pp.mtx)){
    
    cond.pp[fs,] <- combined.pp.mtx[fs,]
    temp.pp.sum <- colSums(cond.pp)
    temp.pp.sum[temp.pp.sum > 1] <- 1
    condit.guide.ll <- compute_perGuide_fs_ll(temp.pp.sum, input.data.list$guide_to_seg_lst)
    
    temp.cond.fs.ll <- 0
    for(i in 1:length(input.data.list$data)){
      temp.hypers <- list(alpha0 = dirichlet.hyper$alpha0[[i]], alpha1 = dirichlet.hyper$alpha1[[i]])
      temp.cond.sgRNa.ll <- estimate_relics_sgrna_log_like(temp.hypers,
                                                           input.data.list$data[[i]],
                                                           condit.guide.ll,
                                                           guide.efficiency)
      temp.cond.fs.ll <- temp.cond.fs.ll +  sum(temp.cond.sgRNa.ll$total_guide_ll)
    }
    
    conditional.fs.ll <- c(conditional.fs.ll, temp.cond.fs.ll)
  }
  
  conditional.fs.ll.diff <- diff(conditional.fs.ll)
  
  out.deltas <- input.deltas
  non.signif.idx <- which(2 * conditional.fs.ll.diff < fs.ll.signif)
  if(length(non.signif.idx) > 0){
    out.deltas[non.signif.idx + 1,] <- 0
  }
  
  out.list <- list(total_model_ll_w_prior = total.model.ll,
                   conditional_fs_ll_w_prior = conditional.fs.ll.diff,
                   posteriors_w_prior = combined.pp.mtx,
                   deltas_w_prior = out.deltas)
  
  return(out.list)
  
}


#' @title compute the model lls fit
#' @param input.pp, matrix of posterior probabilities
#' @param dirichlet.hyper: hyper parameters
#' @param input.data.list, list: $guide_to_seg_lst, $data, $true_pos_seg, $overlapMat, $guide_to_seg_lst, $seg_to_guide_lst
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @param input.model.ll: vector of ll for a given model
#' @return list: total_model_ll, per_fs_ll (the model.ll increase for each FS, taking into account all other FS), model_ll_increase (the step-by-step increase in model ll with each FS)
#' @export compute_model_fit()

compute_model_fit <- function(input.pp, input.data.list, guide.efficiency, 
                              dirichlet.hyper, input.model.ll){
  
  total.model.ll <- input.model.ll
  
  # take the sum across all weighted FS placements
  # compute the model ll
  combined.pp <- colSums(input.pp)
  combined.pp[combined.pp > 1] <- 1
  dirichlet.guide.ll <- compute_perGuide_fs_ll(combined.pp, input.data.list$guide_to_seg_lst)
  
  temp.total.model.ll <- 0
  
  for(i in 1:length(input.data.list$data)){
    temp.hypers <- list(alpha0 = dirichlet.hyper$alpha0[[i]], alpha1 = dirichlet.hyper$alpha1[[i]])
    temp.sgRNa.ll <- estimate_relics_sgrna_log_like(temp.hypers,
                                                    input.data.list$data[[i]],
                                                    dirichlet.guide.ll,
                                                    guide.efficiency)
    temp.dirichlet.ll <- sum(temp.sgRNa.ll$total_guide_ll)
    temp.total.model.ll <- temp.total.model.ll + temp.dirichlet.ll
  }
  
  total.model.ll <- c(total.model.ll, temp.total.model.ll)
  
  # compute the conditional contribution of each FS to the model (but not FS0)
  conditional.fs.ll <- c()
  cond.pp <- matrix(0, ncol = ncol(input.pp), nrow = nrow(input.pp))
  for(fs in 1:nrow(input.pp)){
    
    cond.pp[fs,] <- input.pp[fs,]
    temp.pp.sum <- colSums(cond.pp)
    temp.pp.sum[temp.pp.sum > 1] <- 1
    condit.guide.ll <- compute_perGuide_fs_ll(temp.pp.sum, input.data.list$guide_to_seg_lst)
    
    temp.cond.fs.ll <- 0
    for(i in 1:length(input.data.list$data)){
      temp.hypers <- list(alpha0 = dirichlet.hyper$alpha0[[i]], alpha1 = dirichlet.hyper$alpha1[[i]])
      temp.cond.sgRNa.ll <- estimate_relics_sgrna_log_like(temp.hypers,
                                                           input.data.list$data[[i]],
                                                           condit.guide.ll,
                                                           guide.efficiency)
      temp.cond.fs.ll <- temp.cond.fs.ll +  sum(temp.cond.sgRNa.ll$total_guide_ll)
    }
    
    conditional.fs.ll <- c(conditional.fs.ll, temp.cond.fs.ll)
  }
  
  conditional.fs.ll.diff <- diff(conditional.fs.ll)
  
  out.list <- list(total_model_ll = total.model.ll,
                   conditional_fs_ll = conditional.fs.ll.diff)
  
  return(out.list)
  
}



#' @title Order the FS according to their model ll contribution
#' @param input.fs.list: list, most recent FS calculation
#' @param input.model.lls: list, conatins previous model lls
#' @return index of optimal CS
#' @export order_FS()

order_FS <- function(input.fs.list, input.model.lls){
  
  fs.order.idx <- order(input.model.lls$conditional_fs_ll[1:length(input.model.lls$conditional_fs_ll)], decreasing = T)
  
  ordered.pp <- input.fs.list$posteriors[c(1,fs.order.idx+1),]
  ordered.cs.pp <- input.fs.list$cs_posteriors[c(1,fs.order.idx+1),]
  ordered.delta <- input.fs.list$deltas[c(1,fs.order.idx+1),]
  ordered.lls <- input.fs.list$segment_lls[fs.order.idx,, drop = F]
  ordered.cs.pp.total <- input.fs.list$cs_total_pp[c(1,fs.order.idx+1)]
  
  ordered.model.ll <- input.model.lls$total_model_ll # don't order this
  ordered.model.cond.ll <- input.model.lls$conditional_fs_ll[fs.order.idx]
  
  out.list <- list(posteriors = ordered.pp,
                   deltas = ordered.delta,
                   cs_posteriors = ordered.cs.pp,
                   total_model_ll = ordered.model.ll,
                   conditional_fs_ll = ordered.model.cond.ll,
                   segment_lls = ordered.lls,
                   cs_total_pp = ordered.cs.pp.total)
  
  return(out.list)
  
}


#' @title use the PP to compute the credible set and place the FS
#' @param input.posteriors: posterior probabilies of each FS
#' @param hyper: hyperparameters
#' @param data: list, each element is a replicate
#' @param guide.to.seg.lst: mapping of guides to segments
#' @param seg.to.guide.lst: mapping of segments to guides, both $guide_idx, $nonGuide_idx
#' @param next.guide.lst: list that contains the index of $next_guide_idx and next_nonGuide_idx
#' @param nr.segs: max number of segemnts that can be combined
#' @param geom.p: probability of the geometic distribution used to calculate the probability of there being x segments in the regulatory element
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @param area.of.effect: type of area of effect. Currently either 'normal' or 'uniform'
#' @param guide.dist.to.seg, list of lists, pre-computed distances of guides to FS
#' @param fix.hypers, logical, whether the hyperparameters are supposed to be fixed
#' @param fs.result.lists, list, requires $fs_placement_ll, $fs_total_ll
#' @param input.delta: delta matrix, containing the placement of FS
#' @param cs.params: list of parameters for computing the CS
#' @param input.cs.posteriors: matrix of cs posterior probabilities
#' @return log likelihood for each region
#' @export estimate_FS_CS()

estimate_FS_CS <- function(input.posteriors, hyper, data,
                           guide.to.seg.lst, seg.to.guide.lst,
                           next.guide.lst, nr.segs = 10,
                           geom.p = 0.1, guide.efficiency,
                           area.of.effect,
                           guide.dist.to.seg, fix.hypers, fs.result.lists,
                           input.delta, cs.params, input.cs.posteriors) {
  
  posteriors <- input.posteriors
  delta.mtx <- input.delta
  cs.posteriors <- input.cs.posteriors
  
  segment.lls <- fs.result.lists$segment_lls
  if(is.null(segment.lls)){
    segment.lls <- matrix(0, ncol = ncol(delta.mtx))
  } else {
    segment.lls <- rbind(segment.lls, rep(0, ncol(delta.mtx)))
  }
  
  n.sgrna <- length(guide.to.seg.lst)
  n.region <- length(seg.to.guide.lst)
  cs.total.pp <- c(fs.result.lists$cs_total_pp, 0)
  
  k.idx <- nrow(posteriors)
  if(! fix.hypers){
    k.idx <- 2
    segment.lls <- matrix(0, ncol = ncol(delta.mtx), nrow = nrow(posteriors)-1)
  }
  
  new.delta.mtx <- delta.mtx
  new.cs.posteriors <- cs.posteriors
  for(fs in k.idx:nrow(new.delta.mtx)){
    new.delta.mtx[fs, ] <- 0
    new.cs.posteriors[fs, ] <- 0
    fs.config <- colSums(new.cs.posteriors)
    
    # Compute the ll for each RE length for both replicates, combine the lls and then add to the ll total
    # set up the data structues as lists to access the different replicates
    layer.guide.ll <- compute_perGuide_fs_ll(fs.config, guide.to.seg.lst) # processed.pp
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
    
    fs.pps <- c()
    if(area.of.effect == 'normal'){
      fs.pps.list <- estimate_fs_AoE_pp_norm_RelMax(data.mat.list,
                                                    seg.to.guide.lst, next.guide.lst, nr.segs, geom.p, 
                                                    sgrna.log.like.list, seg.adjust = 0, 
                                                    total.segs = length(seg.to.guide.lst),
                                                    guide.dist.to.seg, fs.config)
      fs.pps <- fs.pps.list$seg_pp
      fs.lls <- fs.pps.list$segment_lls
    } else {
      fs.pps.list <- estimate_fs_pp_norm_RelMax(seg.to.guide.lst, next.guide.lst, 
                                                nr.segs, geom.p, sgrna.log.like.list,
                                                fs.config)
      fs.pps <- fs.pps.list$seg_pp
      fs.lls <- fs.pps.list$segment_lls
    }

    posteriors[fs, ] <- fs.pps
    segment.lls[fs-1,] <- fs.lls
    
    # compute CS
    cs.list <- compute_CS(fs.pps, cs.params, fs.config)
    fs.delta.idx <- cs.list$cs_idx
    new.delta.mtx[fs, fs.delta.idx] <- 1
    
    new.cs.posteriors[fs,fs.delta.idx] <- cs.list$cs_pp
    cs.total.pp[fs] <- cs.list$cs_total_pp
    
  }
  
  out.list <- list(posteriors = posteriors,
                   deltas = new.delta.mtx,
                   segment_lls = segment.lls,
                   cs_posteriors = new.cs.posteriors,
                   cs_total_pp = cs.total.pp)
  
  return(out.list)
  
}


#' @title Compute the credible set from the PP to update the delta matrix. Use a sliding window to determine optimal CS idx
#' @param input.pp: vector of pp
#' @param cs.params: list, $cs_sw_size, $cs_threshold, $cs_sw_continous, $min_cs_sum
#' @param fs.config: vector of posteriors, wherever it's not 0 means an FS has already been placed
#' @return index of optimal CS
#' @export compute_CS()

compute_CS <- function(input.pp, cs.params,fs.config){
  
  window.size <- cs.params$cs_sw_size
  sw.idx.list <- list()
  max.sw.pp <- 0
  max.sw.idx <- 1
  max.idx.total.pp <- 0
  sw.pp <- vector('numeric', length = length(input.pp) - window.size + 1)
  median.pp <- median(input.pp)
  
  for(i in 1:(length(input.pp) - window.size + 1)){
    
    temp.sw.idx <- c(1:window.size) + (i-1)
    if(sum(fs.config[temp.sw.idx]) > 0){
      temp.sw.idx <- temp.sw.idx[-which(fs.config[temp.sw.idx] > 0)]
    }
    
    if(length(temp.sw.idx) > 0){
      temp.pp <- input.pp[temp.sw.idx]
      temp.pp.sum <- sum(temp.pp)
      
      if(temp.pp.sum > 0){
        temp.sw.idx.ordered <- temp.sw.idx[order(temp.pp, decreasing = TRUE)]
        temp.pp.ordered <- temp.pp[order(temp.pp, decreasing = TRUE)]
        
        # compute the relative proportion of the CS
        temp.sw.idx.pass <- c()
        for(j in 1:length(temp.sw.idx.ordered)){
          temp.sw.pp.sum <- sum(temp.pp.ordered[1:j])
          temp.fract.pp <- temp.sw.pp.sum / temp.pp.sum
          temp.sw.idx.pass <- c(temp.sw.idx.pass, temp.sw.idx.ordered[j])
          if(temp.fract.pp > cs.params$cs_threshold){
            if(temp.sw.pp.sum > max.sw.pp){
              max.sw.pp <- temp.sw.pp.sum
              max.sw.idx <- temp.sw.idx.pass
              max.idx.total.pp <- temp.pp.sum
            }
            sw.idx.list[[i]] <- temp.sw.idx.pass
            sw.pp[i] <- temp.sw.pp.sum
            break()
          }
        }
      }
    }
  }
  
  # the the FS output is to be continuous, set the idx to the range from min to max idx
  if(cs.params$cs_sw_continous){
    min.idx <- min(max.sw.idx)
    max.idx <- max(max.sw.idx)
    idx.range <- min.idx:max.idx
    if(sum(fs.config[idx.range]) > 0){
      idx.range <- idx.range[-which(fs.config[idx.range] > 0)]
    }
  } else {
    idx.range <- max.sw.idx
  }
  max.sw.pp.idx <- idx.range
  max.sw.cs.pp <- input.pp[max.sw.pp.idx] / max.idx.total.pp
  
  if(sum(input.pp[max.sw.pp.idx]) < cs.params$min_cs_sum){
    out.list <- list(cs_idx = NULL,
                     rel_cs_pp = 0,
                     cs_total_pp = 0)
    return(out.list)
  } else {
    out.list <- list(cs_idx = max.sw.pp.idx,
                     cs_pp = max.sw.cs.pp,
                     cs_total_pp = max.sw.pp)
    return(out.list)
  }
  
}


#' @title Compute the log likelihood of each possible configuration of the placement of a fs. Take into account distance between segment and the guide target site
#' @param in.data.list: list, each element contains a replicate, the total column has already been removed
#' @param seg.to.guide.lst: list, each element is a segment, contains a list with guide_idx, nonGuide_idx
#' @param next.guide.lst: list that contains the index of $next_guide_idx and next_nonGuide_idx
#' @param nr.segs: max number of segemnts that can be combined
#' @param geom.p: probability of the geometic distribution used to calculate the probability of there being x segments in the regulatory element
#' @param sg.ll.list: list, each element is a data frame per.guide log likelihood given the posteriors for a replicate and alt. ll (total_guide_ll, alt_only_ll)
#' @param seg.adjust: numeric, offset used when calculating lmax
#' @param total.segs: how many segments to calculate PP for, adjusted for lmax, otherwise = length(seg.to.guide.lst)
#' @param guide.dist.to.seg, list of lists, pre-computed distances of guides to FS
#' @param fs.config: previous placements of FS to ensure that no repeated placement happens
#' @return list: segment_lls, total_ll, seg_pp
#' @export estimate_fs_AoE_pp_norm_RelMax()

estimate_fs_AoE_pp_norm_RelMax <- function(in.data.list, seg.to.guide.lst, next.guide.lst, nr.segs = 10, geom.p = 0.1,
                                       sg.ll.list, seg.adjust = 0, total.segs, guide.dist.to.seg, fs.config) {
  
  # set up a list, the length of the number of segments
  # each element contains a vector, which contains the ll of that segment being overlapped by a regulatory element
  segment.ll.list <- list()
  
  region.priors <- rep(1/total.segs, total.segs)
  geom.norm.factr <- pgeom(nr.segs, geom.p)
  
  total.ll <- log(0)
  fs.length <- 1
  
  null.model.ll.list <- list()
  for(repl in 1:length(sg.ll.list)){
    repl.guide.null <- sum(sg.ll.list[[repl]][seg.to.guide.lst[[1]]$guide_idx,1])
    repl.non.guide.null <- sum(sg.ll.list[[repl]][seg.to.guide.lst[[1]]$nonGuide_idx,1])
    null.model.ll.list[[repl]] <- repl.guide.null + repl.non.guide.null
  }
  
  # make an initial pass with one-segment length RE placements
  for(seg in 1:total.segs){
    curr.seg <- seg + seg.adjust
    temp.segment.ll <- 0
    temp.guide.idx<- seg.to.guide.lst[[curr.seg]]$guide_idx
    temp.nonGuide.idx<- seg.to.guide.lst[[curr.seg]]$nonGuide_idx
    
    temp.seg.poi.bkg <- guide.dist.to.seg[[curr.seg]][[fs.length]]
    
    for(repl in 1:length(sg.ll.list)){
      
      temp.guide.ll <- vector('numeric', length = length(temp.seg.poi.bkg))
      
      for(i in 1:length(temp.seg.poi.bkg)){
        temp.guide.ll[i] <- addlogs(sg.ll.list[[repl]][temp.guide.idx[i],2] + log(1 - temp.seg.poi.bkg[i]),
                                    sg.ll.list[[repl]][temp.guide.idx[i],1] + log(temp.seg.poi.bkg[i]))
      }
      
      temp.nonGuide.ll <- sum(sg.ll.list[[repl]][temp.nonGuide.idx,1])
      
      if(fs.config[seg] > 0){
        temp.segment.ll <- temp.segment.ll + null.model.ll.list[[repl]] + log(dgeom(1, geom.p) / geom.norm.factr)
      } else if(sum(temp.guide.ll) + temp.nonGuide.ll < null.model.ll.list[[repl]]){
        temp.segment.ll <- temp.segment.ll + null.model.ll.list[[repl]] + log(dgeom(1, geom.p) / geom.norm.factr)
      } else {
        temp.segment.ll <- temp.segment.ll + sum(temp.guide.ll) + temp.nonGuide.ll + log(dgeom(1, geom.p) / geom.norm.factr)
      }
    }
    segment.ll.list[[seg]] <- temp.segment.ll
    total.ll <- addlogs(total.ll, temp.segment.ll)
  }
  
  all.guide.idx <- c(1:nrow(in.data.list[[1]])) # helps to differentitate between guides for null vs alternative
  
  # now make a pass through all 2 to nr.segs length segments
  if(nr.segs > 1){
    for(ns in 1:(nr.segs-1) ){
      fs.length <- fs.length + 1
      for(seg in 1:(total.segs - ns)){
        temp.stretch.ll <- 0
        temp.stretch.segs <- c(seg:(seg + ns)) # index of segments to be included as part of the stretch and be updated in this iteration
        curr.seg <- seg + seg.adjust
        curr.stretch.segs <- temp.stretch.segs + seg.adjust
        temp.guide.idx <- c(seg.to.guide.lst[[curr.seg]]$guide_idx, unlist(next.guide.lst[(curr.seg + 1):(curr.seg + ns)]))
        temp.nonGuide.idx <- all.guide.idx[-temp.guide.idx]
        
        temp.seg.poi.bkg <- guide.dist.to.seg[[curr.seg]][[fs.length]]
        
        for(repl in 1:length(sg.ll.list)){
          
          temp.guide.ll <- vector('numeric', length = length(temp.seg.poi.bkg))
          
          for(i in 1:length(temp.seg.poi.bkg)){
            temp.guide.ll[i] <- addlogs(sg.ll.list[[repl]][temp.guide.idx[i],2] + log(1 - temp.seg.poi.bkg[i]),
                                        sg.ll.list[[repl]][temp.guide.idx[i],1] + log(temp.seg.poi.bkg[i]))
          }
          
          temp.nonGuide.ll <- sum(sg.ll.list[[repl]][temp.nonGuide.idx,1])
          
          if(sum(fs.config[temp.stretch.segs]) > 0 ){
            temp.stretch.ll <- temp.stretch.ll + null.model.ll.list[[repl]] + log(dgeom(1 + ns, geom.p) / geom.norm.factr)
          } else if(sum(temp.guide.ll) + temp.nonGuide.ll < null.model.ll.list[[repl]]){
            temp.stretch.ll <- temp.stretch.ll + null.model.ll.list[[repl]] + log(dgeom(1 + ns, geom.p) / geom.norm.factr)
          } else {
            temp.stretch.ll <- temp.stretch.ll + sum(temp.guide.ll) + temp.nonGuide.ll + log(dgeom(1 + ns, geom.p) / geom.norm.factr)
          }
        }
        
        total.ll <- addlogs(total.ll, temp.stretch.ll)
        
        for(span in 1:length(temp.stretch.segs)){
          segment.ll.list[[temp.stretch.segs[span]]] <- c(segment.ll.list[[temp.stretch.segs[span]]], temp.stretch.ll)
        }
      }
    }
  }
  
  final.total.ll <- addlogs_vectorized(unlist(segment.ll.list))
  
  # normalize the list by the total, and add the logs across the segments
  segment.lls <- unlist(lapply(segment.ll.list, function(x){
    x.norm <- x - final.total.ll # total.ll
    
    x.norm.sum <- combine_segment_ll(x.norm)
    x.norm.sum
  }))
  segment.posteriors <- exp(segment.lls)
  
  raw.segment.lls <- unlist(lapply(segment.ll.list, function(x){
    x.norm.sum <- combine_segment_ll(x)
    x.norm.sum
  }))
  
  out.list <- list(segment_lls = segment.lls,
                   total_ll = final.total.ll,
                   seg_pp = segment.posteriors,
                   raw_segment_lls = raw.segment.lls)
  
  return(out.list)
  
}


#' @title Compute the log likelihood of each possible configuration of the placement of a fs amonst a set of segments
#' @param seg.to.guide.lst: list, each element is a segment, contains a list with guide_idx, nonGuide_idx
#' @param next.guide.lst: list that contains the index of $next_guide_idx and next_nonGuide_idx
#' @param nr.segs: max number of segemnts that can be combined
#' @param geom.p: probability of the geometic distribution used to calculate the probability of there being x segments in the regulatory element
#' @param sg.ll.list: list, each element is a data frame per.guide log likelihood given the posteriors for a replicate and alt. ll (total_guide_ll, alt_only_ll)
#' @param fs.config: previous placements of FS to ensure that no repeated placement happens
#' @return log likelihood for each segment
#' @export estimate_fs_pp_norm_RelMax()

estimate_fs_pp_norm_RelMax <- function(seg.to.guide.lst, next.guide.lst, nr.segs = 10, geom.p = 0.1, sg.ll.list, fs.config) {
  
  # set up a list, the length of the number of segments
  # each element contains a vector, which contains the ll of that segment being overlapped by a regulatory element
  segment.ll.list <- list()
  
  total.segs <- length(seg.to.guide.lst)
  region.priors <- rep(1/total.segs, total.segs)
  geom.norm.factr <- pgeom(nr.segs, geom.p)
  
  total.ll <- log(0)
  null.model.ll.list <- list()
  for(repl in 1:length(sg.ll.list)){
    repl.guide.null <- sum(sg.ll.list[[repl]][seg.to.guide.lst[[1]]$guide_idx,1])
    repl.non.guide.null <- sum(sg.ll.list[[repl]][seg.to.guide.lst[[1]]$nonGuide_idx,1])
    null.model.ll.list[[repl]] <- repl.guide.null + repl.non.guide.null
  }
  
  # make an initial pass with one-segment length RE placements
  for(seg in 1:total.segs){
    temp.segment.ll <- 0
    temp.guide.idx<- seg.to.guide.lst[[seg]]$guide_idx
    temp.nonGuide.idx<- seg.to.guide.lst[[seg]]$nonGuide_idx
    
    for(repl in 1:length(sg.ll.list)){ 
      temp.guide.ll <- sum(sg.ll.list[[repl]][temp.guide.idx,2])
      temp.nonGuide.ll <- sum(sg.ll.list[[repl]][temp.nonGuide.idx,1])
      if(fs.config[seg] > 0){
        temp.segment.ll <- temp.segment.ll + null.model.ll.list[[repl]] + log(dgeom(1, geom.p) / geom.norm.factr)
      } else if(temp.guide.ll + temp.nonGuide.ll < null.model.ll.list[[repl]]){
        temp.segment.ll <- temp.segment.ll + null.model.ll.list[[repl]] + log(dgeom(1, geom.p) / geom.norm.factr)
      } else {
        temp.segment.ll <- temp.segment.ll + temp.guide.ll + temp.nonGuide.ll + log(dgeom(1, geom.p) / geom.norm.factr)
      }
    }
    segment.ll.list[[seg]] <- temp.segment.ll
    total.ll <- addlogs(total.ll, temp.segment.ll)

  }
  
  all.guide.idx <- c(1:length(sg.ll.list[[1]][,1])) # helps to differentiate between guides for null vs alternative
  
  # now make a pass through all 2 to nr.segs length segments
  if(nr.segs > 1){
    for(ns in 1:(nr.segs-1) ){
      
      for(seg in 1:(total.segs - ns)){
        
        temp.stretch.ll <- 0
        temp.stretch.segs <- c(seg:(seg + ns)) # index of segments to be included as part of the stretch and be updated in this iteration
        temp.guide.idx <- c(seg.to.guide.lst[[seg]]$guide_idx, unlist(next.guide.lst[(seg + 1):(seg + ns)]))
        temp.nonGuide.idx <- all.guide.idx[-temp.guide.idx]
        
        for(repl in 1:length(sg.ll.list)){
          temp.guide.ll <- sum(sg.ll.list[[repl]][temp.guide.idx,2])
          temp.nonGuide.ll <- sum(sg.ll.list[[repl]][temp.nonGuide.idx,1])
          
          if(sum(fs.config[temp.stretch.segs]) > 0){
            temp.stretch.ll <- temp.stretch.ll + null.model.ll.list[[repl]] + log(dgeom(1 + ns, geom.p) / geom.norm.factr)
          } else if(temp.guide.ll + temp.nonGuide.ll < null.model.ll.list[[repl]]){
            temp.stretch.ll <- temp.stretch.ll + null.model.ll.list[[repl]] + log(dgeom(1 + ns, geom.p) / geom.norm.factr)
          } else {
            temp.stretch.ll <- temp.stretch.ll + temp.guide.ll + temp.nonGuide.ll + log(dgeom(1 + ns, geom.p) / geom.norm.factr)
          }
        }
        
        total.ll <- addlogs(total.ll, temp.stretch.ll)
        
        for(span in 1:length(temp.stretch.segs)){
          segment.ll.list[[temp.stretch.segs[span]]] <- c(segment.ll.list[[temp.stretch.segs[span]]], temp.stretch.ll)
        }
      }
    }
    
  }
  
  final.total.ll <- addlogs_vectorized(unlist(segment.ll.list)) # total.ll
  
  # normalize the list by the total, and add the logs across the segments
  segment.lls <- unlist(lapply(segment.ll.list, function(x){
    x.norm <- x - final.total.ll #total.ll
    x.norm.sum <- combine_segment_ll(x.norm)
    x.norm.sum
  }))
  segment.posteriors <- exp(segment.lls)
  
  raw.segment.lls <- unlist(lapply(segment.ll.list, function(x){
    x.norm.sum <- combine_segment_ll(x)
    x.norm.sum
  }))
  
  out.list <- list(segment_lls = segment.lls,
                   total_ll = final.total.ll,
                   seg_pp = segment.posteriors,
                   raw_segment_lls = raw.segment.lls)
  
  return(out.list)
  
}


#' @title wrapper which returns the updated guide efficiency and the corresponding coefficients used to calculate it
#' @param hyper.components: hyperparameters components, 
#' @param param: matrix of all posterior probs
#' @param data: list, each element is a replicate
#' @param guide.seg.idx.lst: list: each element is a guide, contaiing the indexes of the segments it overlaps
#' @param guide.efficiency.scores: matrix, each column containing a different set of scores per guide
#' @param ge.coeff vector, containing the guide efficiency coefficients
#' @return list: $guide_efficiency, $ge_coeff
#' @export recompute_ge_coefficients_w_deltas()

recompute_ge_coefficients_w_deltas <- function(param, hyper.components, data, guide.seg.idx.lst, guide.efficiency.scores, ge.coeff) {
  cumulativedeltas <- colSums(param$deltas) #apply(param$deltadeltas, 2, sum)
  cumulativedeltas[cumulativedeltas > 1] <- 1
  
  guide.lls.list <- compute_perGuide_fs_ll(cumulativedeltas, guide.seg.idx.lst)
  
  ge.coeff.param <- ge.coeff
  disp.type <- c()
  if(is.null(hyper.components$dispersion)){
    disp.type <- hyper.components$bkg_dispersion
    model.disp <- FALSE
  } else {
    disp.type <- hyper.components$dispersion
    model.disp <- TRUE
  }
  
  res <- optim(ge.coeff.param, guide_coeff_ll, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
               data=data, region.ll.list = guide.lls.list,
               bkg.alpha = hyper.components$bkg_alpha, fs.alpha = hyper.components$FS_alpha,
               guide.efficiency.scores = guide.efficiency.scores, input.dispersion = disp.type, model.disp = model.disp)
  
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
    out.list <- list()
    out.list$guide_efficiency <- guide.efficiency.scores
    out.list$ge_coeff <- ge.coeff
  }
  
  return(out.list)
}


#' @title re-estimate the dirichlet sorting proportions (dispersions are already provided)
#' @param param: matrix of all posterior probs
#' @param guide.seg.idx.lst: list: each element is a guide, contaiing the indexes of the segments it overlaps
#' @param hyper.components: individual components of the hyper parameters
#' @param data.par.list: list, contains the elements of the processed and formatted data
#' @param guide.efficiency: either a vector of guide efficiency or NULL
#' @param hyper: hyperparameters, already divided into alpha0 and alpha1
#' @param analysis.parameters: list of analysis paramters
#' @return list of the re-estimated dirichlet hyperparameters / sorting proportions
#' @export recompute_dirichlet_hypers_w_deltas()

recompute_dirichlet_hypers_w_deltas <- function(param, guide.seg.idx.lst, hyper.components, data.par.list, guide.efficiency, hyper, analysis.parameters){
  
  cumulative.deltas <- colSums(param$deltas)
  cumulative.deltas[cumulative.deltas > 1] <- 1
  
  guide.lls.list <- compute_perGuide_fs_ll(cumulative.deltas, guide.seg.idx.lst)
  
  for(i in 1:length(data.par.list)){
    
    temp.hyper.params <- c(sqrt(hyper.components$bkg_alpha[[i]]), sqrt(hyper.components$FS_alpha[[i]]))
    temp.bkg.idx <- c(1:length(hyper.components$bkg_alpha[[i]]))
    temp.fs.idx <- c(1:length(hyper.components$FS_alpha[[i]])) + max(temp.bkg.idx)
    temp.repl.disp <- hyper.components$dispersion[[i]]
    
    temp.res.drch <- c()
    temp.res.drch <- optim(temp.hyper.params, prior_dirichlet_proportions, method= 'L-BFGS-B',
                           data = data.par.list[[i]], 
                           region.ll.list = guide.lls.list,
                           bkg.idx = temp.bkg.idx, 
                           fs.idx = temp.fs.idx, 
                           guide.efficiency = guide.efficiency, 
                           repl.disp = temp.repl.disp,
                           model.disp = analysis.parameters$model_dispersion)
    
    if(temp.res.drch$convergence %in% c(0, 52) ) {
      
      temp.bkg.alpha <- c(1, temp.res.drch$par[temp.bkg.idx]**2) / sum(c(1, temp.res.drch$par[temp.bkg.idx]**2))
      temp.fs.alpha <- c(1, temp.res.drch$par[temp.fs.idx]**2) / sum(c(1, temp.res.drch$par[temp.fs.idx]**2))
      
      hyper.components$bkg_alpha[[i]] <- temp.res.drch$par[temp.bkg.idx]**2
      hyper.components$FS_alpha[[i]] <- temp.res.drch$par[temp.fs.idx]**2
      
      temp.hyper <- reparameterize_hypers(temp.bkg.alpha, temp.fs.alpha, temp.repl.disp, guide.efficiency, analysis.parameters$model_dispersion)
      hyper$alpha0[[i]] <- temp.hyper$alpha0
      hyper$alpha1[[i]] <- temp.hyper$alpha1
    } else {
      warning("estimation of hyperparameters failed to converge")
    }
    
  }
  
  out.list <- list(hyper_pars = hyper,
                   hyper_par_components = hyper.components)
  
  return(out.list)
}


#' @title re-estimate the dirichlet sorting proportions if the dispersion is not modeled
#' @param hyper: hyperparameters, already divided into alpha0 and alpha1
#' @param hyper.components: individual components of the hyper parameters
#' @param param: matrix of all posterior probs
#' @param data: data, consists of: $y1, $y2 and $n
#' @param guide.seg.idx.lst: list: each element is a guide, containing the indexes of the segments it overlaps
#' @param guide.efficiency: either a vector of guide efficiency or NULL
#' @param analysis.parameters: list of analysis paramters
#' @return list of the re-estimated hyperparameters
#' @export recompute_hyper_parameters_w_deltas()

recompute_hyper_parameters_w_deltas <- function(param, hyper, hyper.components, data, guide.seg.idx.lst, guide.efficiency, analysis.parameters) {
  cumulative.deltas <- colSums(param$deltas) #adeltasly(param$delta.deltas, 2, sum)
  cumulative.deltas[cumulative.deltas > 1] <- 1
  
  guide.lls.list <- compute_perGuide_fs_ll(cumulative.deltas, guide.seg.idx.lst)
  
  for(i in 1:length(data)){
    
    temp.hyper.params <- c(sqrt(hyper.components$bkg_alpha[[i]]), sqrt(hyper.components$FS_alpha[[i]]), hyper.components$bkg_dispersion[[i]])
    temp.bkg.idx <- c(1:length(hyper.components$bkg_alpha[[i]]))
    temp.fs.idx <- c(1:length(hyper.components$FS_alpha[[i]])) + max(temp.bkg.idx)
    temp.disp.idx <- c(1:length(hyper.components$bkg_dispersion[[i]])) + max(temp.fs.idx)
    
    temp.res.drch <- optim(temp.hyper.params, prior_dirichlet_parameters, method= 'L-BFGS-B',
                           data = data[[i]], 
                           region.ll.list = guide.lls.list,
                           bkg.idx = temp.bkg.idx, 
                           fs.idx = temp.fs.idx, 
                           disp.idx = temp.disp.idx, 
                           guide.efficiency = guide.efficiency, 
                           mean.var.type = analysis.parameters$mean_var_type,
                           model.disp = analysis.parameters$model_dispersion)
    
    # account for flat surface!
    if(temp.res.drch$convergence %in% c(0, 52) ) {
      # return new estimates of hyperparamers
      
      temp.bkg.alpha <- c(1, temp.res.drch$par[temp.bkg.idx]**2) / sum(c(1, temp.res.drch$par[temp.bkg.idx]**2))
      temp.fs.alpha <- c(1, temp.res.drch$par[temp.fs.idx]**2) / sum(c(1, temp.res.drch$par[temp.fs.idx]**2))
      temp.bkg.disp <- temp.res.drch$par[temp.disp.idx]
      
      hyper.components$bkg_alpha[[i]] <- temp.res.drch$par[temp.bkg.idx]**2
      hyper.components$FS_alpha[[i]] <- temp.res.drch$par[temp.fs.idx]**2
      hyper.components$bkg_dispersion[[i]] <- temp.bkg.disp
      
      if(analysis.parameters$mean_var_type == 'independent'){
        temp.disp <- rep(temp.bkg.disp[1]^2, nrow(data[[i]]))
        temp.hyper <- reparameterize_hypers(temp.bkg.alpha, temp.fs.alpha, temp.disp, guide.efficiency, TRUE)
        hyper$alpha0[[i]] <- temp.hyper$alpha0
        hyper$alpha1[[i]] <- temp.hyper$alpha1
      }
      
    } else {
      warning("estimation of hyperparameters failed to converge")
    }
  }
  
  out.list <- list(hyper_pars = hyper,
                   hyper_par_components = hyper.components)
  
  # return(hyper)
  return(out.list)
}


#' @title Display the credible sets along with cumulative posteriors. Formerly 'display_mrvr_norm_tiff'
#' @param input.L matrix of posteriors
#' @param input.labels: labelling of all segments
#' @param tiff.name: name for plot
#' @param fs.threshold: threshold for defining functional sequences
#' @param seg.info: contains info for each segment, specicially chromosome information
#' @return One plot for each selection step, credible sets colored in purple
#' @export display_relics_segLLs_as_tiff()

display_relics_segLLs_as_tiff <- function(input.seg.lls, input.labels, tiff.name, fs.threshold, seg.info){
  
  cs.col.orig <- rep("darkgrey", ncol(input.seg.lls))
  
  if(length(unique(seg.info$chrom)) > 1){
    even.chroms <- unique(seg.info$chrom)[c(FALSE,TRUE)]
    even.chroms.idx <- which(seg.info$chrom %in% even.chroms)
    cs.col.orig[even.chroms.idx] <- 'lightgrey'
  }
  
  # attempt to put 8 per page
  if(nrow(input.seg.lls) <= 9){
    tiff(paste0(tiff.name, '.tiff'), width = 11520 / 2, height = 5760 / 2, units = "px", res = 400)
    par(mfrow = c(3,3))
    
    for(i in 1:nrow(input.seg.lls)){
      cs.col <- cs.col.orig
      temp.lls <- input.seg.lls[i,]
      max.idx <- which(temp.lls == max(temp.lls)[1])

      plot(input.seg.lls[i,], pch=21,
           main=paste("FS", i, ", Nr. top segments: ", length(max.idx), sep=""), col=cs.col,
           ylab = 'Segment ll', xlab = 'Genome Segment')
      abline(v=max.idx, col="red")
    }
    dev.off()
  } else {
    nr.pages <- ceiling(nrow(input.seg.lls) / 9)
    
    for(pg in 1:nr.pages){
      tiff(paste0(tiff.name, '_', pg, '.tiff'), width = 11520 / 2, height = 5760 / 2, units = "px", res = 400)
      par(mfrow = c(3,3))
      
      for(i in 1:9){
        
        row.index <- (pg - 1)*9 + i
        
        if(row.index > nrow(input.seg.lls)){
          break()
        }
        
        cs.col <- cs.col.orig
        temp.lls <- input.seg.lls[row.index,]
        max.idx <- which(temp.lls == max(temp.lls)[1])
        
        plot(input.seg.lls[row.index,], pch=21,
             main=paste("FS", row.index, ", Nr. top segments: ", length(max.idx), sep=""), col=cs.col,
             ylab = 'Segment ll', xlab = 'Genome Segment')
        abline(v=max.idx, col="red")
      }
      dev.off()
    }
    
  }
  
}


#' @title obtain the splines for each replicate to generate the dispersions
#' @param repl.spline.df: number of dgrees of freedom used for the spline, per replicate
#' @param input.data: list: each element is a replicate
#' @param out.dir: directory to which the file is written
#' @param nr.bins: default 20, number of bins to divide the data into for estimating the hyper parameters
#' @param fs0.label: label of the training guides
#' @param guide.info: data frame containing info about each guide
#' @return sum of the -log likelihood across all guides
#' @export disp_from_spline()

disp_from_spline <- function(repl.spline.df, input.data, out.dir, nr.bins = 20, fs0.label, guide.info) {
  
  repl.guide.disp <- list()
  
  # identify guides used for training and not part of background calculation, remove them
  fs0.label.idx <- which(guide.info$label %in% fs0.label)
  
  for(i in 1:length(input.data)){
    # counts from a replicate. The last column has to be the total of all counts
    temp.all.counts <- input.data[[i]]

    # remove the FS0 labeled data.
    temp.counts <- temp.all.counts[-fs0.label.idx,]
    temp.total.counts <- temp.counts[,ncol(temp.counts)]
    
    # estimate the sorting proportions
    # but only for the background. 
    # Due to identifiability we only need to specify n-1 sorting proportions, where n is the number of pools
    temp.bkg.alpha <- rep(1, ncol(temp.counts) - 2)
    temp.bkg.disp <- c(1)
    temp.hyper.params <- c(temp.bkg.alpha, temp.bkg.disp)
    temp.bkg.idx <- c(1:length(temp.bkg.alpha))
    temp.disp.idx <- length(temp.hyper.params)
    temp.res.drch <- optim(temp.hyper.params, prior_bkg_dirichlet_parameters, method= 'L-BFGS-B',
                           data = temp.counts, 
                           bkg.idx = temp.bkg.idx, 
                           disp.idx = temp.disp.idx)		
    
    # Extract the obtained hyperparameters
    temp.bkg.freq <- c(1, temp.res.drch$par[temp.bkg.idx]**2) / sum(c(1, temp.res.drch$par[temp.bkg.idx]**2))
    temp.bkg.disp.optim <- temp.res.drch$par[temp.disp.idx]^2
    
    # now estimate the dispersion for each of the bins
    # the guides are ordered and then placed into bins
    # the assumption is that with increasing counts in the bins the dispersion will change too
    temp.counts.sort <- temp.counts[order(temp.counts[,ncol(temp.counts)], decreasing = T),]
    temp.bg.disp <- temp.bkg.disp.optim # start the bkg dispersion where the optimizer above left off
    disp.groups <- nr.bins # number of dispersion groups
    guides.per.grp <- round(nrow(temp.counts.sort)/disp.groups)
    guide.grouping <- rep(c(1:disp.groups), each = guides.per.grp)[1:nrow(temp.counts.sort)]
    temp.repl.counts.groups <- split(temp.counts.sort, guide.grouping)
    
    # for each bin, find the optimal dispersion
    temp.repl.group.disp <- unlist(lapply(temp.repl.counts.groups, function(x){
      temp.grp.disp <- sqrt(temp.bg.disp)
      temp.res <- optim(temp.grp.disp, optim_guide_disp,
                        input.hyper = temp.bkg.freq, 
                        input.data = x, method= 'L-BFGS-B')
      temp.res$par^2
    }))
    temp.repl.group.disp[which(temp.repl.group.disp < 0.1)] <- 0.1
    
    # extract the mean counts of each bin
    temp.repl.group.means <- unlist(lapply(temp.repl.counts.groups, function(x){
      mean(x[,ncol(x)])
    }))
    
    # now we can plot the estimated dispersions for each bin
    pdf(paste0(out.dir, '_countsVSdispersion_repl', i, '.pdf'), useDingbats = FALSE)
    plot(x = temp.repl.group.means, y = 1/temp.repl.group.disp, pch = 19,
         xlab = 'Counts', ylab = 'Dispersion', main = 'Dispersion estimates')
    
    # the df parameter stands for degrees of freedom
    # it refers to the number of data points that is spaces around the range of all data points to fit the spline
    # the two corner cases (max and min) are already given. And then it fits another df-1 knots to fitting the data
    data.df <- data.frame(disp = temp.repl.group.disp, counts = temp.repl.group.means, stringsAsFactors = F)
    spline.mdl.1 <- suppressWarnings(lm(disp ~ bs(counts, df = repl.spline.df[[i]]), data = data.df))

    # Now that you have the spline models established, we have to see how well it actually does
    # The predict model uses the models established above to compute what the expected 
    # dispersion is given the range of all possible counts we have in the data.
    # Because the spline models use the same variable names we first have to set up another data frame with the name 'counts' as the column
    vals.df <- data.frame(counts = temp.all.counts[,ncol(temp.counts)], stringsAsFactors = F)
    spline.mdl.1.predict <- suppressWarnings(predict(spline.mdl.1, vals.df))
    
    # With the predictions we can see how well they actually get modeled
    points(x = vals.df$counts[order(vals.df$counts)], y = 1/spline.mdl.1.predict[order(vals.df$counts)], 
           pch = 19, col = 'red', cex = 0.5) 
    lgndTxt<-c('black dot = per-bin estimated dispersion','red dot = per-guide dispersion')
    if(1/temp.repl.group.disp[1] < 
       mean(1/temp.repl.group.disp[(length(temp.repl.group.disp)-min(10, length(temp.repl.group.disp)) ):length(temp.repl.group.disp)])){
      legend('topright',legend=lgndTxt, text.col=c('black','red'), pch = c(19,19), col = c('black','red'))
    } else {
      legend('bottomright',legend=lgndTxt, text.col=c('black','red'), pch = c(19,19), col = c('black','red'))
    }
    dev.off()
    
    # cap the dispersion predictions at the highest mean value
    max.mean.pos <- which(temp.repl.group.means == max(temp.repl.group.means))
    spline.mdl.1.predict.capped <- spline.mdl.1.predict
    spline.mdl.1.predict.capped[which(temp.all.counts[,ncol(temp.counts)] > max(temp.repl.group.means))] <- temp.repl.group.disp[max.mean.pos]
    spline.mdl.1.predict.capped[which(spline.mdl.1.predict.capped < 0.1)] <- 0.1
    
    repl.guide.disp[[i]] <- list(repl_disp = spline.mdl.1.predict.capped,
                                 repl_spline_mdl = spline.mdl.1)
    
  }
 
  return(repl.guide.disp)
}


#' @title optimize the hyper parameters, only one dispersion across the two distributions. But only for the background guides
#' @param hyper.param: hyperparameters
#' @param data: data, consists of: $pool1, pool2... and $n
#' @param bkg.idx: positions in the par vector of the background alphas
#' @param disp.idx: positions in the par vector of the dispersion parameters
#' @return sum of the -log likelihood across all guides
#' @export prior_bkg_dirichlet_parameters()

prior_bkg_dirichlet_parameters <- function(hyper.param, data, bkg.idx, disp.idx) {
  
  bkg.alpha <- c(1, hyper.param[bkg.idx]**2) / sum(c(1, hyper.param[bkg.idx]**2))
  disp.alpha <- hyper.param[disp.idx]
  
  hyper <- list()

  temp.disp <- disp.alpha[1]^2
  hyper$alpha0 <- bkg.alpha * temp.disp

  sgrna.null.log.like <- ddirmnom(data[,1:(ncol(data) - 1)], size = data[,ncol(data)], alpha = hyper$alpha0, log = T)
  
  -sum(sgrna.null.log.like)
}


#' @title optimize the dispersion for a set of dirichlet probabilities and data
#' @param input.hyper: dirichlet probabilities
#' @param input.data: data frame with guide counts
#' @param par: parameter to optimize (dispersion)
#' @return -log likelihood
#' @export optim_guide_disp()

optim_guide_disp <- function(par, input.hyper, input.data){
  
  temp.dir.hyper <- input.hyper * par^2
  
  out.ll <- ddirmnom(input.data[,c(1:(ncol(input.data) - 1))], 
                     size = input.data[,ncol(input.data)], alpha = temp.dir.hyper, log = T)
  
  -sum(out.ll)
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
      
      # adjust in case of l2fc
      if(! is.null(input.parameter.list$l2fc_groups)){
        input.parameter.list$l2fc_groups[[i]] <- repl_pools[[i]][match(input.parameter.list$l2fc_groups[[i]], repl_pools.original[[i]])]
      }
      
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
    sim.counts <- sim.counts[-which(! sim.info$label %in% labelHierarchy),]
    sim.info <- sim.info[-which(! sim.info$label %in% labelHierarchy),]
  }

  if(length(which(is.na(sim.info$start))) > 0){
    print('Removing rows due to NAs in the guide target sites')
    sim.counts <- sim.counts[-which(is.na(sim.info$start)),]
    sim.info <- sim.info[-which(is.na(sim.info$start)),]
    #save.files <- TRUE
  }

  sim.counts <- sim.counts[order(sim.info$chrom, sim.info$start),]
  sim.info <- sim.info[order(sim.info$chrom, sim.info$start),]

  # if GE scores exist, extract them from the info file
  filtered.ges <- NULL
  if(! is.null(input.parameter.list$guide_efficiency_scores)){
    processed.ge.scores <- as.matrix(sim.info[, ges.cols, drop = F])
    filtered.ges <- processed.ge.scores
  }

  # initialize varaible
  sim.guide.seg.list <- c()
  
  print('Formatting data...')
  
  guide.to.seg.lst <- c()

  if(input.parameter.list$crisprSystem %in% c('CRISPRi', 'CRISPRa', 'CRISPRcas9', 'Cas9')){
    # single guide setup
    if(! 'sgTarget' %in% names(sim.info)){
      sim.info$sgTarget <- round( (sim.info$start + sim.info$end) / 2)
    }
    
    # adjust the target positions according to CRISPRi
    sim.info$start <- sim.info$sgTarget - guide.offset
    sim.info$end <- sim.info$sgTarget + guide.offset

    # generate the guide-segment matrix and the per-segment labels
    # seg_info, guide_to_seg_lst, seg_to_guide_lst, counts
    sim.guide.seg.list <- adapt_data_to_regionFormat(sim.counts, repl_pools, sim.info, labelHierarchy, min.seg.dist)
    
    guide.to.seg.lst <- sim.guide.seg.list$guide_to_seg_lst
    
    if(input.parameter.list$areaOfEffect_type == 'uniform'){
      guide.to.seg.lst <- lapply(guide.to.seg.lst, function(x){
        x$dist_to_seg <- rep(1, length(x$dist_to_seg))
        x
      })
    } else if(input.parameter.list$areaOfEffect_type == 'normal'){
      guide.to.seg.lst <- lapply(guide.to.seg.lst, function(x){
        x$dist_to_seg <- dnorm(x$dist_to_seg, mean = 0, sd = input.parameter.list$normal_areaOfEffect_sd) / 
          dnorm(0, mean = 0, sd = input.parameter.list$normal_areaOfEffect_sd)
        
        x$dist_to_seg[which(x$dist_to_seg == max(x$dist_to_seg))] <- 1
        x
      })
    } else if(input.parameter.list$areaOfEffect_type == 'slab_and_spike'){
      
      # create variables to hold flanking distance and slab area of effect
      flanking.distance <- input.parameter.list$flanking.distance
      slab.aoe <- input.parameter.list$slab.aoe
      
      # iterate through every guide in guide.to.seg.lst
      guide.to.seg.lst <- lapply(guide.to.seg.lst, function(x){
        
        # set all segments to have a default AoE of 0
        x$dist_to_seg <- rep(0, length(x$dist_to_seg))
        
        # set AoE to slab value if within flanking distance
        within.flank <- x$dist_to_seg < flanking.distance
        x$dist_to_seg[within.flank] <- input.parameter.list$slab.aoe
        
        # set closest segment to have an AoE of 1 (spike)
        x$dist_to_seg[which(x$dist_to_seg == min(x$dist_to_seg))] <- 1
        
        x
      })
    }
    
  } else if(input.parameter.list$crisprSystem %in% c('dualCRISPR', 'dualCRISPRi', 'dualCRISPRa') ){
    
    if(! 'sg1_Target' %in% names(sim.info)){
      sim.info$sg1_Target <- sim.info$start + input.parameter.list$crisprEffectRange
    }
    if(! 'sg2_Target' %in% names(sim.info)){
      sim.info$sg2_Target <- sim.info$end - input.parameter.list$crisprEffectRange
    }
    
    # no guide offset
    
    sim.guide.seg.list <- adapt_data_to_regionFormat(sim.counts, repl_pools, sim.info, labelHierarchy, min.seg.dist, is.dual.guide = T)
    
    guide.to.seg.lst <- sim.guide.seg.list$guide_to_seg_lst
    
    if(input.parameter.list$areaOfEffect_type == 'uniform'){
      guide.to.seg.lst <- lapply(guide.to.seg.lst, function(x){
        x$dist_to_seg <- rep(1, length(x$dist_to_seg))
        x
      })
    } else if(input.parameter.list$areaOfEffect_type == 'normal'){
      guide.to.seg.lst <- lapply(guide.to.seg.lst, function(x){
        
        min.dist <- dnorm(x$dist_to_seg, mean = 0, sd = input.parameter.list$normal_areaOfEffect_sd) / 
          dnorm(0, mean = 0, sd = input.parameter.list$normal_areaOfEffect_sd)
        
        min.dist[x$sg_overl_idx] <- 1
        
        unif.prob.idx <- min.dist[x$sg_overl_idx[1]:x$sg_overl_idx[2]][which(min.dist[x$sg_overl_idx[1]:x$sg_overl_idx[2]] < input.parameter.list$deletionProb)] <- input.parameter.list$deletionProb
        
        x$dist_to_seg <- min.dist
        
        x
      })
    } else if (input.parameter.list$areaOfEffect_type == 'slab_and_spike'){

      # create variables to hold flanking distance and slab area of effect
      flanking.distance <- input.parameter.list$flanking.distance
      slab.aoe <- input.parameter.list$slab.aoe
      
      # apply slab and spike area of effect for each guide in guide.to.seg.lst
      guide.to.seg.lst <- lapply(guide.to.seg.lst, function(x){
        
        # set baseline AoE of 0
        min.dist <- rep(0, length(x$dist_to_seg))

        # set segments that guide overlaps to 1
        min.dist[x$sg_overl_idx] <- 1
        
        # seg segments within flanking distance of 500bp to slab AoE
        min.dist[x$dist_to_seg < flanking.distance] <- slab.aoe
        
        # return modified dist_to_seg vector accounting for AoE
        x$dist_to_seg <- min.dist
        
        x
      })
    }
  }

  sim.seg.info <- sim.guide.seg.list$seg_info

  # set up the delta vector and the vector containing positions of positive controls
  segment.labels <- as.numeric(as.factor(sim.seg.info$label))

  # next guie index list
  next.guide.list <- generate_next_guide_list(sim.guide.seg.list$seg_to_guide_lst)
  
  guide.dist.to.seg <- NULL
  if(! input.parameter.list$areaOfEffect_type == 'uniform'){
    guide.dist.to.seg <- compute_guide_dist_to_seg(guide.to.seg.lst, sim.guide.seg.list$seg_to_guide_lst, next.guide.list, input.parameter.list$nr_segs)
  }

  format.data.beta <- list(data = sim.guide.seg.list$counts,
                           fs0_idx = which(sim.seg.info$label %in% fs0.label),
                           segLabels = segment.labels,
                           seg_info = sim.seg.info,
                           guide_to_seg_lst = guide.to.seg.lst, #sim.guide.seg.list$guide_to_seg_lst,
                           seg_to_guide_lst = sim.guide.seg.list$seg_to_guide_lst,
                           next_guide_lst = next.guide.list,
                           guide_efficiency_scores = filtered.ges,
                           guide_info = sim.info,
                           guide_dist_to_seg = guide.dist.to.seg)

  write.csv(sim.seg.info, file = paste0(file.save.dir, '_segmentInfo.csv'), row.names = F)
  
  # To Do: adding LFC
  if(! is.null(input.parameter.list$l2fc_groups)){
    
    # plot combined lfc
    # plot per-repl lfc
    
    # plot combined lfc sliding window
    # plot per-repl lfc slinding window
    
    # l2fc_groups, 
    compute_l2fc(input.parameter.list, file.save.dir, sim.counts, sim.info)
  }

  if(save.files){
    write.csv(sim.seg.info, file = paste0(file.save.dir, '_segmentInfo.csv'), row.names = F)
    write.csv(sim.info, file = paste0(file.save.dir, '_guideInfo.csv'), row.names = F)
    write.csv(sim.counts, file = paste0(file.save.dir, '_guideCounts.csv'), row.names = F)
  }

  return(format.data.beta)

}


#' @title pre-ompute distance of the guides to all possible sizes of FS
#' @param guide.to.seg.lst: list, mapping of guides to segments
#' @param seg.to.guide.lst: list, each element is a segment, contains a list with guide_idx, nonGuide_idx
#' @param next.guide.lst: list that contains the index of $next_guide_idx and next_nonGuide_idx
#' @param nr.segs: max number of segemnts that can be combined
#' @return list of lists, each element is a segment, each sub-element a vector containing the distance of all guides overlapping it
#' @export compute_guide_dist_to_seg()

compute_guide_dist_to_seg <- function(guide.to.seg.lst, seg.to.guide.lst, next.guide.lst, nr.segs) {
  
  total.segs <- length(seg.to.guide.lst)
  dist.to.seg.list <- list()
  
  fs.length <- 1
  
  for(seg in 1:total.segs){
    
    temp.guide.idx<- seg.to.guide.lst[[seg]]$guide_idx

    temp.seg.poi.bkg <- unlist(lapply(guide.to.seg.lst[temp.guide.idx], function(x){
      adj.seg.idx <- which(x$seg_overlapped %in% seg )
      temp.dist.to.seg <- x$dist_to_seg[adj.seg.idx]
      dpoibin(0, temp.dist.to.seg)
    }))
    
    dist.to.seg.list[[seg]] <- list()
    dist.to.seg.list[[seg]][[fs.length]] <- temp.seg.poi.bkg
    
  }
  
  if(nr.segs > 1){
    for(ns in 1:(nr.segs-1) ){
      fs.length <- fs.length + 1
      for(seg in 1:(total.segs - ns)){
        
        temp.stretch.segs <- c(seg:(seg + ns)) # index of segments to be included as part of the stretch and be updated in this iteration
        temp.guide.idx <- c(seg.to.guide.lst[[seg]]$guide_idx, unlist(next.guide.lst[(seg + 1):(seg + ns)]))
        
        temp.seg.poi.bkg <- unlist(lapply(guide.to.seg.lst[temp.guide.idx], function(x){
          adj.seg.idx <- which(x$seg_overlapped %in% temp.stretch.segs )
          temp.dist.to.seg <- x$dist_to_seg[adj.seg.idx]
          dpoibin(0, temp.dist.to.seg)
        }))
        
        dist.to.seg.list[[seg]][[fs.length]] <- temp.seg.poi.bkg
        
      }
      
    } 
  }

  return(dist.to.seg.list)
  
}



#' @title compute the log2 fold change of the input data
#' @param par.list: list containing: $l2fc_groups, $l2fc_slidWind_guideNr, $l2fc_slidWind_maxGap, $l2fc_slidWind, and a score column
#' @param out.dir: output directory
#' @param input.counts: counts to compute log2 fold change from
#' @param input.info: info data frame
#' @return log2 flod change files (per replicate) written to specified locations.
#' @export compute_l2fc()

compute_l2fc <- function(par.list, out.dir, input.counts, input.info){
  
  out.bg <- list()
  
  out.l2fc.rs <- c()
  out.l2fc.sw.rs <- c()
  
  out.df.template <- input.info[,c('chrom', 'start', 'end', 'label')]
  
  for(i in 1:length(par.list$l2fc_groups)){
    
    temp.denom <- input.counts[, par.list$l2fc_groups[[i]][1]]
    temp.num <- input.counts[, par.list$l2fc_groups[[i]][2]]
    
    temp.l2fc <- log2( ((temp.num + 1) / sum(temp.num + 1)) / ( (temp.denom + 1) / sum(temp.denom + 1) ) )
    out.l2fc.rs <- rbind(out.l2fc.rs, temp.l2fc)
    
    temp.out.df <- out.df.template
    temp.out.df$l2fc <- temp.l2fc
    
    write.csv(temp.out.df, file = paste0(out.dir, '_l2fc_r', i,'.csv'), row.names = F)
    
    out.bg[[paste0('l2fc_r', i)]] <- temp.out.df
    
    if(par.list$l2fc_slidWind){
      
      # c('chrom', 'windowStart', 'windowEnd', 'label', 'windowScore')
      temp.out.sw <- slide_window(temp.out.df, window.size = par.list$l2fc_slidWind_guideNr, 
                                  max.window.range = par.list$l2fc_slidWind_maxGap, input.score.name = 'l2fc')
      
      out.l2fc.sw.rs <- rbind(out.l2fc.sw.rs, temp.out.sw$windowScore)
      
      write.csv(temp.out.sw, file = paste0(out.dir, '_l2fc_sW_r', i,'.csv'), row.names = F)
      
      out.bg[[paste0('l2fc_sW_r', i)]] <- temp.out.sw
    }
  }
  
  if(length(par.list$l2fc_groups) > 1){
    
    l2fc.rs.means <- colMeans(out.l2fc.rs)
    temp.out.df <- out.df.template
    temp.out.df$l2fc_avg <- l2fc.rs.means
    write.csv(temp.out.df, file = paste0(out.dir, '_l2fc_avg.csv'), row.names = F)
    out.bg[['l2fc_avg']] <- temp.out.df
    
    if(par.list$l2fc_slidWind){
      l2fc.sw.rs.means <- colMeans(out.l2fc.sw.rs)
      temp.out.sw <- temp.out.sw
      temp.out.sw$l2fc_sw_avg <- l2fc.sw.rs.means
      write.csv(temp.out.sw, file = paste0(out.dir, '_l2fc_sW_avg.csv'), row.names = F)
      out.bg[['l2fc_sW_avg']] <- temp.out.sw
    }
    
  }
  
  create_bedgraphs(out.bg, out.dir)
  
}

#' @title merge overlapping information using a sliding window
#' @param input.df: data frame containing: $chrom, $start, $end, $label, and a score column
#' @param window.size: numeric, number of rows to be included in sliding
#' @param max.window.range: max. distance between target sites
#' @param input.label.hierarchy: order by which multiple labels are combined. By default priority is given to label with the least occurances
#' @param input.score.name: name of the column containing the scores
#' @return date frame containing $chrom, $start, $end, $label, and a score column
#' @export slide_window()

slide_window <- function(input.df, window.size, max.window.range = 3000,input.label.hierarchy = NULL, input.score.name = NULL){
  
  input.positions <- input.df
  
  input.positions.chroms <- split(input.positions, input.positions$chrom)
  
  out.df <- c()
  
  for(chr in 1:length(input.positions.chroms)){
    
    input.positions.ordered <- input.positions.chroms[[chr]][order(input.positions.chroms[[chr]]$start),]
    
    label.hierarchy <- c()
    if(is.null(input.label.hierarchy)){
      counted.lables <- table(input.positions.ordered$label)
      
      while(length(counted.lables) > 1){
        max.lab.pos <- which(counted.lables == max(counted.lables)[1])[1]
        max.lab.name <- names(counted.lables)[max.lab.pos]
        
        label.hierarchy<- c(label.hierarchy, max.lab.name)
        
        counted.lables <- counted.lables[-max.lab.pos]
      }
      
      # then last part
      label.hierarchy <- c(label.hierarchy, names(counted.lables))
      
    } else {
      label.hierarchy <- input.label.hierarchy
    }
    
    score.name <- 'guideScore'
    if(! is.null(input.score.name)){
      score.name <- input.score.name
    }
    
    window.scores <- vector('numeric', length = nrow(input.positions.ordered))
    window.start <- vector('numeric', length = nrow(input.positions.ordered))
    window.end <- vector('numeric', length = nrow(input.positions.ordered))
    window.label <- vector('character', length = nrow(input.positions.ordered))
    
    w <- 0 # i is window index
    for(i in 1:(nrow(input.positions.ordered)-1) ){
      temp.window.size <- min(window.size, nrow(input.positions.ordered) - i + 1)
      window.range <- input.positions.ordered$end[i+temp.window.size-1] - input.positions.ordered$start[i]
      if(window.range < max.window.range){
        w <- w + 1
        window.scores[w] <- mean(input.positions.ordered[[score.name]][i:(i+temp.window.size-1)])
        window.start[w] <- input.positions.ordered$start[i]
        window.end[w] <- input.positions.ordered$end[i+temp.window.size-1]
        window.label[w] <- label.hierarchy[max(which(label.hierarchy %in% input.positions.ordered$label[i:(i+temp.window.size-1)]))[1]]
        
      } else {
        for(j in 1:(temp.window.size-1)){
          smaller.window.range <- input.positions.ordered$end[i+temp.window.size-1 - j] - input.positions.ordered$start[i]
          if(smaller.window.range < max.window.range){
            w <- w + 1
            window.scores[w] <- mean(input.positions.ordered[[score.name]][i:(i+temp.window.size-1 - j)])
            window.start[w] <- input.positions.ordered$start[i]
            window.end[w] <- input.positions.ordered$end[i+temp.window.size-1 - j]
            window.label[w] <- label.hierarchy[max(which(label.hierarchy %in% input.positions.ordered$label[i:(i+temp.window.size-1 - j)]))[1]]
          }
        }
      }
    }
    
    temp.out.df <- cbind.data.frame(rep(input.positions.ordered$chrom[1],w),
                               window.start[1:w], window.end[1:w], window.label[1:w], window.scores[1:w], stringsAsFactors = F)
    
    out.df <- rbind(out.df, temp.out.df)
    
  }

  colnames(out.df) <- c('chrom', 'start', 'end', 'label', 'windowScore')
  
  return(out.df)
}


#' @title Adapt an info file to a overlap matrix and corresponding data structures and filter the counts, return by replicate
#' @param input.counts: data.frame, rows are guides, columns are pools
#' @param replicate.list: each element is a replicate
#' @param input.info: data frame: each row contains chrom, start, end, a label
#' @param input.label.hierarchy: ordering lof guide labels, left to right signifies ordering hor hierarchy to override in case of multiple labels overlapping
#' @param min.seg.dist: minimum distance of a segment. Default = 100
#' @param is.dual.guide, logical, whether screen type is dual CRISPR
#' @return list: seg_info, guide_to_seg_lst, seg_to_guide_lst, counts (list, per replicate)
#' @export adapt_data_to_regionFormat()

adapt_data_to_regionFormat <- function(input.counts, replicate.list, input.info, input.label.hierarchy, min.seg.dist = 100, is.dual.guide = FALSE){

  info.targeting.df <- c()
  counts.targeting.df <- c()
  if(length(which(is.na(input.info$start))) > 0){
    print('Error: data was not properly filtered!')
    break()
  }

  info.targeting.df <- input.info
  counts.targeting.df <- input.counts

  # first create guide-segment matrix for targeting guides: seg_info, guide_to_seg_lst, seg_to_guide_lst
  targeting.gs.list <- create_targeting_guide_segment_matrix(info.targeting.df, input.label.hierarchy, min.seg.dist, is.dual.guide)

  targeting.gs.list$counts <- list()
  for(i in 1:length(replicate.list)){
    temp.counts <- counts.targeting.df[, replicate.list[[i]]]
    # colnames(temp.counts) <- paste('pool', c(1:ncol(temp.counts)), sep = '_')
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
#' @param is.dual.guide, logical, whether screen type is dual CRISPR
#' @return list: seg_info, guide_to_seg_lst, seg_to_guide_lst
#' @export create_targeting_guide_segment_matrix()

create_targeting_guide_segment_matrix <- function(input.targeting.info, input.label.hierarchy, 
                                                  min.seg.dist = 100, is.dual.guide = FALSE){


  info.chroms <- split(input.targeting.info, input.targeting.info$chrom)
  region.df <- c()

  for(i in 1:length(info.chroms)){
    temp.info <- info.chroms[[i]]
    all.breaks <- unique(c(temp.info$start, temp.info$end))
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

    temp.region.df <- data.frame(chrom = rep(temp.info$chrom[1], length(bin.start)),
                            start = bin.start, end = bin.ends, stringsAsFactors = F)

    region.df <- rbind(region.df, temp.region.df)
  }

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
  guide.to.seg.list <- generate_guide_to_seg_list(guide.overlap.list, input.targeting.info, region.df.filtered, is.dual.guide)

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
#' @return list. each element corresponds to a segment, contains a list: guide_idx, nonGuide_idx
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
#' @param input.info: info data frame. has to contain $start, $end, $sgTarget
#' @param input.segments: data frame of segment coordinates
#' @param is.dual.guide, logical; whether the screen is single or dual guide
#' @return list of lists. each element is a list: $seg_overlapped = containing a vector of indexes, indexing the segments overlapped by it. $dist_to_seg, distance of sgTarget to segment
#' @export generate_guide_to_seg_list()

generate_guide_to_seg_list <- function(input.guide.overlaps, input.info, input.segments, is.dual.guide = FALSE){
  
  out.guide.to.seg.lst <- c()
  
  if(! is.dual.guide){
    out.guide.to.seg.lst <- lapply(input.guide.overlaps, function(x){
      #x$subjectHits
      temp.guide <- input.info[x$queryHits[1],]
      
      temp.segments <- input.segments[x$subjectHits,,drop = F]
      
      dists.to.seg <- vector('numeric', length = length(x$subjectHits))
      
      for(i in 1:nrow(temp.segments)){
        dists.to.seg[i] <- min(c(abs(temp.guide$sgTarget - temp.segments$start[i]), 
                                 abs(temp.guide$sgTarget - temp.segments$end[i])))
      }
      
      list(seg_overlapped = x$subjectHits,
           dist_to_seg = dists.to.seg)
    })
    
  } else {
    
    out.guide.to.seg.lst <- lapply(input.guide.overlaps, function(x){
      temp.guide <- input.info[x$queryHits[1],]
      
      temp.segments <- input.segments[x$subjectHits,,drop = F]
      
      dists.to.seg <- vector('numeric', length = length(x$subjectHits))
      g1.min.dist <- 1e6
      g2.min.dist <- 1e6
      sg.overl.idx <- vector('numeric', length = 2)
      
      for(i in 1:nrow(temp.segments)){
        temp.g1.dist <- min(c(abs(temp.guide$sg1_Target - temp.segments$start[i]), 
                                 abs(temp.guide$sg1_Target - temp.segments$end[i])))
        temp.g2.dist <- min(c(abs(temp.guide$sg2_Target - temp.segments$start[i]), 
                                 abs(temp.guide$sg2_Target - temp.segments$end[i])))
        
        dists.to.seg[i] <- min(temp.g1.dist, temp.g2.dist)
        
        if(temp.g1.dist < g1.min.dist){
          g1.min.dist <- temp.g1.dist
          sg.overl.idx[1] <- i
        } 
        if(temp.g2.dist < g2.min.dist){
          g2.min.dist <- temp.g2.dist
          sg.overl.idx[2] <- i
        }
      }
      
      list(seg_overlapped = x$subjectHits,
           dist_to_seg = dists.to.seg,
           sg_overl_idx = sg.overl.idx)
    })
    
    
  }

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


#' @title estimate the dirichlet sorting proportions (dispersions are already provided)
#' @param data.par.list: list, contains the elements of the processed and formatted data
#' @param analysis.par.list: list, contains all analysis flags
#' @param input.repl.pools: list, each element is a set of columns which correspond to a replicate
#' @param fs0.label: label used for generating FS0
#' @return list of lists: $hyper_pars $hyper_par_components: hyper parameter estimates for each replicate per list element, $bkg_alpha, $FS_alpha, $bkg_dispersion
#' @export estimate_dirichlet_proportions()

estimate_dirichlet_proportions <- function(data.par.list, analysis.par.list, input.repl.pools, fs0.label){
  
  # most basic implementation of GE
  if(! is.null(data.par.list$guide_efficiency_scores)){
    logit.ge <- apply(data.par.list$guide_efficiency_scores, 2, function(x){
      log(x / (1 - x))
    })
    data.par.list$guide_efficiency <- 1 / (1 + exp(-(logit.ge %*% rep(1, ncol(data.par.list$guide_efficiency_scores)))))
  }
  
  # make sure all guides are considered to be the same category
  fs.assignment <- rep(0, nrow(data.par.list$seg_info))
  fs.assignment[which(data.par.list$seg_info$label %in% fs0.label)] <- 1
  
  dirichlet.guide.ll <- compute_perGuide_fs_ll(fs.assignment, data.par.list$guide_to_seg_lst, hyper.setup = TRUE)
  
  final.dirichlet.pars <- list()
  final.dirichlet.pars$bkg_alpha <- list()
  final.dirichlet.pars$FS_alpha <- list()
  final.dirichlet.pars$dispersion <- list()

  final.alpha <- list()
  final.alpha$alpha0 <- list()
  final.alpha$alpha1 <- list()
  
  orig.fs.par <- list()
  
  model.ll <- 0
  
  for(i in 1:length(data.par.list$data)){

    temp.bkg.alpha <- rep(1, length(input.repl.pools[[i]]) - 1)
    temp.fs.alpha <- rep(1, length(input.repl.pools[[i]]) - 1)

    temp.hyper.params <- c(temp.bkg.alpha, temp.fs.alpha) #, temp.bkg.disp)
    temp.bkg.idx <- c(1:length(temp.bkg.alpha))
    temp.fs.idx <- c(1:length(temp.fs.alpha)) + max(temp.bkg.idx)

    temp.res.drch <- c()
    repl.data <- data.par.list$data[[i]]
    temp.repl.disp <- analysis.par.list$repl_disp[[i]]$repl_disp
    
    temp.res.drch <- optim(temp.hyper.params, prior_dirichlet_proportions, method= 'L-BFGS-B',
                           data = repl.data, 
                           region.ll.list = dirichlet.guide.ll,
                           bkg.idx = temp.bkg.idx, 
                           fs.idx = temp.fs.idx, 
                           guide.efficiency = data.par.list$guide_efficiency, 
                           repl.disp = temp.repl.disp,
                           model.disp = analysis.par.list$model_dispersion)
    
    bkg.par <- temp.res.drch$par[temp.bkg.idx]**2
    fs.par <- temp.res.drch$par[temp.fs.idx]**2
    orig.fs.par[[i]] <- fs.par

    if(analysis.par.list$hyper_adj < 1){
      fs.par <- adjust_hypers(bkg.par, fs.par, analysis.par.list$hyper_adj)
    }
    
    temp.bkg.alpha <- c(1, bkg.par) / sum(c(1, bkg.par))
    temp.fs.alpha <- c(1, fs.par) / sum(c(1, fs.par))
    
    final.dirichlet.pars$bkg_alpha[[i]] <- bkg.par
    final.dirichlet.pars$FS_alpha[[i]] <- fs.par
    final.dirichlet.pars$dispersion[[i]] <- temp.repl.disp
    
    temp.hyper <- reparameterize_hypers(temp.bkg.alpha, temp.fs.alpha, temp.repl.disp, data.par.list$guide_efficiency, analysis.par.list$model_dispersion)
    final.alpha$alpha0[[i]] <- temp.hyper$alpha0
    final.alpha$alpha1[[i]] <- temp.hyper$alpha1
    # final.alpha$alpha0[[i]] <- t(temp.bkg.alpha %*% t(temp.repl.disp) )
    # final.alpha$alpha1[[i]] <- t(temp.fs.alpha %*% t(temp.repl.disp) )

    model.ll <- model.ll + temp.res.drch$value
  }
  
  final.alpha$L <- 1

  out.list <- list(hyper_pars = final.alpha,
                   hyper_par_components = final.dirichlet.pars,
                   init_model_ll = model.ll,
                   orig_hyper_par_components = orig.fs.par)
  return(out.list)
}


#' @title adjust the hyper parameters to detect elements which are only part of the strength of the original signal
#' @param input.bkg: background hyperparameter
#' @param input.fs: fs hyperparameter
#' @param prct.adj: by what percent the hyperparameters is similar to FS0 (at 0.9 the FS0 signal is at 0.9*FS0)
#' @return adjusted FS hyperparameter
#' @export adjust_hypers()

adjust_hypers <- function(input.bkg, input.fs, prct.adj){
  
  fs.adj.out <- c()
  bkg.prop <- c(1, input.bkg) / sum(c(1, input.bkg))
  fs.prop <- c(1, input.fs) / sum(c(1, input.fs))
  
  fs.prop.adj <- (bkg.prop - fs.prop) * (1-prct.adj) + fs.prop
  
  fs.adj.out <- (fs.prop.adj / fs.prop.adj[1])[2:length(fs.prop.adj)]
  
  return(fs.adj.out)
}


#' @title optimize the dirichlet hyper parameter proportions, dispersion is given
#' @param hyper.param: hyperparameters
#' @param data: data, consists of: $pool1, pool2... and $n
#' @param region.ll.list: list containing the ll and the indexes of the null alternative and both
#' @param bkg.idx: positions in the par vector of the background alphas
#' @param fs.idx: positions in the par vector of the FS alphas
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @param repl.disp: per-guide dispersion based on total counts
#' @param model.disp: whether the dispersion is modeled
#' @return sum of the -log likelihood across all guides
#' @export prior_dirichlet_proportions()

prior_dirichlet_proportions <- function(hyper.param, data, region.ll.list, bkg.idx, fs.idx, 
                                       guide.efficiency, repl.disp, model.disp) {
  
  bkg.alpha <- c(1, hyper.param[bkg.idx]**2) / sum(c(1, hyper.param[bkg.idx]**2))
  fs.alpha <- c(1, hyper.param[fs.idx]**2) / sum(c(1, hyper.param[fs.idx]**2))

  # hyper <- list()
  # hyper$alpha0 <- t(bkg.alpha %*% t(repl.disp) )
  # hyper$alpha1 <- t(fs.alpha %*% t(repl.disp) )
  
  hyper <- reparameterize_hypers(bkg.alpha, fs.alpha, repl.disp, guide.efficiency, model.disp)

  out.sg.ll <- estimate_relics_sgrna_log_like(hyper, data, region.ll.list, guide.efficiency)
  
  -sum(out.sg.ll$total_guide_ll)
}


#' @title estimate the hyperparamters using guide efficiency
#' @param input.bkg.alpha: dirichlet proportions for the background
#' @param input.fs.alpha: dirichlet proportions for FS
#' @param input.disp: vector, dispersion
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @param model.disp: logical, whether the dispersion is modeled
#' @return list of hyperparameters
#' @export reparameterize_hypers()

reparameterize_hypers <- function(input.bkg.alpha, input.fs.alpha, input.disp, guide.efficiency, model.disp){
  
  alpha0 <- c()
  if(model.disp){
    alpha0 <- t(input.bkg.alpha %*% t(input.disp) )
  } else {
    alpha0 <- input.bkg.alpha * input.disp
  }
  
  alpha1 <- c()
  
  if(is.null(guide.efficiency)){
    if(model.disp){
      alpha1 <- t(input.fs.alpha %*% t(input.disp) )
    } else {
      alpha1 <- input.fs.alpha * input.disp
    }
  } else {
    alpha0.matrix <- t(apply(matrix(0, ncol = length(input.bkg.alpha), nrow = length(guide.efficiency)), 1, function(x){x + input.bkg.alpha}))
    alpha.diffs <- input.bkg.alpha - input.fs.alpha
    
    alpha1.matrix <- matrix(0, ncol = length(input.bkg.alpha), nrow = length(guide.efficiency))
    
    for(i in 1:length(guide.efficiency)){
      # this cound be an issue if the dispersion is no longer identical...
      temp.alpha1 <- alpha0.matrix[i,] - alpha.diffs * guide.efficiency[i]
      
      alpha1.matrix[i,] <- temp.alpha1
    }
    
    alpha1 <- alpha1.matrix * input.disp
  }
  
  out.list <- list(alpha0 = alpha0, alpha1 = alpha1)
  
  return(out.list)
  
}


#' @title estimate the hyperparamters for the dirichlet proportions and dispersion (dispersion assumed to be independent of the counts)
#' @param data.par.list: list, contains the elements of the processed and formatted data
#' @param analysis.par.list: list, contains all analysis flags
#' @param input.repl.pools: list, each element is a set of columns which correspond to a replicate
#' @param fs0.label: label used for generating FS0
#' @return list of lists: $hyper_pars $hyper_par_components: hyper parameter estimates for each replicate per list element, $bkg_alpha, $FS_alpha, $bkg_dispersion
#' @export estimate_hyper_parameters()

estimate_hyper_parameters <- function(data.par.list, analysis.par.list, input.repl.pools,  fs0.label){

  # most basic implementation of GE
  if(! is.null(data.par.list$guide_efficiency_scores)){
    logit.ge <- apply(data.par.list$guide_efficiency_scores, 2, function(x){
      log(x / (1 - x))
    })
    data.par.list$guide_efficiency <- 1 / (1 + exp(-(logit.ge %*% rep(1, ncol(data.par.list$guide_efficiency_scores)))))
  }

  # make sure all guides are considered to be the same category
  fs.assignment <- rep(0, nrow(data.par.list$seg_info))
  fs.assignment[which(data.par.list$seg_info$label %in% fs0.label)] <- 1

  dirichlet.guide.ll <- compute_perGuide_fs_ll(fs.assignment, data.par.list$guide_to_seg_lst, hyper.setup = TRUE)

  final.dirichlet.pars <- list()
  final.dirichlet.pars$bkg_alpha <- list()
  final.dirichlet.pars$FS_alpha <- list()
  final.dirichlet.pars$bkg_dispersion <- list()

  final.alpha <- list()
  final.alpha$alpha0 <- list()
  final.alpha$alpha1 <- list()
  
  orig.fs.par <- list()
  
  model.ll <- 0

  for(i in 1:length(data.par.list$data)){

    temp.bkg.disp <- c() 
    if(analysis.par.list$mean_var_type == 'independent'){
      temp.bkg.disp <- c(1)
    }
    
    temp.bkg.alpha <- rep(1, length(input.repl.pools[[i]]) - 1)
    temp.fs.alpha <- rep(1, length(input.repl.pools[[i]]) - 1)
 
    temp.hyper.params <- c(temp.bkg.alpha, temp.fs.alpha, temp.bkg.disp)
    temp.bkg.idx <- c(1:length(temp.bkg.alpha))
    temp.fs.idx <- c(1:length(temp.fs.alpha)) + max(temp.bkg.idx)
    temp.disp.idx <- c(1:length(temp.bkg.disp)) + max(temp.fs.idx)

    temp.res.drch <- c()
    repl.data <- data.par.list$data[[i]]
    
    temp.res.drch <- optim(temp.hyper.params, prior_dirichlet_parameters, method= 'L-BFGS-B',
                           data = repl.data, 
                           region.ll.list = dirichlet.guide.ll,
                           bkg.idx = temp.bkg.idx, 
                           fs.idx = temp.fs.idx, 
                           disp.idx = temp.disp.idx, 
                           guide.efficiency = data.par.list$guide_efficiency, 
                           mean.var.type = analysis.par.list$mean_var_type,
                           model.disp = analysis.par.list$model_dispersion)
    
    bkg.par <- temp.res.drch$par[temp.bkg.idx]**2
    fs.par <- temp.res.drch$par[temp.fs.idx]**2
    orig.fs.par[[i]] <- fs.par
    
    if(analysis.par.list$hyper_adj < 1){
      fs.par <- adjust_hypers(bkg.par, fs.par, analysis.par.list$hyper_adj)
    }
    
    temp.bkg.alpha <- c(1, bkg.par) / sum(c(1, bkg.par))
    temp.fs.alpha <- c(1, fs.par) / sum(c(1, fs.par))
    temp.bkg.disp <- temp.res.drch$par[temp.disp.idx]
    
    final.dirichlet.pars$bkg_alpha[[i]] <- bkg.par
    final.dirichlet.pars$FS_alpha[[i]] <- fs.par
    final.dirichlet.pars$bkg_dispersion[[i]] <- temp.bkg.disp
    
    if(analysis.par.list$mean_var_type == 'independent'){
      temp.disp <- rep(temp.bkg.disp[1]^2, nrow(repl.data))
      temp.hyper <- reparameterize_hypers(temp.bkg.alpha, temp.fs.alpha, temp.disp, data.par.list$guide_efficiency, TRUE)
      final.alpha$alpha0[[i]] <- temp.hyper$alpha0
      final.alpha$alpha1[[i]] <- temp.hyper$alpha1
      # final.alpha$alpha0[[i]] <- t(temp.bkg.alpha %*% t(temp.disp) )
      # final.alpha$alpha1[[i]] <- t(temp.fs.alpha %*% t(temp.disp) )
    }
    
    model.ll <- model.ll + temp.res.drch$value

  }

  final.alpha$L <- 1

  out.list <- list(hyper_pars = final.alpha,
                   hyper_par_components = final.dirichlet.pars,
                   init_model_ll = model.ll,
                   orig_hyper_par_components = orig.fs.par)
  return(out.list)
}


#' @title optimize the hyper parameters, only one dispersion across the two distributions
#' @param hyper.param: hyperparameters
#' @param data: data, consists of: $pool1, pool2... and $n
#' @param region.ll.list: list containing the ll and the indexes of the null alternative and both
#' @param bkg.idx: positions in the par vector of the background alphas
#' @param fs.idx: positions in the par vector of the FS alphas
#' @param disp.idx: positions in the par vector of the dispersion parameters
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @param mean.var.type: type of mean-variance relationship
#' @param model.disp: whether the dispersion is modeled
#' @return sum of the -log likelihood across all guides
#' @export prior_dirichlet_parameters()

prior_dirichlet_parameters <- function(hyper.param, data, region.ll.list, bkg.idx, fs.idx, disp.idx, 
                                       guide.efficiency, mean.var.type, model.disp) {
  
  bkg.alpha <- c(1, hyper.param[bkg.idx]**2) / sum(c(1, hyper.param[bkg.idx]**2))
  fs.alpha <- c(1, hyper.param[fs.idx]**2) / sum(c(1, hyper.param[fs.idx]**2))
  disp.alpha <- hyper.param[disp.idx]

  temp.disp <- c()
  # hyper <- list()
  if(mean.var.type == 'independent'){
    temp.disp <- disp.alpha[1]^2
    # hyper$alpha0 <- bkg.alpha * temp.disp
    # hyper$alpha1 <- fs.alpha * temp.disp
  }
  
  hyper <- reparameterize_hypers(bkg.alpha, fs.alpha, temp.disp, guide.efficiency, FALSE)

  out.sg.ll <- estimate_relics_sgrna_log_like(hyper, data, region.ll.list, guide.efficiency)

  -sum(out.sg.ll$total_guide_ll)
}


#' @title Compute log likelihoods of observed counts for each sgRNA, given. Formerly 'compute_dirichlet_sgrna_log_like'
#' @param hyper: hyperparameters, $alpha0, alpha1
#' @param data: data, consists of: pool1, pool2... and $n
#' @param region.ll.list: list containing the ll and the indexes of the null alternative and both
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @param return.model.ll: logical, whether the function should terminate early and return only the null and alternative model ll
#' @return data frame: total_guide_ll, alt_only_ll
#' @export estimate_relics_sgrna_log_like()

estimate_relics_sgrna_log_like <- function(hyper, data, region.ll.list, guide.efficiency, return.model.ll = FALSE, hyper.components){

  pool.cols <- c(1:(ncol(data) - 1))

  sgrna.null.log.like <- ddirmnom(data[,pool.cols], size = data[,ncol(data)], alpha = hyper$alpha0, log = T)
  sgrna.alt.log.like <- ddirmnom(data[,pool.cols], size = data[,ncol(data)], alpha = hyper$alpha1, log = T)
  

  # sgrna.alt.log.like <- c()
  # if(is.null(guide.efficiency)){
  #   sgrna.alt.log.like <- ddirmnom(data[,pool.cols], size = data[,ncol(data)], alpha = hyper$alpha1, log = T)
  # } else {
  #   # toDo, difference in alphas or in proportions? 
  #   alpha0.matrix <- t(apply(matrix(0, ncol = length(hyper$alpha0), nrow = length(guide.efficiency)), 1, function(x){x + hyper$alpha0}))
  #   alpha.diffs <- hyper$alpha0 - hyper$alpha1
  # 
  #   alpha1.matrix <- matrix(0, ncol = length(hyper$alpha0), nrow = length(guide.efficiency))
  # 
  #   for(i in 1:length(guide.efficiency)){
  #     # this cound be an issue if the dispersion is no longer identical...
  #     alpha1.matrix[i,] <- alpha0.matrix[i,] - alpha.diffs * guide.efficiency[i]
  #   }
  # 
  #   sgrna.alt.log.like <- ddirmnom(data[,pool.cols], size = data[,ncol(data)], alpha = alpha1.matrix, log = T)
  # }

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
    region.pp <- cumulative.pp[guide.seg.idx.lst[[j]]$seg_overlapped]
    
    region.pp.distance.adj <- region.pp * guide.seg.idx.lst[[j]]$dist_to_seg
    
    p.k.eq.0 <- 0
    
    if(max(region.pp.distance.adj) < 1){
      # compute probability that number of regulatory regions overlapped by this sgRNA is 0 or >0, given posterior
      # probabilities, using poisson binomial probability mass function
      p.k.eq.0 <- dpoibin(0, region.pp.distance.adj)
    } else if(max(region.pp.distance.adj) == 0){
      p.k.eq.0 <- 1
    }

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

  # if there are overlaps between the indexes that's a problem!
  if(sum(c((null.only.counter - 1), alt.only.counter - 1, both.counter - 1)) != length(guide.seg.idx.lst) & !hyper.setup){
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
#' @param pool.names, if not NULL, the names of the pools to be used when recording output
#' @param mean.var.type: type of mean-variance relationship
#' @param pp.calculation: what version to use when calculating PP. 'v2'is for computing with normal AoE
#' @param analysis.parameters: list, contains various elements used for analysis
#' @param guide.dist.to.seg, list of lists, pre-computed distances of guides to FS
#' @return list of final per-layer posteriors
#' @export run_RELICS_2()

run_RELICS_2 <- function(input.data, final.layer.nr, out.dir = NULL,
                     fix.hypers = FALSE,
                     iterative.hyper.est = FALSE,
                     input.hypers = NULL,
                     input.hyper.components,
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
                     local.max, local.max.range,
                     pool.names,
                     mean.var.type,
                     pp.calculation,
                     analysis.parameters,
                     guide.dist.to.seg){
  
  #posteriors, alpha0, alpha1, bkg_hyper, fs_hyper, bkg_disp, fs_ll, per_fs_ll, correlation, ge_coeff, model_ll, ll_ratio
  fs.result.lists <- initialize_results_list()

  relics.hyper <- input.hypers #c()
  hyper.components <- input.hyper.components

  relics.param <- init_relics_param(relics.hyper, input.data, local.max, recompute.fs0)

  coverg.tol <- rep(input.convergence.tol, final.layer.nr)

  for(i in 1:final.layer.nr){
    print(paste0('Computing FS: ', i))
    
    fs.data <- relics_compute_FS_k(input.param = relics.param,
                                   input.hyper = relics.hyper,
                                   input.data.list = input.data,
                                   input.tol = coverg.tol[i],
                                   nr.segs, geom.p,
                                   min.pp = input.min.rs.pp, input.data$guide_efficiency,
                                   one.dispersion,
                                   local.max, local.max.range, analysis.parameters$areaOfEffect_type,
                                   guide.dist.to.seg)
    
    # fs.result.lists <- process_fs_pp(fs.data, fs.result.lists, analysis.parameters, input.data, i, relics.hyper, hyper.components)
    current.result.lists <- process_fs_pp_output(fs.data, fs.result.lists, analysis.parameters, input.data, i, relics.hyper, hyper.components)
   
    if(record.all.fs){
      
      record_results(current.result.lists, i, input.data, analysis.parameters, '', hyper.components, relics.hyper)
      
    }

    if(i > 1){

      fs.correlation.cutoff.list <- determine_FS_nr_cutoff(current.result.lists$correlation, analysis.parameters$fs_correlation_cutoff)

      if(fs.correlation.cutoff.list$need_to_stop){
        print(fs.correlation.cutoff.list$why_to_stop)
        
        if(auto.stop){
          record_results(fs.result.lists, i - 1, input.data, analysis.parameters, '_final', hyper.components, relics.hyper)
          break()
        } else {
          print(paste0('Recommend stopping RELICS 2'))
          print('Continuing to run because auto.stop was set to FALSE')
          record_results(fs.result.lists, i - 1, input.data, analysis.parameters, '_recommendedFinal', hyper.components, relics.hyper)
        }

      }

    }
    
    # update the results
    fs.result.lists <- current.result.lists
    
    # prep the parameters for the next iteration (both the posteriors and hypers)
    relics.hyper$L <- relics.hyper$L + 1
    relics.param$delta.pp <- rbind(fs.result.lists$posteriors, rep(0, ncol(relics.param$delta.pp)))  # fs.result.lists$posteriors[[i]]
    if(local.max){
      relics.param$ll_rt <- rbind(fs.result.lists$ll_ratio, rep(0, ncol(relics.param$delta.pp)))  # fs.result.lists$ll_ratio[[i]]
    }

    # if the hyper parameters are reestimated after convergence of the posteriors
    if(! fix.hypers & !iterative.hyper.est){
      
      relics.hyper.list <- c()
      if(analysis.parameters$model_dispersion){
        relics.hyper.list <- recompute_dirichlet_hypers(relics.param, 
                                                        input.data$guide_to_seg_lst, 
                                                        hyper.components, 
                                                        input.data$data, 
                                                        input.data$guide_efficiency, 
                                                        relics.hyper,
                                                        analysis.parameters)
      } else {
        relics.hyper.list <- recompute_hyper_parameters(relics.param,
                                                        relics.hyper,
                                                        hyper.components,
                                                        input.data$data,
                                                        input.data$guide_to_seg_lst,
                                                        input.data$guide_efficiency,
                                                        analysis.parameters)
      }

      relics.hyper <- relics.hyper.list$hyper_pars
      hyper.components <- relics.hyper.list$hyper_par_components
    }

    if(! is.null(input.data$guide_efficiency_scores)){
      if(! input.data$fixed_ge_coeff){
        ge.list <- recompute_ge_coefficients(relics.param,
                                             hyper.components, #relics.hyper,
                                             input.data$data,
                                             input.data$guide_to_seg_lst,
                                             input.data$guide_efficiency_scores,
                                             input.data$ge_coeff)

        input.data$guide_efficiency <- ge.list$guide_efficiency
        input.data$ge_coeff <- ge.list$ge_coeff
      }

    }

  }

}


#' @title wrapper for reporting the per-segment ll ratio while accounting for the area of effect
#' @param relics.hyper: list: hyper parameters
#' @param input.data: list: $seg_to_guide_lst, $guide_efficiency_scores (and $guide_efficiency if former not NULL)
#' @param analysis.parameters: list containing all the analysis parameters
#' @param fs.iter: current FS interation
#' @param file.extension: extended label for the file
#' @return list with parameters for running RELICS or logical FALSE if required parameter is missing
#' @export ll_ratio_recording()

ll_ratio_recording <- function(input.data, analysis.parameters, relics.hyper, fs.iter, file.extension = ''){
  
  out.dir <- paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName)
  
  if(!analysis.parameters$fix_hypers | (!is.null(input.data$guide_efficiency_scores) &! input.data$fixed_ge_coeff) ){
    out.pars <- list(out_dir = out.dir,
                     iter = fs.iter)
    # record_AoE_ll_ratio(relics.hyper, input.data, out.pars, file.extension) # record_ll_ratio
    record_ll_ratio(relics.hyper, input.data, out.pars, analysis.parameters, file.extension)
    record_guide_ll_ratio(relics.hyper, input.data, out.pars, analysis.parameters, file.extension)
  }
  
}


#' @title given hyper parameters and guide efficiency (optional), report the per-guide ll ratio to extract efficient guide (pairs)
#' @param input.hypers: list: hyper parameters
#' @param input.data: list: $seg_to_guide_lst, $guide_efficiency_scores (and $guide_efficiency if former not NULL)
#' @param input.ge: guide efficiency
#' @param input.seg.to.guide.lst: list: hyper parameters
#' @return list with parameters for running RELICS or logical FALSE if required parameter is missing
#' @export record_guide_ll_ratio()

record_guide_ll_ratio <- function(input.hypers, input.data, out.pars, input.parameters, file.extension){
  
  not.used <- c()
  guide.efficiency <- NULL
  if(! is.null(input.data$guide_efficiency_scores)){
    guide.efficiency <- input.data$guide_efficiency
  }
  
  guide.ll.rt <- vector('numeric', length = nrow(input.data$guide_info))
  for(i in 1:length(input.data$data)){
    temp.data <- input.data$data[[i]]
    temp.hyper <- list(alpha0 = input.hypers$alpha0[[i]],
                       alpha1 = input.hypers$alpha1[[i]])
    
    guide.model.ll <- estimate_relics_sgrna_log_like(temp.hyper, temp.data, not.used, guide.efficiency, return.model.ll = TRUE)
    temp.guide.ll.rt <- -2*(guide.model.ll$null_only_ll - guide.model.ll$alt_only_ll)
    guide.ll.rt <- guide.ll.rt + temp.guide.ll.rt
  }
  
  # combine with segment info
  guide.info <- input.data$guide_info
  guide.info$score <- round(guide.ll.rt, 3)
  
  write.csv(guide.info, file = paste0(out.pars$out_dir, file.extension, '_FS', out.pars$iter, '_guideLLR.csv'), row.names = F)
  
  # # create bedgraph
  # to.bg.list <- list(seg_llRt = segment.info)
  # 
  # # write bedgraph to output
  # # out.dir <- paste0(out.pars$out_dir, '_FS', out.pars$iter)
  # out.dir <- paste0(out.pars$out_dir, file.extension, '_FS', out.pars$iter)
  # create_bedgraphs(to.bg.list, out.dir)
}



#' @title given hyper parameters and guide efficiency (optional), report the per-segment ll ratio while accounting for the area of effect
#' @param input.hypers: list: hyper parameters
#' @param input.data: list: $seg_to_guide_lst, $guide_efficiency_scores (and $guide_efficiency if former not NULL)
#' @param input.parameters: list containing all the analysis parameters
#' @param file.extension: extended label for the file
#' @return list with parameters for running RELICS or logical FALSE if required parameter is missing
#' @export record_AoE_ll_ratio()

record_AoE_ll_ratio <- function(input.hypers, input.data, input.parameters, file.extension = ''){
  
  not.used <- c()
  guide.efficiency <- NULL
  if(! is.null(input.data$guide_efficiency_scores)){
    guide.efficiency <- input.data$guide_efficiency
  }
  
  ll.rt <- vector('numeric', length = length(input.data$seg_to_guide_lst))
  for(i in 1:length(input.data$data)){
    temp.data <- input.data$data[[i]]
    temp.hyper <- list(alpha0 = input.hypers$alpha0[[i]],
                       alpha1 = input.hypers$alpha1[[i]])
    
    guide.model.ll <- estimate_relics_sgrna_log_like(temp.hyper, temp.data, not.used, guide.efficiency, return.model.ll = TRUE)
    temp.ll.rt <- compute_local_AoE_ll_ratio(guide.model.ll, input.data$seg_to_guide_lst, input.data$guide_to_seg_lst)
    ll.rt <- ll.rt + temp.ll.rt
  }
  
  # combine with segment info
  segment.info <- input.data$seg_info
  segment.info$score <- ll.rt
  
  # create bedgraph
  to.bg.list <- list(seg_llRt = segment.info)
  
  # write bedgraph to output
  out.dir <- paste0(input.parameters$out_dir, file.extension, '_FS', input.parameters$iter)
  create_bedgraphs(to.bg.list, out.dir)
}


#' @title comput log-lik. ratio for the selected idx while accounting for the area of effect
#' @param local.seg.to.guide.lst: list, $guide_idx, $nonGuide_idx
#' @param guide.ll.df: data.frame, null_only_ll, alt_only_ll
#' @param guide.to.seg.lst: list, $seg_overlapped, $dist_to_seg
#' @return vector with idx of region to look for signal
#' @export compute_local_AoE_ll_ratio()

compute_local_AoE_ll_ratio <- function(guide.ll.df, local.seg.to.guide.lst, guide.to.seg.lst){
  
  total.segs <- length(local.seg.to.guide.lst)
  fs.ll <- vector('numeric', length = total.segs)
  
  for(seg in 1:total.segs){
    curr.seg <- seg
    temp.guide.idx<- local.seg.to.guide.lst[[curr.seg]]$guide_idx

    temp.seg.poi.bkg <- unlist(lapply(guide.to.seg.lst[temp.guide.idx], function(x){
      adj.seg.idx <- which(x$seg_overlapped %in% curr.seg )
      temp.dist.to.seg <- x$dist_to_seg[adj.seg.idx]
      dpoibin(0, temp.dist.to.seg)
    }))
    
    temp.guide.ll <- vector('numeric', length = length(temp.seg.poi.bkg))
    
    for(i in 1:length(temp.seg.poi.bkg)){
      temp.guide.ll[i] <- addlogs(guide.ll.df[temp.guide.idx[i],2] + log(1 - temp.seg.poi.bkg[i]),
                                  guide.ll.df[temp.guide.idx[i],1] + log(temp.seg.poi.bkg[i]))
    }
    
    fs.ll[seg] <- sum(temp.guide.ll)

  }
  
  background.ll <- unlist(lapply(local.seg.to.guide.lst, function(x){
    sum(guide.ll.df$null_only_ll[x$guide_idx])
  }))
  
  out.ll.rt <- -2 * (background.ll - fs.ll)
  
  return(out.ll.rt)
  
}


#' @title Wrapper function for recording hyperparameters
#' @param input.results.list: lest, each element keeps track of results from a number of FS identified
#' @param analysis.parameters: list containing all the analysis parameters
#' @param hyper.components: hyperpparameter components (alpha0 and alpha1 proportions and dispersions)
#' @param file.extension: extension for files, either '', '_final', or '_recommendedFinal'
#' @param fs.iter: current FS interation
#' @return csv for all current hyperparameters
#' @export hyperparameter_recording()

hyperparameter_recording <- function(input.results.list, fs.iter, hyper.components, analysis.parameters, file.extension){
  
  record.bkg.hyper <- input.results.list$bkg_hyper#[[fs.iter]]
  record.fs.hyper <- input.results.list$fs_hyper#[[fs.iter]]
  out.dir <- paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName)
  
  if(analysis.parameters$model_dispersion){
    record_hyperparameters(input.bkg.alpha = record.bkg.hyper, 
                           input.fs.alpha = record.fs.hyper, 
                           input.bkg.disp = hyper.components$dispersion, 
                           paste0(out.dir, file.extension), fs.iter, analysis.parameters$pool_names)
  } else {
    record.bkg.disp <- input.results.list$bkg_disp#[[fs.iter]]
    record_hyper_components(input.bkg.alpha = record.bkg.hyper,
                            input.fs.alpha = record.fs.hyper,
                            input.bkg.disp = record.bkg.disp,
                            paste0(out.dir, file.extension), fs.iter, analysis.parameters$pool_names)
  }
  
}


#' @title Compute log likelihoods of observed counts for each sgRNA, given. Formerly 'compute_dirichlet_sgrna_log_like'
#' @param hyper: hyperparameters, $alpha0, alpha1
#' @param data: data, consists of: pool1, pool2... and $n
#' @param region.ll.list: list containing the ll and the indexes of the null alternative and both
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @return data frame: total_guide_ll, alt_only_ll
#' @export estimate_relics_sgrna_log_like_eSize()

estimate_relics_sgrna_log_like_eSize <- function(hyper, data, region.ll.list, guide.efficiency){

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


#' @title wrapper which returns the updated guide efficiency and the corresponding coefficients used to calculate it
#' @param hyper.components: hyperparameters components, 
#' @param param: matrix of all posterior probs
#' @param data: list, each element is a replicate
#' @param guide.seg.idx.lst: list: each element is a guide, contaiing the indexes of the segments it overlaps
#' @param guide.efficiency.scores: matrix, each column containing a different set of scores per guide
#' @param ge.coeff vector, containing the guide efficiency coefficients
#' @return list: $guide_efficiency, $ge_coeff
#' @export recompute_ge_coefficients()

recompute_ge_coefficients <- function(param, hyper.components, data, guide.seg.idx.lst, guide.efficiency.scores, ge.coeff) {
  cumulative.pp <- colSums(param$delta.pp) #apply(param$delta.pp, 2, sum)
  cumulative.pp[cumulative.pp > 1] <- 1

  guide.lls.list <- compute_perGuide_fs_ll(cumulative.pp, guide.seg.idx.lst)

  ge.coeff.param <- ge.coeff
  disp.type <- c()
  if(is.null(hyper.components$dispersion)){
    disp.type <- hyper.components$bkg_dispersion
    model.disp <- FALSE
  } else {
    disp.type <- hyper.components$dispersion
    model.disp <- TRUE
  }

  res <- optim(ge.coeff.param, guide_coeff_ll, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
               data=data, region.ll.list = guide.lls.list,
               bkg.alpha = hyper.components$bkg_alpha, fs.alpha = hyper.components$FS_alpha,
               guide.efficiency.scores = guide.efficiency.scores, input.dispersion = disp.type, model.disp = model.disp)

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
    out.list <- list()
    out.list$guide_efficiency <- guide.efficiency.scores
    out.list$ge_coeff <- ge.coeff
  }

  return(out.list)
}


#' @title optimize the guide efficiency coefficients
#' @param ge.coeff.param: guide efficiency coefficients (beta0, beta1, ...)
#' @param data: data, consists of: $pool1, pool2... and $n
#' @param region.ll.list: list containing the ll and the indexes of the null alternative and both
#' @param alpha0.input: hyper parameters for background
#' @param alpha1.input: hyper parameters for FS
#' @param guide.efficiency.scores: matrix, each column containing a different set of scores per guide
#' @param model.disp: whether the dispersion is modeled
#' @return sum of the -log likelihood across all guides
#' @export guide_coeff_ll()

guide_coeff_ll <- function(ge.coeff.param, data, region.ll.list, bkg.alpha, fs.alpha, 
                           guide.efficiency.scores, input.dispersion, model.disp) {

  guide.efficiency.scores.logit <- apply(guide.efficiency.scores, 2, function(x){
    log(x / (1-x))
  })
  guide.efficiency <- 1 / (1 + exp(-(ge.coeff.param[1] + guide.efficiency.scores.logit %*% ge.coeff.param[2:length(ge.coeff.param)])))

  total.neg.ll <- 0

  for(i in 1:length(data)){
    # hyper <- list(alpha0 = alpha0.input[[i]],
    #               alpha1 = alpha1.input[[i]])
    
    temp.bkg.alpha <- c(1,bkg.alpha[[i]])/ sum(c(1, bkg.alpha[[i]]))
    temp.fs.alpha <- c(1,fs.alpha[[i]])/ sum(c(1, fs.alpha[[i]]))
    
    hyper <- reparameterize_hypers(temp.bkg.alpha, temp.fs.alpha, input.dispersion[[i]], guide.efficiency, model.disp)
    
    temp.neg.ll <- estimate_relics_sgrna_log_like(hyper, data[[i]], region.ll.list, guide.efficiency)

    total.neg.ll <- total.neg.ll + -sum(temp.neg.ll$total_guide_ll)

  }

  total.neg.ll
}


#' @title Based fs.threshold, determine record the location of each FS
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


#' @title Record the hyper parameters for the background and the functional sorting probabilities, dispersions separately
#' @param input.bkg.alpha: list of the background alphas resulting the max log-lik. for that FS
#' @param input.fs.alpha: list of the FS alphas resulting the max log-lik. for that FS
#' @param input.bkg.disp: list of the background dispersion parameters resulting the max log-lik. for that FS
#' @param input.alpha.outDir: directory to which the alphas are to be written
#' @param layer.nr: number of functional sequences recorded
#' @param pool.names: names of the pools given, includes the total column ($n)
#' @return .csv file
#' @export record_hyperparameters()

record_hyperparameters <- function(input.bkg.alpha, input.fs.alpha, input.bkg.disp,
                                    input.alpha.outDir, layer.nr, pool.names, print.out = FALSE){
  
  out.alpha.df <- c()
  
  total.rows <- length(input.bkg.alpha) + length(input.fs.alpha)
  total.cols <- max(c( unlist(lapply(input.bkg.alpha, function(x){length(x) + 1})),
                       unlist(lapply(input.fs.alpha, function(x){length(x) + 1}))))
  
  alpha.names <- c(paste(rep('background', length(input.bkg.alpha)), 'r', c(1:length(input.bkg.alpha)), sep = '_'),
                   paste(rep('FS', length(input.fs.alpha)), 'r', c(1:length(input.fs.alpha)), sep = '_'))
  
  alpha.matrix <- matrix(0, nrow = total.rows, ncol = total.cols)
  
  alpha.list <- list()
  alpha.list$bkg_prop <- list()
  alpha.list$fs_prop <- list()
  
  if(! is.null(pool.names)){
    temp.pool.names <- unique(unlist(pool.names))
    
    for(i in 1:length(input.bkg.alpha)){
      alpha0.scores <- rep(0, length(length(temp.pool.names)))
      alpha1.scores <- rep(0, length(length(temp.pool.names)))
      
      temp.bkg.alpha <- c(1,input.bkg.alpha[[i]]) / sum(c(1, input.bkg.alpha[[i]]))
      temp.fs.alpha <- c(1,input.fs.alpha[[i]]) / sum(c(1, input.fs.alpha[[i]]))
      
      alpha.list$bkg_prop[[i]] <- temp.bkg.alpha
      alpha.list$fs_prop[[i]] <- temp.fs.alpha
      
      alpha0.scores[match(pool.names[[i]], temp.pool.names)] <- round(temp.bkg.alpha, 3)
      alpha1.scores[match(pool.names[[i]], temp.pool.names)] <- round(temp.fs.alpha, 3)
      
      alpha.matrix[i,] <- alpha0.scores
      alpha.matrix[i + length(input.bkg.alpha),] <- alpha1.scores

      
      out.alpha.df <- cbind(alpha.names, alpha.matrix)
      
      colnames(out.alpha.df) <- c('hyperPar_type', temp.pool.names)
      
    }
    
  } else {
    
    for(i in 1:length(input.bkg.alpha)){
      temp.bkg.alpha <- c(1,input.bkg.alpha[[i]]) / sum(c(1, input.bkg.alpha[[i]]))
      alpha.list$bkg_prop[[i]] <- temp.bkg.alpha
      alpha.matrix[i, c(1:length(temp.bkg.alpha))] <- round(temp.bkg.alpha, 3)
    }
    
    for(j in (length(input.bkg.alpha) + 1):total.rows){
      temp.fs.alpha <- c(1,input.fs.alpha[[j - length(input.bkg.alpha)]]) / sum(c(1, input.fs.alpha[[j - length(input.bkg.alpha)]]))
      alpha.list$fs_prop[[j - length(input.bkg.alpha)]] <- temp.fs.alpha
      alpha.matrix[j, c(1:length(temp.fs.alpha))] <- round(temp.fs.alpha, 3)
    }
    
    out.alpha.df <- cbind(alpha.names, alpha.matrix)
    
    colnames(out.alpha.df) <- c('hyperPar_type', paste0('pool', c(1:(ncol(out.alpha.df) - 1))))
  }
  
  write.csv(out.alpha.df, file = paste0(input.alpha.outDir, '_k', layer.nr, '_hyperPars.csv'), row.names = F, quote = F)
  
  disp.df <- do.call(cbind, input.bkg.disp)
  colnames(disp.df) <- paste0('repl_', c(1:length(input.bkg.disp)))
  write.csv(disp.df, file = paste0(input.alpha.outDir, '_k', layer.nr, '_replDispersions.csv'), row.names = F, quote = F)
  
  if(print.out){
    print(out.alpha.df)
    return(alpha.list)
  }
}

#' @title Record the hyper parameters for the background and the functional sorting probabilities
#' @param input.bkg.alpha: list of the background alphas resulting the max log-lik. for that FS
#' @param input.fs.alpha: list of the FS alphas resulting the max log-lik. for that FS
#' @param input.bkg.disp: list of the background dispersion parameters resulting the max log-lik. for that FS
#' @param input.alpha.outDir: directory to which the alphas are to be written
#' @param layer.nr: number of functional sequences recorded
#' @param pool.names: names of the pools given, includes the total column ($n)
#' @return .csv file
#' @export record_hyper_components()

record_hyper_components <- function(input.bkg.alpha, input.fs.alpha, input.bkg.disp,
                                    input.alpha.outDir, layer.nr, pool.names){
  
  out.alpha.df <- c()
  
  total.rows <- length(input.bkg.alpha) + length(input.fs.alpha)
  total.cols <- max(c( unlist(lapply(input.bkg.alpha, function(x){length(x) + 1})),
                       unlist(lapply(input.fs.alpha, function(x){length(x) + 1}))))
  
  alpha.names <- c(paste(rep('background', length(input.bkg.alpha)), 'r', c(1:length(input.bkg.alpha)), sep = '_'),
                   paste(rep('FS', length(input.fs.alpha)), 'r', c(1:length(input.fs.alpha)), sep = '_'))
  
  alpha.matrix <- matrix(0, nrow = total.rows, ncol = total.cols + length(input.bkg.disp[[1]]))
  
  if(! is.null(pool.names)){
    temp.pool.names <- unique(unlist(pool.names))
    
    for(i in 1:length(input.bkg.alpha)){
      alpha0.scores <- rep(0, length(length(temp.pool.names)))
      alpha1.scores <- rep(0, length(length(temp.pool.names)))
      
      temp.bkg.alpha <- c(1,input.bkg.alpha[[i]]) / sum(c(1, input.bkg.alpha[[i]]))
      temp.fs.alpha <- c(1,input.fs.alpha[[i]]) / sum(c(1, input.fs.alpha[[i]]))
      
      alpha0.scores[match(pool.names[[i]], temp.pool.names)] <- round(temp.bkg.alpha, 3)
      alpha1.scores[match(pool.names[[i]], temp.pool.names)] <- round(temp.fs.alpha, 3)
      
      alpha0.scores.wDisp <- c(alpha0.scores, round(input.bkg.disp[[i]], 3))
      alpha1.scores.wDisp <- c(alpha1.scores, round(input.bkg.disp[[i]], 3))
      
      alpha.matrix[i,] <- alpha0.scores.wDisp
      alpha.matrix[i + length(input.bkg.alpha),] <- alpha1.scores.wDisp
      
      out.alpha.df <- cbind(alpha.names, alpha.matrix)
      
      colnames(out.alpha.df) <- c('hyperPar_type', temp.pool.names, paste0(rep('dispersion', length(input.bkg.disp[[1]])), '_', c(1:length(input.bkg.disp[[1]]))) )
      
    }
    
  } else {
    
    for(i in 1:length(input.bkg.alpha)){
      temp.bkg.alpha <- c(1,input.bkg.alpha[[i]]) / sum(c(1, input.bkg.alpha[[i]]))
      alpha.matrix[i, c(1:length(temp.bkg.alpha))] <- round(temp.bkg.alpha, 3)
      alpha.matrix[i, (length(temp.bkg.alpha) + 1):ncol(alpha.matrix)] <- round(input.bkg.disp[[i]], 3)
    }
    
    for(j in (length(input.bkg.alpha) + 1):total.rows){
      temp.fs.alpha <- c(1,input.fs.alpha[[j - length(input.bkg.alpha)]]) / sum(c(1, input.fs.alpha[[j - length(input.bkg.alpha)]]))
      alpha.matrix[j, c(1:length(temp.fs.alpha))] <- round(temp.fs.alpha, 3)
      alpha.matrix[j, (length(temp.fs.alpha) + 1):ncol(alpha.matrix)] <- round(input.bkg.disp[[j - length(input.bkg.alpha)]], 3)
    }
    
    out.alpha.df <- cbind(alpha.names, alpha.matrix)
    
    colnames(out.alpha.df) <- c('hyperPar_type', paste0('pool', c(1:(nrow(out.alpha.df) - 2))), paste0(rep('dispersion', length(input.bkg.disp[[1]])), '_', c(1:length(input.bkg.disp[[1]]))))
  }
  
  write.csv(out.alpha.df, file = paste0(input.alpha.outDir, '_k', layer.nr, '_hyperPars.csv'), row.names = F, quote = F)
  
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
      } else if('l2fc' %in% names(temp.score.df)){
        temp.score.df$formatScores <- temp.score.df$l2fc
      } else if('windowScore' %in% names(temp.score.df)){
        temp.score.df$formatScores <- temp.score.df$windowScore
      } else if('l2fc_avg' %in% names(temp.score.df)){
        temp.score.df$formatScores <- temp.score.df$l2fc_avg
      } else if('l2fc_sw_avg' %in% names(temp.score.df)){
        temp.score.df$formatScores <- temp.score.df$l2fc_sw_avg
      } else {
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


#' @title Set initial hyper parameters, generate the fs posterior matrix
#' @param hyper: hyperparameters
#' @param in.data.list: list, $seg_info, $fs0_idx
#' @param local.max: whether local.max is computed
#' @param recompute.fs0: logical, whether FS0 should be recomputed
#' @return list
#' @export init_relics_param()

init_relics_param <- function(hyper, in.data.list, local.max, recompute.fs0) {
  param <- list()

  # delta posterior probabilities
  param$delta.pp <- matrix(0, nrow=hyper$L+1, ncol= nrow(in.data.list$seg_info)) #nrow(in.data.list$overlapMat))

  # first row is special and contains posterior probs of known regulatory elements
  param$delta.pp[1, in.data.list$fs0_idx] <- 1

  param$ll_rt <- list()

  if(local.max){
    param$ll_rt <- matrix(0, nrow=hyper$L+1, ncol= nrow(in.data.list$seg_info))
  }
  
  if(recompute.fs0){
    param$delta.pp[1,] <- 0
  }

  return(param)
}


#' @title comput log-lik. ratio for the selected idx
#' @param local.seg.to.guide.lst: list, $guide_idx, $nonGuide_idx
#' @param guide.ll.df: data.frame, null_only_ll, alt_only_ll
#' @return vector with idx of region to look for signal
#' @export compute_seg_ll_ratio()

compute_seg_ll_ratio <- function(guide.ll.df, seg.to.guide.lst, guide.dist.to.seg, input.aoe){
  
  background.ll <- unlist(lapply(seg.to.guide.lst, function(x){
    sum(guide.ll.df$null_only_ll[x$guide_idx])
  }))
  
  fs.ll <- c()
  
  if(input.aoe != 'normal'){
    fs.ll <- unlist(lapply(seg.to.guide.lst, function(x){
      sum(guide.ll.df$alt_only_ll[x$guide_idx])
    }))
  } else {
    fs.ll <- rep(0, length(seg.to.guide.lst))
    for(seg in 1:length(seg.to.guide.lst)){
      temp.guide.idx<- seg.to.guide.lst[[seg]]$guide_idx
      temp.seg.poi.bkg <- guide.dist.to.seg[[seg]][[1]]
      temp.guide.ll <- vector('numeric', length = length(temp.seg.poi.bkg))
      
      # issue with sg.ll.list and temp.guide.idx
      for(i in 1:length(temp.seg.poi.bkg)){
        temp.guide.ll[i] <- addlogs(guide.ll.df[temp.guide.idx[i],2] + log(1 - temp.seg.poi.bkg[i]),
                                  guide.ll.df[temp.guide.idx[i],1] + log(temp.seg.poi.bkg[i]))
      }
      fs.ll[seg] <- sum(temp.guide.ll)
    }
  }

  out.ll.rt <- -2 * (background.ll - fs.ll)
  
  return(out.ll.rt)
  
}

compute_seg_ll_ratio_old <- function(guide.ll.df, local.seg.to.guide.lst){

  fs.ll <- unlist(lapply(local.seg.to.guide.lst, function(x){
    sum(guide.ll.df$alt_only_ll[x$guide_idx])
  }))

  background.ll <- unlist(lapply(local.seg.to.guide.lst, function(x){
    sum(guide.ll.df$null_only_ll[x$guide_idx])
  }))

  out.ll.rt <- -2 * (background.ll - fs.ll)

  return(out.ll.rt)

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


#' @title Numerically stable addition of logs but vectorized. from https://leimao.github.io/blog/LogSumExp/
#' @param loga: first log value
#' @param logb: second log value
#' @return log addition
#' @export addlogs_vectorized()

addlogs_vectorized <- function(logs) {
  max.log <- max(logs)
  max.log + log(sum(exp(logs - max.log)))
}
