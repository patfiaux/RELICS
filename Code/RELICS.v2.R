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

  
  # # if no spline has yet been specified, identify the optimal one
  # # To Do 5-fold crossvalidation
  # repl.splines <- c() # spline df for each replicate
  # if(analysis.parameters$estimateSpline){
  #   print('not implemented yet')
  #   # repl.splines <- identify_splines()
  #   break()
  # } else {
  #   repl.spline.df <- analysis.parameters$repl_spline_df
    # repl.disp <- disp_from_spline(repl.spline.df, data.setup$data, 
    #                               paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName), 
    #                               analysis.parameters$nr_disp_bins, 
    #                               analysis.parameters$FS0_label)
    # analysis.parameters$repl_disp <- repl.disp 
  # }
  
  # 
  
  
  # check if hyper parameters are provided, otherwise estimate from data
  if(! 'hyper_pars' %in% names(analysis.parameters)){

    print('Estimating initial hyperparameters ...')
    #labels for background
    background.labels <- c()
    if(analysis.parameters$background_label_specified){
      background.labels <- analysis.parameters$background_label
    } else {
      background.labels <- analysis.parameters$labelHierarchy[-which(analysis.parameters$labelHierarchy %in% analysis.parameters$FS0_label)]
    }
    
    fs0.alphas <-
    if(analysis.parameters$model_dispersion){
      fs0.alphas <- estimate_dirichlet_proportions(data.setup,
                                                   analysis.parameters,
                                                   input.repl.pools = analysis.parameters$repl_groups,
                                                   fs0.label = analysis.parameters$FS0_label)
    } else {
      fs0.alphas <- estimate_hyper_parameters(data.setup,
                                              analysis.parameters,
                                              input.repl.pools = analysis.parameters$repl_groups,
                                              fs0.label = analysis.parameters$FS0_label,
                                              analysis.parameters$one_dispersion)
    }

    analysis.parameters$hyper_pars <- fs0.alphas$hyper_pars
    analysis.parameters$hyper_par_components <- fs0.alphas$hyper_par_components
  }

  # To Do: need to find better way to deal with estimated hyperparameters
  # if(return.init.hypers){
  #   # record alphas used for posterior calculation
  #   alpha.out.dir <- paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName)
  #   # record_alphas(fs0.alphas$alpha0, fs0.alphas$alpha1, alpha.out.dir, 'init', analysis.parameters$pool_names)
  #   record_hyper_components(input.bkg.alpha = fs0.alphas$hyper_par_components$bkg_alpha, 
  #                           input.fs.alpha = fs0.alphas$hyper_par_components$FS_alpha, 
  #                           input.bkg.disp = fs0.alphas$hyper_par_components$bkg_dispersion,
  #                           alpha.out.dir, 'init', analysis.parameters$pool_names)
  #   
  #   break()
  # }
  
  # plot_counts_vs_dispersion(data.setup$data, analysis.parameters$hyper_par_components, 
  #                           paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName), 
  #                           analysis.parameters$nr_disp_bins, analysis.parameters$mean_var_type,
  #                           fs0.label = analysis.parameters$FS0_label, data.setup$guide_info)

  data.setup$fixed_ge_coeff <- analysis.parameters$fixed_ge_coeff

  # To Do: adjust the guide efficiency to work with 
  # # if guide efficiency scores are provided, calculate guide efficiency and include in the model
  # if(! is.null(data.setup$guide_efficiency_scores)){
  #   
  #   if(! analysis.parameters$fixed_ge_coeff){
  # 
  #     temp.relics.params <- init_relics_param(analysis.parameters$hyper_pars, data.setup, analysis.parameters$local_max)
  # 
  #     ge.list <- recompute_ge_coefficients(temp.relics.params,
  #                                          analysis.parameters$hyper_pars,
  #                                          data.setup$data,
  #                                          data.setup$guide_to_seg_lst,
  #                                          data.setup$guide_efficiency_scores,
  #                                          c(0, rep(1, ncol(data.setup$guide_efficiency_scores))))
  # 
  #     data.setup$guide_efficiency <- ge.list$guide_efficiency
  #     data.setup$ge_coeff <- ge.list$ge_coeff
  #   }
  # 
  # }

  # plot the per-segment ll-ratio
  out.pars <- list(out_dir = paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName),
                   iter = 0)
  record_ll_ratio(analysis.parameters$hyper_pars, data.setup, out.pars)

  run_RELICS_2(input.data = data.setup,
               final.layer.nr = analysis.parameters$min_FS_nr,
               out.dir = paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName),
               input.hypers = analysis.parameters$hyper_pars,
               input.hyper.components = analysis.parameters$hyper_par_components,
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
               local.max.range = analysis.parameters$local_max_range,
               pool.names = analysis.parameters$pool_names,
               mean.var.type = analysis.parameters$mean_var_type,
               pp.calculation = analysis.parameters$pp_calculation)


}


#' @title obtain the splines for each replicate to generate the dispersions
#' @param repl.spline.df: number of dgrees of freedom used for the spline, per replicate
#' @param input.data: list: each element is a replicate
#' @param out.dir: directory to which the file is written
#' @param nr.bins: default 20, number of bins to divide the data into for estimating the hyper parameters
#' @param fs0.label: label of the training guides
#' @return sum of the -log likelihood across all guides
#' @export disp_from_spline()

disp_from_spline <- function(repl.spline.df, input.data, out.dir, nr.bins = 20, fs0.label) {
  
  repl.guide.disp <- list()
  
  for(i in 1:length(input.data)){
    # counts from a replicate. The last column has to be the total of all counts
    temp.counts <- input.data[[i]]
    
    # remove the FS0 labeled data.
    temp.counts <- temp.counts[-fs0.label.idx,]
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
    disp.groups <- 100 # number of dispersion groups
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
    
    # extract the mean counts of each bin
    temp.repl.group.means <- unlist(lapply(temp.repl.counts.groups, function(x){
      mean(x[,ncol(x)])
    }))
    
    # now we can plot the estimated dispersions for each bin
    pdf(paste0(out.dir, '_countsVSdispersion_repl', repl, '.pdf'), height = 3, width = 8)
    plot(x = temp.repl.group.means, y = temp.repl.group.disp, xlab = 'mean bin count', ylab = 'bin dispersion', main = 'Bin dispersion estimate')
    
    # the df parameter stands for degrees of freedom
    # it refers to the number of data points that is spaces around the range of all data points to fit the spline
    # the two corner cases (max and min) are already given. And then it fits another df-1 knots to fitting the data
    data.df <- data.frame(disp = temp.repl.group.disp, counts = temp.repl.group.means, stringsAsFactors = F)
    spline.mdl.1 <- lm(disp ~ bs(counts, df = repl.spline.df[[i]]), data = data.df)
    repl.spline.mdls[[i]] <- spline.mdl.1
    
    # Now that you have the spline models established, we have to see how well it actually does
    # The predict model uses the models established above to compute what the expected 
    # dispersion is given the range of all possible counts we have in the data.
    # Because the spline models use the same variable names we first have to set up another data frame with the name 'counts' as the column
    vals.df <- data.frame(counts = temp.total.counts, stringsAsFactors = F)
    spline.mdl.1.predict <- predict(spline.mdl.1, vals.df) 
    
    # With the predictions we can see how well they actually get modeled
    points(x = temp.total.counts, y = spline.mdl.1.predict, col = 'red') 
    dev.off()
    
    # cap the dispersion predictions at the highest mean value
    max.mean.pos <- which(temp.repl.group.means == max(temp.repl.group.means))
    spline.mdl.1.predict.capped <- spline.mdl.1.predict
    spline.mdl.1.predict.capped[which(spline.mdl.1.predict > max(temp.total.counts))] <- temp.repl.group.disp[max.mean.pos]
    
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


#' @title plot the relationship between the total counts and the dispersion
#' @param input.hypers: list: hyper parameters
#' @param input.data: list: each element is a replicate
#' @param out.dir: directory to which the file is written
#' @param mean.var.type: type of count-varaince relationship
#' @param nr.bins: default 20, number of bins to divide the data into for estimating the hyper parameters
#' @param fs0.label: label of the training guides
#' @param guide.info: info file
#' @return plot
#' @export plot_counts_vs_dispersion()

plot_counts_vs_dispersion <- function(input.data, input.hypers, out.dir, nr.bins = 20, mean.var.type,
                                      fs0.label, guide.info){
  
  # identify guides used for training and not part of background calculation, remove them
  fs0.label.idx <- which(guide.info$label %in% fs0.label)
  
  out.est.hypers <- c()
  
  for(repl in 1:length(input.data)){
    
    temp.counts <- input.data[[repl]]
    temp.counts <- temp.counts[-fs0.label.idx,]
    
    temp.total.counts <- temp.counts[,ncol(temp.counts)]
    
    ############
    # compute the hyper parameters with just proportion and one dispersion
    # only do so for background
    temp.bkg.alpha <- rep(1, ncol(temp.counts) - 2)
    temp.bkg.disp <- c(1)

    temp.hyper.params <- c(temp.bkg.alpha, temp.bkg.disp)
    temp.bkg.idx <- c(1:length(temp.bkg.alpha))
    temp.disp.idx <- length(temp.hyper.params)
    
    temp.res.drch <- optim(temp.hyper.params, prior_bkg_dirichlet_parameters, method= 'L-BFGS-B',
                           data = temp.counts, 
                           bkg.idx = temp.bkg.idx, 
                           disp.idx = temp.disp.idx)

    temp.bkg.freq <- c(1, temp.res.drch$par[temp.bkg.idx]**2) / sum(c(1, temp.res.drch$par[temp.bkg.idx]**2))
    temp.bkg.disp.optim <- temp.res.drch$par[temp.disp.idx]^2
    
    out.est.hypers <- rbind(out.est.hypers, c(temp.bkg.freq, temp.bkg.disp.optim))
    
    # temp.bkg.freq <- c(1, input.hypers$bkg_alpha[[repl]]) / sum(c(1,input.hypers$bkg_alpha[[repl]]) )
    # temp.fs.freq <- input.hypers$hyper_par_components$FS_alpha / sum(input.hypers$hyper_par_components$FS_alpha)
    
    temp.bkg.disp <- input.hypers$bkg_dispersion[[repl]]
    temp.disp <- c()
    
    if(mean.var.type == 'radical'){
      temp.disp <- temp.bkg.disp[1] + temp.bkg.disp[2] * temp.total.counts + temp.bkg.disp[3] * sqrt(temp.total.counts)
      temp.disp[which(temp.disp < 0.1)] <- 0.1
    } else if(mean.var.type == 'exponential'){
      temp.disp <- temp.bkg.disp[1] + temp.bkg.disp[2] * log(temp.total.counts)
      # temp.disp <- temp.bkg.disp[1] + temp.bkg.disp[2] * temp.total.counts + temp.bkg.disp[3] * log(temp.total.counts + 1)
      temp.disp[which(temp.disp < 0.1)] <- 0.1
      temp.disp[which(is.na(temp.disp))] <- 0.1
      temp.disp[which(temp.disp == Inf)] <- max(temp.disp[which(temp.disp < Inf)])
    } else if(mean.var.type == 'independent'){
      temp.disp <- rep(temp.bkg.disp[1]^2, length(temp.total.counts))
    } else if(mean.var.type == 'quadratic'){
      temp.disp <- temp.bkg.disp[1] + temp.bkg.disp[2] * temp.total.counts + temp.bkg.disp[3] * temp.total.counts^2
      temp.disp[which(temp.disp < 0.1)] <- 0.1
      temp.disp[which(is.na(temp.disp))] <- 0.1
      temp.disp[which(temp.disp == Inf)] <- max(temp.disp[which(temp.disp < Inf)])
    }
    
    # sort counts and set up bins
    temp.counts.sort <- temp.counts[order(temp.counts[,ncol(temp.counts)], decreasing = T),]
    temp.bg.disp <- 20 # randomly start the bkg dispersion
    
    disp.groups <- nr.bins
    guides.per.grp <- round(nrow(temp.counts.sort)/disp.groups)
    guide.grouping <- rep(c(1:disp.groups), each = guides.per.grp)[1:nrow(temp.counts.sort)]

    # for each bin of guides estimate a dispersion parameter
    temp.repl.counts.groups <- split(temp.counts.sort, guide.grouping)
    
    temp.repl.group.disp <- unlist(lapply(temp.repl.counts.groups, function(x){
      
      temp.grp.disp <- sqrt(temp.bg.disp)
      temp.res <- optim(temp.grp.disp, optim_guide_disp,
                        input.hyper = temp.bkg.freq, 
                        input.data = x, method= 'L-BFGS-B')
      temp.res$par^2
      
    }))
    
    temp.repl.group.means <- unlist(lapply(temp.repl.counts.groups, function(x){
      mean(x[,ncol(x)])
    }))
    
    # compute the fit
    vals <- seq(min(temp.total.counts):max(temp.total.counts))
    if(mean.var.type == 'radical'){
      mdl <- lm(temp.disp ~ temp.total.counts + sqrt(temp.total.counts))
      fit <- mdl$coefficients[1] + mdl$coefficients[2]*vals + mdl$coefficients[3]*sqrt(vals)
    } else if(mean.var.type == 'exponential'){
      mdl <- lm(temp.disp ~ log(temp.total.counts + 1))
      fit <- mdl$coefficients[1] + mdl$coefficients[2]*log(vals + 1)
      
      # mdl <- lm(temp.disp ~ temp.total.counts + log(temp.total.counts))
      # fit <- mdl$coefficients[1] + mdl$coefficients[2] * vals + mdl$coefficients[3]*log(vals)
    } else if(mean.var.type == 'independent'){
      mdl <- lm(temp.disp ~ temp.total.counts)
      fit <- rep(mdl$coefficients[1], length(vals))
    } else if(mean.var.type == 'quadratic'){
      squared.counts <- temp.total.counts^2
      mdl <- lm(temp.disp ~ temp.total.counts + squared.counts)
      fit <- mdl$coefficients[1] + mdl$coefficients[2]*vals + mdl$coefficients[3]*vals^2
    }
    
    fit[fit < 0] <- 0
    
    pdf(paste0(out.dir, '_countsVSdispersion_repl', repl, '.pdf'), height = 3, width = 8)
    par(mfrow = c(1,3))
    plot(x = temp.repl.group.means, y = temp.repl.group.disp, xlab = 'mean bin count', ylab = 'bin dispersion', 
         main = 'Bin dispersion estimate')
    lines(x = vals, y = fit)
    plot(x = log2(temp.repl.group.means + 1), y = log2(temp.repl.group.disp + 1), xlab = 'log2 mean bin count', ylab = 'log2 bin dispersion', 
         main = 'log2 Bin dispersion estimate')
    lines(log2(vals + 1), log2(fit + 1))
    plot(x = log2(temp.total.counts + 1), y = log2(temp.disp + 1), xlab = 'log2 total guide counts', ylab = 'log2 guide dispersion', 
         main = 'Per-guide dispersion estimate')
    dev.off()
    
  }
  
  write.csv(out.est.hypers, file = paste0(out.dir, '_countsVSdispersion_hyperEst.csv'), row.names = F)
  
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


#' @title given hyper parameters and guide efficiency (optional), report the per-segment ll ratio
#' @param input.hypers: list: hyper parameters
#' @param input.data: list: $seg_to_guide_lst, $guide_efficiency_scores (and $guide_efficiency if former not NULL)
#' @param input.ge: guide efficiency
#' @param input.seg.to.guide.lst: list: hyper parameters
#' @return list with parameters for running RELICS or logical FALSE if required parameter is missing
#' @export record_ll_ratio()

record_ll_ratio <- function(input.hypers, input.data, input.parameters){

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
    temp.ll.rt <- compute_local_ll_ratio(guide.model.ll, input.data$seg_to_guide_lst)
    ll.rt <- ll.rt + temp.ll.rt
  }

  # combine with segment info
  segment.info <- input.data$seg_info
  segment.info$score <- ll.rt

  # create bedgraph
  to.bg.list <- list(seg_llRt = segment.info)

  # write bedgraph to output
  out.dir <- paste0(input.parameters$out_dir, '_FS', input.parameters$iter)
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
  
  if(! 'model_dispersion' %in% par.given){
    out.parameter.list$model_dispersion <- TRUE
  }
  if(! 'estimateSpline' %in% par.given){
    out.parameter.list$estimateSpline <- TRUE
  }
  if(! 'nr_disp_bins' %in% par.given){
    out.parameter.list$nr_disp_bins <- 20
  }
  if(! 'pp_calculation' %in% par.given){
    out.parameter.list$pp_calculation <- 'v2' #TRUE
  }
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

    if(! 'ge_beta_estimation' %in% par.given){
      out.parameter.list$ge_beta_estimation <- FALSE
    } else {
      out.parameter.list$ge_beta_estimation <- input.parameter.list$ge_beta_estimation
    }

  } else {
    out.parameter.list$guide_efficiency_scores <- NULL
    out.parameter.list$fixed_ge_coeff <- NULL
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
  
  # type of area of effect. options are 'normal' (default), 'uniform', 'logistic' (need to implement)
  if(! 'areaOfEffect_type' %in% par.given){
    out.parameter.list$areaOfEffect_type <- 'normal'
  }

  if(! 'mean_var_type' %in% par.given){
    out.parameter.list$mean_var_type <- 'radical'
  }
  
  minimum.parameters <- c()
  if(data.file.split){
    minimum.parameters <- c('dataName','repl_groups', 'CountFileLoc', 'sgRNAInfoFileLoc', 'min_FS_nr', 'crisprSystem', 'FS0_label')
  } else {
    minimum.parameters <- c('dataName','repl_groups', 'DataInputFileLoc', 'min_FS_nr', 'crisprSystem', 'FS0_label')
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
      if(! out.parameter.list$crisprSystem %in% c('CRISPRi', 'CRISPRa', 'Cas9', 'CRISPRcas9', 'dualCRISPR')){
          print(paste0('Please provide a valid CRISPR system: CRISPRi, CRISPRa, Cas9 (or CRISPRcas9) or dualCRISPR'))
          missing.parameters <- TRUE
        }
    #   if(! 'crisprEffectRange' %in% par.given){
    #     if(input.parameter.list$crisprSystem == 'CRISPRi'){
    #       out.parameter.list$crisprEffectRange <- 200
    #     } else if(input.parameter.list$crisprSystem == 'CRISPRa'){
    #       out.parameter.list$crisprEffectRange <- 200
    #     } else if(input.parameter.list$crisprSystem %in% c('CRISPRcas9', 'Cas9')){
    #       out.parameter.list$crisprEffectRange <- 20
    #     } else if(input.parameter.list$crisprSystem == 'dualCRISPR'){
    #       out.parameter.list$crisprEffectRange <- 0
    #     } else{
    #       print(paste0('Please provide a valid CRISPR system: CRISPRi, CRISPRa, Cas9 (or CRISPRcas9), or dualCRISPR'))
    #       missing.parameters <- TRUE
    #     }
    #   } else if(! out.parameter.list$crisprSystem %in% c('CRISPRi', 'CRISPRa', 'Cas9', 'CRISPRcas9', 'dualCRISPR')){
    #     print(paste0('Please provide a valid CRISPR system: CRISPRi, CRISPRa, Cas9 (or CRISPRcas9) or dualCRISPR'))
    #     missing.parameters <- TRUE
    #   }
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
  
  if(out.parameter.list$areaOfEffect_type == 'normal'){
    if(! 'normal_areaOfEffect_sd' %in% par.given){
      if(out.parameter.list$crisprSystem %in% c('CRISPRi', 'CRISPRa', 'dualCRISPRi', 'dualCRISPRa') ){
        out.parameter.list$normal_areaOfEffect_sd <- 170
        if(! 'crisprEffectRange' %in% par.given){
          out.parameter.list$crisprEffectRange <- 415
        }
      } else if(out.parameter.list$crisprSystem %in% c('Cas9','CRISPRcas9', 'dualCRISPR') ){
        out.parameter.list$normal_areaOfEffect_sd <- 8.5
        if(! 'crisprEffectRange' %in% par.given){
          out.parameter.list$crisprEffectRange <- 21
        }
      } else {
        print("Error: please specify a valid CRISPR system for the 'uniform' area of effect (Cas9 (or CRISPRcas9), CRISPRi, CRISPRa, dualCRISPR)")
        break()
      }
    }
  }
  
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
    
    guide.to.seg.lst <- sim.guide.seg.list$guide_to_seg_lst

  } else if(input.parameter.list$crisprSystem %in% c('CRISPRi', 'CRISPRa', 'CRISPRcas9', 'Cas9')){
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
    }
    
  } else if(input.parameter.list$crisprSystem == 'dualCRISPR'){
    
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
    }
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
                           guide_to_seg_lst = guide.to.seg.lst, #sim.guide.seg.list$guide_to_seg_lst,
                           seg_to_guide_lst = sim.guide.seg.list$seg_to_guide_lst,
                           next_guide_lst = next.guide.list,
                           guide_efficiency_scores = filtered.ges,
                           guide_info = sim.info)

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
    # colnames(temp.counts) <- paste('pool', c(1:ncol(temp.counts)), sep = '_')
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
  guide.to.seg.list <- generate_guide_to_seg_list(guide.overlap.list) # generate_guide_to_seg_list(guide.overlap.list, input.targeting.info, region.df.filtered)

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

estimate_dirichlet_proportions <- function(data.par.list, analysis.par.list, input.repl.pools,  fs0.label){
  
  # most basic implementation of GE
  if(! is.null(data.par.list$guide_efficiency_scores)){
    if(! analysis.par.list$fixed_ge_coeff){
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
  
  final.dirichlet.pars <- list()
  final.dirichlet.pars$bkg_alpha <- list()
  final.dirichlet.pars$FS_alpha <- list()
  final.dirichlet.pars$dispersion <- list()

  final.alpha <- list()
  final.alpha$alpha0 <- list()
  final.alpha$alpha1 <- list()
  
  for(i in 1:length(data.par.list$data)){

    temp.bkg.alpha <- rep(1, length(input.repl.pools[[i]]) - 1)
    temp.fs.alpha <- rep(1, length(input.repl.pools[[i]]) - 1)

    temp.hyper.params <- c(temp.bkg.alpha, temp.fs.alpha) #, temp.bkg.disp)
    temp.bkg.idx <- c(1:length(temp.bkg.alpha))
    temp.fs.idx <- c(1:length(temp.fs.alpha)) + max(temp.bkg.idx)

    temp.res.drch <- c()
    repl.data <- data.par.list$data[[i]]
    temp.repl.disp <- analysis.par.list[[i]]$repl_disp
    
    temp.res.drch <- optim(temp.hyper.params, prior_dirichlet_proportions, method= 'L-BFGS-B',
                           data = repl.data, 
                           region.ll.list = dirichlet.guide.ll,
                           bkg.idx = temp.bkg.idx, 
                           fs.idx = temp.fs.idx, 
                           guide.efficiency = data.par.list$guide_efficiency, 
                           repl.disp = temp.repl.disp)
    
    temp.bkg.alpha <- c(1, temp.res.drch$par[temp.bkg.idx]**2) / sum(c(1, temp.res.drch$par[temp.bkg.idx]**2))
    temp.fs.alpha <- c(1, temp.res.drch$par[temp.fs.idx]**2) / sum(c(1, temp.res.drch$par[temp.fs.idx]**2))
    
    final.dirichlet.pars$bkg_alpha[[i]] <- temp.res.drch$par[temp.bkg.idx]**2
    final.dirichlet.pars$FS_alpha[[i]] <- temp.res.drch$par[temp.fs.idx]**2
    final.dirichlet.pars$dispersion[[i]] <- temp.repl.disp

    final.alpha$alpha0[[i]] <- t(temp.bkg.alpha %*% t(temp.repl.disp) )
    final.alpha$alpha1[[i]] <- t(temp.fs.alpha %*% t(temp.repl.disp) )

  }
  
  final.alpha$L <- 1

  out.list <- list(hyper_pars = final.alpha,
                   hyper_par_components = final.dirichlet.pars)
  return(out.list)
}


#' @title optimize the dirichlet hyper parameter proportions, dispersion is given
#' @param hyper.param: hyperparameters
#' @param data: data, consists of: $pool1, pool2... and $n
#' @param region.ll.list: list containing the ll and the indexes of the null alternative and both
#' @param bkg.idx: positions in the par vector of the background alphas
#' @param fs.idx: positions in the par vector of the FS alphas
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @param repl.disp: per-guide dispersion based on total counts
#' @return sum of the -log likelihood across all guides
#' @export prior_dirichlet_proportions()

prior_dirichlet_proportions <- function(hyper.param, data, region.ll.list, bkg.idx, fs.idx, 
                                       guide.efficiency, repl.disp) {
  
  bkg.alpha <- c(1, hyper.param[bkg.idx]**2) / sum(c(1, hyper.param[bkg.idx]**2))
  fs.alpha <- c(1, hyper.param[fs.idx]**2) / sum(c(1, hyper.param[fs.idx]**2))

  hyper <- list()
  hyper$alpha0 <- t(bkg.alpha %*% t(repl.disp) )
  hyper$alpha1 <- t(fs.alpha %*% t(repl.disp) )

  out.sg.ll <- estimate_relics_sgrna_log_like(hyper, data, region.ll.list, guide.efficiency)
  
  -sum(out.sg.ll$total_guide_ll)
}



#' @title extract hyper parameters given a label for multiple replicates
#' @param data.par.list: list, contains the elements of the processed and formatted data
#' @param analysis.par.list: list, contains all analysis flags
#' @param input.repl.pools: list, each element is a set of columns which correspond to a replicate
#' @param fs0.label: label used for generating FS0
#' @param one.dispersion: logical, whether analysis is happening with one ro two dispersions
#' @return list of lists: $hyper_pars $hyper_par_components: hyper parameter estimates for each replicate per list element, $bkg_alpha, $FS_alpha, $bkg_dispersion
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

  final.dirichlet.pars <- list()
  final.dirichlet.pars$bkg_alpha <- list()
  final.dirichlet.pars$FS_alpha <- list()
  final.dirichlet.pars$bkg_dispersion <- list()

  final.alpha <- list()
  final.alpha$alpha0 <- list()
  final.alpha$alpha1 <- list()

  for(i in 1:length(data.par.list$data)){
    # temp.drch.hypers <- list(bkg_alpha = rep(1 / length(input.repl.pools[[i]]), length(input.repl.pools[[i]])),
    #                          FS_alpha = rep(1 / length(input.repl.pools[[i]]), length(input.repl.pools[[i]])))
    
    temp.bkg.disp <- c() 
    if(analysis.par.list$mean_var_type == 'radical' | analysis.par.list$mean_var_type == 'quadratic' ){
      temp.bkg.disp <- c(1, 0, 0)
    } else if(analysis.par.list$mean_var_type == 'exponential'){
      temp.bkg.disp <- c(1, 0)
    } else if(analysis.par.list$mean_var_type == 'independent'){
      temp.bkg.disp <- c(1)
    }
    
    temp.bkg.alpha <- rep(1, length(input.repl.pools[[i]]) - 1)
    temp.fs.alpha <- rep(1, length(input.repl.pools[[i]]) - 1)

    
    # temp.hyper.params <- c(sqrt(temp.drch.hypers$alpha0), sqrt(temp.drch.hypers$alpha1))
    # temp.alpha0.idx <- c(1:length(temp.drch.hypers$alpha0))
    # temp.alpha1.idx <- c(1:length(temp.drch.hypers$alpha1)) + max(temp.alpha0.idx)
     
    temp.hyper.params <- c(temp.bkg.alpha, temp.fs.alpha, temp.bkg.disp)
    temp.bkg.idx <- c(1:length(temp.bkg.alpha))
    temp.fs.idx <- c(1:length(temp.fs.alpha)) + max(temp.bkg.idx)
    temp.disp.idx <- c(1:length(temp.bkg.disp)) + max(temp.fs.idx)

    temp.res.drch <- c()
    repl.data <- data.par.list$data[[i]]

    if(one.dispersion){
      temp.res.drch <- optim(temp.hyper.params, prior_dirichlet_parameters, method= 'L-BFGS-B',
                             data = repl.data, 
                             region.ll.list = dirichlet.guide.ll,
                             bkg.idx = temp.bkg.idx, 
                             fs.idx = temp.fs.idx, 
                             disp.idx = temp.disp.idx, 
                             guide.efficiency = data.par.list$guide_efficiency, 
                             mean.var.type = analysis.par.list$mean_var_type)
      
      # temp.res.drch <- optim(temp.hyper.params, prior_dirichlet_ll_singleDisp, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
      #                        data = data.par.list$data[[i]], region.ll.list = dirichlet.guide.ll,
      #                        alpha0.idx = temp.alpha0.idx, alpha1.idx = temp.alpha1.idx,
      #                        guide.efficiency = data.par.list$guide_efficiency)
      
      temp.bkg.alpha <- c(1, temp.res.drch$par[temp.bkg.idx]**2) / sum(c(1, temp.res.drch$par[temp.bkg.idx]**2))
      temp.fs.alpha <- c(1, temp.res.drch$par[temp.fs.idx]**2) / sum(c(1, temp.res.drch$par[temp.fs.idx]**2))
      temp.bkg.disp <- temp.res.drch$par[temp.disp.idx]
      
      final.dirichlet.pars$bkg_alpha[[i]] <- temp.res.drch$par[temp.bkg.idx]**2
      final.dirichlet.pars$FS_alpha[[i]] <- temp.res.drch$par[temp.fs.idx]**2
      final.dirichlet.pars$bkg_dispersion[[i]] <- temp.bkg.disp
      
      if(analysis.par.list$mean_var_type == 'radical'){
        temp.disp <- temp.bkg.disp[1] + temp.bkg.disp[2] * repl.data[,ncol(repl.data)] + temp.bkg.disp[3] * sqrt(repl.data[,ncol(repl.data)])
        temp.disp[which(temp.disp < 0.1)] <- 0.1
        final.alpha$alpha0[[i]] <- t(temp.bkg.alpha %*% t(temp.disp) )
        final.alpha$alpha1[[i]] <- t(temp.fs.alpha %*% t(temp.disp) )
      } else if(analysis.par.list$mean_var_type == 'exponential'){
        temp.disp <- temp.bkg.disp[1] + temp.bkg.disp[2] * log(repl.data[,ncol(repl.data)])
        # temp.disp <- temp.bkg.disp[1] + temp.bkg.disp[2] * data.par.list$data[[i]][,ncol(data.par.list$data[[i]])] + temp.bkg.disp[3] * log(data.par.list$data[[i]][,ncol(data.par.list$data[[i]])] + 1)
        temp.disp[which(temp.disp < 0.1)] <- 0.1
        temp.disp[which(is.na(temp.disp))] <- 0.1
        temp.disp[which(temp.disp == Inf)] <- max(temp.disp[which(temp.disp < Inf)])
        final.alpha$alpha0[[i]] <- t(temp.bkg.alpha %*% t(temp.disp) )
        final.alpha$alpha1[[i]] <- t(temp.fs.alpha %*% t(temp.disp) )
      } else if(analysis.par.list$mean_var_type == 'independent'){
        temp.disp <- rep(temp.bkg.disp[1]^2, nrow(repl.data))
        final.alpha$alpha0[[i]] <- t(temp.bkg.alpha %*% t(temp.disp) )
        final.alpha$alpha1[[i]] <- t(temp.fs.alpha %*% t(temp.disp) )
      } else if(analysis.par.list$mean_var_type == 'quadratic'){
        temp.disp <- temp.bkg.disp[1] + temp.bkg.disp[2] * repl.data[,ncol(repl.data)] + temp.bkg.disp[3] * repl.data[,ncol(repl.data)]^2
        temp.disp[which(temp.disp < 0.1)] <- 0.1
        temp.disp[which(is.na(temp.disp))] <- 0.1
        temp.disp[which(temp.disp == Inf)] <- max(temp.disp[which(temp.disp < Inf)])
        final.alpha$alpha0[[i]] <- t(temp.bkg.alpha %*% t(temp.disp) )
        final.alpha$alpha1[[i]] <- t(temp.fs.alpha %*% t(temp.disp) )
      }

    } else {
      print('two-dispersion option currently disabled!')
      break()
      # temp.res.drch <- optim(temp.hyper.params, prior_dirichlet_ll, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
      #                        data = data.par.list$data[[i]], region.ll.list = dirichlet.guide.ll,
      #                        alpha0.idx = temp.alpha0.idx, alpha1.idx = temp.alpha1.idx,
      #                        guide.efficiency = data.par.list$guide_efficiency)
      # 
      # temp.alpha0s <- temp.res.drch$par[temp.alpha0.idx]**2
      # temp.alpha1s <- temp.res.drch$par[temp.alpha1.idx]**2
      # 
      # final.alpha$alpha0[[i]] <- temp.alpha0s
      # final.alpha$alpha1[[i]] <- temp.alpha1s

    }

  }

  final.alpha$L <- 1
  # return(final.alpha)

  out.list <- list(hyper_pars = final.alpha,
                   hyper_par_components = final.dirichlet.pars)
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
#' @return sum of the -log likelihood across all guides
#' @export prior_dirichlet_parameters()

prior_dirichlet_parameters <- function(hyper.param, data, region.ll.list, bkg.idx, fs.idx, disp.idx, 
                                       guide.efficiency, mean.var.type) {
  
  bkg.alpha <- c(1, hyper.param[bkg.idx]**2) / sum(c(1, hyper.param[bkg.idx]**2))
  fs.alpha <- c(1, hyper.param[fs.idx]**2) / sum(c(1, hyper.param[fs.idx]**2))
  disp.alpha <- hyper.param[disp.idx]

  hyper <- list()
  
  if(mean.var.type == 'radical'){
    temp.disp <- disp.alpha[1] + disp.alpha[2] * data[,ncol(data)] + disp.alpha[3] * sqrt(data[,ncol(data)])
    temp.disp[which(temp.disp < 0.1)] <- 0.1
    hyper$alpha0 <- t(bkg.alpha %*% t(temp.disp) )
    hyper$alpha1 <- t(fs.alpha %*% t(temp.disp) )
  } else if(mean.var.type == 'exponential'){
    temp.disp <- disp.alpha[1] + disp.alpha[2] * log(data[,ncol(data)])
    # temp.disp <- disp.alpha[1] + disp.alpha[2] * data[,ncol(data)] + disp.alpha[3] * log(data[,ncol(data)] + 1)
    temp.disp[which(temp.disp < 0.1)] <- 0.1
    temp.disp[which(is.na(temp.disp))] <- 0.1
    temp.disp[which(temp.disp == Inf)] <- max(temp.disp[which(temp.disp < Inf)])
    hyper$alpha0 <- t(bkg.alpha %*% t(temp.disp) )
    hyper$alpha1 <- t(fs.alpha %*% t(temp.disp) )
  } else if(mean.var.type == 'independent'){
    temp.disp <- disp.alpha[1]^2
    hyper$alpha0 <- bkg.alpha * temp.disp
    hyper$alpha1 <- fs.alpha * temp.disp
  } else if(mean.var.type == 'quadratic'){
    temp.disp <- disp.alpha[1] + disp.alpha[2] * data[,ncol(data)] + disp.alpha[3] * data[,ncol(data)]^2
    temp.disp[which(temp.disp < 0.1)] <- 0.1
    temp.disp[which(is.na(temp.disp))] <- 0.1
    temp.disp[which(temp.disp == Inf)] <- max(temp.disp[which(temp.disp < Inf)])
    hyper$alpha0 <- t(bkg.alpha %*% t(temp.disp) )
    hyper$alpha1 <- t(fs.alpha %*% t(temp.disp) )
  }

  out.sg.ll <- estimate_relics_sgrna_log_like(hyper, data, region.ll.list, guide.efficiency)

  -sum(out.sg.ll$total_guide_ll)
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
#' @param return.model.ll: logical, whether the function should terminate early and return only the null and alternative model ll
#' @return data frame: total_guide_ll, alt_only_ll
#' @export estimate_relics_sgrna_log_like()

estimate_relics_sgrna_log_like <- function(hyper, data, region.ll.list, guide.efficiency, return.model.ll = FALSE){

  pool.cols <- c(1:(ncol(data) - 1))

  sgrna.null.log.like <- ddirmnom(data[,pool.cols], size = data[,ncol(data)], alpha = hyper$alpha0, log = T)

  sgrna.alt.log.like <- c()
  if(is.null(guide.efficiency)){
    sgrna.alt.log.like <- ddirmnom(data[,pool.cols], size = data[,ncol(data)], alpha = hyper$alpha1, log = T)
  } else {
    # toDo, difference in alphas or in proportions? 
    alpha0.matrix <- t(apply(matrix(0, ncol = length(hyper$alpha0), nrow = length(guide.efficiency)), 1, function(x){x + hyper$alpha0}))
    alpha.diffs <- hyper$alpha0 - hyper$alpha1

    alpha1.matrix <- matrix(0, ncol = length(hyper$alpha0), nrow = length(guide.efficiency))

    for(i in 1:length(guide.efficiency)){
      # this cound be an issue if the dispersion is no longer identical...
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
                     pp.calculation){

  final.layer.posterior <- list()
  final.layer.alpha0 <- list()
  final.layer.alpha1 <- list()
  
  hyper.components <- input.hyper.components
  final.bkg.hyper <- list()
  final.fs.hyper <- list()
  final.bkg.disp <- list()
  
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
    browser()
    print(paste0('Computing FS: ', i))
    layer.time <- proc.time()
    layer.data <- relics_compute_FS_k(input.param = relics.param,
                                      input.hyper = relics.hyper,
                                      input.data.list = input.data,
                                      input.tol = coverg.tol[i],
                                      fix.hypers, iterative.hyper.est, nr.segs, geom.p,
                                      min.pp = input.min.rs.pp, input.data$guide_efficiency,
                                      one.dispersion,
                                      local.max, local.max.range,  pp.calculation)

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
    
    final.bkg.hyper[[i]] <- hyper.components$bkg_alpha
    final.fs.hyper[[i]] <- hyper.components$FS_alpha
    final.bkg.disp[[i]] <- hyper.components$bkg_dispersion
    
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

      dirichlet.guide.ll <- compute_perGuide_fs_ll(order.pps.lst$pp_ordered, input.data$guide_to_seg_lst)
      total.per.guide.ll <- 0
      for(repl in 1:length(input.data$data)){
        temp.hypers <- list(alpha0 = relics.hyper$alpha0[[repl]], alpha1 = relics.hyper$alpha1[[repl]])
        temp.sgRNa.ll <- estimate_relics_sgrna_log_like(temp.hypers,
                                                        input.data$data[[repl]],
                                                        dirichlet.guide.ll,
                                                        input.data$guide_efficiency, return.model.ll = TRUE)

        total.per.guide.ll <- total.per.guide.ll + (-2 * (temp.sgRNa.ll$null_only_ll - temp.sgRNa.ll$alt_only_ll) )
      }

      data.info.names <- colnames(input.data$guide_info)
      out.guide.info <- cbind(input.data$guide_info, total.per.guide.ll)
      colnames(out.guide.info) <- c(data.info.names, 'guide_ll')

      write.csv(out.guide.info, file = paste0(out.dir, '_k',i,'_guide_ll_info.csv'), row.names = F)

      out.guide.info$score <- out.guide.info$guide_ll
      to.bg.list$guide_ll <- out.guide.info
      
      # total_effSize, repl_effSize
      abs.sum.effect.size <- record_sum_effectSizes(order.pps.lst$pp_ordered,
                                                      input.min.rs.pp,
                                                      relics.hyper,
                                                      input.data$data,
                                                      input.data$guide_to_seg_lst,
                                                      input.data$seg_to_guide_lst,
                                                      input.data$guide_efficiency,
                                                      one.dispersion,
                                                      input.data$seg_info,
                                                   abs.es = TRUE, pool.names = pool.names,
                                                   mean.var.type, hyper.components)
      
      sum.effect.size <- record_sum_effectSizes(order.pps.lst$pp_ordered,
                                                input.min.rs.pp,
                                                relics.hyper,
                                                input.data$data,
                                                input.data$guide_to_seg_lst,
                                                input.data$seg_to_guide_lst,
                                                input.data$guide_efficiency,
                                                one.dispersion,
                                                input.data$seg_info,
                                                abs.es = FALSE, pool.names = pool.names,
                                                mean.var.type, hyper.components)

      to.bg.list$abs_sumEffSize <- abs.sum.effect.size$total_effSize
      to.bg.list$sumEffSize <- sum.effect.size$total_effSize

      write.csv(abs.sum.effect.size$fs_effSize_ID, file = paste0(out.dir, '_k', i, '_abs_sumEffSize.csv'), row.names = F)
      write.csv(sum.effect.size$fs_effSize_ID, file = paste0(out.dir, '_k', i, '_sumEffSize.csv'), row.names = F)

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
      # record_alphas(final.layer.alpha0[[i]], final.layer.alpha1[[i]], out.dir, i, pool.names)
      # record_hyper_components(input.bkg.alpha = final.bkg.hyper[[i]], 
      #                         input.fs.alpha = final.fs.hyper[[i]], 
      #                         input.bkg.disp = final.bkg.disp[[i]],
      #                         out.dir, i, pool.names)
      record_hyperparameters(input.bkg.alpha = final.bkg.hyper[[i]], 
                             input.fs.alpha = final.fs.hyper[[i]], 
                             input.bkg.disp = input.hypers$hyper_par_components$dispersion, #final.bkg.disp[[i]],
                             out.dir, i, pool.names)

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
                               input.min.rs.pp,
                               input.data$seg_info)

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

        final.alpha0s <- final.layer.alpha0[[i - 1]]
        final.alpha1s <- final.layer.alpha1[[i - 1]]
        final.hypers <- list(alpha0 = final.alpha0s, alpha1 = final.alpha1s)

        abs.sum.effect.size <- record_sum_effectSizes(final.pp.out,
                                                     input.min.rs.pp,
                                                     final.hypers,
                                                     input.data$data,
                                                     input.data$guide_to_seg_lst,
                                                     input.data$seg_to_guide_lst,
                                                     input.data$guide_efficiency,
                                                     one.dispersion,
                                                     input.data$seg_info,
                                                     abs.es = TRUE, pool.names = pool.names,
                                                     mean.var.type, hyper.components)
        
        sum.effect.size <- record_sum_effectSizes(final.pp.out,
                                                      input.min.rs.pp,
                                                      final.hypers,
                                                      input.data$data,
                                                      input.data$guide_to_seg_lst,
                                                      input.data$seg_to_guide_lst,
                                                      input.data$guide_efficiency,
                                                      one.dispersion,
                                                      input.data$seg_info,
                                                      abs.es = FALSE, pool.names = pool.names,
                                                  mean.var.type, hyper.components)

        to.bg.list$abs_sumEffSize <- abs.sum.effect.size$total_effSize
        to.bg.list$sumEffSize <- sum.effect.size$total_effSize

        write.csv(abs.sum.effect.size$fs_effSize_ID, file = paste0(out.dir, '_final_k', i - 1, '_abs_sumEffSize.csv'), row.names = F)
        write.csv(sum.effect.size$fs_effSize_ID, file = paste0(out.dir, '_final_k', i - 1, '_sumEffSize.csv'), row.names = F)

         dirichlet.guide.ll <- compute_perGuide_fs_ll(final.pp.out, input.data$guide_to_seg_lst)
         total.per.guide.ll <- 0
         for(repl in 1:length(input.data$data)){
           temp.hypers <- list(alpha0 = final.hypers$alpha0[[repl]], alpha1 = final.hypers$alpha1[[repl]])
           temp.sgRNa.ll <- estimate_relics_sgrna_log_like(temp.hypers,
                                                           input.data$data[[repl]],
                                                           dirichlet.guide.ll,
                                                           input.data$guide_efficiency, return.model.ll = TRUE)

           total.per.guide.ll <- total.per.guide.ll + (-2 * (temp.sgRNa.ll$null_only_ll - temp.sgRNa.ll$alt_only_ll) )
         }

         data.info.names <- colnames(input.data$guide_info)
         out.guide.info <- cbind(input.data$guide_info, total.per.guide.ll)
         colnames(out.guide.info) <- c(data.info.names, 'guide_ll')

         out.guide.info$score <- out.guide.info$guide_ll
         to.bg.list$guide_ll <- out.guide.info

        if(auto.stop){

          write.csv(out.guide.info, file = paste0(out.dir, '_final_k', i - 1,'_guide_ll_info.csv'), row.names = F)

          # record bedgraph
          create_bedgraphs(to.bg.list, paste0(out.dir, '_final_k', i - 1) )

          # save FS locations
          write.table(all.seg.fs.df.final, file = paste0(out.dir, '_final_k', i - 1,'_FS_locations.bed'),
                      sep = '\t', quote = F, row.names = F, col.names = F)

          # record alphas used for posterior calculation
          # record_alphas(final.layer.alpha0[[i - 1]], final.layer.alpha1[[i - 1]], paste0(out.dir, '_final'), i - 1, pool.names)
          # record_hyper_components(input.bkg.alpha = final.bkg.hyper[[i - 1]], 
          #                         input.fs.alpha = final.fs.hyper[[i - 1]], 
          #                         input.bkg.disp = final.bkg.disp[[i - 1]],
          #                         paste0(out.dir, '_final'), i - 1, pool.names)
          record_hyperparameters(input.bkg.alpha = final.bkg.hyper[[i - 1]], 
                                 input.fs.alpha = final.fs.hyper[[i - 1]], 
                                 input.bkg.disp = input.hypers$hyper_par_components$dispersion, #final.bkg.disp[[i - 1]],
                                 paste0(out.dir, '_final'), i - 1, pool.names)

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
                                    input.min.rs.pp,
                                    input.data$seg_info)

          display_relics_fs_as_tiff(final.pp.out,
                                    input.data$segLabels,
                                    paste0(out.dir, '_final_k', i-1),
                                    input.min.rs.pp,
                                    input.data$seg_info)

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

          write.csv(out.guide.info, file = paste0(out.dir, '_recommendedFinal_k', i - 1,'_guide_ll_info.csv'), row.names = F)

          # record bedgraph
          create_bedgraphs(to.bg.list, paste0(out.dir, '_recommendedFinal_k', i - 1) )

          # save FS locations
          write.table(all.seg.fs.df.final, file = paste0(out.dir, '_recommendedFinal_k', i - 1,'_FS_locations.bed'),
                      sep = '\t', quote = F, row.names = F, col.names = F)

          # record alphas used for posterior calculation
          # record_alphas(final.layer.alpha0[[i - 1]], final.layer.alpha1[[i - 1]], paste0(out.dir, '_recommendedFinal'), i - 1, pool.names)
          # record_hyper_components(input.bkg.alpha = final.bkg.hyper[[i - 1]], 
          #                         input.fs.alpha = final.fs.hyper[[i - 1]], 
          #                         input.bkg.disp = final.bkg.disp[[i - 1]],
          #                         paste0(out.dir, '_recommendedFinal'), i - 1, pool.names)
          record_hyperparameters(input.bkg.alpha = final.bkg.hyper[[i - 1]], 
                                 input.fs.alpha = final.fs.hyper[[i - 1]], 
                                 input.bkg.disp = input.hypers$hyper_par_components$dispersion, #final.bkg.disp[[i - 1]],
                                 paste0(out.dir, '_recommendedFinal'), i - 1, pool.names)

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
                                    input.min.rs.pp,
                                    input.data$seg_info)

          display_relics_fs_as_tiff(final.pp.out,
                                    input.data$segLabels,
                                    paste0(out.dir, '_recommendedFinal_k', i-1),
                                    input.min.rs.pp,
                                    input.data$seg_info)

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
      relics.hyper.list <- recompute_hyper_parameters(relics.param,
                                                 relics.hyper,
                                                 hyper.components,
                                                 input.data$data,
                                                 input.data$guide_to_seg_lst,
                                                 input.data$guide_efficiency,
                                                 one.dispersion,
                                                 mean.var.type)
      relics.hyper <- relics.hyper.list$hyper_pars
      hyper.components <- relics.hyper.list$hyper_par_components
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

    if(!fix.hypers | !is.null(input.data$guide_efficiency_scores)){
      #out.pars <- list(out_dir = out.dir, dataName = '')
      out.pars <- list(out_dir = out.dir,
                       iter = i)
      record_ll_ratio(relics.hyper, input.data, out.pars)
    }

  }

}


#' @title records the sum of (absolute) effect sizes across all pools (avg across replicates) for each FS
#' @param hyper: hyperparameters
#' @param input.pp: matrix of all posterior probs
#' @param input.min.rs.pp: threshold for FS detection
#' @param data: list, each element is a replicate
#' @param guide.seg.idx.lst: list: each element is a guide, containing the indexes of the segments it overlaps
#' @param seg.guide.idx.lst: list: each element is a segment and all guide that overlap it
#' @param guide.efficiency: either a vector of guide efficiency or NULL
#' @param one.dispersion. logical, if there should be one or two dispersions for the hyper parameters
#' @param seg.info: data frame, coordinates for each segment
#' @param abs.es: logical, whether or not the absolute difference should be returned
#' @param pool.names: if given, the names of the individual pools
#' @param mean.var.type: type of mean-variance relationship
#' @param hyper.components: hyperparameter estimates of $bkg_alpha, $FS_alpha, $bkg_dispersion
#' @return list: total_effSize, repl_effSize, fs_effSize_ID
#' @export record_sum_effectSizes()

record_sum_effectSizes <- function(input.pp, input.min.rs.pp, hyper, data,
                                  guide.seg.idx.lst, seg.guide.idx.lst,
                                  guide.efficiency, one.dispersion, seg.info,
                                  abs.es = FALSE, pool.names, 
                                  mean.var.type, hyper.components) {

  # list of all alpha 1s (one for each FS)
  fs.alpha1.eff.size <- list()
  fs.total.eff.size <- c()
  fs.eff.size.ID <- c()

  for(e in 1:nrow(input.pp)){

    temp.fs.idx <- which(input.pp[e,] > input.min.rs.pp)

    fs.alpha1.eff.size[[e]] <- list()

    if(length(temp.fs.idx) > 0){
      temp.pp <- input.pp[e,]

      temp.seg.info <- seg.info[temp.fs.idx,]

      # 1. extract the index of guides overlapping all the FS segments
      temp.fs.guide.idx <- sort(unique(unlist(lapply(seg.guide.idx.lst[temp.fs.idx], function(x){
        x$guide_idx
      }))))

      did.converge <- TRUE

      # for each replicate
      for(i in 1:length(data)){

        # extract the subset of guides overlapping the FS
        fs.data <- data[[i]][temp.fs.guide.idx,]
        fs.guide.efficiency <- guide.efficiency[temp.fs.guide.idx,]

        # 2. use fs.data to get overlapping pp and calculate the per-guide ll given all the FS overlapping that guide
        # Issue here is that I actually include more pp than the FS is long...
        temp.fs.guide.pp.idx <- sort(unique(unlist(lapply(guide.seg.idx.lst[temp.fs.guide.idx], function(x){
          x$seg_overlapped
        }))))

        temp.guide.lls.list <- FS_only_guide_ll(temp.pp, guide.seg.idx.lst, temp.fs.guide.pp.idx, temp.fs.guide.idx)
        
        temp.bkg.alpha <- c(1, hyper.components$bkg_alpha[[i]]) / sum(c(1, hyper.components$bkg_alpha[[i]]))
        temp.fs.alpha <- sqrt(hyper.components$FS_alpha[[i]])
        
        temp.disp <- hyper.components$dispersion[[i]]
        temp.alpha0 <- t(temp.bkg.alpha %*% t(temp.disp) )
        
        # temp.bkg.disp <- hyper.components$bkg_dispersion[[i]]
        # 
        # temp.disp <- c()
        # temp.alpha0 <- c()
        # 
        # if(mean.var.type == 'radical'){
        #   temp.disp <- (temp.bkg.disp[1] + temp.bkg.disp[2] * fs.data[,ncol(fs.data)] + temp.bkg.disp[3] * sqrt(fs.data[,ncol(fs.data)]))
        #   temp.disp[which(temp.disp < 0.1)] <- 0.1
        #   temp.alpha0 <- t(temp.bkg.alpha %*% t(temp.disp) )
        # } else if(mean.var.type == 'exponential'){
        #   temp.disp <- (temp.bkg.disp[1] + temp.bkg.disp[2] * log(fs.data[,ncol(fs.data)]))
        #   # temp.disp <- (temp.bkg.disp[1] + temp.bkg.disp[2] * fs.data[,ncol(fs.data)] + temp.bkg.disp[3] * log(fs.data[,ncol(fs.data)] + 1))
        #   temp.disp[which(temp.disp < 0.1)] <- 0.1
        #   temp.disp[which(is.na(temp.disp))] <- 0.1
        #   temp.disp[which(temp.disp == Inf)] <- max(temp.disp[which(temp.disp < Inf)])
        #   temp.alpha0 <- t(temp.bkg.alpha %*% t(temp.disp) )
        # } else if(mean.var.type == 'independent'){
        #   temp.disp <- temp.bkg.disp[1]^2
        #   temp.alpha0 <- temp.bkg.alpha * temp.disp
        # } else if(mean.var.type == 'quadratic'){
        #   temp.disp <- temp.bkg.disp[1] + temp.bkg.disp[2] * fs.data[,ncol(fs.data)] + temp.bkg.disp[3] * fs.data[,ncol(fs.data)]^2
        #   temp.disp[which(temp.disp < 0.1)] <- 0.1
        #   temp.disp[which(is.na(temp.disp))] <- 0.1
        #   temp.disp[which(temp.disp == Inf)] <- max(temp.disp[which(temp.disp < Inf)])
        #   temp.alpha0 <- t(temp.bkg.alpha %*% t(temp.disp) )
        # }
        # 
        # temp.alpha0 <- hyper$alpha0[[i]]
        # temp.alpha1 <- sqrt(hyper$alpha1[[i]])

        res <- c()

        if(one.dispersion){
          # res <- optim(temp.alpha1, prior_dirichlet_ll_singleDisp_eSize, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
          #              data = fs.data, region.ll.list = temp.guide.lls.list,
          #              guide.efficiency = fs.guide.efficiency, alpha.zero = temp.alpha0)
          
          res <- optim(temp.fs.alpha, FS_prop_eSize, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
                       alpha.zero = temp.alpha0, input.disp = temp.disp,
                       data = fs.data, region.ll.list = temp.guide.lls.list,
                       guide.efficiency = fs.guide.efficiency)
        } else {
          print('Two dispersions currently not possible!')
          break()
          # res <- optim(temp.alpha1, prior_dirichlet_ll_eSize, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
          #              data = fs.data, region.ll.list = temp.guide.lls.list,
          #              guide.efficiency = fs.guide.efficiency, alpha.zero = temp.alpha0)
        }

        alpha1s <- res$par**2

        if(one.dispersion){

          # alpha1s.norm <- alpha1s / sum(alpha1s)
          # alpha1s.adj <- alpha1s.norm * sum(temp.alpha0)
          
          alpha1s.adj <- c(1, alpha1s) / sum(c(1, alpha1s))

          if(abs.es){
            # fs.alpha1.eff.size[[e]][[i]] <- abs((alpha1s.adj / sum(alpha1s.adj)) - (temp.alpha0 / sum(temp.alpha0 )))
            fs.alpha1.eff.size[[e]][[i]] <- abs(alpha1s.adj - temp.bkg.alpha)
          } else {
            # fs.alpha1.eff.size[[e]][[i]] <- (alpha1s.adj / sum(alpha1s.adj)) - (temp.alpha0 / sum(temp.alpha0 ))
            fs.alpha1.eff.size[[e]][[i]] <- alpha1s.adj - temp.bkg.alpha
          }
        } else {
          print('Two dispersions currently not possible!')
          break()
          # if(abs.es){
          #   fs.alpha1.eff.size[[e]][[i]] <- abs( (alpha1s / sum(alpha1s)) - (temp.alpha0 / sum(temp.alpha0)) )
          # } else {
          #   fs.alpha1.eff.size[[e]][[i]] <- (alpha1s / sum(alpha1s)) - (temp.alpha0 / sum(temp.alpha0))
          # }
        }

          if(res$convergence != 0) {
            warning("estimation of hyperparameters failed to converge")
            did.converge <- FALSE
          }

      }

      # at this stage, all FS alpha1s should be obtained, present in list of lists at fs.alpha1[[f]]
      # compute the average effect size across all replicates
      
      temp.total.avg.eff <- c()
      # if replicates don't have the same number of pools and no pool names for combining has been arranged
      if(! is.null(pool.names)){
        all.pool.names <- unique(unlist(pool.names))
        temp.eff.list <- list()
        
        # for all functional sequences found
        for(fs in 1:length(fs.alpha1.eff.size)){
          temp.repl.alphas <- fs.alpha1.eff.size[[fs]]
          
          if(length(temp.repl.alphas) > 0){
            # for all replicates
            for(r in 1:length(temp.repl.alphas)){
              # for each pool / alpha
              for(a in 1:length(temp.repl.alphas[[r]]))
                temp.eff.list[[pool.names[[r]][a]]] <- c(temp.eff.list[[pool.names[[r]][a]]], temp.repl.alphas[[r]][a])
            }
          }
        }
        
        temp.total.eff <- unlist(lapply(temp.eff.list, function(x){mean(x)}))
        temp.total.avg.eff <- sum(temp.total.eff)
        
      } else {
        if(length(unique(unlist(lapply(data, function(x){nrow(x)})))) != 1){
          print('Warning! Not sure how to combine per-pool effect sizes!')
        }
        temp.eff.df <- do.call(rbind, fs.alpha1.eff.size[[e]])
        temp.total.avg.eff <- sum(colSums(temp.eff.df) / nrow(temp.eff.df))
      }
      
      temp.eff.size.loc <- temp.seg.info[, c('chrom', 'start', 'end')]
      temp.eff.size.loc$FS <- paste0('FS', e-1)
      temp.eff.size.loc$score <- temp.total.avg.eff # preparing for bg format
      temp.eff.size.loc$converged <- as.character(did.converge)

      fs.eff.size.ID <- rbind(fs.eff.size.ID, data.frame(FS = paste0('FS', e-1), FS_effSize = temp.total.avg.eff, stringsAsFactors = F))

      fs.total.eff.size <- rbind(fs.total.eff.size, temp.eff.size.loc)
    }

  }

  return(list(total_effSize = fs.total.eff.size, repl_effSize = fs.alpha1.eff.size, fs_effSize_ID= fs.eff.size.ID))
}


#' @title optimize the proportion parameters for a functional sequence, only one dispersion across the two distributions
#' @param data: data, consists of: $pool1, pool2... and $n
#' @param region.ll.list: list containing the ll and the indexes of the null alternative and both
#' @param fs.prop: proportion parameters for functional sequence
#' @param alpha.zero: hyper parameters for background
#' @param input.disp: dispersion
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @return sum of the -log likelihood across all guides
#' @export FS_prop_eSize()

FS_prop_eSize <- function(fs.alpha, alpha.zero, input.disp, data, region.ll.list, guide.efficiency) {
  
  adj.alpha <- fs.alpha**2
  adj.alpha[which(adj.alpha == 0)] <- 0.001
  
  fs.alpha.adj <- c(1, adj.alpha) / sum(c(1, adj.alpha))
  alpha.one <- t(fs.alpha.adj %*% t(input.disp) )

  hyper <- list(alpha0 = alpha.zero,
                alpha1 = alpha.one)
  
  out.sg.ll <- estimate_relics_sgrna_log_like_eSize(hyper, data, region.ll.list, guide.efficiency)
  
  -sum(out.sg.ll$total_guide_ll)
}


#' @title Compute log likelihoods of each sgRNA that overlaps a specifc FS.
#' @param cumulative.pp: sum of all previous posteriors
#' @param guide.seg.idx.lst: list, each element is a guide and contains the indexes of the segments it overlaps
#' @param fs.idx: vector of segment indeces of the functional sequence
#' @param guide.fs.idx: index of the guides overlapping the FS
#' @return list: null_only_idx, null_only_lls, alt_only_idx, alt_only_lls, both_idx, both_lls (last one is a df)
#' @export FS_only_guide_ll()

FS_only_guide_ll <- function(cumulative.pp, guide.seg.idx.lst, fs.idx, guide.fs.idx, hyper.setup = FALSE){

  # fs.pp <- cumulative.pp[fs.idx]

  nr.fs.sgrna <- length(guide.fs.idx)

  null.ll <- vector('numeric', length = nr.fs.sgrna)
  alt.ll <- vector('numeric', length = nr.fs.sgrna)
  null.indicator <- vector('numeric', length = nr.fs.sgrna)
  alt.indicator <- vector('numeric', length = nr.fs.sgrna)

  null.only.idx <- vector('numeric', length = nr.fs.sgrna)
  null.only.counter <- 1
  alt.only.idx <- vector('numeric', length = nr.fs.sgrna)
  alt.only.counter <- 1
  both.idx <- vector('numeric', length = nr.fs.sgrna)
  both.counter <- 1

  for(j in 1:nr.fs.sgrna){
    
    region.pp <- cumulative.pp[guide.seg.idx.lst[[guide.fs.idx[j] ]]$seg_overlapped ]
    
    region.pp.distance.adj <- region.pp * guide.seg.idx.lst[[guide.fs.idx[j] ]]$dist_to_seg

    # compute probability that number of regulatory regions overlapped by this sgRNA is 0 or >0, given posterior
    # probabilities, using poisson binomial probability mass function
    p.k.eq.0 <- dpoibin(0, region.pp.distance.adj)
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

  out.list <- list(null_only_idx = null.idxs,
                   null_only_lls = null.ll[null.idxs],
                   alt_only_idx = alt.idxs,
                   alt_only_lls = alt.ll[alt.idxs],
                   both_idx = both.idxs,
                   both_lls = data.frame(null_lls = null.ll[both.idxs], alt_lls = alt.ll[both.idxs]))

  return(out.list)

}


#' @title optimize the hyper parameters, only one dispersion across the two distributions
#' @param scaling.par: scaling parameter for alpha 1 hypers
#' @param data: data, consists of: $pool1, pool2... and $n
#' @param region.ll.list: list containing the ll and the indexes of the null alternative and both
#' @param alpha.one: hyper parameters for functional sequence
#' @param alpha.zero: hyper parameters for background
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @return sum of the -log likelihood across all guides
#' @export prior_dirichlet_ll_singleDisp_eSize()

prior_dirichlet_ll_singleDisp_eSize <- function(alpha.one, alpha.zero, data, region.ll.list, guide.efficiency) {

  alpha1s <- alpha.one**2

  alpha1s.norm <- alpha1s / sum(alpha1s)
  alpha1s.adj <- alpha1s.norm * sum(alpha.zero)

  # cat(u0, u1, v0, v1, "\n")
  hyper <- list(alpha0 = alpha.zero,
                alpha1 = alpha1s.adj)

  out.sg.ll <- estimate_relics_sgrna_log_like_eSize(hyper, data, region.ll.list, guide.efficiency)

  -sum(out.sg.ll$total_guide_ll)
}


#' @title optimize the hyper parameters, each distribution with it's own variance
#' @param scaling.par: scaling parameter for alpha 1 hypers
#' @param data: data, consists of: $pool1, pool2... and $n
#' @param region.ll.list: list containing the ll and the indexes of the null alternative and both
#' @param alpha.one: hyper parameters for functional sequence
#' @param alpha.zero: hyper parameters for background
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @return sum of the -log likelihood across all guides
#' @export prior_dirichlet_ll_eSize()

prior_dirichlet_ll_eSize <- function(alpha.one, alpha.zero, data, region.ll.list, guide.efficiency) {

  alpha1s <- alpha.one**2

  # cat(u0, u1, v0, v1, "\n")
  hyper <- list(alpha0 = alpha.zero,
                alpha1 = alpha1s)

  out.sg.ll <- estimate_relics_sgrna_log_like_eSize(hyper, data, region.ll.list, guide.efficiency)

  -sum(out.sg.ll$total_guide_ll)

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


#' @title records the  effect size as scaling parameter for the hypers across all pools (avg across replicates) for each FS
#' @param hyper: hyperparameters
#' @param input.pp: matrix of all posterior probs
#' @param input.min.rs.pp: threshold for FS detection
#' @param data: list, each element is a replicate
#' @param guide.seg.idx.lst: list: each element is a guide, contaiing the indexes of the segments it overlaps
#' @param seg.guide.idx.lst: list: each elemnt is a segment and all guide that overlap it
#' @param guide.efficiency: either a vector of guide efficiency or NULL
#' @param one.dispersion. logical, if there should be one or two dispersions for the hyper parameters
#' @param seg.info: data frame, coordinates for each segment
#' @return list: total_scaledEffSize, repl_scaledEffSize, fs_effSize_ID
#' @export record_scaling_effectSize()

record_scaling_effectSize <- function(input.pp, input.min.rs.pp, hyper, data,
                                  guide.seg.idx.lst, seg.guide.idx.lst,
                                  guide.efficiency, one.dispersion, seg.info) {

  # list of all alpha 1s (one for each FS)
  fs.alpha1.scaling <- list()
  fs.alpha1.scaling.df <- c()
  fs.eff.size.ID <- c()

  for(f in 1:nrow(input.pp)){

    temp.fs.idx <- which(input.pp[f,] > input.min.rs.pp)

    fs.alpha1.scaling[[f]] <- rep(0, length(data))

    if(length(temp.fs.idx) > 0){
      temp.pp <- input.pp[f,]

      temp.seg.info <- seg.info[temp.fs.idx,]

      # 1. extract the index of guides overlapping all the FS segments
      temp.fs.guide.idx <- sort(unique(unlist(lapply(seg.guide.idx.lst[temp.fs.idx], function(x){
        x$guide_idx
      }))))

      did.converge <- TRUE

      # for each replicate
      for(i in 1:length(data)){

        # extract the subset of guides overlapping the FS
        fs.data <- data[[i]][temp.fs.guide.idx,]
        fs.guide.efficiency <- guide.efficiency[temp.fs.guide.idx]

        # 2. use fs.data to get overpalling pp and calculate the per-guide ll given all the FS overlapping that guide
        # Issue here is that I actually include more pp than the FS is long...
        temp.fs.guide.pp.idx <- sort(unique(unlist(lapply(guide.seg.idx.lst[temp.fs.guide.idx], function(x){
          x
        }))))

        temp.guide.lls.list <- FS_only_guide_ll(temp.pp, guide.seg.idx.lst, temp.fs.guide.pp.idx, temp.fs.guide.idx)

        temp.alpha0 <- hyper$alpha0[[i]]
        temp.alpha1 <- hyper$alpha1[[i]]
        temp.scaling <- 1

        res <- c()

        if(one.dispersion){
          res <- optim(temp.scaling, prior_dirichlet_ll_singleDisp_eSize_scaled, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
                       data = fs.data, region.ll.list = temp.guide.lls.list,
                       guide.efficiency = fs.guide.efficiency, alpha.zero = temp.alpha0, alpha.one = temp.alpha1)
        } else {
          res <- optim(temp.scaling, prior_dirichlet_ll_eSize_scaled, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
                       data = fs.data, region.ll.list = temp.guide.lls.list,
                       guide.efficiency = fs.guide.efficiency, alpha.zero = temp.alpha0, alpha.one = temp.alpha1)
        }

        fs.alpha1.scaling[[f]][i] <- abs(res$par)

        if(res$convergence != 0) {
          warning("estimation of hyperparameters failed to converge")
          did.converge <- FALSE

        }

      }

      # at this stage, all FS alpha1s should be obtained, present in list of lists at fs.alpha1[[f]]
      # compute the average effect size across all replicates
      temp.avg.scaling <- mean(fs.alpha1.scaling[[f]])

      temp.eff.size.loc <- temp.seg.info[, c('chrom', 'start', 'end')]
      temp.eff.size.loc$FS <- paste0('FS', f-1)
      temp.eff.size.loc$score <- temp.avg.scaling
      temp.eff.size.loc$converged <- as.character(did.converge)

      fs.alpha1.scaling.df <- rbind(fs.alpha1.scaling.df, temp.eff.size.loc)

      fs.eff.size.ID <- rbind(fs.eff.size.ID, data.frame(FS = paste0('FS', f-1), FS_effSize = temp.avg.scaling, stringsAsFactors = F))

    } else {
      fs.alpha1.scaling[[f]] <- rep(0, length(data))
    }


  }

  return(list(total_scaledEffSize = fs.alpha1.scaling.df, repl_scaledEffSize = fs.alpha1.scaling, fs_effSize_ID = fs.eff.size.ID))
}


#' @title optimize the hyper parameters, only one dispersion across the two distributions
#' @param scaling.par: scaling parameter for alpha 1 hypers
#' @param data: data, consists of: $pool1, pool2... and $n
#' @param region.ll.list: list containing the ll and the indexes of the null alternative and both
#' @param alpha.one: hyper parameters for functional sequence
#' @param alpha.zero: hyper parameters for background
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @return sum of the -log likelihood across all guides
#' @export prior_dirichlet_ll_singleDisp_eSize_scaled()

prior_dirichlet_ll_singleDisp_eSize_scaled <- function(scaling.par, alpha.one, alpha.zero, data, region.ll.list, guide.efficiency) {

  if(scaling.par == 0){
    scaling.par <- 0.000001
  }

  alpha1s <- alpha.one * abs(scaling.par)

  alpha1s.norm <- alpha1s / sum(alpha1s)
  alpha1s.adj <- alpha1s.norm * sum(alpha.zero)

  # cat(u0, u1, v0, v1, "\n")
  hyper <- list(alpha0 = alpha.zero,
                alpha1 = alpha1s.adj)

  out.sg.ll <- estimate_relics_sgrna_log_like_eSize(hyper, data, region.ll.list, guide.efficiency)

  -sum(out.sg.ll$total_guide_ll)
}


#' @title optimize the hyper parameters, each distribution with it's own variance
#' @param scaling.par: scaling parameter for alpha 1 hypers
#' @param data: data, consists of: $pool1, pool2... and $n
#' @param region.ll.list: list containing the ll and the indexes of the null alternative and both
#' @param alpha.one: hyper parameters for functional sequence
#' @param alpha.zero: hyper parameters for background
#' @param guide.efficiency: data.frame of guide efficiences. Either vector of guide efficiency or NULL
#' @return sum of the -log likelihood across all guides
#' @export prior_dirichlet_ll_eSize_scaled()

prior_dirichlet_ll_eSize_scaled <- function(scaling.par, alpha.one, alpha.zero, data, region.ll.list, guide.efficiency) {

  if(scaling.par == 0){
    scaling.par <- 0.000001
  }

  alpha1s <- alpha.one * abs(scaling.par)

  # cat(u0, u1, v0, v1, "\n")
  hyper <- list(alpha0 = alpha.zero,
                alpha1 = alpha1s)

  out.sg.ll <- estimate_relics_sgrna_log_like_eSize(hyper, data, region.ll.list, guide.efficiency)

  -sum(out.sg.ll$total_guide_ll)

}


#' @title wrapper which returns the updated guide efficiency and the corresponding coefficients used to calculate it
#' @param hyper: hyperparameters
#' @param param: matrix of all posterior probs
#' @param data: list, each element is a replicate
#' @param guide.seg.idx.lst: list: each element is a guide, contaiing the indexes of the segments it overlaps
#' @param guide.efficiency.scores: matrix, each column containing a different set of scores per guide
#' @param ge.coeff vector, containing the guide efficiency coefficients
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
#' @param alpha0.input: hyper parameters for background
#' @param alpha1.input: hyper parameters for FS
#' @param guide.efficiency.scores: matrix, each column containing a different set of scores per guide
#' @return sum of the -log likelihood across all guides
#' @export guide_coeff_ll()

guide_coeff_ll <- function(ge.coeff.param, data, region.ll.list, alpha0.input, alpha1.input, guide.efficiency.scores) {

  # , alpha0.input, alpha1.input
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
#' @param seg.info: contains info for each segment, specicially chromosome information
#' @return One plot for each selection step, credible sets colored in purple
#' @export display_relics_fs_as_tiff()

display_relics_fs_as_tiff <- function(input.L, input.labels, tiff.name, fs.threshold, seg.info){

  cs.col.orig <- rep("darkgrey", ncol(input.L))
  
  if(length(unique(seg.info$chrom)) > 1){
    even.chroms <- unique(seg.info$chrom)[c(FALSE,TRUE)]
    even.chroms.idx <- which(seg.info$chrom %in% even.chroms)
    cs.col.orig[even.chroms.idx] <- 'lightgrey'
  }

  # attempt to put 8 per page
  if(nrow(input.L) <= 9){
    tiff(paste0(tiff.name, '.tiff'), width = 11520 / 2, height = 5760 / 2, units = "px", res = 400)
    par(mfrow = c(3,3))

    for(i in 1:nrow(input.L)){
      cs.col <- cs.col.orig
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
            p.sum <- input.L[1:2,]
          }
          p.sum[p.sum > 1] <- 1
          plot(p.sum, pch=21,
               main="Sum of Posteriors",
               col=cs.col.orig, bg = input.labels, ylab = 'PP', xlab = 'Genome Segment')
          #dev.off()
          break()
        } else{

          cs.col <- cs.col.orig #rep("black", ncol(input.L))
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
                 col=cs.col.orig, bg = input.labels, ylab = 'PP', xlab = 'Genome Segment') #"black"
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
                                    input.alpha.outDir, layer.nr, pool.names){
  
  out.alpha.df <- c()
  
  total.rows <- length(input.bkg.alpha) + length(input.fs.alpha)
  total.cols <- max(c( unlist(lapply(input.bkg.alpha, function(x){length(x) + 1})),
                       unlist(lapply(input.fs.alpha, function(x){length(x) + 1}))))
  
  alpha.names <- c(paste(rep('background', length(input.bkg.alpha)), 'r', c(1:length(input.bkg.alpha)), sep = '_'),
                   paste(rep('FS', length(input.fs.alpha)), 'r', c(1:length(input.fs.alpha)), sep = '_'))
  
  alpha.matrix <- matrix(0, nrow = total.rows, ncol = total.cols) # + length(input.bkg.disp[[1]]))
  
  if(! is.null(pool.names)){
    temp.pool.names <- unique(unlist(pool.names))
    
    for(i in 1:length(input.bkg.alpha)){
      alpha0.scores <- rep(0, length(length(temp.pool.names)))
      alpha1.scores <- rep(0, length(length(temp.pool.names)))
      
      temp.bkg.alpha <- c(1,input.bkg.alpha[[i]]) / sum(c(1, input.bkg.alpha[[i]]))
      temp.fs.alpha <- c(1,input.fs.alpha[[i]]) / sum(c(1, input.fs.alpha[[i]]))
      
      alpha0.scores[match(pool.names[[i]], temp.pool.names)] <- round(temp.bkg.alpha, 3)
      alpha1.scores[match(pool.names[[i]], temp.pool.names)] <- round(temp.fs.alpha, 3)
      
      alpha.matrix[i,] <- alpha0.scores
      alpha.matrix[i + length(input.bkg.alpha),] <- alpha1.scores
      
      # alpha0.scores.wDisp <- c(alpha0.scores, round(input.bkg.disp[[i]], 3))
      # alpha1.scores.wDisp <- c(alpha1.scores, round(input.bkg.disp[[i]], 3))
      
      # alpha.matrix[i,] <- alpha0.scores.wDisp
      # alpha.matrix[i + length(input.bkg.alpha),] <- alpha1.scores.wDisp
      
      out.alpha.df <- cbind(alpha.names, alpha.matrix)
      
      colnames(out.alpha.df) <- c('hyperPar_type', temp.pool.names) #, paste0(rep('dispersion', length(input.bkg.disp[[1]])), '_', c(1:length(input.bkg.disp[[1]]))) )
      
    }
    
  } else {
    
    for(i in 1:length(input.bkg.alpha)){
      temp.bkg.alpha <- c(1,input.bkg.alpha[[i]]) / sum(c(1, input.bkg.alpha[[i]]))
      alpha.matrix[i, c(1:length(temp.bkg.alpha))] <- round(temp.bkg.alpha, 3)
      # alpha.matrix[i, (length(temp.bkg.alpha) + 1):ncol(alpha.matrix)] <- round(input.bkg.disp[[i]], 3)
    }
    
    for(j in (length(input.bkg.alpha) + 1):total.rows){
      temp.fs.alpha <- c(1,input.fs.alpha[[j - length(input.bkg.alpha)]]) / sum(c(1, input.fs.alpha[[j - length(input.bkg.alpha)]]))
      alpha.matrix[j, c(1:length(temp.fs.alpha))] <- round(temp.fs.alpha, 3)
      # alpha.matrix[j, (length(temp.fs.alpha) + 1):ncol(alpha.matrix)] <- round(input.bkg.disp[[j - length(input.bkg.alpha)]], 3)
    }
    
    out.alpha.df <- cbind(alpha.names, alpha.matrix)
    
    colnames(out.alpha.df) <- c('hyperPar_type', paste0('pool', c(1:(nrow(out.alpha.df) - 2)))) #, paste0(rep('dispersion', length(input.bkg.disp[[1]])), '_', c(1:length(input.bkg.disp[[1]]))))
  }
  
  write.csv(out.alpha.df, file = paste0(input.alpha.outDir, '_k', layer.nr, '_hyperPars.csv'), row.names = F, quote = F)
  
  
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


#' @title Record the hyper parameters for the background and the functional sorting probabilities
#' @param input.alpha0.list: list of the null-alphas resulting the max lolg-lik. for that layer
#' @param input.alpha1.list: list of the alternative-alphas resulting the max lolg-lik. for that layer
#' @param input.alpha.outDir: directory to which the alphas are to be written
#' @param layer.nr: number of functional sequences recorded
#' @param pool.names: names of the pools given, includes the total column ($n)
#' @return .csv file
#' @export record_alphas()

record_alphas <- function(input.alpha0.list, input.alpha1.list, input.alpha.outDir, layer.nr, pool.names){
  
  out.alpha.df <- c()
  
  total.rows <- length(input.alpha0.list) + length(input.alpha1.list)
  total.cols <- max(c( unlist(lapply(input.alpha0.list, function(x){length(x)})),
                       unlist(lapply(input.alpha1.list, function(x){length(x)}))))
  
  alpha.names <- c(paste(rep('alpha0', length(input.alpha0.list)), 'r', c(1:length(input.alpha0.list)), sep = '_'),
                   paste(rep('alpha1', length(input.alpha1.list)), 'r', c(1:length(input.alpha1.list)), sep = '_'))
  
  alpha.matrix <- matrix(0, nrow = total.rows, ncol = total.cols + 1)
  
  if(! is.null(pool.names)){
    temp.pool.names <- unique(unlist(pool.names))
    
    for(i in 1:length(input.alpha0.list)){
      alpha0.scores <- rep(0, length(length(temp.pool.names)))
      alpha1.scores <- rep(0, length(length(temp.pool.names)))
      
      alpha0.scores[match(pool.names[[i]], temp.pool.names)] <- round(input.alpha0.list[[i]] / sum(input.alpha0.list[[i]]), 3)
      alpha1.scores[match(pool.names[[i]], temp.pool.names)] <- round(input.alpha1.list[[i]] / sum(input.alpha1.list[[i]]), 3)
      
      alpha0.scores.wDisp <- c(alpha0.scores, round(sum(input.alpha0.list[[i]]), 3))
      alpha1.scores.wDisp <- c(alpha1.scores, round(sum(input.alpha1.list[[i]]), 3))
      
      alpha.matrix[i,] <- alpha0.scores.wDisp
      alpha.matrix[i + length(input.alpha0.list),] <- alpha1.scores.wDisp
      
      out.alpha.df <- cbind(alpha.names, alpha.matrix)
      
      colnames(out.alpha.df) <- c('alpha_type', temp.pool.names, 'dispersion')

    }
    
  } else {
    
    for(i in 1:length(input.alpha0.list)){
      alpha.matrix[i, c(1:length(input.alpha0.list[[i]]))] <- round(input.alpha0.list[[i]] / sum(input.alpha0.list[[i]]), 3)
      alpha.matrix[i, length(input.alpha0.list[[i]]) + 1] <- round(sum(input.alpha0.list[[i]]), 3)
    }
    
    for(j in (length(input.alpha0.list) + 1):total.rows){
      alpha.matrix[j, c(1:length(input.alpha1.list[[j - length(input.alpha0.list)]]))] <- round(input.alpha1.list[[j - length(input.alpha0.list)]] / sum(input.alpha1.list[[j - length(input.alpha0.list)]]), 3)
      alpha.matrix[j, length(input.alpha1.list[[j - length(input.alpha0.list)]]) + 1] <- round(sum(input.alpha1.list[[j - length(input.alpha0.list)]]), 3)
    }
    
    out.alpha.df <- cbind(alpha.names, alpha.matrix)
    
    colnames(out.alpha.df) <- c('alpha_type', paste0('pool', c(1:(nrow(out.alpha.df) - 2))), 'dispersion')
  }

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
#' @param local.max: whether local.max is computed
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
#' @param pp.calculation: type of PP calculation to use. 'v2' is for normal AoE
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
                          local.max, local.max.range,  pp.calculation){


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
                                          local.max, local.max.range, pp.calculation)


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
      print('iterative hyperparameter estimation currently not possible!')
      break()
      # dirichlet.hyper <- recompute_hyper_parameters(dirichlet.param,
      #                                                dirichlet.hyper,
      #                                                input.data.list$data,
      #                                                input.data.list$guide_to_seg_lst,
      #                                               guide.efficiency,
      #                                               one.dispersion)
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
#' @param hyper: hyperparameters, already divided into alpha0 and alpha1
#' @param hyper.components: individual components of the hyper parameters
#' @param param: matrix of all posterior probs
#' @param data: data, consists of: $y1, $y2 and $n
#' @param guide.seg.idx.lst: list: each element is a guide, contaiing the indexes of the segments it overlaps
#' @param guide.efficiency: either a vector of guide efficiency or NULL
#' @param one.dispersion. logical, if there should be one or two dispersions for the hyper parameters
#' @param mean.var.type: type of mean-variance relationship
#' @return list of the re-estimated hyperparameters
#' @export recompute_hyper_parameters()

recompute_hyper_parameters <- function(param, hyper, hyper.components, data, guide.seg.idx.lst, guide.efficiency, one.dispersion, mean.var.type) {
  cumulative.pp <- colSums(param$delta.pp) #apply(param$delta.pp, 2, sum)
  cumulative.pp[cumulative.pp > 1] <- 1

  guide.lls.list <- compute_perGuide_fs_ll(cumulative.pp, guide.seg.idx.lst)

  for(i in 1:length(data)){
    # need alphas to be greater than 0, taking square root now to exponentiate after
    # hyper.param <- c(sqrt(hyper$alpha0[[i]]), sqrt(hyper$alpha1[[i]]))
    # alpha0.idx <- c(1:length(hyper$alpha0[[i]]))
    # alpha1.idx <- c(1:length(hyper$alpha0[[i]])) + max(alpha0.idx)
    
    temp.hyper.params <- c(sqrt(hyper.components$bkg_alpha[[i]]), sqrt(hyper.components$FS_alpha[[i]]), hyper.components$bkg_dispersion[[i]])
    temp.bkg.idx <- c(1:length(hyper.components$bkg_alpha[[i]]))
    temp.fs.idx <- c(1:length(hyper.components$FS_alpha[[i]])) + max(temp.bkg.idx)
    temp.disp.idx <- c(1:length(hyper.components$bkg_dispersion[[i]])) + max(temp.fs.idx)

    temp.res.drch <- c()

    if(one.dispersion){
      # res <- optim(hyper.param, prior_dirichlet_ll_singleDisp, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
      #              data=data[[i]], region.ll.list = guide.lls.list,
      #              alpha0.idx = alpha0.idx, alpha1.idx = alpha1.idx,
      #              guide.efficiency = guide.efficiency)
      temp.res.drch <- optim(temp.hyper.params, prior_dirichlet_parameters, method= 'L-BFGS-B',
                             data = data[[i]], 
                             region.ll.list = guide.lls.list,
                             bkg.idx = temp.bkg.idx, 
                             fs.idx = temp.fs.idx, 
                             disp.idx = temp.disp.idx, 
                             guide.efficiency = guide.efficiency, 
                             mean.var.type = mean.var.type)
      
      
    } else {
      print('currently one one dispersion available')
      break()
      # res <- optim(hyper.param, prior_dirichlet_ll, method= 'L-BFGS-B', #'BFGS', #"Nelder-Mead",
      #              data=data[[i]], region.ll.list = guide.lls.list,
      #              alpha0.idx = alpha0.idx, alpha1.idx = alpha1.idx,
      #              guide.efficiency = guide.efficiency)
    }

    # account for flat surface!
    if(temp.res.drch$convergence %in% c(0, 52) ) {
      # return new estimates of hyperparamers

      if(one.dispersion){
        # alpha0s <- res$par[alpha0.idx]**2
        # alpha1s <- res$par[alpha1.idx]**2
        # 
        # alpha1s.norm <- alpha1s / sum(alpha1s)
        # alpha1s.adj <- alpha1s.norm * sum(alpha0s)
        # 
        # hyper$alpha0[[i]] <- alpha0s
        # hyper$alpha1[[i]] <- alpha1s.adj
        
        temp.bkg.alpha <- c(1, temp.res.drch$par[temp.bkg.idx]**2) / sum(c(1, temp.res.drch$par[temp.bkg.idx]**2))
        temp.fs.alpha <- c(1, temp.res.drch$par[temp.fs.idx]**2) / sum(c(1, temp.res.drch$par[temp.fs.idx]**2))
        temp.bkg.disp <- temp.res.drch$par[temp.disp.idx]
        
        hyper.components$bkg_alpha[[i]] <- temp.res.drch$par[temp.bkg.idx]**2
        hyper.components$FS_alpha[[i]] <- temp.res.drch$par[temp.fs.idx]**2
        hyper.components$bkg_dispersion[[i]] <- temp.bkg.disp
        
        if(mean.var.type == 'radical'){
          temp.disp <- temp.bkg.disp[1] + temp.bkg.disp[2] * data[[i]][,ncol(data[[i]])] + temp.bkg.disp[3] * sqrt(data[[i]][,ncol(data[[i]])])
          temp.disp[which(temp.disp < 0.1)] <- 0.1
          hyper$alpha0[[i]] <- t(temp.bkg.alpha %*% t(temp.disp) )
          hyper$alpha1[[i]] <- t(temp.fs.alpha %*% t(temp.disp) )
        } else if(mean.var.type == 'exponential'){
          temp.disp <- temp.bkg.disp[1] + temp.bkg.disp[2] * log(data[[i]][,ncol(data[[i]])])
          # temp.disp <- temp.bkg.disp[1] + temp.bkg.disp[2] * data[[i]][,ncol(data[[i]])] + temp.bkg.disp[3] * log(data[[i]][,ncol(data[[i]])])
          temp.disp[which(temp.disp < 0.1)] <- 0.1
          temp.disp[which(is.na(temp.disp))] <- 0.1
          temp.disp[which(temp.disp == Inf)] <- max(temp.disp[which(temp.disp < Inf)])
          hyper$alpha0[[i]] <- t(temp.bkg.alpha %*% t(temp.disp) )
          hyper$alpha1[[i]] <- t(temp.fs.alpha %*% t(temp.disp) )
        } else if(mean.var.type == 'independent'){
          temp.disp <- rep(temp.bkg.disp[1]^2, nrow(data[[i]]))
          hyper$alpha0[[i]] <- t(temp.bkg.alpha %*% t(temp.disp) )
          hyper$alpha1[[i]] <- t(temp.fs.alpha %*% t(temp.disp) )
        } else if(mean.var.type == 'quadratic'){
          temp.disp <- temp.bkg.disp[1] + temp.bkg.disp[2] * data[[i]][,ncol(data[[i]])] + temp.bkg.disp[3] * data[[i]][,ncol(data[[i]])]^2
          temp.disp[which(temp.disp < 0.1)] <- 0.1
          temp.disp[which(is.na(temp.disp))] <- 0.1
          temp.disp[which(temp.disp == Inf)] <- max(temp.disp[which(temp.disp < Inf)])
          hyper$alpha0[[i]] <- t(temp.bkg.alpha %*% t(temp.disp) )
          hyper$alpha1[[i]] <- t(temp.fs.alpha %*% t(temp.disp) )
        }
        
      } else {
        print('currently one one dispersion available')
        break()
        # hyper$alpha0[[i]] <- res$par[alpha0.idx]**2
        # hyper$alpha1[[i]] <- res$par[alpha1.idx]**2
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
#' @param pp.calculation: tyep of PP calculation. 'v2' is for normal AoE
#' @return log likelihood for each region
#' @export relics_estimate_pp()

relics_estimate_pp <- function(param, hyper, data, known.reg,
                                guide.to.seg.lst, seg.to.guide.lst,
                                next.guide.lst, nr.segs = 10,
                               geom.p = 0.1, guide.efficiency,
                               local.max, local.max.range, pp.calculation) {
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

    # sgrna.log.like.list <- guide_ll_calculation(pp, data, guide.to.seg.lst, hyper, guide.efficiency)

    delta.pps <- c()
    if(pp.calculation == 'v4'){
      delta.pps <- estimate_fs_AoE_pp_v2(hyper, data.mat.list, data.total.list, guide.to.seg.lst,
                                         seg.to.guide.lst, next.guide.lst, nr.segs, geom.p, sgrna.log.like.list)
    # } else if(pp.calculation == 'v2'){
    #   delta.pps <- estimate_fs_AoE_pp(hyper, data.mat.list, data.total.list, guide.to.seg.lst,
    #                                   seg.to.guide.lst, next.guide.lst, nr.segs, geom.p, sgrna.log.like.list)
    # } else if(pp.calculation == 'v3'){
    #   seg.ll.list <- seg_ll_calculation(pp, data, seg.to.guide.lst, hyper, guide.efficiency, guide.to.seg.lst)
    #   
    #   delta.pps <- estimate_seg_fs_AoE_pp(hyper, data.mat.list, data.total.list, guide.to.seg.lst,
    #                                       seg.to.guide.lst, next.guide.lst, nr.segs, geom.p, 
    #                                       seg.ll.list$out_seg_ll, seg.ll.list$out_seg_pp_ll)
    } else {
      delta.pps <- estimate_fs_pp(#hyper, data.mat.list, data.total.list,
        seg.to.guide.lst, next.guide.lst, nr.segs, geom.p, sgrna.log.like.list)
    }
    

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

seg_pp_calculation <- function(cumulative.pp, seg.to.guide.lst, hyper, guide.efficiency, guide.to.seg.lst,
                               in.data.list, nr.segs = 10, geom.p = 0.1){
  
  # set up a list, the length of the number of segments
  # each element contains a vector, which contains the ll of that segment being overlapped by a regulatory element
  segment.ll.list <- list()
  
  total.segs <- length(seg.to.guide.lst)
  region.priors <- rep(1/total.segs, total.segs)
  geom.norm.factr <- pgeom(nr.segs, geom.p)
  
  total.ll <- log(0)
  
  for(seg in 1:total.segs){
    temp.segment.ll <- 0
    temp.guide.idx<- seg.to.guide.lst[[seg]]$guide_idx
    # temp.nonGuide.idx<- seg.to.guide.lst[[seg]]$nonGuide_idx
    
    # for each guide, identify what segments we're processing (what index matches current segment)
    # return the max distance / probability of all segments covered
    temp.seg.dist <- unlist(lapply(guide.to.seg.lst[temp.guide.idx], function(x){
      adj.seg.idx <- which(x$seg_overlapped %in% seg)
      temp.dist.to.seg <- x$dist_to_seg[adj.seg.idx]
      max(temp.dist.to.seg)
    }))
    
    temp.pp <- cumulative.pp[seg]
    
    for(repl in 1:length(in.data.list)){
      
      # initialize the segment lls, both the new ll (temp.seg.pp.ll) as well as the one using the previous FS
      temp.seg.ll <- 0
      temp.seg.pp.ll <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx,],
                                     size = in.data.totals[[repl]][temp.guide.idx],
                                     alpha = hyper$alpha1[[repl]][temp.guide.idx,], log = T) + log(temp.seg.dist))
      
      if(temp.pp == 0){
        temp.seg.ll <- sum(log(temp.seg.dist) + ddirmnom(in.data.list[[repl]][temp.guide.idx,],
                                                         size = in.data.totals[[repl]][temp.guide.idx],
                                                         alpha = hyper$alpha0[[repl]][temp.guide.idx,], log = T))
      } else if(temp.pp == 1){
        temp.seg.ll <- sum(log(temp.seg.dist) + ddirmnom(in.data.list[[repl]][temp.guide.idx,],
                                                         size = in.data.totals[[repl]][temp.guide.idx],
                                                         alpha = hyper$alpha1[[repl]][temp.guide.idx,], log = T))
      } else {
        temp.seg.bkg.ll <- sum(log(temp.seg.dist) + ddirmnom(in.data.list[[repl]][temp.guide.idx,],
                                                             size = in.data.totals[[repl]][temp.guide.idx],
                                                             alpha = hyper$alpha0[[repl]][temp.guide.idx,], log = T) + log(1 - temp.pp))
        temp.seg.fs.ll <- sum(log(temp.seg.dist) + ddirmnom(in.data.list[[repl]][temp.guide.idx,],
                                                            size = in.data.totals[[repl]][temp.guide.idx],
                                                            alpha = hyper$alpha1[[repl]][temp.guide.idx,], log = T) + log(temp.pp))
        temp.seg.ll <- addlogs(temp.seg.fs.ll, temp.seg.bkg.ll)
      }
      
      temp.segment.ll <- temp.segment.ll + (temp.seg.pp.ll - temp.seg.ll) + log(dgeom(1, geom.p) / geom.norm.factr)

    }
    segment.ll.list[[seg]] <- temp.segment.ll
    total.ll <- addlogs(total.ll, temp.segment.ll)
  }

  # now make a pass through all 2 to nr.segs length segments
  if(nr.segs > 1){
    for(ns in 1:(nr.segs-1) ){
      
      for(seg in 1:(total.segs - nr.segs)){
        
        temp.stretch.ll <- 0
        temp.stretch.segs <- c(seg:(seg + ns)) # index of segments to be included as part of the stretch and be updated in this iteration
        temp.guide.idx <- c(seg.to.guide.lst[[seg]]$guide_idx, unlist(next.guide.lst[(seg + 1):(seg + ns)]))
        # temp.nonGuide.idx <- all.guide.idx[-temp.guide.idx]
        
        temp.seg.dist <- unlist(lapply(guide.to.seg.lst[temp.guide.idx], function(x){
          adj.seg.idx <- which(x$seg_overlapped %in% temp.stretch.segs)
          temp.dist.to.seg <- x$dist_to_seg[adj.seg.idx]
          max(temp.dist.to.seg)
        }))
        
        for(repl in 1:ncol(sg.ll.df)){
          
          temp.seg.ll <- 0
          temp.seg.pp.ll <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx,],
                                         size = in.data.totals[[repl]][temp.guide.idx],
                                         alpha = hyper$alpha1[[repl]][temp.guide.idx,], log = T) + log(temp.seg.dist))
          
          if(temp.pp == 0){
            temp.seg.ll <- sum(log(temp.seg.dist) + ddirmnom(in.data.list[[repl]][temp.guide.idx,],
                                                             size = in.data.totals[[repl]][temp.guide.idx],
                                                             alpha = hyper$alpha0[[repl]][temp.guide.idx,], log = T))
          } else if(temp.pp == 1){
            temp.seg.ll <- sum(log(temp.seg.dist) + ddirmnom(in.data.list[[repl]][temp.guide.idx,],
                                                             size = in.data.totals[[repl]][temp.guide.idx],
                                                             alpha = hyper$alpha1[[repl]][temp.guide.idx,], log = T))
          } else {
            temp.seg.bkg.ll <- sum(log(temp.seg.dist) + ddirmnom(in.data.list[[repl]][temp.guide.idx,],
                                                                 size = in.data.totals[[repl]][temp.guide.idx],
                                                                 alpha = hyper$alpha0[[repl]][temp.guide.idx,], log = T) + log(1 - temp.pp))
            temp.seg.fs.ll <- sum(log(temp.seg.dist) + ddirmnom(in.data.list[[repl]][temp.guide.idx,],
                                                                size = in.data.totals[[repl]][temp.guide.idx],
                                                                alpha = hyper$alpha1[[repl]][temp.guide.idx,], log = T) + log(temp.pp))
            temp.seg.ll <- addlogs(temp.seg.fs.ll, temp.seg.bkg.ll)
          }
          
          temp.segment.ll <- temp.segment.ll + (temp.seg.pp.ll - temp.seg.ll) + log(dgeom(1, geom.p) / geom.norm.factr)
          
        }
        
        total.ll <- addlogs(total.ll, temp.stretch.ll)
        
        for(span in 1:length(temp.stretch.segs)){
          segment.ll.list[[temp.stretch.segs[span]]] <- c(segment.ll.list[[temp.stretch.segs[span]]], temp.stretch.ll)
        }
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

seg_ll_calculation <- function(cumulative.pp, data, seg.to.guide.lst, hyper, guide.efficiency, guide.to.seg.lst){
  
  n.seg <- length(seg.to.guide.lst)
  out.seg.ll <- c()
  out.seg.pp.ll <- c()
  
  for(repl in 1:length(data)){
    
    repl.seg.ll <- vector('numeric', length = n.seg)
    repl.data <- data[[repl]]
    pool.cols <- c(1:(ncol(repl.data) - 1))
    temp.hypers <- list(alpha0 = hyper$alpha0[[repl]], alpha1 = hyper$alpha1[[repl]])
    
    sgrna.null.log.like <- ddirmnom(repl.data[,pool.cols], size = repl.data[,ncol(repl.data)], alpha = temp.hypers$alpha0, log = T)
    
    sgrna.alt.log.like <- c()
    if(is.null(guide.efficiency)){
      sgrna.alt.log.like <- ddirmnom(repl.data[,pool.cols], size = repl.data[,ncol(repl.data)], alpha = temp.hypers$alpha1, log = T)
    } else {
      # toDo, difference in alphas or in proportions? 
      alpha0.matrix <- t(apply(matrix(0, ncol = length(temp.hypers$alpha0), nrow = length(guide.efficiency)), 1, function(x){x + temp.hypers$alpha0}))
      alpha.diffs <- temp.hypers$alpha0 - temp.hypers$alpha1
      
      alpha1.matrix <- matrix(0, ncol = length(temp.hypers$alpha0), nrow = length(guide.efficiency))
      
      for(i in 1:length(guide.efficiency)){
        # this cound be an issue if the dispersion is no longer identical...
        alpha1.matrix[i,] <- alpha0.matrix[i,] - alpha.diffs * guide.efficiency[i]
      }
      
      sgrna.alt.log.like <- ddirmnom(repl.data[,pool.cols], size = repl.data[,ncol(repl.data)], alpha = alpha1.matrix, log = T)
    }
    
    repl.seg.ll <- vector('numeric', length = n.seg)
    repl.seg.pp.ll <- vector('numeric', length = n.seg)
    # for every guide, extract the distances
    # compute the null and the alternative for each of the distances
    for(j in 1:n.seg){
      
      temp.pp <- cumulative.pp[j]
      temp.overl.guide.idx <- seg.to.guide.lst[[j]]$guide_idx
      
      temp.seg.dist <- unlist(lapply(guide.to.seg.lst[temp.overl.guide.idx], function(x){
        adj.seg.idx <- which(x$seg_overlapped %in% j)
        temp.dist.to.seg <- x$dist_to_seg[adj.seg.idx]
        max(temp.dist.to.seg)
      }))
      
      temp.seg.ll <- 0
      temp.seg.pp.ll <- sum(log(temp.seg.dist) + sgrna.alt.log.like[temp.overl.guide.idx])
      
      if(temp.pp == 0){
        temp.seg.ll <- sum(log(temp.seg.dist) + sgrna.null.log.like[temp.overl.guide.idx])
      } else if(temp.pp == 1){
        temp.seg.ll <- sum(log(temp.seg.dist) + sgrna.alt.log.like[temp.overl.guide.idx])
      } else {
        temp.seg.bkg.ll <- sum(log(temp.seg.dist) + sgrna.null.log.like[temp.overl.guide.idx] + log(1 - temp.pp))
        temp.seg.fs.ll <- sum(log(temp.seg.dist) + sgrna.alt.log.like[temp.overl.guide.idx] + log(temp.pp))
        temp.seg.ll <- addlogs(temp.seg.fs.ll, temp.seg.bkg.ll)
      }
      
      repl.seg.ll[j] <- temp.seg.ll
      repl.seg.pp.ll[j] <- temp.seg.pp.ll
    }
    
    out.seg.ll <- cbind(out.seg.ll, repl.seg.ll)
    out.seg.pp.ll <- cbind(out.seg.pp.ll, repl.seg.pp.ll)
    
  }
  
  out.list <- list(out_seg_ll = out.seg.ll, out_seg_pp_ll = out.seg.pp.ll)
  return(out.list)
}

guide_ll_calculation <- function(cumulative.pp, data, guide.seg.idx.lst, hyper, guide.efficiency){
  
  n.sgrna <- length(guide.seg.idx.lst)
  out.guide.ll <- list()
  
  for(repl in 1:length(data)){
    
    repl.guide.ll <- vector('numeric', length = n.sgrna)
    repl.data <- data[[repl]]
    pool.cols <- c(1:(ncol(repl.data) - 1))
    temp.hypers <- list(alpha0 = hyper$alpha0[[repl]], alpha1 = hyper$alpha1[[repl]])
    
    sgrna.null.log.like <- ddirmnom(repl.data[,pool.cols], size = repl.data[,ncol(repl.data)], alpha = temp.hypers$alpha0, log = T)
    
    sgrna.alt.log.like <- c()
    if(is.null(guide.efficiency)){
      sgrna.alt.log.like <- ddirmnom(repl.data[,pool.cols], size = repl.data[,ncol(repl.data)], alpha = temp.hypers$alpha1, log = T)
    } else {
      # toDo, difference in alphas or in proportions? 
      alpha0.matrix <- t(apply(matrix(0, ncol = length(temp.hypers$alpha0), nrow = length(guide.efficiency)), 1, function(x){x + temp.hypers$alpha0}))
      alpha.diffs <- temp.hypers$alpha0 - temp.hypers$alpha1
      
      alpha1.matrix <- matrix(0, ncol = length(temp.hypers$alpha0), nrow = length(guide.efficiency))
      
      for(i in 1:length(guide.efficiency)){
        # this cound be an issue if the dispersion is no longer identical...
        alpha1.matrix[i,] <- alpha0.matrix[i,] - alpha.diffs * guide.efficiency[i]
      }
      
      sgrna.alt.log.like <- ddirmnom(repl.data[,pool.cols], size = repl.data[,ncol(repl.data)], alpha = alpha1.matrix, log = T)
    }
    
    # for every guide, extract the distances
    # compute the null and the alternative for each of the distances
    for(j in 1:n.sgrna){
      
      temp.dist <- guide.seg.idx.lst[[j]]$dist_to_seg
      temp.pp <- cumulative.pp[guide.seg.idx.lst[[j]]$seg_overlapped]
      
      temp.zero.pp.idx <- which(temp.pp == 0)
      temp.one.pp.idx <- which(temp.pp == 1)
      
      temp.bkg.guide.ll.seg <- log(temp.dist) + log(1 - temp.pp) + sgrna.null.log.like[j]
      temp.fs.guide.ll.seg <- log(temp.dist) + log(temp.pp) + sgrna.alt.log.like[j]
      
      # data.frame for each of the segments
      # compute bth the null and the alternative
      # for 0 pp compute only the null
      # for 1 pp compute only the alternative
      
      if(length(temp.zero.pp.idx) > 0){
        temp.fs.guide.ll.seg[temp.zero.pp.idx] <- 0
      }
      if(length(temp.one.pp.idx) > 0){
        temp.bkg.guide.ll.seg[temp.one.pp.idx] <- 0
      }
      
      temp.seg.ll <- temp.bkg.guide.ll.seg + temp.fs.guide.ll.seg
      temp.guide.total.ll <- temp.seg.ll[1]
      
      if(length(temp.seg.ll) > 1){
        for(i in 2:length(temp.seg.ll)){
          temp.guide.total.ll <- addlogs(temp.guide.total.ll, temp.seg.ll[i])
        }
      }
      
      
      # temp.guide.ll <- addlogs(sum(temp.bkg.guide.ll.seg), sum(temp.fs.guide.ll.seg))
      repl.guide.ll[j] <- temp.guide.total.ll
    }
    
    out.guide.ll[[repl]] <- repl.guide.ll
    
  }
  
  return(out.guide.ll)
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


#' @title Compute the log likelihood of each possible configuration of the placement of a fs. Take into account distance between segment and the guide target site
#' @param hyper: list, hyperparameters
#' @param in.data.list: list, each element contains a replicate, the total column has already been removed
#' @param in.data.totals: list, each element contains the per-guide totals for a replicate
#' @param seg.to.guide.lst: list, each element is a segment, contains a list with guide_idx, nonGuide_idx
#' @param next.guide.lst: list that contains the index of $next_guide_idx and next_nonGuide_idx
#' @param nr.segs: max number of segemnts that can be combined
#' @param geom.p: probability of the geometic distribution used to calculate the probability of there being x segments in the regulatory element
#' @param sg.ll.list: list, each element is a data frame per.guide log likelihood given the posteriors for a replicate and alt. ll (total_guide_ll, alt_only_ll)
#' @return log likelihood for each segment
#' @export estimate_fs_AoE_pp_v2()

estimate_fs_AoE_pp_v2 <- function(hyper, in.data.list, in.data.totals, guide.to.seg.lst,
                               seg.to.guide.lst, next.guide.lst, nr.segs = 10, geom.p = 0.1,
                               sg.ll.list) {
  
  # set up a list, the length of the number of segments
  # each element contains a vector, which contains the ll of that segment being overlapped by a regulatory element
  segment.ll.list <- list()
  
  total.segs <- length(seg.to.guide.lst)
  region.priors <- rep(1/total.segs, total.segs)
  geom.norm.factr <- pgeom(nr.segs, geom.p)
  
  # # initialize the total likelihood with the first, one-segment likelihood
  # initial.segment.ll <- 0
  # initial.guide.idx<- seg.to.guide.lst[[1]]$guide_idx
  # initial.nonGuide.idx<- seg.to.guide.lst[[1]]$nonGuide_idx
  # 
  # initial.seg.dist <- guide.to.seg.lst[[initial.guide.idx]]$dist_to_seg[1] #unlist(lapply(guide.to.seg.lst[[initial.guide.idx]]$dist_to_seg, function(x){x[1]}))
  # 
  # for(repl in 1:length(sg.ll.list)){ # in.data.list
  #   
  #   initial.guide.ll <- sum(ddirmnom(in.data.list[[repl]][initial.guide.idx,],
  #                                    size = in.data.totals[[repl]][initial.guide.idx],
  #                                    alpha = hyper$alpha1[[repl]][initial.guide.idx,], log = T) + 
  #                             log(1 - dpoibin(0, initial.seg.dist) ))
  #   
  #   initial.guide.ll.bkg <- sum(ddirmnom(in.data.list[[repl]][initial.guide.idx[init.with.bkg.idx],],
  #                                        size = in.data.totals[[repl]][initial.guide.idx[init.with.bkg.idx]],
  #                                        alpha = hyper$alpha0[[repl]][initial.guide.idx[init.with.bkg.idx],], log = T) + 
  #                                 log(dpoibin(0, initial.seg.dist)))
  #   
  #   init.with.bkg.idx <- which(initial.seg.dist < 1)
  #   
  #   if(length(init.with.bkg.idx) > 0){
  #     initial.guide.ll.bkg <- sum(ddirmnom(in.data.list[[repl]][initial.guide.idx[init.with.bkg.idx],],
  #                                          size = in.data.totals[[repl]][initial.guide.idx[init.with.bkg.idx]],
  #                                          alpha = hyper$alpha0[[repl]][initial.guide.idx[init.with.bkg.idx],], log = T) + 
  #                                   log(1 - initial.seg.dist[init.with.bkg.idx])) # 1 - 
  #   }
  #   
  #   # initial.guide.ll <- sum(sg.ll.list[[repl]][initial.guide.idx,2])
  #   
  #   
  #   initial.nonGuide.ll <- sum(sg.ll.list[[repl]][initial.nonGuide.idx,1])
  #   # initial.nonGuide.ll <- sum(sg.ll.list[[repl]][initial.nonGuide.idx])
  #   initial.segment.ll <- initial.segment.ll + addlogs(initial.guide.ll, initial.guide.ll.bkg) + initial.nonGuide.ll + log(dgeom(1, geom.p) / geom.norm.factr)
  # }
  # 
  # segment.ll.list[[1]] <- initial.segment.ll
  
  total.ll <- log(0)
  
  # make an initial pass with one-segment length RE placements
  for(seg in 1:total.segs){
    temp.segment.ll <- 0
    temp.guide.idx<- seg.to.guide.lst[[seg]]$guide_idx
    temp.nonGuide.idx<- seg.to.guide.lst[[seg]]$nonGuide_idx
    # 
    # region.pp <- cumulative.pp[guide.seg.idx.lst[[j]]$seg_overlapped]
    # 
    # region.pp.distance.adj <- region.pp * guide.seg.idx.lst[[j]]$dist_to_seg
    # 
    
    # for each guide, identify what segments we're processing (what index matches current segment)
    # return the max distance / probability of all segments covered
    # temp.seg.dist <- unlist(lapply(guide.to.seg.lst[temp.guide.idx], function(x){
    #   adj.seg.idx <- which(x$seg_overlapped %in% seg)
    #   temp.dist.to.seg <- x$dist_to_seg[adj.seg.idx]
    #   max(temp.dist.to.seg)
    # }))
    
    temp.seg.poi.bkg <- unlist(lapply(guide.to.seg.lst[temp.guide.idx], function(x){
      adj.seg.idx <- which(x$seg_overlapped %in% seg)
      temp.dist.to.seg <- x$dist_to_seg[adj.seg.idx]
      dpoibin(0, temp.dist.to.seg)
    }))
    
    # p.fs.idx <- which(temp.seg.poi.bkg < 1)
    # p.bkg.idx <- which(temp.seg.poi.bkg > 0)

    for(repl in 1:length(sg.ll.list)){ # in.data.list
      
      temp.guide.ll <- vector('numeric', length = length(temp.seg.poi.bkg))
      
      for(i in 1:length(temp.seg.poi.bkg)){
        temp.guide.ll[i] <- addlogs(sg.ll.list[[repl]][temp.guide.idx[i],2] + log(1 - temp.seg.poi.bkg[i]),
                                    sg.ll.list[[repl]][temp.guide.idx[i],1] + log(temp.seg.poi.bkg[i]))
        # temp.guide.ll <- sum(temp.guide.ll,
        #                      addlogs(ddirmnom(in.data.list[[repl]][temp.guide.idx[i],],
        #                                       size = in.data.totals[[repl]][temp.guide.idx[i]],
        #                                       alpha = hyper$alpha1[[repl]][temp.guide.idx[i],], log = T) +
        #                                log(1 - temp.seg.poi.bkg[i]),
        #                              ddirmnom(in.data.list[[repl]][temp.guide.idx[i],],
        #                                       size = in.data.totals[[repl]][temp.guide.idx[i]],
        #                                       alpha = hyper$alpha0[[repl]][temp.guide.idx[i],], log = T) +
        #                                log(temp.seg.poi.bkg[i]) ) )
      }

      # temp.fs.ll <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx[p.fs.idx],],
      #                            size = in.data.totals[[repl]][temp.guide.idx[p.fs.idx]],
      #                            alpha = hyper$alpha1[[repl]][temp.guide.idx[p.fs.idx],], log = T) + 
      #                     log(1 - temp.seg.poi.bkg[p.fs.idx]) )
      # 
      # temp.bkg.ll <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx[p.bkg.idx],],
      #                             size = in.data.totals[[repl]][temp.guide.idx[p.bkg.idx]],
      #                             alpha = hyper$alpha0[[repl]][temp.guide.idx[p.bkg.idx],], log = T) + 
      #                      log(temp.seg.poi.bkg[p.bkg.idx]) )

      
      # temp.guide.ll <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx,],
      #                               size = in.data.totals[[repl]][temp.guide.idx],
      #                               alpha = hyper$alpha1[[repl]][temp.guide.idx,], log = T) +
      #                        log(1 - temp.seg.poi.bkg) ) # dpoibin(0, temp.seg.dist)
      # 
      # 
      # 
      # temp.guide.ll.bkg <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx,],
      #                                   size = in.data.totals[[repl]][temp.guide.idx],
      #                                   alpha = hyper$alpha0[[repl]][temp.guide.idx,], log = T) +
      #                            log(temp.seg.poi.bkg) ) # temp.seg.poi.bkg
      
      # temp.with.bkg.idx <- which(temp.seg.dist < 1)
      # 
      # if(length(temp.with.bkg.idx) > 0){
      #   temp.guide.ll.bkg <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx[temp.with.bkg.idx],],
      #                                     size = in.data.totals[[repl]][temp.guide.idx[temp.with.bkg.idx]],
      #                                     alpha = hyper$alpha1[[repl]][temp.guide.idx[temp.with.bkg.idx],], log = T) + 
      #                              log(1 - temp.seg.dist[temp.with.bkg.idx])) # 1 - 
      # }
      # 
      # temp.guide.ll <- sum(sg.ll.list[[repl]][temp.guide.idx,2])
      temp.nonGuide.ll <- sum(sg.ll.list[[repl]][temp.nonGuide.idx,1])
      # temp.nonGuide.ll <- sum(sg.ll.list[[repl]][temp.nonGuide.idx])
      # temp.segment.ll <- temp.segment.ll + addlogs(temp.guide.ll, temp.guide.ll.bkg) + temp.nonGuide.ll + log(dgeom(1, geom.p) / geom.norm.factr)
      # temp.segment.ll <- temp.segment.ll + addlogs(temp.fs.ll, temp.bkg.ll) + temp.nonGuide.ll + log(dgeom(1, geom.p) / geom.norm.factr)
      temp.segment.ll <- temp.segment.ll + sum(temp.guide.ll) + temp.nonGuide.ll + log(dgeom(1, geom.p) / geom.norm.factr)
    }
    segment.ll.list[[seg]] <- temp.segment.ll
    total.ll <- addlogs(total.ll, temp.segment.ll)
  }
  
  # all.guide.idx <- c(1:length(sg.ll.list[[1]][,1])) # helps to differentitate between guides for null vs alternative
  all.guide.idx <- c(1:nrow(in.data.list[[1]])) # helps to differentitate between guides for null vs alternative
  
  # now make a pass through all 2 to nr.segs length segments
  for(ns in 1:(nr.segs-1) ){
    
    for(seg in 1:(total.segs - ns)){
      
      temp.stretch.ll <- 0
      temp.stretch.segs <- c(seg:(seg + ns)) # index of segments to be included as part of the stretch and be updated in this iteration
      temp.guide.idx <- c(seg.to.guide.lst[[seg]]$guide_idx, unlist(next.guide.lst[(seg + 1):(seg + ns)]))
      temp.nonGuide.idx <- all.guide.idx[-temp.guide.idx]
      
      # temp.seg.dist <- unlist(lapply(guide.to.seg.lst[temp.guide.idx], function(x){
      #   adj.seg.idx <- which(x$seg_overlapped %in% temp.stretch.segs)
      #   temp.dist.to.seg <- x$dist_to_seg[adj.seg.idx]
      #   max(temp.dist.to.seg)
      # }))
      
      temp.seg.poi.bkg <- unlist(lapply(guide.to.seg.lst[temp.guide.idx], function(x){
        adj.seg.idx <- which(x$seg_overlapped %in% temp.stretch.segs)
        temp.dist.to.seg <- x$dist_to_seg[adj.seg.idx]
        dpoibin(0, temp.dist.to.seg)
      }))
      
      # p.fs.idx <- which(temp.seg.poi.bkg < 1)
      # p.bkg.idx <- which(temp.seg.poi.bkg > 0)
      
      for(repl in 1:length(sg.ll.list)){ # in.data.list
        
        temp.guide.ll <- vector('numeric', length = length(temp.seg.poi.bkg))

        for(i in 1:length(temp.seg.poi.bkg)){
          temp.guide.ll[i] <- addlogs(sg.ll.list[[repl]][temp.guide.idx[i],2] + log(1 - temp.seg.poi.bkg[i]),
                                      sg.ll.list[[repl]][temp.guide.idx[i],1] + log(temp.seg.poi.bkg[i]))
          # temp.guide.ll <- sum(temp.guide.ll,
          #                      addlogs(ddirmnom(in.data.list[[repl]][temp.guide.idx[i],],
          #                                       size = in.data.totals[[repl]][temp.guide.idx[i]],
          #                                       alpha = hyper$alpha1[[repl]][temp.guide.idx[i],], log = T) +
          #                                log(1 - temp.seg.poi.bkg[i]),
          #                              ddirmnom(in.data.list[[repl]][temp.guide.idx[i],],
          #                                       size = in.data.totals[[repl]][temp.guide.idx[i]],
          #                                       alpha = hyper$alpha0[[repl]][temp.guide.idx[i],], log = T) +
          #                                log(temp.seg.poi.bkg[i]) ) )
        }
        
        # 
        # temp.fs.ll <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx[p.fs.idx],],
        #                            size = in.data.totals[[repl]][temp.guide.idx[p.fs.idx]],
        #                            alpha = hyper$alpha1[[repl]][temp.guide.idx[p.fs.idx],], log = T) + 
        #                     log(1 - temp.seg.poi.bkg[p.fs.idx]) )
        # 
        # temp.bkg.ll <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx[p.bkg.idx],],
        #                             size = in.data.totals[[repl]][temp.guide.idx[p.bkg.idx]],
        #                             alpha = hyper$alpha0[[repl]][temp.guide.idx[p.bkg.idx],], log = T) + 
        #                      log(temp.seg.poi.bkg[p.bkg.idx]) )
        
        # temp.guide.ll <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx,],
        #                               size = in.data.totals[[repl]][temp.guide.idx],
        #                               alpha = hyper$alpha1[[repl]][temp.guide.idx,], log = T) +
        #                        log(1 - temp.seg.poi.bkg) ) # temp.seg.poi.bkg
        # 
        # temp.guide.ll.bkg <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx,],
        #                                   size = in.data.totals[[repl]][temp.guide.idx],
        #                                   alpha = hyper$alpha0[[repl]][temp.guide.idx,], log = T) +
        #                            log(temp.seg.poi.bkg) ) # temp.seg.poi.bkg
        
        # temp.with.bkg.idx <- which(temp.seg.dist < 1)
        # 
        # if(length(temp.with.bkg.idx) > 0){
        #   temp.guide.ll.bkg <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx[temp.with.bkg.idx],],
        #                                     size = in.data.totals[[repl]][temp.guide.idx[temp.with.bkg.idx]],
        #                                     alpha = hyper$alpha1[[repl]][temp.guide.idx[temp.with.bkg.idx],], log = T) + 
        #                              log(1 - temp.seg.dist[temp.with.bkg.idx])) # 1 - 
        # }
        # 
        # temp.guide.ll <- sum(sg.ll.list[[repl]][temp.guide.idx,2])
        temp.nonGuide.ll <- sum(sg.ll.list[[repl]][temp.nonGuide.idx,1])
        # temp.nonGuide.ll <- sum(sg.ll.list[[repl]][temp.nonGuide.idx])
        
        # likelihood of this continous stretch of bins to contain a regulatory element
        # temp.stretch.ll <- temp.stretch.ll + addlogs(temp.guide.ll, temp.guide.ll.bkg) + temp.nonGuide.ll + log(dgeom(1 + ns, geom.p) / geom.norm.factr)
        # temp.stretch.ll <- temp.stretch.ll + addlogs(temp.fs.ll, temp.bkg.ll) + temp.nonGuide.ll + log(dgeom(1 + ns, geom.p) / geom.norm.factr)
        temp.stretch.ll <- temp.stretch.ll + sum(temp.guide.ll) + temp.nonGuide.ll + log(dgeom(1 + ns, geom.p) / geom.norm.factr)
        
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



#' @title Compute the log likelihood of each possible configuration of the placement of a fs. Take into account distance between segment and the guide target site
#' @param hyper: list, hyperparameters
#' @param in.data.list: list, each element contains a replicate, the total column has already been removed
#' @param in.data.totals: list, each element contains the per-guide totals for a replicate
#' @param seg.to.guide.lst: list, each element is a segment, contains a list with guide_idx, nonGuide_idx
#' @param next.guide.lst: list that contains the index of $next_guide_idx and next_nonGuide_idx
#' @param nr.segs: max number of segemnts that can be combined
#' @param geom.p: probability of the geometic distribution used to calculate the probability of there being x segments in the regulatory element
#' @param sg.ll.df: data frame, each column in the segment log-likelihood based on the PP
#' @return log likelihood for each segment
#' @export estimate_seg_fs_AoE_pp()

estimate_seg_fs_AoE_pp <- function(hyper, in.data.list, in.data.totals, guide.to.seg.lst,
                               seg.to.guide.lst, next.guide.lst, nr.segs = 10, geom.p = 0.1,
                               sg.ll.df, sg.pp.ll.df) {
  
  # set up a list, the length of the number of segments
  # each element contains a vector, which contains the ll of that segment being overlapped by a regulatory element
  segment.ll.list <- list()
  
  total.segs <- length(seg.to.guide.lst)
  region.priors <- rep(1/total.segs, total.segs)
  geom.norm.factr <- pgeom(nr.segs, geom.p)

  total.ll <- log(0)
  
  # make an initial pass with one-segment length RE placements
  for(seg in 1:total.segs){
    temp.segment.ll <- 0
    # temp.guide.idx<- seg.to.guide.lst[[seg]]$guide_idx
    # temp.nonGuide.idx<- seg.to.guide.lst[[seg]]$nonGuide_idx
    
    # for each guide, identify what segments we're processing (what index matches current segment)
    # return the max distance / probability of all segments covered
    # temp.seg.dist <- unlist(lapply(guide.to.seg.lst[temp.guide.idx], function(x){
    #   adj.seg.idx <- which(x$seg_overlapped %in% seg)
    #   temp.dist.to.seg <- x$dist_to_seg[adj.seg.idx]
    #   max(temp.dist.to.seg)
    # }))
    
    for(repl in 1:ncol(sg.ll.df)){
      # temp.seg.ll <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx,],
      #                               size = in.data.totals[[repl]][temp.guide.idx],
      #                               alpha = hyper$alpha1[[repl]][temp.guide.idx,], log = T) + log(temp.seg.dist))
      temp.seg.ll <- sg.pp.ll.df[seg, repl]
      temp.seg.noFS.ll <- sg.ll.df[seg, repl]
      temp.segment.ll <- temp.segment.ll + (temp.seg.ll - temp.seg.noFS.ll) + log(dgeom(1, geom.p) / geom.norm.factr)

      # temp.nonSeg.ll <- sg.ll.df[-seg, repl]
      # temp.segment.ll <- temp.segment.ll + addlogs_vectorized(c(temp.seg.ll, temp.nonSeg.ll)) + log(dgeom(1, geom.p) / geom.norm.factr)
    }
    segment.ll.list[[seg]] <- temp.segment.ll
    total.ll <- addlogs(total.ll, temp.segment.ll)
  }
  
  # all.guide.idx <- c(1:nrow(sg.ll.df)) # helps to differentitate between guides for null vs alternative
  
  # now make a pass through all 2 to nr.segs length segments
  if(nr.segs > 1){
    for(ns in 1:(nr.segs-1) ){
      
      for(seg in 1:(total.segs - nr.segs)){
        
        temp.stretch.ll <- 0
        temp.stretch.segs <- c(seg:(seg + ns)) # index of segments to be included as part of the stretch and be updated in this iteration
        # temp.guide.idx <- c(seg.to.guide.lst[[seg]]$guide_idx, unlist(next.guide.lst[(seg + 1):(seg + ns)]))
        # temp.nonGuide.idx <- all.guide.idx[-temp.guide.idx]
        
        # temp.seg.dist <- unlist(lapply(guide.to.seg.lst[temp.guide.idx], function(x){
        #   adj.seg.idx <- which(x$seg_overlapped %in% temp.stretch.segs)
        #   temp.dist.to.seg <- x$dist_to_seg[adj.seg.idx]
        #   max(temp.dist.to.seg)
        # }))
        
        for(repl in 1:ncol(sg.ll.df)){
          # temp.seg.ll <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx,],
          #                               size = in.data.totals[[repl]][temp.guide.idx],
          #                               alpha = hyper$alpha1[[repl]][temp.guide.idx,], log = T) + log(temp.seg.dist))
          temp.seg.ll <- addlogs_vectorized(sg.pp.ll.df[temp.stretch.segs,repl])
          temp.seg.noFS.ll <- addlogs_vectorized(sg.ll.df[temp.stretch.segs,repl])
          temp.stretch.ll <- temp.stretch.ll + (temp.seg.ll - temp.seg.noFS.ll) + log(dgeom(1 + ns, geom.p) / geom.norm.factr)
  
          # temp.nonSeg.ll <- sg.ll.df[-temp.stretch.segs,repl]
          # temp.stretch.ll <- temp.stretch.ll + addlogs_vectorized(c(temp.seg.ll, temp.nonSeg.ll)) + log(dgeom(1 + ns, geom.p) / geom.norm.factr)
          
        }
        
        total.ll <- addlogs(total.ll, temp.stretch.ll)
        
        for(span in 1:length(temp.stretch.segs)){
          segment.ll.list[[temp.stretch.segs[span]]] <- c(segment.ll.list[[temp.stretch.segs[span]]], temp.stretch.ll)
        }
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


#' @title Compute the log likelihood of each possible configuration of the placement of a fs. Take into account distance between segment and the guide target site
#' @param hyper: list, hyperparameters
#' @param in.data.list: list, each element contains a replicate, the total column has already been removed
#' @param in.data.totals: list, each element contains the per-guide totals for a replicate
#' @param seg.to.guide.lst: list, each element is a segment, contains a list with guide_idx, nonGuide_idx
#' @param next.guide.lst: list that contains the index of $next_guide_idx and next_nonGuide_idx
#' @param nr.segs: max number of segemnts that can be combined
#' @param geom.p: probability of the geometic distribution used to calculate the probability of there being x segments in the regulatory element
#' @param sg.ll.list: list, each element is a data frame per.guide log likelihood given the posteriors for a replicate and alt. ll (total_guide_ll, alt_only_ll)
#' @return log likelihood for each segment
#' @export estimate_fs_AoE_pp()

estimate_fs_AoE_pp <- function(hyper, in.data.list, in.data.totals, guide.to.seg.lst,
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
  
  initial.seg.dist <- guide.to.seg.lst[[initial.guide.idx]]$dist_to_seg[1] #unlist(lapply(guide.to.seg.lst[[initial.guide.idx]]$dist_to_seg, function(x){x[1]}))
  
  for(repl in 1:length(sg.ll.list)){ # in.data.list
    
    initial.guide.ll <- sum(ddirmnom(in.data.list[[repl]][initial.guide.idx,],
                                     size = in.data.totals[[repl]][initial.guide.idx],
                                     alpha = hyper$alpha1[[repl]][initial.guide.idx,], log = T) + log(initial.seg.dist))
    initial.guide.ll.bkg <- 0
    init.with.bkg.idx <- which(initial.seg.dist < 1)
    
    if(length(init.with.bkg.idx) > 0){
      initial.guide.ll.bkg <- sum(ddirmnom(in.data.list[[repl]][initial.guide.idx[init.with.bkg.idx],],
                                           size = in.data.totals[[repl]][initial.guide.idx[init.with.bkg.idx]],
                                           alpha = hyper$alpha0[[repl]][initial.guide.idx[init.with.bkg.idx],], log = T) + 
                                    log(1 - initial.seg.dist[init.with.bkg.idx])) # 1 - 
    }
    
    # initial.guide.ll <- sum(sg.ll.list[[repl]][initial.guide.idx,2])
    
    
    initial.nonGuide.ll <- sum(sg.ll.list[[repl]][initial.nonGuide.idx,1])
    # initial.nonGuide.ll <- sum(sg.ll.list[[repl]][initial.nonGuide.idx])
    initial.segment.ll <- initial.segment.ll + initial.guide.ll + initial.nonGuide.ll + log(dgeom(1, geom.p) / geom.norm.factr) + initial.guide.ll.bkg
  }
  
  segment.ll.list[[1]] <- initial.segment.ll
  total.ll <- initial.segment.ll
  
  # make an initial pass with one-segment length RE placements
  for(seg in 2:total.segs){
    temp.segment.ll <- 0
    temp.guide.idx<- seg.to.guide.lst[[seg]]$guide_idx
    temp.nonGuide.idx<- seg.to.guide.lst[[seg]]$nonGuide_idx
    # 
    # region.pp <- cumulative.pp[guide.seg.idx.lst[[j]]$seg_overlapped]
    # 
    # region.pp.distance.adj <- region.pp * guide.seg.idx.lst[[j]]$dist_to_seg
    # 
    
    # for each guide, identify what segments we're processing (what index matches current segment)
    # return the max distance / probability of all segments covered
    temp.seg.dist <- unlist(lapply(guide.to.seg.lst[temp.guide.idx], function(x){
      adj.seg.idx <- which(x$seg_overlapped %in% seg)
      temp.dist.to.seg <- x$dist_to_seg[adj.seg.idx]
      max(temp.dist.to.seg)
      }))
    
    for(repl in 1:length(sg.ll.list)){ # in.data.list
      
      temp.guide.ll <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx,],
                                    size = in.data.totals[[repl]][temp.guide.idx],
                                    alpha = hyper$alpha1[[repl]][temp.guide.idx,], log = T) + log(temp.seg.dist))
      temp.guide.ll.bkg <- 0
      
      temp.with.bkg.idx <- which(temp.seg.dist < 1)
      
      if(length(temp.with.bkg.idx) > 0){
        temp.guide.ll.bkg <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx[temp.with.bkg.idx],],
                                          size = in.data.totals[[repl]][temp.guide.idx[temp.with.bkg.idx]],
                                          alpha = hyper$alpha1[[repl]][temp.guide.idx[temp.with.bkg.idx],], log = T) + 
                                   log(1 - temp.seg.dist[temp.with.bkg.idx])) # 1 - 
      }
      
      # temp.guide.ll <- sum(sg.ll.list[[repl]][temp.guide.idx,2])
      temp.nonGuide.ll <- sum(sg.ll.list[[repl]][temp.nonGuide.idx,1])
      # temp.nonGuide.ll <- sum(sg.ll.list[[repl]][temp.nonGuide.idx])
      temp.segment.ll <- temp.segment.ll + temp.guide.ll + temp.nonGuide.ll + log(dgeom(1, geom.p) / geom.norm.factr) + temp.guide.ll.bkg
    }
    segment.ll.list[[seg]] <- temp.segment.ll
    total.ll <- addlogs(total.ll, temp.segment.ll)
  }
  
  # all.guide.idx <- c(1:length(sg.ll.list[[1]][,1])) # helps to differentitate between guides for null vs alternative
  all.guide.idx <- c(1:length(sg.ll.list[[1]])) # helps to differentitate between guides for null vs alternative
  
  # now make a pass through all 2 to nr.segs length segments
  for(ns in 1:(nr.segs-1) ){
    
    for(seg in 1:(total.segs - nr.segs)){
      
      temp.stretch.ll <- 0
      temp.stretch.segs <- c(seg:(seg + ns)) # index of segments to be included as part of the stretch and be updated in this iteration
      temp.guide.idx <- c(seg.to.guide.lst[[seg]]$guide_idx, unlist(next.guide.lst[(seg + 1):(seg + ns)]))
      temp.nonGuide.idx <- all.guide.idx[-temp.guide.idx]
      
      temp.seg.dist <- unlist(lapply(guide.to.seg.lst[temp.guide.idx], function(x){
        adj.seg.idx <- which(x$seg_overlapped %in% temp.stretch.segs)
        temp.dist.to.seg <- x$dist_to_seg[adj.seg.idx]
        max(temp.dist.to.seg)
      }))
      
      for(repl in 1:length(sg.ll.list)){ # in.data.list
        temp.guide.ll <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx,],
                                      size = in.data.totals[[repl]][temp.guide.idx],
                                      alpha = hyper$alpha1[[repl]][temp.guide.idx,], log = T) + log(temp.seg.dist))
        
        temp.guide.ll.bkg <- 0
        
        temp.with.bkg.idx <- which(temp.seg.dist < 1)
        
        if(length(temp.with.bkg.idx) > 0){
          temp.guide.ll.bkg <- sum(ddirmnom(in.data.list[[repl]][temp.guide.idx[temp.with.bkg.idx],],
                                            size = in.data.totals[[repl]][temp.guide.idx[temp.with.bkg.idx]],
                                            alpha = hyper$alpha1[[repl]][temp.guide.idx[temp.with.bkg.idx],], log = T) + 
                                     log(1 - temp.seg.dist[temp.with.bkg.idx])) # 1 - 
        }
        
        # temp.guide.ll <- sum(sg.ll.list[[repl]][temp.guide.idx,2])
        temp.nonGuide.ll <- sum(sg.ll.list[[repl]][temp.nonGuide.idx,1])
        # temp.nonGuide.ll <- sum(sg.ll.list[[repl]][temp.nonGuide.idx])
        
        # likelihood of this continous stretch of bins to contain a regulatory element
        temp.stretch.ll <- temp.stretch.ll + temp.guide.ll + temp.nonGuide.ll + log(dgeom(1 + ns, geom.p) / geom.norm.factr) + temp.guide.ll.bkg
        
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

    for(seg in 1:(total.segs - nr.segs)){ # pretty sure nr.segs should be 'ns'....

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

#' @title Numerically stable addition of logs but vectorized. from https://leimao.github.io/blog/LogSumExp/
#' @param loga: first log value
#' @param logb: second log value
#' @return log addition
#' @export addlogs_vectorized()

addlogs_vectorized <- function(logs) {
  max.log <- max(logs)
  max.log + log(sum(exp(logs - max.log)))
}
