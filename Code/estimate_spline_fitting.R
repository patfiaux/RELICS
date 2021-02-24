library(splines)

#' @title Function to help estimate the correct spline parameters for RELICS analysis
#' @param input.parameter.file: location of file containing all analysis parameters
#' @param input.parameter.list: RELICS analysis parameters in list format. Default = NULL
#' @param data.file.split: logical, whether the input data is split acorss info and count file or part of the same file
#' @return list of final per-layer posteriors
#' @export spline_fitting()

spline_fitting <- function(input.parameter.file, input.parameter.list = NULL, data.file.split = FALSE, repl, deg.freed = c(3,5,10,15)){
  
  if(length(deg.freed) > 6){
    print('Please select at most 6 different degrees of freedom')
    break()
  } 
  analysis.parameters <- read_parameters(input.parameter.list, input.parameter.file, data.file.split)
  
  # if parameters are missing, then do not run the analysis
  if(is.logical(analysis.parameters)){
    print('Please fill in the missing parameters before running spline_fitting (or RELICS)')
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
  
  # plot_splines_dispersion_f1(data.setup$data, analysis.parameters$hyper_par_components,
  #                            paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName),
  #                            analysis.parameters$nr_disp_bins, analysis.parameters$mean_var_type,
  #                            fs0.label = analysis.parameters$FS0_label, data.setup$guide_info)
  
  plot_splines_dispersion(data.setup$data, analysis.parameters$hyper_par_components,
                          paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName),
                          analysis.parameters$nr_disp_bins, analysis.parameters$mean_var_type,
                          fs0.label = analysis.parameters$FS0_label, data.setup$guide_info, repl, deg.freed)
  
  # plot_counts_vs_dispersion(data.setup$data, analysis.parameters$hyper_par_components,
  #                           paste0(analysis.parameters$out_dir, '/', analysis.parameters$dataName),
  #                           analysis.parameters$nr_disp_bins, analysis.parameters$mean_var_type,
  #                           fs0.label = analysis.parameters$FS0_label, data.setup$guide_info)
  
  
}


#' @title plot the spline fits
#' @param input.parameter.file: location of file containing all analysis parameters
#' @param input.parameter.list: RELICS analysis parameters in list format. Default = NULL
#' @param data.file.split: logical, whether the input data is split acorss info and count file or part of the same file
#' @return list of final per-layer posteriors
#' @export spline_fitting()

plot_splines_dispersion <- function(input.data, input.hypers, out.dir, nr.bins = 20, mean.var.type,
                                    fs0.label, guide.info, input.repl, deg.freed){
  
  # identify guides used for training and not part of background calculation, remove them
  fs0.label.idx <- which(guide.info$label %in% fs0.label)
  
  out.est.hypers <- c()
  if(length(deg.freed) <= 2){
    par(mfrow = c(1,length(deg.freed)))
  } else if(length(deg.freed) <= 4){
    par(mfrow = c(2,2))
  } else if(length(deg.freed) <= 6){
    par(mfrow = c(2,3))
  }
  # par(mfrow = c(1,length(input.data)))
  
  # for(repl in 1:length(input.data)){
    
  temp.counts <- input.data[[input.repl]]
  
  # remove the FS0 labeled data.
  temp.counts <- temp.counts[-fs0.label.idx,]
  temp.total.counts <- temp.counts[,ncol(temp.counts)]
  # estimate the sorting proportions
  # but only for the background. 
  # Due to identifiabiity we only need to specify n-1 sorting proportions, where n is the number of pools
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
  # plot(x = temp.repl.group.means, y = temp.repl.group.disp, xlab = 'mean bin count', ylab = 'bin dispersion', main = 'Bin dispersion estimate')
  
  data.df <- data.frame(disp = temp.repl.group.disp, counts = temp.repl.group.means, stringsAsFactors = F)
  # spline.mdl.0 <- lm(disp ~ bs(counts, df = 2), data = data.df)
  spline.mdl.list <- list()
  for(degFreed in 1:length(deg.freed)){
    spline.mdl.list[[degFreed]] <- lm(disp ~ bs(counts, df = deg.freed[degFreed]), data = data.df)
  }
  
  # spline.mdl.1 <- lm(disp ~ bs(counts, df = 3), data = data.df)
  # spline.mdl.2 <- lm(disp ~ bs(counts, df = 5), data = data.df)
  # spline.mdl.3 <- lm(disp ~ bs(counts, df = 10), data = data.df)
  # spline.mdl.4 <- lm(disp ~ bs(counts, df = 15), data = data.df)
  
  vals.df <- data.frame(counts = temp.total.counts, stringsAsFactors = F)
  spline.mdl.predict.lst <- list()
  spline.mdl.cols <- c('red', 'orange','blue','purple','grey','green','cyan')
  for(degFreed in 1:length(deg.freed)){
    spline.mdl.predict.lst[[degFreed]] <- predict(spline.mdl.list[[degFreed]], vals.df) 
    plot(x = temp.repl.group.means, y = temp.repl.group.disp, pch = 19, 
         xlab = 'Counts', ylab = 'Dispersion', main = paste0('Dispersion estimates for df = ', deg.freed[degFreed]))

    points(x = temp.total.counts, y = spline.mdl.predict.lst[[degFreed]], 
           col = spline.mdl.cols[degFreed], pch = 19, cex = 0.5) 

    lgndTxt<-c('black dots = per-bin estimated dispersion', paste0(spline.mdl.cols[degFreed],' dots = per-guide dispersion'))
    if(temp.repl.group.disp[1] < 
       mean(temp.repl.group.disp[(length(temp.repl.group.disp)-min(10, length(temp.repl.group.disp)) ):length(temp.repl.group.disp)])){
      legend('topright',legend=lgndTxt, text.col=c('black',spline.mdl.cols[degFreed]), pch = c(19,19), col = c('black',spline.mdl.cols[degFreed]))
    } else {
      legend('bottomright',legend=lgndTxt, text.col=c('black',spline.mdl.cols[degFreed]), pch = c(19,19), col = c('black',spline.mdl.cols[degFreed]))
    }
  }
  # spline.mdl.0.predict <- predict(spline.mdl.0, vals.df)
  # spline.mdl.1.predict <- predict(spline.mdl.1, vals.df) 
  # spline.mdl.2.predict <- predict(spline.mdl.2, vals.df)
  # spline.mdl.3.predict <- predict(spline.mdl.3, vals.df)
  # spline.mdl.4.predict <- predict(spline.mdl.4, vals.df)
  
  # With the predictions we can see how well they actually get modeled
  # points(x = temp.total.counts, y = spline.mdl.0.predict, col = 'darkgreen') 
  # points(x = temp.total.counts, y = spline.mdl.1.predict, col = 'red') 
  # points(x = temp.total.counts, y = spline.mdl.2.predict, col = 'orange')
  # points(x = temp.total.counts, y = spline.mdl.3.predict, col = 'blue')
  # points(x = temp.total.counts, y = spline.mdl.4.predict, col = 'purple')
    
  # }
  
}
