plot_vgm <- function(values, data, cutoff = NULL, title = NULL, fit.kappa = FALSE, ylim = NULL){
  
  vgm <- variogram(values ~ 1,
                   data, 
                   # cutoff = cutoff,
                   # width = width,
                   cressie = TRUE)

  fit.vgm <- fit.variogram(vgm, vgm("Mat"), fit.kappa = fit.kappa)
  print(fit.vgm)

  if (!is.null(ylim)){
    plot(vgm,
         fit.vgm,
         ylim = ylim,
         xlab = "Distance", main = title) -> plot
  }else{
    plot(vgm,
         fit.vgm,
         xlab = "Distance", main = title) -> plot
  }

  
  return(plot)

}
