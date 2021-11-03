plot_vgm <- function(values, loc, title = NULL){
  
  vgm <- variogram(values ~ 1,
                   loc, 
                   cressie = TRUE)

  fit.vgm <- fit.variogram(vgm, vgm("Mat"), fit.kappa = TRUE)
  print(fit.vgm)

  plot(vgm,
       fit.vgm,
       xlab = "Distance", main = title) -> plot
  
  return(plot)

}
