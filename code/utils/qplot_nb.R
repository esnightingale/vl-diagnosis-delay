qplot.nb <- function(fit, obs){
  
  x <- qnbinom(ppoints(obs), 
               size=fit$summary.hyperpar[1,"mean"], 
               mu=fit$summary.fitted.values[1:nrow(dat),"mean"])
  y <- sort(dat$delay_mth)
  
  plot(x, y)
  abline(0, 1)
  
}
