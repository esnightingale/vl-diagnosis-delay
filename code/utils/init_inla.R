init_inla <- function(f, data){
  fit <- inla(f,
              family = "poisson",
              data = data,
              control.predictor = list(compute = TRUE, 
                                       link = 1),
              control.compute = list(dic = TRUE, 
                                     waic = TRUE, 
                                     config = TRUE),
              control.fixed = list(mean = 0, 
                                   prec = 0.1, 
                                   mean.intercept = 0, 
                                   prec.intercept = 0.001))
  
  out <- list(fit = fit, f = f)
  return(out)
}