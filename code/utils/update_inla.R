update_inla <- function(fit.init, stk, interval = c(0.1,0.9)){
  
  fit <- inla(fit.init$f,
              family = "poisson",
              data = stk,
              control.predictor = list(compute = TRUE, 
                                       link = 1),
              control.compute = list(dic = TRUE, 
                                     waic = TRUE, 
                                     config = TRUE),
              control.mode = list(result = fit.init$fit, restart = TRUE),
              control.fixed = list(mean = 0, 
                                   prec = 0.1, 
                                   mean.intercept = 0, 
                                   prec.intercept = 0.001))
  return(out)
} 
