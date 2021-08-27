init_inla <- function(f, data.stack, family){
  fit <- inla(f,
              family = family,
              data = inla.stack.data(data.stack),
              control.predictor = list(
                compute = TRUE, link = 1,
                A = inla.stack.A(data.stack)),
              control.compute = list(dic = TRUE, 
                                     waic = TRUE, 
                                     config = TRUE),
              control.fixed = list(mean = 0, 
                                   prec = 0.1, 
                                   mean.intercept = 0, 
                                   prec.intercept = 0.001),
              verbose = TRUE)
  
  out <- list(fit = fit, f = f)
  return(out)
}
