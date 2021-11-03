init_inla <- function(f, data = NULL, data.stack = NULL, family, cpo = FALSE){
  
  if (!is.null(data.stack)){
    fit <- inla(f,
                family = family,
                data = inla.stack.data(data.stack),
                control.predictor = list(
                  compute = TRUE, link = 1,
                  A = inla.stack.A(data.stack)),
                control.compute = list(dic = TRUE, 
                                       waic = TRUE, 
                                       config = TRUE,
                                       cpo = cpo),
                control.fixed = list(mean = 0, 
                                     prec = 0.1, 
                                     mean.intercept = 0, 
                                     prec.intercept = 0.1),
                verbose = TRUE)
  }else if (!is.null(data)){
    fit <- inla(f,
                family = family,
                data = data,
                control.predictor = list(
                  compute = TRUE, link = 1),
                control.compute = list(dic = TRUE, 
                                       waic = TRUE, 
                                       config = TRUE,
                                       cpo = cpo),
                control.fixed = list(mean = 0, 
                                     prec = 0.1, 
                                     mean.intercept = 0, 
                                     prec.intercept = 0.1),
                verbose = TRUE)
  }else{
    print("Provide data.frame or inla.stack")
    return()
  }

  out <- list(fit = fit, f = f)
  return(out)
}
