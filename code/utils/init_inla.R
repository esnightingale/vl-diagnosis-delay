init_inla <- function(f, data = NULL, data.stack = NULL, family, cpo = FALSE, marg.fitted = FALSE){
  
  if (!is.null(data.stack)){
    fit <- inla(f,
                family = family,
                data = inla.stack.data(data.stack),
                Ntrials = 1,
                control.predictor = list(
                  compute = TRUE, link = 1,
                  A = inla.stack.A(data.stack)),
                control.compute = list(dic = FALSE, 
                                       waic = TRUE, 
                                       config = TRUE,
                                       cpo = cpo,
                                       return.marginals.predictor = marg.fitted),
                control.fixed = list(mean = 0, 
                                     prec = 0.1, 
                                     mean.intercept = 0, 
                                     prec.intercept = 0.1),
                # control.family(hyper = list(size = list(prior="pc.mgamma", param=7))),
                verbose = TRUE)
  }else if (!is.null(data)){
    fit <- inla(f,
                family = family,
                data = data,
                Ntrials = 1,
                control.predictor = list(
                  compute = TRUE, link = 1),
                control.compute = list(dic = FALSE, 
                                       waic = TRUE, 
                                       config = TRUE,
                                       cpo = cpo,
                                       return.marginals.predictor = marg.fitted),
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
