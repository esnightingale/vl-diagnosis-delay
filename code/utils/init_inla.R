init_inlabru <- function(cmp, formula = NULL, data = NULL, family, cpo = FALSE){

    fit <- inlabru::bru(f,
                        family = family,
                        data = data,
                        Ntrials = 1,
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
    
    
    cmp <- delay ~ Intercept + field(main = coordinates, model = spde)
    fit <- bru(cmp,
               family = family, 
        data = data)

  out <- list(fit = fit, f = f)
  return(out)
  
}
