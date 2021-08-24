fit_inla <- function(f, data, cpo = FALSE, predictor = FALSE, verbose = TRUE){
  
  if(cpo){
  res <- inla(f,
              family = "poisson",
              control.family = list(link = "log"),
              data = data,
              control.predictor = list(compute = predictor, link = 1),
              control.compute = list(dic = TRUE, 
                                     waic = TRUE, 
                                     cpo = cpo,
                                     config = TRUE),
              control.fixed = list(mean = 0, prec = 0.1, 
                                   mean.intercept = 0, prec.intercept = 0.001),
              verbose = verbose
  )
  }else{
  
  res <- inla(f,
       family = "poisson",
       control.family = list(link = "log"),
       data = data,
       control.predictor = list(compute = predictor, link = 1),
       control.compute = list(dic = TRUE, 
                              waic = TRUE, 
                              cpo = cpo,
                              config = TRUE),
       control.fixed = list(mean = 0, prec = 0.1, 
                            mean.intercept = 0, prec.intercept = 0.001),
       verbose = verbose
  )
  }
  
  return(res)
  
}
