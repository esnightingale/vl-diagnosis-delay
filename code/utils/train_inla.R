train_inla <- function(formula, fit.init, test.percent = 0.2){
  
  # Split training and tuning data
  test.idx <- sample(1:nrow(dat.fit), floor(nrow(dat.fit)*test.percent))
  
  dat.train <- dat.fit
  dat.train$days_fever[test.idx] <- NA
  
  fit <- inla(formula,
              family = "poisson",
              data = dat.train,
              control.predictor = list(compute = TRUE, 
                                       link = 1),
              control.compute = list(dic = TRUE, 
                                     waic = TRUE, 
                                     config = TRUE),
              control.fixed = list(mean = 0, 
                                   prec = 0.1, 
                                   mean.intercept = 0, 
                                   prec.intercept = 0.001))
  
  out <- list(fit = fit, dat.train = dat.train, test.idx = test.idx)
  return(out)
} 
