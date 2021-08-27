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
              control.mode = list(result = fit.init, restart = TRUE),
              control.fixed = list(mean = 0, 
                                   prec = 0.1, 
                                   mean.intercept = 0, 
                                   prec.intercept = 0.001))
  
  # Prediction at test obs (for later plotting and RMSE)
  dat.pred <- dat.fit %>%
                slice(test.idx) %>%
                mutate(sq.err = (days_fever - fit$summary.fitted.values[test.idx])) 

  # Fitted SD of IID effects, SD and range of SPDE
  marg.hyperpars <- list(IID.sd = inla.tmarginal(function(x) sqrt(1/x), 
                                                 fit$marginals.hyperpar$`Precision for id`),
                         SPDE.sd = fit$marginals.hyperpar$`Stdev for v`,
                         SPDE.range = fit$marginals.hyperpar$`Range for v`)
          
  out <- list(RMSE = RMSE,
              marg.hyperpars = marg.hyperpars,
              dat.pred = dat.pred, test.idx = test.idx)
  return(out)
} 
