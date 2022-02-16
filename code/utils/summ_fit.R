# Tabulate fitted values alongside observed to calculate model fit summary measures

summ_fit <- function(fit, data.stack, C = 30){
  
  # Identify training indices
  # index <- data.stack$data$index$train
  index <- 1:data.stack$data$nrow
  
  # Extract summary stats of fitted values at these indices
  pred <- data.frame(ll = fit$summary.fitted.values[index, "0.025quant"],
                     pred = fit$summary.fitted.values[index, "mean"],
                     sd = fit$summary.fitted.values[index, "sd"],
                     med = fit$summary.fitted.values[index, "0.5quant"],
                     ul = fit$summary.fitted.values[index, "0.975quant"],
                     exc.prob = sapply(fit$marginals.fitted.values[index],
                                       FUN = function(marg){1-inla.pmarginal(q = C, marginal = marg)}),
                     obs = data.stack$data$data$y[index]) %>%
    dplyr::mutate(exc.obs = (obs > C))
  
  # Adjust observed for binomial model
  if (fit$.args$family == "binomial"){
    pred$obs <- as.numeric(pred$exc.obs)
  }
  
  return(pred)
  
}


get_mse <- function(pred){
  mean((pred$pred - pred$obs)^2)
}


get_stdmse <- function(pred){
  mean(((pred$pred - pred$obs)^2)/(pred$sd^2))
}


get_mae <- function(pred){
  mean(abs(pred$pred - pred$obs))
}
