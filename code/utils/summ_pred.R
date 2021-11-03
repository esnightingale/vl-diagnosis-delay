# Tabulate fitted values alongside observed to calculate model fit summary measures

summ_pred <- function(fit, data.stack){
  
  # Identify training indices
  index <- data.stack$data$index$train
  
  # Extract summary stats of fitted values at these indices
  pred <- data.frame(ll = fit$summary.fitted.values[index, "0.025quant"],
                     med = fit$summary.fitted.values[index, "0.5quant"],
                     ul = fit$summary.fitted.values[index, "0.975quant"],
                     # exc.prob = sapply(fit$marginals.fitted.values,
                     #                   FUN = function(marg){1-inla.pmarginal(q = 30, marginal = marg)}),
                     obs = data.stack$data$data$y[index])
  
  return(pred)
  
}


get_rmse <- function(pred){
  sqrt(mean((pred$med - pred$obs)^2))
}

get_mae <- function(pred){
  mean(abs(pred$med - pred$obs))
}