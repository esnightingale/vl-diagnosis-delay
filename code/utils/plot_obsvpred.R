plot_obsvpred <- function(fit, stk = NULL, idx = NULL, obs = NULL, title = c("formula","family"), trans = FALSE, smooth = FALSE){
 
  if (title == "formula"){t <- fit$.args$formula
  }else if(title == "family"){t <- fit$.args$family}
  
  if (!is.null(stk)){
    x <- stk$data$data$y[idx]

  }else{
    x <- obs
    idx <- 1:length(obs)
  }
  
  df <- data.frame(y = fit$summary.fitted.values$mean[idx],
                   ylo = fit$summary.fitted.values$`0.025quant`[idx],
                   yhi = fit$summary.fitted.values$`0.975quant`[idx],
                   x = x)
  
  if(fit$.args$family == "nbinomial2"){
    # df <- df %>%
    #   mutate(across(starts_with("y"), function(p) (1-p)/p))
    
    tmargs <- lapply(fit$marginals.fitted.values[idx],
                     function(marg) inla.tmarginal(function(p) (1-p)/p, marg))
    exp_delay <- unlist(lapply(tmargs, 
                               function(marg) inla.emarginal(function(m) m, marg)))
    qi_delay <- lapply(tmargs, 
                       function(marg) inla.qmarginal(c(0.025,0.975), marg))
    
    df <- data.frame(y = exp_delay,
                     ylo = unlist(lapply(qi_delay, function(q) q[1])),
                     yhi =  unlist(lapply(qi_delay, function(q) q[2])),
                     x = x)
  }
  

  r <- round(cor(df$x, df$y, method = "pearson"),4)
  # r <- NULL
  
  if(fit$.args$family == "binomial"){
    df$x <- as.factor(df$x)
    ggplot(data = df,
           aes(x, y)) +
      geom_boxplot() +
      labs(y = "Estimated", x = "Observed",
           subtitle = t) -> p
  } else{
    
    ggplot(data = df,
           aes(x, y, ymin = ylo, ymax = yhi)) +
      geom_abline(col = "grey") +
      geom_errorbar(alpha = 0.1, col = "grey") +
      geom_point(alpha = 0.1) +
      labs(y = "Estimated", x = "Observed",
           subtitle = t,
           caption = paste0("Correlation = ", r)) -> p
  }

  if (trans == TRUE){
    p <- p +
      scale_x_continuous(trans = "sqrt") +
      scale_y_continuous(trans = "sqrt") 
      
  }
  
  if (smooth == TRUE){
    p <- p + 
      geom_smooth()
  }
  
  return(p)
}
