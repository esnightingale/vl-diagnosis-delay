plot_obsvpred <- function(fit, stk){
  
  idx <- stk$data$index$train
  
  ggplot(data = data.frame(y = fit$summary.fitted.values$mean[idx],
                           ylo = fit$summary.fitted.values$`0.025`[idx],
                           yhi = fit$summary.fitted.values$`0.975`[idx],
                           x = stk$data$data$y[idx]),
         aes(x, y, ymin = ylo, ymax = yhi)) +
    geom_abline(col = "grey") +
    geom_errorbar(alpha = 0.1, col = "grey") +
    geom_point(alpha = 0.1) +
    geom_smooth() +
    scale_x_continuous(trans = "sqrt") +
    scale_y_continuous(trans = "sqrt") +
    labs(y = "Estimated", x = "Observed") -> p
  
  return(p)
}
