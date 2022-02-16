# Plot observed values versus predicted - raw and exceedance

plot_preds <- function(x, name){
  
  ggplot(x, aes(obs, pred, ymin = ll, ymax = ul)) +
    # geom_errorbar(width = 0.5) +
    geom_jitter(alpha = 0.3) + 
    geom_abline(slope = 1, intercept = 0, lty = "dashed", col = "grey") +
    # scale_x_continuous(trans = "log10") +
    # scale_y_continuous(trans = "log10") +
    labs(x = "Observed delay (days)", y = "Posterior predicted delay",
         title = name)
  
}

plot_exc <- function(x, name){
  
  ggplot(x, aes(exc.obs, exc.prob, col = exc.obs)) +
    geom_boxplot() +
    # geom_hline(yintercept = 0.5, lty = "dashed", col = "grey") + 
    guides(col = "none") + 
    # scale_y_continuous(trans = "logit") +
    # ylim(c(0,1)) +
    labs(x = "Observed delay > 30 days", y = "Fitted exceedance probability",
         title = name) 

}

plot_resids <- function(x, name){
  
  ggplot(x, aes(pred, pred-obs)) +
    geom_jitter(alpha = 0.3) + 
    geom_hline(yintercept = 0, lty = "dashed", col = "grey") +
    labs(x = "Predicted", y = "Error",
         title = name)
  
}
