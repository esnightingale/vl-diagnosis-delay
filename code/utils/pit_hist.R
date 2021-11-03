## PIT HISTOGRAM ##
pit_hist <- function(fit, bins = 30, title = ""){
  print(summary(fit$cpo$failure))
  plotdata <- data.frame(pit = fit$cpo$pit)
  return(
    ggplot(plotdata, aes(pit, after_stat(density))) + 
      geom_histogram(fill = "white", col = "black", bins = bins) + 
      geom_hline(yintercept = 1, lty = "dashed", col = "red") + 
      ylim(0,1.5) +
      theme_minimal()+ 
      labs(x = "Probability Integral Transform (PIT)", y = "Density",
           title = title)
  )
}