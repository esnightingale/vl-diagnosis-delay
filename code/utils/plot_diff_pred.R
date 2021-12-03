plot_diff_pred <- function(fit, name, 
                               limit.mean = c(0,100)
                               # breaks = c(-Inf, -2.5, 0, 2.5, 5, 7.5, 10, Inf)
                           ){
  
  
  # Extract fitted values from this model at prediction points
  mean_i <- fit$summary.fitted.values[pidx, "mean"]
  sd_i <- fit$summary.fitted.values[pidx, "sd"]
  
  # Add estimates to data frame and calculate the differences
  df <- df %>%
    dplyr::mutate(mean_i = as.vector(mean_i),
                  sd_i = as.vector(sd_i),
                  diff_mean = (mean_i - mean_base),
                  diff_sd = (sd_i - sd_base))  %>%
    dplyr::mutate(across(starts_with("diff"), function(x) cut(x, 5))) 
  
  # Map out the values themselves
  pal <- viridis::viridis(2)
  gmean <- ggplot() + 
    geom_raster(data = df, aes(x = x, y = y, fill = mean_i)) +
    scale_fill_viridis_c(na.value = "transparent", option = "magma", limits = limit.mean) +
    gg(boundary.spdf, fill = "transparent") +
    labs(subtitle = "Mean",
         fill = "",
         x = "", y = "") +
    coord_fixed(ratio = 1) + 
    theme_bw()
  
  # Map out the differences
  gdiff <- ggplot() + 
    geom_raster(data = df, aes(x = x, y = y, fill = diff_mean)) +
    scale_fill_viridis_d(na.value = "transparent", na.translate = FALSE)  +
    gg(boundary.spdf, fill = "transparent") +
    labs(subtitle = "Difference from null",
         fill = "",
         x = "", y = "") +
    coord_fixed(ratio = 1) + 
    theme_bw() 

  plots <- cowplot::plot_grid(gmean, gdiff, rel_heights = c(1,1))
  
  title <- cowplot::ggdraw() + 
    cowplot::draw_label(
      name,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  
  cowplot::plot_grid(
    title, plots,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  ) -> p_final
  
  return(p_final)
  
}
