# Plot difference in SPDE (wrt mean and sd) between a given fit and the null model

plot_diff_spde <- function(fit, name, absolute = FALSE, 
                           limit.mean = c(-0.1,0.1), limit.sd = c(-0.02,0.02)){
  
  # Evaluate fitted SPDE from given model over this projection
  mean_i <- inla.mesh.project(proj, fit$summary.random$s$mean)
  sd_i <- inla.mesh.project(proj, fit$summary.random$s$sd)
  
  # Add estimates to data frame and calculate the differences
  df <- df %>%
    dplyr::mutate(mean_i = as.vector(mean_i),
                  sd_i = as.vector(sd_i),
                  diff_mean = (mean_i - mean_base),
                  diff_sd = (sd_i - sd_base))
  
  if(absolute == TRUE){
    df <- df %>%
      mutate(across(starts_with("diff"),abs))
  }
  
  # Map out the differences as a raster
  pal <- viridis::viridis(2)
  gmean <- ggplot() + 
    geom_raster(data = df, aes(x = x, y = y, fill = diff_mean)) +
    gg(boundary.spdf, fill = "transparent") +
    scale_fill_viridis_c(na.value = "transparent", limits = limit.mean) +
    labs(subtitle = "Difference in mean",
         fill = "",
         x = "", y = "") +
    coord_fixed(ratio = 1) + 
    theme_bw()
  
  gsd <- ggplot() + 
    geom_raster(data = df, aes(x = x, y = y, fill = diff_sd)) +
    gg(boundary.spdf, fill = "transparent") +
    scale_fill_viridis_c(na.value = "transparent", limits = limit.sd) +
    labs(subtitle = "Difference in stdev",
         fill = "",
         x = "", y = "") +
    coord_fixed(ratio = 1) + 
    theme_bw()
  
  plots <- cowplot::plot_grid(gmean, gsd)
  
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