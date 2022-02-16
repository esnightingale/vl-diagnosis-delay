# Plot difference in SPDE (wrt mean and sd) between a given fit and the null model

plot_diff_spde <- function(fit, name, absolute = FALSE, 
                           limit.mean = c(-1,1), 
                           limit.diff = c(-0.5,0.5),
                           limit.sd = c(-0.02,0.02)){
  
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
      dplyr::mutate(across(starts_with("diff"), abs))
      # dplyr::mutate(mean_i = abs(as.vector(mean_i)),
      #               sd_i = abs(as.vector(sd_i)),
      #               diff_mean = (mean_i - abs(mean_base)),
      #               diff_sd = (sd_i - abs(sd_base)))
  }

  
  df <- st_as_sf(df, coords = c("x","y"), remove = FALSE) %>%
    st_intersection(boundary)
  
  # Map out the values themselves
  pal <- viridis::viridis(2)
  gmean <- ggplot() + 
    geom_raster(data = df, aes(x = x, y = y, fill = mean_i)) +
    scale_fill_viridis_c(na.value = "transparent", limits = limit.mean, direction = -1, end = 0.9) +
    gg(boundary.spdf, fill = "transparent") +
    labs(subtitle = "Fitted mean",
         fill = "",
         x = "", y = "") +
    coord_fixed(ratio = 1) + 
    theme_bw()
  
  # Map out the differences
  gdiff <- ggplot() + 
    geom_raster(data = df, aes(x = x, y = y, fill = diff_mean)) +
    scale_fill_viridis_c(option = "plasma", na.value = "transparent", limits = limit.diff) +
    gg(boundary.spdf, fill = "transparent") +
    labs(subtitle = "Difference in mean",
         fill = "",
         x = "", y = "") +
    coord_fixed(ratio = 1) + 
    theme_bw()
  
  gsd <- ggplot() + 
    geom_raster(data = df, aes(x = x, y = y, fill = diff_sd))+
    scale_fill_viridis_c(option = "plasma", na.value = "transparent", limits = limit.sd)  +
    gg(boundary.spdf, fill = "transparent") +
    labs(subtitle = "Difference in stdev",
         fill = "",
         x = "", y = "") +
    coord_fixed(ratio = 1) + 
    theme_bw()

  plots <- cowplot::plot_grid(gmean, gdiff)
  
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
