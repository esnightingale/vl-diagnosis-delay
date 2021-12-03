
plot_diff_spde_cat <- function(fit, name, 
                               limit.mean = c(-1,1),
                               breaks = c(-Inf, -0.1,-0.05, 0, 0.05, 0.1, Inf),
                               break.labs = c("Reduced","Increased")){
  
  # Evaluate fitted SPDE from given model over this projection
  mean_i <- inla.mesh.project(proj, fit$summary.random$s$mean)
  sd_i <- inla.mesh.project(proj, fit$summary.random$s$sd)
  
  # Add estimates to data frame and calculate the differences
  df <- df %>%
    dplyr::mutate(mean_i = as.vector(mean_i),
                  sd_i = as.vector(sd_i),
                  diff_mean = (mean_i - mean_base),
                  diff_sd = (sd_i - sd_base)) %>%
    dplyr::mutate(diff_mean_cat = cut(diff_mean, breaks = breaks),
                  abs_diff_mean = cut(abs(diff_mean), 6),
                  diff_abs_mean = cut(abs(mean_i) - abs(mean_base), breaks = breaks)
                    # across(starts_with("diff"), function(x) cut(x, breaks = breaks)
                    )
   
  # Map out the values themselves
  pal <- viridis::viridis(2)
  gmean <- ggplot() + 
    geom_raster(data = df, aes(x = x, y = y, fill = mean_i)) +
    scale_fill_viridis_c(na.value = "transparent", limits = limit.mean) +
    gg(boundary.spdf, fill = "transparent") +
    labs(subtitle = "Mean",
         fill = "",
         x = "", y = "") +
    coord_fixed(ratio = 1) + 
    theme_bw()
  
  # Map out the differences
  gdiff <- ggplot() + 
    geom_raster(data = df, aes(x = x, y = y, fill = diff_mean_cat)) +
    scale_fill_viridis_d(na.value = "transparent", option = "magma", na.translate = FALSE)  +
    gg(boundary.spdf, fill = "transparent") +
    labs(subtitle = "Difference from null",
         fill = "",
         x = "", y = "") +
    coord_fixed(ratio = 1) + 
    theme_bw() 

  # gsd <- ggplot() + 
  #   geom_raster(data = df, aes(x = x, y = y, fill = diff_sd))+
  #   scale_fill_viridis_d(na.value = "transparent", na.translate = FALSE)  +
  #   gg(boundary.spdf, fill = "transparent") +
  #   labs(subtitle = "Difference in stdev",
  #        fill = "",
  #        x = "", y = "") +
  #   coord_fixed(ratio = 1) + 
  #   theme_bw()

  plots <- cowplot::plot_grid(gmean, gdiff, rel_heights = c(1,1.1))
  
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
