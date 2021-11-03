plot_spde <- function(res, title = "Fitted SPDE", limit.mean = NULL, limit.sd = NULL) {
  rang <- apply(mesh$loc[, c(1, 2)], 2, range)
  
  proj <- inla.mesh.projector(mesh, 
                              xlim = rang[, 1], 
                              ylim = rang[, 2], 
                              dims = c(300, 300))
  
  mean_i <- inla.mesh.project(proj, res$summary.random$s$mean)
  sd_i <- inla.mesh.project(proj, res$summary.random$s$sd)
  
  df <- expand.grid(x = proj$x, y = proj$y)
  df$mean_i <- as.vector(mean_i)
  df$sd_i <- as.vector(sd_i)
  
  pal <- viridis::viridis(2)
  gmean <- ggplot() + 
    geom_raster(data = df, aes(x = x, y = y, fill = mean_i)) +
    gg(boundary.spdf, fill = "transparent") +
    scale_fill_viridis_c(na.value = "transparent", limits = limit.mean) +
    labs(subtitle = "Mean",
         fill = "",
         x = "", y = "") +
    # scale_fill_gradient2(na.value = "transparent", low = pal[1], mid = "white", high = pal[2]) +
    coord_fixed(ratio = 1) + 
    theme_bw()
  
  gsd <- ggplot() + 
    geom_raster(data = df, aes(x = x, y = y, fill = sd_i)) +
    gg(boundary.spdf, fill = "transparent") +
    scale_fill_viridis_c(na.value = "transparent", limits = limit.sd) +
    labs(subtitle = "Stdev",
         fill = "",
         x = "", y = "") +
    coord_fixed(ratio = 1) + 
    theme_bw()
  
  plots <- cowplot::plot_grid(gmean, gsd)
  
  title <- cowplot::ggdraw() + 
    cowplot::draw_label(
      title,
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
