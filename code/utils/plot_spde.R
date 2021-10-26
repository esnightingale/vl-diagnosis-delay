plot_spde <- function(res) {
  rang <- apply(mesh$loc[, c(1, 2)], 2, range)
  
  proj <- inla.mesh.projector(mesh, 
                              xlim = rang[, 1], 
                              ylim = rang[, 2], 
                              dims = c(300, 300))
  
  mean_i <- inla.mesh.project(proj, res$fit$summary.random$v$mean)
  sd_i <- inla.mesh.project(proj, res$fit$summary.random$v$sd)
  
  df <- expand.grid(x = proj$x, y = proj$y)
  df$mean_i <- as.vector(mean_i)
  df$sd_i <- as.vector(sd_i)
  
  pal <- viridis::viridis(2)
  gmean <- ggplot() + 
    geom_raster(data = df, aes(x = x, y = y, fill = mean_i)) +
    gg(boundary.spdf, fill = "transparent") +
    scale_fill_viridis_c(na.value = "transparent") +
    # scale_fill_gradient2(na.value = "transparent", low = pal[1], mid = "white", high = pal[2]) +
    coord_fixed(ratio = 1) + 
    theme_bw()
  
  gsd <- ggplot() + 
    geom_raster(data = df, aes(x = x, y = y, fill = sd_i)) +
    gg(boundary.spdf, fill = "transparent") +
    scale_fill_viridis_c(na.value = "transparent") +
    coord_fixed(ratio = 1) + 
    theme_bw()
  
  return(cowplot::plot_grid(gmean, gsd))
  
}
