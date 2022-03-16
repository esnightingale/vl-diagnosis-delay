plot_spde_mean <- function(res, title = NULL, limits = NULL, trans = FALSE, palopt = "viridis") {
  
  rang <- apply(mesh$loc[, c(1, 2)], 2, range)
  
  proj <- inla.mesh.projector(mesh, 
                              xlim = rang[, 1], 
                              ylim = rang[, 2], 
                              dims = c(300, 300))
  
  df <- expand.grid(x = proj$x, y = proj$y) %>%
    mutate(v = as.vector(inla.mesh.project(proj, res$summary.random$s$mean))) %>%
    st_as_sf(coords = c("x","y"), remove = FALSE, crs = 7759) %>%
    st_intersection(boundary)
  
  if (trans == TRUE){
    df$v <- INLA::inla.link.invlogit(df$v)
  }
  
  
  if (is.null(title)){
    title <- res$.args$formula
  }
  
  
  pal <- viridis::viridis(2)
  p <- ggplot() + 
    geom_raster(data = df, aes(x = x, y = y, fill = v)) +
    gg(boundary.spdf, fill = "transparent") +
    scale_fill_viridis_c(na.value = "transparent", option = palopt, limits = limits, direction = -1, end = 0.9) +
    labs(title = title,
         fill = "",
         x = "", y = "") +
    coord_fixed(ratio = 1) + 
    theme_bw()
  
  return(p)
  
}
