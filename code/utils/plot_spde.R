plot_spde <- function(res, title = NULL,
                      stat = "mean",
                      limit1 = NULL, limit2 = NULL) {
  
  rang <- apply(mesh$loc[, c(1, 2)], 2, range)
  
  proj <- inla.mesh.projector(mesh, 
                              xlim = rang[, 1], 
                              ylim = rang[, 2], 
                              dims = c(300, 300))
  
  df <- expand.grid(x = proj$x, y = proj$y) 
  
  if (stat == "mean"){
  v1 <- inla.mesh.project(proj, res$summary.random$s$mean)
  v2 <- inla.mesh.project(proj, res$summary.random$s$sd)

  df <- df %>%
    mutate(v1 = as.vector(v1),
           v2 = as.vector(v2))
  
  stt <- c("Mean","Stdev")
  
  }else if (stat == "quantile"){
    v1 <- inla.mesh.project(proj, res$summary.random$s$`0.025quant`)
    v2 <- inla.mesh.project(proj, res$summary.random$s$`0.5quant`)
    v3 <- inla.mesh.project(proj, res$summary.random$s$`0.975quant`)

    df <- df %>%
      mutate(v1 = as.vector(v1),
             v2 = as.vector(v2),
             v3 = as.vector(v3))

    stt <- c("2.5%","Median","97.5%")
  }
  
  df <- st_as_sf(df, coords = c("x","y"), remove = FALSE) %>%
    st_intersection(boundary)
  
  pal <- viridis::viridis(2)
  g1 <- ggplot() + 
    geom_raster(data = df, aes(x = x, y = y, fill = v1)) +
    gg(boundary.spdf, fill = "transparent") +
    scale_fill_viridis_c(na.value = "transparent", limits = limit1, direction = -1, end = 0.9) +
    labs(subtitle = stt[1],
         fill = "",
         x = "", y = "") +
    # scale_fill_gradient2(na.value = "transparent", low = pal[1], mid = "white", high = pal[2]) +
    coord_fixed(ratio = 1) + 
    # scale_x_continuous(limits = c(700,1250)) +
    # scale_y_continuous(limits = c(800,1200)) +
    theme_bw()
  
  g2 <- ggplot() + 
    geom_raster(data = df, aes(x = x, y = y, fill = v2)) +
    gg(boundary.spdf, fill = "transparent") +
    scale_fill_viridis_c(na.value = "transparent", limits = limit2, direction = -1, end = 0.9) +
    labs(subtitle = stt[2],
         fill = "",
         x = "", y = "") +
    coord_fixed(ratio = 1) +  
    theme_bw()
  
  plots <- cowplot::plot_grid(g1, g2)
  
  if (stat == "quantile"){
    g3 <- ggplot() +
      geom_raster(data = df, aes(x = x, y = y, fill = v3)) +
      gg(boundary.spdf, fill = "transparent") +
      scale_fill_viridis_c(na.value = "transparent", limits = limit1, direction = -1, end = 0.9) +
      labs(subtitle = stt[3],
           fill = "",
           x = "", y = "") +
      coord_fixed(ratio = 1) + 
      theme_bw()

    plots <- cowplot::plot_grid(g1,g2,g3) #, nrow = 1

  }

  if (is.null(title)){
    title <- res$.args$formula
  }
  
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
