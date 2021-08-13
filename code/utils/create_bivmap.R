
create_bivmap <- function(map, data, xvar, yvar, xlab, ylab,  pal = "DkBlue", style = "quantile", title = ""){
  
  plotdat <- biscale::bi_class(data,
                               x = !!sym(xvar), y = !!sym(yvar), 
                               style = style, dim = 3) 
  
  ggplot() +
    geom_sf(data = map, fill = "white") +
    geom_sf(data = plotdat, aes(fill = bi_class), show.legend = FALSE) +
    biscale::bi_scale_fill(pal = pal, dim = 3, na.value = "white") +
    biscale::bi_theme() +
    ggtitle(title) -> biv_map
  
  legend <- biscale::bi_legend(pal = pal,
                               dim = 3,
                               xlab = xlab,
                               ylab = ylab,
                               size = 10)
  
  final <- cowplot::ggdraw() +
    cowplot::draw_plot(biv_map, 0, 0, 1, 1) +
    cowplot::draw_plot(legend, 0.6, 0.7, 0.3, 0.3) 
  
  return(final)
  
}
