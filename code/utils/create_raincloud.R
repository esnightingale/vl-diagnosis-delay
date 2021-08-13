create_raincloud <- function(data, xvar, yvar, xlab, ylab, y_trans = "identity",
                             col_by = NA, drop_na = FALSE, adj = 1) {
  
  if (drop_na == TRUE) {
    data <- filter(data, !is.na(!!sym(xvar)))
  }
  
  data %>%
    ggplot(aes(!!sym(xvar), !!sym(yvar))) +
    geom_jitter(alpha = 0.5, cex = 0.8, height = 0, width = 0.2, aes(col = !!sym(col_by))) +
    geom_boxplot(alpha = 0, width = 0.4, outlier.shape = NA) + 
    ggdist::stat_halfeye(
      aes(fill = !!sym(col_by)),
      adjust = adj, 
      width = .6, 
      justification = -0.4, 
      .width = 0, 
      point_colour = NA) +
  scale_y_continuous(trans = y_trans) + 
  guides(col = "none", fill = "none") +
  # scale_colour_viridis_d(option = "plasma") +
  # scale_fill_viridis_d(option = "plasma") +
  labs(x = xlab, y = ylab)
  
}

