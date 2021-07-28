
# Organise the data, projection matrices and fixed/random effects into a "stack".
# Need separate ones for estimation at observed locations and prediction across
# defined grid

# Fixed effects are intercept
# Random effect is spatial GRF

# List A contains 1 to indicate fixed effects are mapped one to one, and a projection
# matrix for the random effects

# Prediction stack has the response set as NA

# stack for estimation stk.e
stk.e <- inla.stack(
  tag = "est", 
  data = list(y = dat$excess_delay),
  A = list(1, A),
  effects = list(data.frame(b0 = rep(1, nrow(coo))), s = indexs)
)

# stack for prediction stk.p
stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, Ap),
  effects = list(data.frame(b0 = rep(1, nrow(coop))), s = indexs)
)

# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.p)

#------------------------------------------------------------------------------#

# Model formula
formula <- y ~ 0 + b0 + f(s, model = spde) # + f(s2, model = "iid") # w/ nugget effect

#------------------------------------------------------------------------------#

# Fitting
res <- inla(formula,
            family = "binomial",
            control.family = list(link = "logit"),
            data = inla.stack.data(stk.full),
            control.predictor = list(
              compute = TRUE, link = 1,
              A = inla.stack.A(stk.full)
            )
)

summary(res)

autoplot(res)

INLAutils::ggplot_inla_residuals(res, dat$excess_delay)

saveRDS(res,here::here("output","fit_base_bin90.rds"))

#------------------------------------------------------------------------------#

# Plot the results

# Identify indices of summary.fitted.values which correspond to the predictions
# via inla.stack.index tagged with "pred".
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

# Extract summary stats of fitted values at these indices
med <- res$summary.fitted.values[index, "0.5quant"]
ll <- res$summary.fitted.values[index, "0.025quant"]
ul <- res$summary.fitted.values[index, "0.975quant"]

pred <- data.frame(x = coop[,1], y = coop[,2], med = med, ll = ll, ul = ul, iqr = ul - ll) 
# %>%
#   tidyr::pivot_longer(-x:-y)

ggmap(bh_lines, 
      base_layer = ggplot(data = pred,
                          aes(x = x, y = y, fill = med, alpha = 1/iqr))) + 
  geom_tile() +
  scale_alpha_continuous(trans = "log10") +
  scale_fill_viridis_c(trans = "log10", direction = -1, option = "plasma") +
  labs(x = "", y = "", fill = "Median", alpha = "1/IQR") +
  coord_fixed(ratio = 1)

ggsave(here::here(figdir,"base_fit_bin.png"), height = 6, width = 9, units = "in")

#------------------------------------------------------------------------------#

# Define a raster of these values to plot
r_mean <- raster::rasterize(
  x = coop, y = access, field = mean, 
)

ggmap(bh_lines, 
      base_layer = ggplot(data = r_mean, aes(x = x, y = y, fill = fct_elevation_2))) +
  geom_raster() +
  labs(x = "", y = "", col = "Mean") 


pal <- colorNumeric("viridis", c(0,530), na.color = "transparent")

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_mean, colors = pal) %>%
  addLegend("bottomright",
            pal = pal,
            values = values(r_mean), title = "Mean"
  ) %>%
  addScaleBar(position = c("bottomleft"))

#------------------------------------------------------------------------------#

# Exceedance probabilities

# Again at prediction locations
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

# Extract fitted marginals
# Calculate probability of exceeding 0.2 from this marginal distribution
pred$excprob <- sapply(res$marginals.fitted.values[index],
                       FUN = function(marg){1 - inla.pmarginal(q = 0.10, marginal = marg)})

ggmap(bh_lines, 
      base_layer = ggplot(data = pred,
                          aes(x = x, y = y, fill = excprob))) +
  geom_tile() +
  scale_fill_viridis_c(direction = 1, option = "plasma") +
  labs(x = "", y = "", fill = "P(risk > 10%)") +
  coord_fixed(ratio = 1) -> map_exc

ggsave(here::here(figdir,"base_bin_exc25.png"), map_exc, height = 6, width = 9, units = "in")


# Again define a raster and plot
r_excprob <- raster::rasterize(
  x = coop, y = access, field = excprob,
  fun = mean
)

################################################################################
################################################################################
