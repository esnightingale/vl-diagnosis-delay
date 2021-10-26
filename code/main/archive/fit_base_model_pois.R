################################################################################
# Description: Setup and fit geostatistical model to diagnosis delay data from
# case details forms.
# Define spatial and temporal effect structure only, without covariate effects.
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit"
outdir <- "output"

#------------------------------------------------------------------------------#

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
  data = list(y = dat$days_fever),
  A = list(1, 1, A),
  effects = list(data.frame(b0 = rep(1, nrow(coo))), i = dat$i, s = indexs)
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
formula <- y ~ 0 + b0 + f(s, model = spde) + f(i, model = "iid") # w/ nugget effect

#------------------------------------------------------------------------------#

# Fitting (~12mins for boundary model, ~8mins without boundary)
res <- inla(formula,
            family = "poisson",
            control.family = list(link = "log"),
            data = inla.stack.data(stk.full),
            control.predictor = list(
              compute = TRUE, link = 1,
              A = inla.stack.A(stk.full)
            ),
            control.compute = list(dic = TRUE, waic = TRUE),
            verbose = TRUE
)

summary(res)

autoplot(res)

INLAutils::ggplot_inla_residuals(res, dat$days_fever)

saveRDS(res,here::here(outdir,"fit_base_pois.rds"))

#------------------------------------------------------------------------------#

# Plot the results

# Identify indices of summary.fitted.values which correspond to the predictions
# via inla.stack.index tagged with "pred".
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

# Extract summary stats of fitted values at these indices
mean <- res$summary.fitted.values[index, "mean"]
ll <- res$summary.fitted.values[index, "0.025quant"]
ul <- res$summary.fitted.values[index, "0.975quant"]

pred <- data.frame(x = coop[,1], y = coop[,2], Mean = mean, Lower = ll, Upper = ul, iqr = ul - ll) 

ggmap(bh_lines, 
      base_layer = ggplot(data = pred,
                          aes(x = x, y = y, fill = Mean))) + # , alpha = 1/iqr
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = -1) +
  labs(title = "Fitted mean of days fever prior to diagnosis", 
       x = "", y = "", fill = "Delay (days)") +
  coord_fixed(ratio = 1)

ggsave(here::here(figdir,"base_fit_pois.png"), height = 6, width = 9, units = "in")


pred.long <- pred %>%
  tidyr::pivot_longer(-c("x","y","iqr"))

ggmap(bh_lines, 
      base_layer = ggplot(data = pred.long,
                          aes(x = x, y = y, fill = value))) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = -1) +
  labs(title = "Fitted mean and 2.5-97.5 quantiles for days fever prior to diagnosis", 
       x = "", y = "", fill = "Delay (days)") +
  facet_wrap(~name) +
  coord_fixed(ratio = 1)

ggsave(here::here(figdir,"base_fit_pois_quants.png"), height = 6, width = 9, units = "in")

#------------------------------------------------------------------------------#

# Exceedance probabilities

# Again at prediction locations
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

# Extract fitted marginals
# Calculate probability of exceeding 30 days from this marginal distribution
pred$excprob <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){1 - inla.pmarginal(q = 30, marginal = marg)})

ggmap(bh_lines, 
      base_layer = ggplot(data = pred,
                          aes(x = x, y = y, fill = excprob))) +
  geom_tile() +
  scale_fill_viridis_c(direction = 1, option = "plasma") + #trans = "log10", 
  labs(title = "Predicted probability of delay exceeding 30 days",
       x = "", y = "", fill = "P(delay > 30)") +
  coord_fixed(ratio = 1) -> map_exc
map_exc

ggsave(here::here(figdir,"base_pois_exc30.png"), map_exc, height = 6, width = 9, units = "in")


#------------------------------------------------------------------------------#
# Define a raster of these values to plot
# r_mean <- raster::rasterize(
#   x = coop, y = access, field = mean, 
# )
# 
# ggmap(bh_lines, 
#       base_layer = ggplot(data = r_mean, aes(x = x, y = y, fill = fct_elevation_2))) +
#   geom_raster() +
#   labs(x = "", y = "", col = "Mean") 
# 
# 
# pal <- colorNumeric("viridis", c(0,530), na.color = "transparent")
# 
# leaflet() %>%
#   addProviderTiles(providers$CartoDB.Positron) %>%
#   addRasterImage(r_mean, colors = pal) %>%
#   addLegend("bottomright",
#             pal = pal,
#             values = values(r_mean), title = "Mean"
#   ) %>%
#   addScaleBar(position = c("bottomleft"))
# 
# 
# # Again define a raster and plot
# r_excprob <- raster::rasterize(
#   x = coop, y = access, field = excprob,
#   fun = mean
# )
################################################################################
################################################################################
