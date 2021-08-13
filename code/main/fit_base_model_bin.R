################################################################################
# Description: Setup and fit geostatistical model to diagnosis delay data from
# case details forms.
# Define spatial and temporal effect structure only, without covariate effects.
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit"

# Need:
# + Village points
# + Bihar boundary
# + 2D Mesh

# Define the model
# Poisson likelihood for the number of days onset-diagnosis, conditional on true rate, lambda.
# OR
# Binomial likelihood for the number of cases with delay < 30 days out of N per village, conditional
# on true risk of excess delay, r.

# Consider only a spatial random field.
# Define the spatial random field as a zero-mean GP with matern covariance function.
# This has a smoothness parameter nu.

#------------------------------------------------------------------------------#

# Organise the data, projection matrices and fixed/random effects into a "stack".
# Need separate ones for estimation at observed locations and prediction across
# defined grid

# Fixed effects are intercept and altitude
# Random effect is spatial GRF

# List A contains 1 to indicate fixed effects are mapped one to one, and a projection
# matrix for the random effects

# Prediction stack has the response set as NA

# stack for estimation stk.e
stk.e <- inla.stack(
  tag = "est", 
  data = list(y = dat$le30),
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

# Fitting (~13 mins, 95 mins with DIC/WAIC)
res <- inla(formula,
            family = "binomial",
            control.family = list(link = "logit"),
            data = inla.stack.data(stk.full),
            control.predictor = list(
              compute = TRUE, link = 1,
              A = inla.stack.A(stk.full)
            ),
            control.compute = list(dic = TRUE, waic = TRUE),
            verbose = TRUE
)

summary(res)

saveRDS(res,here::here("output","fit_base_bin.rds"))

#------------------------------------------------------------------------------#

# Plot the results

# Identify indices of summary.fitted.values which correspond to the predictions
# via inla.stack.index tagged with "pred".
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

# Extract summary stats of fitted values at these indices
mean <- res$summary.fitted.values[index, "mean"]
ll <- res$summary.fitted.values[index, "0.025quant"]
ul <- res$summary.fitted.values[index, "0.975quant"]

pred <- data.frame(x = coop[,1], y = coop[,2], Mean = mean, Lower = ll, Upper = ul) 

ggmap(bh_lines, 
      base_layer = ggplot(data = pred,
                          aes(x = x, y = y, fill = mean))) + 
  geom_tile() +
  scale_fill_viridis_c(trans = "identity", direction = 1) +
  labs(x = "", y = "", fill = "Mean", alpha = "Precision",
       title = "Fitted probability of diagnosis within 30 days of fever onset",
       subtitle = "Posterior mean") +
  coord_fixed(ratio = 1)

ggsave(here::here(figdir,"base_fit_bin.png"), height = 6, width = 9, units = "in")

pred.long <- pred %>%
  tidyr::pivot_longer(-c("x","y"))

ggmap(bh_lines, 
      base_layer = ggplot(data = pred.long,
                          aes(x = x, y = y, fill = value))) + 
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = 1) +
  labs(title = "Fitted probability of diagnosis within 30 days of fever onset",
       subtitle = "Posterior mean and 2.5-97.5% quantiles",
       x = "", y = "", fill = "P(delay < 30)") +
  facet_wrap(~name) +
  coord_fixed(ratio = 1)

ggsave(here::here(figdir,"base_fit_bin_quants.png"), height = 6, width = 9, units = "in")

#------------------------------------------------------------------------------#

# Exceedance probabilities

# Again at prediction locations
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

# Extract fitted marginals
# Calculate probability of more than 50% chance of early diagnosis from this marginal distribution
pred$excprob50 <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){inla.pmarginal(q = 0.5, marginal = marg)})

ggmap(bh_lines, 
      base_layer = ggplot(data = pred,
                          aes(x = x, y = y, fill = excprob50))) +
  geom_tile() +
  scale_fill_viridis_c(direction = 1, option = "plasma") +
  labs(x = "", y = "", fill = "P(p < 0.5)") +
  coord_fixed(ratio = 1) +
  labs(title = "Predicted probability of long delay (> 30 days) being more likely than short (<= 30 days)")

ggsave(here::here(figdir,"base_bin_exc50.png"), height = 6, width = 9, units = "in")

# Define a raster 
# r_excprob <- raster::rasterize(
#   x = coop, y = access, field = excprob,
#   fun = mean
# )

################################################################################
################################################################################
