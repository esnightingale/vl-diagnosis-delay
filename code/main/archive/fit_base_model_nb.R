################################################################################
# Description: Setup and fit geostatistical model to diagnosis delay data from
# case details forms.
# Define spatial and temporal effect structure only, without covariate effects.
################################################################################
################################################################################

library(tidyverse)
library(INLA)
library(inlabru)
library(raster)
library(sf)
library(ggmap)
library(mapr)

theme_set(theme_minimal())
figdir <- "figures/fit"

dat <- readRDS(here::here("data/analysis","dat_nona.rds")) %>%
  st_transform(7759)
spde <- readRDS(here::here("data/analysis","spde.rds"))
mesh <- readRDS(here::here("data/analysis","mesh.rds"))
coop <- readRDS(here::here("data/analysis","coop.rds"))

# # Setup map context
blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) 

# boundary <- sf::st_union(blockmap)
# boundary.spdf <- as_Spatial(boundary)

extent <- setNames(st_bbox(blockmap), c("left","bottom","right","top"))
bh_lines <- get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)

# Consider only an intercept and a spatial random field.
# Define the spatial random field as a zero-mean GP with matern covariance function.
# This has a smoothness parameter nu.

################################################################################
# Set up stack
################################################################################

# Generate the index set for the SPDE
indexv <- inla.spde.make.index("v", spde$n.spde)
lengths(indexv)

#------------------------------------------------------------------------------#

# Generate a projection matrix A - projects the spatially continuous GRF from
# the observations to the mesh nodes
coo <- st_coordinates(dat)
A <- inla.spde.make.A(mesh = mesh, loc = coo)

dim(A)
nrow(dat)

# Repeat for nugget/IID effect
# A2 <- inla.spde.make.A(mesh = mesh, loc = coo)

#------------------------------------------------------------------------------#

# Construct another matrix that projects the continuous GRF to the prediction
# points
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)

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
  A = list(1, A),
  effects = list(data.frame(b0 = rep(1, nrow(coo))), 
                 v = indexv)
)

# stack for prediction stk.p
stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, Ap),
  effects = list(data.frame(b0 = rep(1, nrow(coop))), 
                 v = indexv)
)

# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.p)

#------------------------------------------------------------------------------#

# Model formula
formula <- y ~ 0 + b0 + f(v, model = spde) 

#------------------------------------------------------------------------------#

# Fitting
res <- inla(formula,
            family = "nbinomial",
            # control.family = list(link = "log"),
            data = inla.stack.data(stk.full),
            control.predictor = list(
              compute = TRUE, link = 1,
              A = inla.stack.A(stk.full)),
            control.fixed = list(mean = 0, prec = 0.1, 
                                 mean.intercept = 0, prec.intercept = 0.1),
            verbose = TRUE)

summary(res)

autoplot(res)
INLAutils::ggplot_inla_residuals(res, dat$days_fever)

saveRDS(res,here::here("output","fit_base.rds"))

#------------------------------------------------------------------------------#

# Plot the results

# Identify indices of summary.fitted.values which correspond to the predictions
# via inla.stack.index tagged with "pred".
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

# Extract summary stats of fitted values at these indices
med <- res$summary.fitted.values[index, "0.5quant"]
ll <- res$summary.fitted.values[index, "0.025quant"]
ul <- res$summary.fitted.values[index, "0.975quant"]

coop_proj <- coop %>%
  st_multipoint() %>%
  st_sfc(crs = 7759) %>%
  st_transform(4326) %>%
  st_coordinates()

pred <- data.frame(x = coop_proj[,1], y = coop_proj[,2], 
                   med = med, 
                   ll = ll, ul = ul) %>%
  pivot_longer(-x:-y)
  

# Map fitted value
ggmap(bh_lines, 
      base_layer = ggplot(data = pred,
                          aes(x = x, y = y, col = value))) +
  geom_point(alpha = 0.5) +
  scale_colour_viridis_c(direction = -1, option = "plasma") +
  labs(x = "", y = "", col = "Median") +
  coord_fixed(ratio = 1) +
  facet_wrap(~name)

ggsave(here::here(figdir,"base_fit.png"), height = 3, width = 10, units = "in")

#------------------------------------------------------------------------------#

# Exceedance probabilities

# Again at prediction locations
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

get_excprob <- function(marg){
  # exponentiate marginal distribution
  tmarg <- inla.tmarginal(exp, marg)
  prob <- 1 - inla.pmarginal(q = 30, marginal = tmarg)
  return(prob)
}

# Extract fitted marginals
# Calculate probability of exceeding 30 days from this marginal distribution
pred$excprob <- sapply(res$marginals.fitted.values[index],
                  FUN = get_excprob)

ggmap(bh_lines, 
      base_layer = ggplot(data = pred,
                          aes(x = x, y = y, fill = excprob))) +
  geom_tile() +
  scale_fill_viridis_c(direction = 1, option = "plasma") +
  labs(x = "", y = "", fill = "P(delay > 30)") +
  coord_fixed(ratio = 1) -> map_exc
map_exc

ggsave(here::here(figdir,"base_nb_exc30.png"), map_exc, height = 6, width = 9, units = "in")


# Again define a raster and plot
r_excprob <- raster::rasterize(
  x = coop, y = access, field = excprob,
  fun = mean
)

################################################################################
################################################################################
