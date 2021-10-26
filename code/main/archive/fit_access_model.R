################################################################################
# Description: Setup and fit geostatistical model to diagnosis delay data from
# case details forms.
# Define spatial and temporal effect structure only, without covariate effects.
################################################################################
################################################################################

library(INLA)
library(inlabru)
library(raster)
library(sf)
library(ggmap)
library(mapr)

theme_set(theme_minimal())
figdir <- "figures/fit"

dat <- readRDS(here::here("data","analysisdata_individual.rds"))

dat.spdf <- dat
coordinates(dat.spdf) <- ~ longitude + latitude

dat <- sf::st_as_sf(dat.spdf) %>%
  st_set_crs(4326)  %>%
  mutate(excess_delay = as.numeric(gt90_cdf))

# st_transform(crs = st_crs(24380))

# Setup map context
blockmap <- readRDS(here::here("data","geography","bihar_block.rds"))
boundary <- sf::st_union(blockmap)  %>%
  st_transform(24380) 

buff <- st_buffer(boundary, 1e3) %>%
  st_transform(4326) 
plot(buff)
boundary.spdf <- as_Spatial(buff)

extent <- setNames(st_bbox(blockmap), c("left","bottom","right","top"))
bh_lines <- get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)

# All Bihar villages
# villages <- readRDS(here::here("data","geography","village_centroids_2011.rds")) %>%
#   st_transform(crs = st_crs(24380))

# Travel time raster
access <- raster::raster(here::here("data","covariates","diag_facility_travel_time.tif"))

# Need:
# + Village points
# + Bihar boundary
# + 2D Mesh

# Define the model
# Poisson likelihood for the number of days onset-diagnosis, conditional on true rate, lambda.
# OR
# Binomial likelihood for the number of cases with delay > 90 days out of N per village, conditional
# on true risk of excess delay, r.

# Consider only a spatial random field.
# Define the spatial random field as a zero-mean GP with matern covariance function.
# This has a smoothness parameter nu.


################################################################################
# Set up mesh
################################################################################

# Keep only points within boundary
# dat <- st_intersection(dat, boundary.buf)

# Pull coordinates
coo <- st_coordinates(dat)

# Define a mesh across the points

mesh <- inla.mesh.2d(
  loc = coo,
  offset = c(0.2, 1),
  max.edge = c(0.5, 2), # small triangles within the region and large around the edge
  cutoff = 0.01) # avoid constructing many small triangles where points are close together

# Number of vertices
mesh$n

ggplot() +
  gg(mesh) +
  gg(boundary.spdf) +
  gg(as_Spatial(dat), col = "red", cex = 0.5) +
  coord_fixed()
#------------------------------------------------------------------------------#

# Build the SPDE which corresponds to the GP over this mesh, with an integrate
# to zero constraint.
# Alpha is related to the smoothness parameter nu of the Matern covariance
# function:
#   - alpha = nu + d/2, where d = dimensions
#   - Setting nu = 1 gives alpha = 1 + 2/2 = 2
spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)

#------------------------------------------------------------------------------#

# Generate the index set for this SPDE
indexs <- inla.spde.make.index("s", spde$n.spde)
lengths(indexs)

#------------------------------------------------------------------------------#

# Generate a projection matrix A - projects the spatially continuous GRF from
# the observations to the mesh nodes
A <- inla.spde.make.A(mesh = mesh, loc = coo)

dim(A)
nrow(dat)

#------------------------------------------------------------------------------#

# Define a grid of points at which to make predictions. For simplicity, use the
# points from the altitude raster
dp <- rasterToPoints(access)
dim(dp)
# 
# # Lower the resolution for speed
ra <- raster::aggregate(access, fact = 5, fun = median)
# # 
dp <- rasterToPoints(ra)
dim(dp)
# # 
coop <- dp[, c("x", "y")]

# Construct another matrix that projects the continuous GRF to the prediction
# points
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)

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
  data = list(y = dat$excess_delay),
  A = list(1, A),
  effects = list(data.frame(b0 = 1, 
                            traveltime = dat$traveltime), 
                 s = indexs, 
                 s2 = indexs)
)

# stack for prediction stk.p
stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, Ap),
  effects = list(data.frame(b0 = 1, 
                            traveltime = dp[,3]), 
                 s = indexs, 
                 s2 = indexs)
)

# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.p)

#------------------------------------------------------------------------------#

# Model formula
formula <- y ~ 0 + b0 + traveltime + f(s, model = spde) + f(s2, model = "iid") # w/ nugget effect

#------------------------------------------------------------------------------#

# Fitting
res_cov <- inla(formula,
            family = "binomial",
            control.family = list(link = "logit"),
            data = inla.stack.data(stk.full),
            control.predictor = list(
              compute = TRUE, link = 1,
              A = inla.stack.A(stk.full)
            )
)

summary(res)

saveRDS(res,here::here("output","fit_access_bin90.rds"))

#------------------------------------------------------------------------------#

# Plot the results

# Identify indices of summary.fitted.values which correspond to the predictions
# via inla.stack.index tagged with "pred".
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

# Extract summary stats of fitted values at these indices
mean <- res$summary.fitted.values[index, "mean"]
sd <- res$summary.fitted.values[index, "sd"]
ll <- res$summary.fitted.values[index, "0.025quant"]
ul <- res$summary.fitted.values[index, "0.975quant"]

pred <- data.frame(x = coop[,1], y = coop[,2], mean = mean, sd = sd, ll = ll, ul = ul) 
# %>%
#   tidyr::pivot_longer(-x:-y)

ggmap(bh_lines, 
      base_layer = ggplot(data = pred,
                          aes(x = x, y = y, fill = mean, alpha = 1/(sd^2)))) +
  geom_tile() +
  scale_alpha_continuous(trans = "log10") +
  scale_fill_viridis_c(trans = "log10", direction = -1, option = "plasma") +
  labs(x = "", y = "", fill = "Mean", alpha = "Precision") +
  coord_fixed(ratio = 1)

ggsave(here::here(figdir,"access_fit_bin.png"), height = 6, width = 9, units = "in")

# +
#   facet_wrap(~name)

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
pred$excprob30 <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){1 - inla.pmarginal(q = 30, marginal = marg)})
pred$excprob90 <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){1 - inla.pmarginal(q = 90, marginal = marg)})

pred.long <- tidyr::pivot_longer(pred, -x:-ul)
ggmap(bh_lines, 
      base_layer = ggplot(data = pred.long,
                          aes(x = x, y = y, fill = value))) +
  geom_tile() +
  scale_fill_viridis_c(direction = -1, option = "plasma") +
  labs(x = "", y = "", fill = "P(mu > X)") +
  facet_wrap(~name) +
  coord_fixed(ratio = 1)

ggmap(bh_lines, 
      base_layer = ggplot(data = pred,
                          aes(x = x, y = y, fill = excprob30))) +
  geom_tile() +
  scale_fill_viridis_c(direction = -1, option = "plasma") +
  labs(x = "", y = "", fill = "P(mu > 30)") +
  coord_fixed(ratio = 1)
ggsave(here::here(figdir,"base_pois_exc30.png"), height = 6, width = 9, units = "in")

ggmap(bh_lines, 
      base_layer = ggplot(data = pred,
                          aes(x = x, y = y, fill = excprob90))) +
  geom_tile() +
  scale_fill_viridis_c(direction = -1, option = "plasma") +
  labs(x = "", y = "", fill = "P(mu > 90)") +
  coord_fixed(ratio = 1)
ggsave(here::here(figdir,"base_pois_exc90.png"), height = 6, width = 9, units = "in")


# Again define a raster and plot
r_excprob <- raster::rasterize(
  x = coop, y = access, field = excprob,
  fun = mean
)

################################################################################
################################################################################
