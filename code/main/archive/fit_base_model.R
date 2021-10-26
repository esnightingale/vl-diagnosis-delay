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

dat <- readRDS(here::here("data","analysisdata_individual.rds")) %>%
  st_set_crs(4326) %>%
  mutate(excess_delay = as.numeric(gt30))
# 
# dat.spdf <- dat
# coordinates(dat.spdf) <- ~ longitude + latitude
# 
# dat <- sf::st_as_sf(dat.spdf) 

# st_transform(crs = st_crs(24380))

# Setup map context
blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
  st_transform(4326)
boundary <- sf::st_union(blockmap)  
boundary.spdf <- as_Spatial(boundary)

extent <- setNames(st_bbox(blockmap), c("left","bottom","right","top"))
bh_lines <- get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)

# Consider only a spatial random field.
# Define the spatial random field as a zero-mean GP with matern covariance function.
# This has a smoothness parameter nu.


################################################################################
# Set up mesh
################################################################################

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

# Define a grid of points at which to make predictions. 

bb <- st_bbox(boundary)
x <- seq(bb[1] - 1, bb[3] + 1, length.out = 200)
y <- seq(bb[2] - 1, bb[4] + 1, length.out = 200)

grid <- st_multipoint(as.matrix(expand.grid(x, y))) %>%
  st_sfc() 

coop <- st_sf(geometry = grid) %>%
  st_set_crs(4326) %>%
  st_intersection(boundary) %>%
  st_coordinates() 

# Remove L1 var
coop <- coop[,1:2]

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
  data = list(y = dat$excess_delay),
  A = list(1, 1, A),
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

saveRDS(res,here::here("output","fit_base_bin30.rds"))

#------------------------------------------------------------------------------#

# Plot the results

# Identify indices of summary.fitted.values which correspond to the predictions
# via inla.stack.index tagged with "pred".
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

# Extract summary stats of fitted values at these indices
mean <- res$summary.fitted.values[index, "mean"]
ll <- res$summary.fitted.values[index, "0.025quant"]
ul <- res$summary.fitted.values[index, "0.975quant"]

pred <- data.frame(x = coop[,1], y = coop[,2], mean = mean, ll = ll, ul = ul)#, iqr = ul - ll

pred.long <- pred %>%
 tidyr::pivot_longer(-x:-y)

ggmap(bh_lines, 
      base_layer = ggplot(data = pred.long,
                          aes(x = x, y = y, fill = value))) +  #, alpha = 1/iqr
  geom_tile() +
  facet_wrap(~name) +
  # scale_alpha_continuous(trans = "log10") +
  scale_fill_viridis_c(direction = -1, option = "plasma") + #trans = "log10", 
  labs(x = "", y = "", fill = "Mean", alpha = "1/IQR",
       title = "Fitted mean and 2.5-97.5 quantiles for the risk of delay greater than 30 days") +
  coord_fixed(ratio = 1) 

ggsave(here::here(figdir,"base_fit_bin30.png"), height = 4, width = 12, units = "in")

#------------------------------------------------------------------------------#
# 
# # Define a raster of these values to plot
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

#------------------------------------------------------------------------------#

# Exceedance probabilities

# Again at prediction locations
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

# Extract fitted marginals
# Calculate probability of exceeding 0.2 from this marginal distribution
pred$excprob <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){1 - inla.pmarginal(q = 0.50, marginal = marg)})

ggmap(bh_lines, 
      base_layer = ggplot(data = pred,
                          aes(x = x, y = y, fill = excprob))) +
  geom_tile() +
  scale_fill_viridis_c(direction = 1, option = "plasma") +
  labs(x = "", y = "", fill = "P(risk > 50%)") +
  coord_fixed(ratio = 1) -> map_exc
map_exc

ggsave(here::here(figdir,"base_bin30_exc50.png"), map_exc, height = 6, width = 9, units = "in")


# Again define a raster and plot
r_excprob <- raster::rasterize(
  x = coop, y = access, field = excprob,
  fun = mean
)

################################################################################
################################################################################
