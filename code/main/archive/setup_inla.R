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

dat <- readRDS(here::here("data","analysisdata_individual.rds")) 

dat.spdf <- dat
coordinates(dat.spdf) <- ~ longitude + latitude

dat <- sf::st_as_sf(dat.spdf) %>%
  st_set_crs(4326) %>%
  mutate(excess_delay = as.numeric(gt90_cdf))

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

################################################################################
################################################################################