################################################################################
# Set up SPDE and projection matrices
################################################################################

# Build the SPDE which corresponds to the GP over this mesh, with an integrate
# to zero constraint.
# Alpha is related to the smoothness parameter nu of the Matern covariance
# function:
#   - alpha = nu + d/2, where d = dimensions
#   - Setting nu = 1 gives alpha = 1 + 2/2 = 2

# mesh <- readRDS(here::here("output","mesh.rds"))

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
nrow(dat.train)

#------------------------------------------------------------------------------#

# Projection matrix for validation points
coot <- st_coordinates(dat.test)
At <- inla.spde.make.A(mesh = mesh, loc = coot)

dim(At)
nrow(dat.test)

#------------------------------------------------------------------------------#

# Define a grid of points at which to make predictions 
# Extract travel time at these points and form a data frame 

# bnd.sfc <- st_multipoint(bnd$loc) %>%
#   st_sfc() %>%
#   st_cast("POLYGON")
# 
# bnd.sf <- st_sf(geometry = bnd.sfc) %>%
#   st_set_crs(4326)
# 
# bb <- st_bbox(bnd.sf)
# x <- seq(bb[1] - 1, bb[3] + 1, length.out = 200)
# y <- seq(bb[2] - 1, bb[4] + 1, length.out = 200)
# 
# grid <- st_multipoint(as.matrix(expand.grid(x, y))) %>%
#   st_sfc()
# 
# coop <- st_sf(geometry = grid) %>%
#   st_set_crs(4326) %>%
#   # Bihar state boundary
#   st_intersection(boundary) %>%
#   # Model estimation boundary
#   st_intersection(bnd.sf) %>%
#   st_coordinates()
# 
# # Remove L1 var
# coop <- coop[,1:2]
# 
# # Projection matrix for prediction points
# Ap <- inla.spde.make.A(mesh = mesh, loc = coop)
# 
# dat.pred <- data.frame(coop, 
#                        traveltime = raster::extract(access, coop)) %>%
#   dplyr::mutate(traveltime_s = as.numeric(scale(traveltime, center = T)))

################################################################################
################################################################################
