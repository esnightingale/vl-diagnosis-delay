################################################################################
# Set up SPDE, projection matrices and stacks
################################################################################

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
nrow(dat.fit)

#------------------------------------------------------------------------------#

# Projection matrix for validation points
coov <- st_coordinates(dat.val)
Av <- inla.spde.make.A(mesh = mesh, loc = coov)

dim(Av)
nrow(dat.val)

#------------------------------------------------------------------------------#

# Define a grid of points at which to make predictions. 

bnd.sfc <- st_multipoint(bnd$loc) %>%
  st_sfc() %>%
  st_cast("POLYGON")

bnd.sf <- st_sf(geometry = bnd.sfc) %>%
  st_set_crs(4326)

bb <- st_bbox(bnd.sf)
x <- seq(bb[1] - 1, bb[3] + 1, length.out = 200)
y <- seq(bb[2] - 1, bb[4] + 1, length.out = 200)

grid <- st_multipoint(as.matrix(expand.grid(x, y))) %>%
  st_sfc()

coop <- st_sf(geometry = grid) %>%
  st_set_crs(4326) %>%
  st_intersection(bnd.sf) %>%
  st_coordinates()

# Remove L1 var
coop <- coop[,1:2]

# Construct another matrix that projects the continuous GRF to the prediction
# points
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)

################################################################################
################################################################################
