################################################################################
# Description: Define mesh on case village locations
################################################################################
################################################################################

# Split training and tuning data
# test.idx <- sample(1:nrow(dat.fit), floor(nrow(dat.fit)*0.25))
# 
# dat.train <- dat.fit[-test.idx,]
# dat.test <- dat.fit[test.idx,]
# 
# dat.train.df <- st_drop_geometry(dat.train)
# dat.test.df <- st_drop_geometry(dat.test)

dat.fit <- readRDS(here::here("data/analysis","dat_fit.rds"))

# Pull coordinates
coo <- st_coordinates(dat.fit)

bnd <- inla.nonconvex.hull(coo, convex = 0.3)

# Define a mesh across the points within boundary
mesh <- inla.mesh.2d(
  boundary = bnd,
  # loc = coo,
  offset = c(0.5, 0.3),
  max.edge = c(0.1, 2), # small triangles within the region and large around the edge
  cutoff = 0.05) # avoid constructing many small triangles where points are close together

# Number of vertices
mesh$n

# plot(mesh)
# points(coo, col = "red")

ggplot() +
  gg(mesh) +
  # gg(boundary.spdf) +
  gg(as_Spatial(dat.fit), col = "red", cex = 0.5) +
  coord_fixed()
ggsave(here::here("figures/fit","mesh.png"), height = 5, width = 7, unit = "in")

# saveRDS(mesh, here::here("output","mesh.rds"))

# Define SPDE model from this mesh
spde <- inla.spde2.pcmatern(mesh = mesh, 
                            prior.range = c(0.05, 0.01), # P(practic.range < 0.05) = 0.01
                            prior.sigma = c(1, 0.01), # P(sigma > 1) = 0.01
                            constr = TRUE)
saveRDS(spde, here::here("data/analysis","spde.rds"))

################################################################################
################################################################################