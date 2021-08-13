################################################################################
# Description: Define mesh on case village locations
################################################################################
################################################################################

# Pull coordinates
coo <- st_coordinates(dat.fit)

bnd <- inla.nonconvex.hull(coo)

# Define a mesh across the points
# 
# mesh <- inla.mesh.2d( #mesh_nobnd
#   loc = coo,
#   offset = c(0.1, 1),
#   max.edge = c(0.5, 5), # small triangles within the region and large around the edge
#   cutoff = 0.01) # avoid constructing many small triangles where points are close together
# 
# # Number of vertices
# mesh$n
# 
# ggplot() +
#   gg(mesh) +
#   gg(boundary.spdf) +
#   gg(as_Spatial(dat), col = "red", cex = 0.5) +
#   coord_fixed()

# With boundary
mesh <- inla.mesh.2d(
  boundary = bnd,
  loc = coo,
  offset = c(0.001, 0.5),
  max.edge = c(0.2, 2), # small triangles within the region and large around the edge
  cutoff = 0.01) # avoid constructing many small triangles where points are close together

# Number of vertices
mesh$n

plot(mesh)
points(coo, col = "red")

################################################################################
################################################################################