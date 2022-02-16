################################################################################
# Description: 
# + Define mesh on case village locations (all data) 
# + Specify SPDE on this mesh
# + Generate grid of points across mesh range to predict at
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

# Bihar state boundary
boundary <- readRDS(here::here("data","geography","boundary.rds")) %>%
  st_set_crs(7759)

# Read and merge cleaned analysis data 
dat <- read_data()

# Pull village coordinates
coo <- sf::st_coordinates(dat)

# Calculate average distance between neighbouring villages
nn_dist <- nndist(unique(coo), k = 1)
summary(nn_dist)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.02822  0.76986  1.51595  2.12726  2.72769 50.48795 

# Bihar state is ~ 500km at its widest.

# Define smooth boundary around points
bnd <- inla.nonconvex.hull(coo, convex = -0.1)
plot(bnd)

mesh <- inla.mesh.2d(
  # loc = coo,
  boundary = list(bnd),
  max.edge = c(10,50),
  offset = c(10, 100),
  min.angle = 25,
  cutoff = 2) # ~2km mean distance to nearest neighbouring village

plot(mesh)
# plot_mesh <- autoplot(mesh)
# plot_mesh 

# Number of vertices
mesh$n
# 3674

out <- inla.mesh.assessment(mesh,
                            spatial.range = 30, 
                            alpha = 2,
                            dims = c(500, 500))

ggplot() + 
  gg(out, aes(color = sd.dev)) + 
  coord_equal() +
  # scale_color_gradient(limits = range(out$sd.dev, na.rm = TRUE)) +
  scale_colour_viridis_c() -> assess_mesh #limits = c(0.9,1.1)
assess_mesh  

png(here::here("figures/fit","mesh.png"), height = 9, width = 12, unit = "in", res = 320)
plot(mesh)
dev.off()
ggsave(here::here("figures/fit","mesh_sd.png"), assess_mesh, height = 9, width = 12, unit = "in", dpi = 320)

saveRDS(mesh, here::here("data/analysis","mesh.rds"))

#------------------------------------------------------------------------------#
# Define a grid of points at which to make predictions 

bnd.sfc <- st_multipoint(bnd$loc) %>%
  st_sfc() %>%
  st_cast("POLYGON")

bnd.sf <- st_sf(geometry = bnd.sfc) %>%
  st_set_crs(7759)

bb <- st_bbox(bnd.sf) 
x <- seq(bb[1] - 1, bb[3] + 1, length.out = 200)
y <- seq(bb[2] - 1, bb[4] + 1, length.out = 200)

grid <- st_multipoint(as.matrix(expand.grid(x, y))) %>%
  st_sfc()

coop <- st_sf(geometry = grid) %>%
  st_set_crs(7759) %>%
  # Model estimation boundary
  st_intersection(bnd.sf) %>%
  # Bihar state boundary
  st_intersection(boundary) %>%
  st_coordinates()

# Remove L1 var
coop <- coop[,1:2]

saveRDS(coop, here::here("data/analysis","coop.rds"))

################################################################################
################################################################################