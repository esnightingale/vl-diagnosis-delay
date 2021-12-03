################################################################################
# Description: 
# + Define mesh on case village locations (all data) 
# + Specify SPDE on this mesh
# + Generate grid of points across mesh range to predict at
################################################################################
################################################################################

# Bihar state boundary
boundary <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
  sf::st_transform(7759) %>%
  sf::st_union() 
# %>% st_transform(4326)

# Define mesh using all data
dat.train <- readRDS(here::here("data/analysis","dat_nona.rds")) %>%
  sf::st_transform(7759)

# Pull village coordinates
coo <- sf::st_coordinates(dat.train)

# Calculate average distance between neighbouring villages to inform prior range
nn_dist <- nndist(unique(coo), k = 1)
summary(nn_dist)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 28.22   768.03  1510.48  2131.17  2727.69 50487.95

# Distance to second nearest neighbour:
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 111.5  1530.6  2547.5  3271.5  4059.4 52934.0 

# Distance to 4th nearest neighbour
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 533.4  2641.1  3995.8  4884.1  5774.8 55516.7 

# Average of 1.5-2km to nearest neighbour, minimum of 30m and 
# maximum 50km. 2.5-3.5km to second nearest.

  
# Define smooth boundary around points
bnd <- inla.nonconvex.hull(coo, convex = -0.1)
# plot(bnd)

mesh <- inla.mesh.create(
  loc = coo,
  boundary = list(bnd),
  refine = list(max.edge = 2e4, # 20km
                min.angle = 28,
                cutoff = 2e3)) # 2km mean distance to nearest neighbouring village

# Number of vertices
mesh$n

ggplot() +
  gg(mesh) +
  gg(as_Spatial(boundary)) +
  # gg(as_Spatial(dat.fit), col = "red", cex = 0.5) +
  coord_fixed() +
  labs(x = "", y = "") +
  theme(panel.grid = element_blank(),
        axis.text = element_blank()) -> plot_mesh
plot_mesh

out <- inla.mesh.assessment(mesh,
                            spatial.range = 5e5, #50km
                            alpha = 2,
                            dims = c(200, 200))

ggplot() + 
  gg(out, aes(color = sd.dev)) + 
  coord_equal() +
  scale_color_gradient(limits = range(out$sd.dev, na.rm = TRUE)) -> assess_mesh
assess_mesh 

ggsave(here::here("figures/fit","mesh.png"), plot_mesh, height = 7, width = 9, unit = "in", dpi = 320)
ggsave(here::here("figures/fit","mesh_sd.png"), assess_mesh, height = 7, width = 9, unit = "in", dpi = 320)

saveRDS(mesh, here::here("data/analysis","mesh.rds"))

# ---------------------------------------------------------------------------- #
# Define SPDE model from this mesh 

# Range of similarity between villages in terms of diagnosis delay, i.e. due to 
# targeted surveillance, unlikely to be smaller than between nearest couple of 
# neighbours => prior P(range < 5km) = 0.01

# Prior sigma on log scale, so multiplicative - relative size of peaks and troughs
# to average in SPDE.
# If average delay is around 30-40 days, unlikely to be greater than 50% variation above
# or below, after accounting for covariates

spde <- inla.spde2.pcmatern(mesh = mesh, 
                            prior.range = c(5000, 0.01), # P(range < U) = a
                            prior.sigma = c(1.5, 0.1), # P(sigma > U) = a
                            constr = TRUE)
saveRDS(spde, here::here("data/analysis","spde.rds"))

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
  # Bihar state boundary
  st_intersection(boundary) %>%
  # Model estimation boundary
  st_intersection(bnd.sf) %>%
  st_coordinates()

# Remove L1 var
coop <- coop[,1:2]

saveRDS(coop, here::here("data/analysis","coop.rds"))

################################################################################
################################################################################