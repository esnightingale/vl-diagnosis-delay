################################################################################
# Description: 
# + Define mesh on case village locations (all data) 
# + Specify SPDE on this mesh
# + Generate grid of points across mesh range to predict at
################################################################################
################################################################################

# Define mesh using all data
dat.train <- readRDS(here::here("data/analysis","dat_nona.rds")) 

# Pull coordinates
coo <- st_coordinates(dat.train)

# Define smooth boundary around points
bnd <- inla.nonconvex.hull(coo, convex = 0.3)

# Define a mesh across the points within boundary
mesh <- inla.mesh.2d(
  boundary = bnd,
  loc = coo,
  offset = c(0.2, 0.2),
  max.edge = c(0.1, 2),
  cutoff = 0.01) 

# Number of vertices
mesh$n

ggplot() +
  gg(mesh) +
  # gg(boundary.spdf) +
  # gg(as_Spatial(dat.fit), col = "red", cex = 0.5) +
  coord_fixed()

ggsave(here::here("figures/fit","mesh.png"), height = 7, width = 9, unit = "in", dpi = 320)

saveRDS(mesh, here::here("data/analysis","mesh.rds"))

# Define SPDE model from this mesh 
# Raw variogram suggests range of ~75km - what scale is prior range on?

dat.train %>% 
  st_drop_geometry() %>%
  group_by(v) %>% 
  summarise(vil.mean = mean(log(days_fever))) %>%
  ungroup() %>%
  summarise(mean.all = mean(vil.mean),
            sd.all = sd(vil.mean))
# mean.all sd.all
#     3.63  0.565

spde <- inla.spde2.pcmatern(mesh = mesh, 
                            prior.range = c(10, 0.01), # P(range < 10) = 0.01
                            prior.sigma = c(1, 0.01), # P(sigma > 1) = 0.01
                            constr = TRUE)
saveRDS(spde, here::here("data/analysis","spde.rds"))

#------------------------------------------------------------------------------#
# Define a grid of points at which to make predictions 

# State boundary
boundary <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
  sf::st_transform(7759) %>%
  sf::st_union() %>% 
  st_transform(4326)

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