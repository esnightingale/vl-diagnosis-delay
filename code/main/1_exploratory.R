################################################################################
# Description: Univariate regression on covariates and assessment of empirical 
# variogram for spatial dependence
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures"

# dat <- readRDS(here::here("data","analysisdata_individual.rds")) %>%
#   dplyr::mutate(vill_inc_2017_gt0 = factor((replace_na(vill_inc_2017,0) > 0), labels = c("No","Yes")),
#                 vill_inc_2017_t = vill_inc_2017*1e3 + 1e-4,
#                 IRS_2017_1 = factor((IRS_2017 != 0), labels = c("No","Yes")),
#                 age_s = as.numeric(scale(age, center = T)),
#                 traveltime_s = as.numeric(scale(traveltime, center = T)),
#                 i = row_number())
# # dat.df <- st_drop_geometry(dat) 
# 
# dat.fit <- dat %>%
#   st_drop_geometry() %>%
#   dplyr::select(days_fever, vil_code, i,
#                 age_s, sex, hiv, marg_caste, occ4_cat, detection, 
#                 prv_tx, conslt_cat, 
#                 latitude, longitude, traveltime_s, 
#                 vill_inc_2017, vill_inc_2017_gt0, vill_inc_2017_t, 
#                 IRS_2017_1, block_endm_2017) %>%
#   drop_na() 

# village <- readRDS(here::here("data","analysisdata_village.rds")) 

# Setup map context
blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
  st_transform(4326)
boundary <- sf::st_union(blockmap)  
boundary.spdf <- as_Spatial(boundary)

extent <- setNames(st_bbox(blockmap), c("left","bottom","right","top"))
bh_lines <- get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)

dat <- readRDS(here::here("data/analysis","dat_nona.rds"))

# ---------------------------------------------------------------------------- #
# Visualise the spatial range of the data and pattern of delays

ggplot() +
  geom_sf(data = blockmap, fill = NA) +
  geom_sf(data = st_jitter(dat), aes(col = days_fever), alpha = 0.8, cex = 0.5) +
  scale_colour_viridis_c(trans = "log2") +
  labs(col = "Delay (days)") + 
  theme(axis.text = element_blank(), panel.grid = element_blank()) +
  annotation_scale(location = "br") +
  annotation_north_arrow(location = "tr", pad_x = unit(2, "cm"), pad_y = unit(1, "cm"))

ggsave(here::here(figdir, "descriptive/fig1.png"), height = 7, width = 9, units = "in")

# ---------------------------------------------------------------------------- #
# Visualise the spatial range of the data and pattern of delays

# Variogram and fit variogram
vgm <- variogram(log(days_fever) ~ 1, dat)
fit.vgm <- fit.variogram(vgm, vgm("Sph"))

plot(vgm, fit.vgm)
# model      psill    range
# 1   Nug 0.29107975  0.00000
# 2   Sph 0.08658465 78.62025

# Kriging

# Define grid of points to interpolate over, based on available locations and state boundary 
bnd <- inla.nonconvex.hull(st_coordinates(dat), convex = 0.3)
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

pred.grid <- st_sf(geometry = grid) %>%
  st_set_crs(4326) %>%
  # Bihar state boundary
  st_intersection(boundary) %>%
  # Model estimation boundary
  st_intersection(bnd.sf) 

krg <- krige(log(days_fever) ~ 1, dat, pred.grid, model = fit.vgm)

#Add estimates to pred.grid
pred.grid$delay.krg <- krg$var1.pred
pred.grid$delay.krg.sd <- sqrt(krg$var1.var)

ggplot() +
  geom_tile(data = pred.grid, aes(fill = delay.krg)) +
  geom_sf(data = blockmap, fill = NA) +
  scale_colour_viridis_c() +
  labs(Fill = "Mean")

ggsave(here::here(figdir, "descriptive/delay_krige.png"), height = 7, width = 9, units = "in")

# ---------------------------------------------------------------------------- #
# Correlation between all covariates

vars <- dat %>%
  dplyr::select(-longitude, -latitude, -vil_code, -patient_id, -vill_inc_2017_t) %>%
  dplyr::mutate(across(everything(), as.numeric))

M <- cor(vars)
corrplot::corrplot(M, method = "color", type = "lower", order = "hclust")

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# matrix of the p-value of the correlation
p.mat <- cor.mtest(vars)

png(here::here(figdir,"descriptive/covariate effects", "corrplot.png"), 
    height = 8, width = 8, units = "in", res = 300)
corrplot::corrplot(M, method = "color", type = "lower", #order = "hclust", 
                   p.mat = p.mat, sig.level = 0.01)
dev.off()

################################################################################