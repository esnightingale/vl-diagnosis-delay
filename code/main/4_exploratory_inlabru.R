################################################################################
# Description: Exploratory fits using inlabru with no covariates, to compare IID
# vs IID+SPDE models for spatial variation
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit/exploratory"
outdir <- "output/exploratory"

dat <- readRDS(here::here("data/analysis","dat_nona.rds"))  %>%
  dplyr::filter(delay >= 0) %>%
  dplyr::mutate(days_fever_cat = gtools::quantcut(days_fever, 5)) %>%
  sf::st_transform(7759)

# %>%
#   filter(detection == "PCD")
mesh <- readRDS(here::here("data/analysis","mesh.rds")) 
spde <- readRDS(here::here("data/analysis","spde.rds"))

# Setup map context
blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
  sf::st_transform(7759)

boundary <- blockmap %>%
  sf::st_union()

# ---------------------------------------------------------------------------- #

cmp1 <- delay ~ Intercept(1, mean.linear = 0, prec.linear = 0.1) +
  rand(v, model = "iid", prior = "pc.prec", param = c(1, 0.01)) +
  field(main = st_coordinates, model = spde)

fit1 <- inlabru::bru(components = cmp1, 
                     data = dat,
                     family = "nbinomial",
                     options = bru_options(bru_verbose = TRUE,
                                           control.compute = list(dic = TRUE, 
                                                                  waic = TRUE, 
                                                                  config = TRUE,
                                                                  cpo = TRUE)))

cmp2 <- delay ~ Intercept(1, mean.linear = 0, prec.linear = 0.1) + 
  field(main = st_coordinates, model = spde)

fit2 <- inlabru::bru(components = cmp2, 
                    data = dat,
                    family = "nbinomial",
                    options = bru_options(bru_verbose = TRUE,
                                          control.compute = list(dic = TRUE, 
                                                                 waic = TRUE, 
                                                                 config = TRUE,
                                                                 cpo = TRUE)))

# Add spatial covariate
library(raster)
access.diag <- raster(here::here("data","covariates","diag_facility_travel_time.tif")) 
# access.trt <- raster(here::here("data","covariates","trt_facility_travel_time.tif"))

access.diag <- raster::projectRaster(access.diag, crs = 7759) 
access.crop <- raster::mask(access.diag, boundary.spdf)

access.spdf <- as(access.crop, "SpatialPixelsDataFrame") 

ggplot() +
  gg(access.spdf) +
  gg(boundary.spdf, alpha = 0) +
  scale_fill_viridis_c(trans = "sqrt") +
  coord_fixed()

# Centre values for more stable estimation
access.spdf$diag_facility_travel_time <- access.spdf$diag_facility_travel_time - mean(access.spdf$diag_facility_travel_time, na.rm = TRUE)

# Define a function to extract travel time at any particular point
f.travel <- function(x, y) {
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(access.spdf))
  proj4string(spp) <- fm_sp_get_crs(access.spdf)
  # Extract travel time values at spp coords, from our access.spdf SpatialGridDataFrame
  v <- over(spp, access.spdf)
  if (any(is.na(v$diag_facility_travel_time))) {
    v$diag_facility_travel_time <- inlabru:::bru_fill_missing(access.spdf, spp, v$diag_facility_travel_time)
  }
  return(v$diag_facility_travel_time)
}

cmp3a <- delay ~ Intercept(1, mean.linear = 0, prec.linear = 0.1) + 
  field(main = st_coordinates, model = spde) + 
  travel(f.travel(x, y), model = "linear", mean.linear = 0, prec.linear = 0.1)

fit3a <- inlabru::bru(components = cmp3a, 
                     data = dat,
                     family = "nbinomial",
                     options = bru_options(bru_verbose = TRUE,
                                           control.compute = list(dic = TRUE, 
                                                                  waic = TRUE, 
                                                                  config = TRUE,
                                                                  cpo = TRUE)))


cmp3b <- delay ~ Intercept(1, mean.linear = 0, prec.linear = 0.1) + 
  field(main = st_coordinates, model = spde) + 
  travel(traveltime_s, model = "linear", mean.linear = 0, prec.linear = 0.1)

fit3b <- inlabru::bru(components = cmp3b, 
                     data = dat,
                     family = "nbinomial",
                     options = bru_options(bru_verbose = TRUE,
                                           control.compute = list(dic = TRUE, 
                                                                  waic = TRUE, 
                                                                  config = TRUE,
                                                                  cpo = TRUE)))


#  Try traditional specification of NB process - no. trials before r = 1 successes
#  => Assumes equal and independent probability of diagnosis each day

fit4 <- inlabru::bru(components = cmp2, 
                      data = dat,
                      family = "nbinomial2",
                      Ntrials = 1,
                      options = bru_options(bru_verbose = TRUE,
                                            control.compute = list(dic = TRUE, 
                                                                   waic = TRUE, 
                                                                   config = TRUE,
                                                                   cpo = TRUE)))


fits <- list(IID_SPDE = fit1, SPDE = fit2,
             SPDE_cov = fit3a, SPDE_cov2 = fit3b,
             SPDE_nbmod = fit4)

lapply(fits, function(x) -mean(log(x$cpo$cpo)))

hist(fit4$cpo$pit, breaks = 10, prob = TRUE)

summary(fit4$cpo$failure)

pwdiff_cpo <- log(fit1$cpo$cpo) - log(fit2$cpo$cpo)
hist(pwdiff_cpo, breaks = 50)
# more negative than positive => scores from fit 2 generally higher

pwdiff_cpo <- log(fit2$cpo$cpo) - log(fit3$cpo$cpo)
hist(pwdiff_cpo, breaks = 50)
# more positive than negative => scores from fit 2 generally higher

pwdiff_cpo <- log(fit2$cpo$cpo) - log(fit3b$cpo$cpo)
hist(pwdiff_cpo, breaks = 50)
# more positive than negative => scores from fit 2 generally higher

pwdiff_cpo <- log(fit2$cpo$cpo) - log(fit4$cpo$cpo)
hist(pwdiff_cpo, breaks = 50)
# more negative than positive => scores from fit 4 generally higher

# Overall seems that fit comparison is determined by scores for a small minority of observations

fit.int <- predict(fit2, pixels(mesh, mask = boundary.spdf), ~ exp(field + Intercept))

ggplot() +
  gg(e.int) +
  gg(boundary, alpha = 0) +
  gg(nests, shape = "+") +
  coord_equal()