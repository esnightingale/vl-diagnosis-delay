################################################################################
# Description: Read geotagged case data, map context and access raster
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

# # Setup map context
# blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
#   st_transform(4326)
# boundary <- sf::st_union(blockmap)
# boundary.spdf <- as_Spatial(boundary)
# 
# extent <- setNames(st_bbox(blockmap), c("left","bottom","right","top"))
# bh_lines <- get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)

dat <- readRDS(here::here("data","analysisdata_individual.rds")) %>%
  mutate(le30 = as.numeric(days_fever <= 30),
         id = row_number(),
         v = as.numeric(as.factor(vil_code)),
         rain = as.numeric(diag_rainseason),
         age_s = as.numeric(scale(age, center = T)),
         consult_gt1 = (num_conslt > 1),
         traveltime_s = as.numeric(scale(traveltime, center = T)),
         IRS_2017_1 = factor((IRS_2017 != 0), labels = c("No","Yes")),
         vill_inc_2017_t = vill_inc_2017*1e3 + 1e-4,
         vill_inc_2017_s = as.numeric(scale(vill_inc_2017, center = T)),
         vill_inc_2017_gt0 = factor((replace_na(vill_inc_2017,0) > 0), labels = c("No","Yes")),
  )

# Travel time raster
# access <- raster::raster(here::here("data","covariates","diag_facility_travel_time.tif"))
#
# # All Bihar villages
# villages <- readRDS(here::here("data","geography","village_centroids_2011.rds")) 

# ---------------------------------------------------------------------------- #
# Split fitting and validation data

n_all <- nrow(dat)
dat <- dat %>%
  dplyr::select(days_fever, age_s, sex, hiv, marg_caste, occ4_cat, detection, prv_tx, 
                consult_gt1, latitude, longitude, traveltime, traveltime_s, vill_inc_2017_t,
                vill_inc_2017_gt0, IRS_2017_1, block_endm_2017, id, v, rain) %>%
  drop_na() 

paste(n_all - nrow(dat),"observations deleted due to missingness")
saveRDS(dat, here::here("data/analysis","dat_nona.rds"))

# Split out final validation data - not to be used for fitting or tuning
v.idx <- sample(1:nrow(dat), floor(nrow(dat)*0.25))

dat.fit <- dat[-v.idx,]
dat.val <- dat[v.idx,]

dat.fit.df <- st_drop_geometry(dat.fit)
dat.val.df <- st_drop_geometry(dat.val)

# Full dataset with outcome set to NA for validation points
dat.fit.val <- dat
dat.fit.val$days_fever[v.idx] <- NA

saveRDS(dat.fit, here::here("data/analysis","dat_fit.rds"))
saveRDS(dat.val, here::here("data/analysis","dat_val.rds"))
saveRDS(dat.fit.val, here::here("data/analysis","dat_fit_val.rds"))


################################################################################
################################################################################