################################################################################
# Description: Read geotagged case data, map context and access raster
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

dat <- readRDS(here::here("data","analysisdata_individual.rds")) %>%
  mutate(le30 = as.numeric(days_fever <= 30),
         i = row_number(),
         age_s = as.numeric(scale(age, center = T)),
         num_conslt_cat = as.factor(cut(num_conslt, breaks = c(0,1,3,5,8), include.lowest = TRUE)),
         traveltime_s = as.numeric(scale(traveltime, center = T)),
         IRS_2017_1 = factor((IRS_2017 != 0), labels = c("No","Yes")),
         vill_inc_2017_t = vill_inc_2017*1e3 + 1e-4,
         vill_inc_2017_s = as.numeric(scale(vill_inc_2017, center = T)),
         vill_inc_2017_gt0 = factor((replace_na(vill_inc_2017,0) > 0), labels = c("No","Yes")),
  )

# Setup map context
blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
  st_transform(4326)
boundary <- sf::st_union(blockmap)
boundary.spdf <- as_Spatial(boundary)

extent <- setNames(st_bbox(blockmap), c("left","bottom","right","top"))
bh_lines <- get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)
# 
# # Travel time raster
# access <- raster::raster(here::here("data","covariates","diag_facility_travel_time.tif"))
#
# # All Bihar villages
# villages <- readRDS(here::here("data","geography","village_centroids_2011.rds")) 

# ---------------------------------------------------------------------------- #
# Split fitting and validation data

dat <- dat %>%
  dplyr::select(days_fever, age_s, sex, hiv, marg_caste, detection, prv_tx_ka, 
                num_conslt_cat, latitude, longitude, traveltime_s, vill_inc_2017_gt0, 
                IRS_2017_1, block_endm_2017, vil_code, i) %>%
  drop_na() 

idx <- sample(1:nrow(dat), floor(nrow(dat)*0.75))

dat.fit <- dat[idx,]
dat.val <- dat[-idx,]

################################################################################
################################################################################