################################################################################
# Description: Setup data set of CDF data from individuals included in Dubey et 
# al. matched to village GPS locations, for primary analysis
################################################################################
################################################################################

dat <- readRDS(here::here("data","kamis_cdf.rds")) %>%
  dplyr::filter(with_gps) %>% 
  dplyr::select(country:patient_id, diag_month_cdf, diag_year_cdf, days_fever_cdf,
                gt30_cdf, gt90_cdf, poss_acd) %>%
  dplyr::mutate(poss_acd = (poss_acd == 1))
  
  # dplyr::mutate(marg_caste_kamis = (caste_category %in% c("SC","ST"))) %>%
  # dplyr::select(c(country:patient_id, ends_with("_kamis"), `treated_for_kala-azar_earlier`,
  #               caste_category, special_caste, basis_of_diagnosis,
  #               hiv_tests_done:tb,whether_other_place, ends_with("_cdf"),
  #               caste4, caste4_r, mahadalit4, prv_tx_ka, prvtx_pkdl, occ4_cat, adv_asha, 
  #               poss_acd, num_conslt)) 

# ---------------------------------------------------------------------------- #
# Exclude points beyond a small buffer of Bihar boundary

# Setup map context
blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) 
boundary <- sf::st_union(blockmap)  %>%
  st_transform(crs = st_crs(24380))
extent <- setNames(st_bbox(blockmap), c("left","bottom","right","top"))
bh_lines <- ggmap::get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)

# Transform to India projection so that buffer is calculated accurately
dat.sf <- st_as_sf(dat, coords = c("longitude","latitude")) %>%
  st_set_crs(4326) %>%
  st_transform(24380)

buffer <- st_buffer(boundary, 1e4)

ggplot() +
  geom_sf(data = boundary) +
  geom_sf(data = buffer, fill = NA) +
  geom_sf(data = dat.sf)

# Exclude likely erroneous points beyond a 10km buffer of the state border
dat_incl <- st_intersection(dat.sf, buffer)

# # Identify points within a small margin of the border
# dat.out <- st_difference(dat.sf, boundary)
# 
# ggplot() +
#   geom_sf(data = boundary) +
#   geom_sf(data = dat.out)
# 
# # Snap points into state boundary
# dat.snap <- st_snap(dat.sf, boundary, tolerance = 600)
# 
# # Check no remaining points outside boundary
# # dat.out2 <- st_difference(dat.snap, boundary) # 0
# 
# ggplot() +
#   geom_sf(data = boundary) +
#   geom_sf(data = dat.out, 
#           col = "red", alpha = 0.5) +
#   geom_sf(data = dat.snap[dat.snap$vil_code %in% dat.out$vil_code,], 
#           col = "green", alpha = 0.5)


# Return final included locations to lat/long projection and map out
dat_clean <- dat_incl %>%
  st_transform(4326)

ggmap(bh_lines,
      base_layer = ggplot(data = dat_clean)) +
  geom_sf(cex = 0.3) +
  labs(x = "", y = "", title = "Bihar village locations - VL affected 2013-18")

ggsave(here::here("figures", "affected_villages.png"), height = 6, width = 8, units = "in")

# ---------------------------------------------------------------------------- #
# Add village IRS targeting and block endemicity

village_chars <- readRDS(here::here("data","village_irs_endemicity.rds")) %>% 
  dplyr::mutate(vil_code = as.factor(vil_code),
         IRS_2016 = factor(IRS_2016_R1 + IRS_2016_R2),
         IRS_2017 = factor(IRS_2017_R1 + IRS_2017_R2)) %>%
  dplyr::select(vil_code, vl_affected_village, IRS_2016, IRS_2017, block_endm_2016, block_endm_2017)

dat_clean %>%
  dplyr::left_join(village_chars) -> dat_wchars
  
# ---------------------------------------------------------------------------- #
# Add travel time values for each village

access.raster <- raster::raster(here::here("data","covariates","diag_facility_travel_time.tif"))
# malariaAtlas::autoplot_MAPraster(access.raster)

# Extract raster values at village points
dat_wchars$traveltime <- raster::extract(access.raster, dat_wchars)

# waccess <- st_as_sf(dat_wchars, coords = c("longitude","latitude")) %>%
#   st_set_crs(4326)

# Check all points have an extracted raster value
ggplot() +
  geom_sf(data = boundary) +
  geom_sf(data = dat_wchars, aes(col = is.na(dat_wchars$traveltime)))

# Distribution of travel time
summary(dat_wchars$traveltime)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   6.534  11.690  14.448  18.419 123.846 

ggplot(dat_wchars, aes(traveltime)) + 
  geom_histogram(bins = 50) +
  labs(title = "Travel time to most accessible health facility",
       x = "Time (minutes)",
       y = "Count")

ggsave(here::here("figures", "covariates", "traveltime_hist.png"), height = 6, width = 8, units = "in")

# ---------------------------------------------------------------------------- #
# Aggregate to village level

# Transform back to data frame with coordinate variables
# coords <- setNames(as.data.frame(st_coordinates(dat_clean)), c("longitude","latitude"))
# dat_clean <- st_drop_geometry(dat_clean) %>%
#   bind_cols(coords)

dat_wchars %>%
  group_by(district, block, village, vil_code, population,
           vl_affected_village, IRS_2017, block_endm_2017, traveltime, 
           diag_year_cdf, diag_month_cdf) %>%
  summarise(n_cases = n(),
            n_acd = sum(poss_acd),
            p_acd = n_acd/n_cases,
            gt90 = sum(gt90_cdf),
            gt30 = sum(gt30_cdf)) -> dat_village

# ---------------------------------------------------------------------------- #
# Save analysis datasets

saveRDS(dat_wchars, here::here("data","analysisdata_individual.rds"))

saveRDS(dat_village, here::here("data","analysisdata_village.rds"))

################################################################################
################################################################################
