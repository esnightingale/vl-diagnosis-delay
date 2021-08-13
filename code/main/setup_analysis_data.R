################################################################################
# Description: Setup data set of CDF data from individuals included in Dubey et 
# al. matched to village GPS locations, for primary analysis
################################################################################
################################################################################

dat <- readRDS(here::here("data","kamis_cdf.rds")) %>%
  dplyr::filter(with_gps) %>% 
  dplyr::mutate(age_child = factor(as.numeric(age_cdf < 16), levels = c(1, 0), labels = c("Under 16", "16+")),
                detection = recode(poss_acd, "0" = "PCD", "1" = "ACD"))  %>%
  dplyr::rename(marg_caste = caste4_r) %>%
  dplyr::select(country:patient_id, age_cdf, age_child, sex_cdf, hiv_cdf, marg_caste,
                prv_tx_ka, num_conslt, poss_acd, detection,
                diag_date_cdf, diag_month_cdf, diag_year_cdf, days_fever_cdf,
                gt30_cdf, gt90_cdf) 

names(dat) <- gsub("_cdf", "",names(dat)) 
  
summary(dat)
  
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
dat.sf <- st_as_sf(dat, coords = c("longitude","latitude"), remove = FALSE) %>%
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

# ---------------------------------------------------------------------------- #
# Add annual village incidence

# Note: one patient in study data (IN/BI/ARR/BHA/2018/05/003) from BAIJUPATTI village in 
# DHANESHARI HSC which has missing population size - excluded.

vil_inc <- readRDS(here::here("data","kamis.rds")) %>%
  dplyr::filter(!is.na(vil_code) & !is.na(population) & diag_year != 2021) %>%
  dplyr::group_by(diag_year, vil_code, latitude, longitude, population) %>%
  dplyr::count() %>%
  dplyr::mutate(inc = n/population) %>%
  dplyr::select(-n) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from = "diag_year", values_from = inc, names_prefix = "vill_inc_") %>%
  dplyr::mutate(across(starts_with("vill_inc"), function(x) replace_na(x, replace = 0)))

ggmap(bh_lines,
      base_layer = ggplot(data = vil_inc, aes(x = longitude, y = latitude))) +
  geom_point(cex = 0.3) +
  labs(x = "", y = "", title = "Bihar geo-tagged villages with VL incidence 2013-2020")

ggsave(here::here("figures", "affected_villages.png"), height = 6, width = 8, units = "in")

dat_clean %>%
  dplyr::inner_join(vil_inc) -> dat_winc

# ---------------------------------------------------------------------------- #
# Add village IRS targeting and block endemicity

village_chars <- readRDS(here::here("data","village_irs_endemicity.rds")) %>% 
  dplyr::mutate(vil_code = as.factor(vil_code),
                IRS_2015 = as.factor(IRS_2015_R1 + IRS_2015_R2),
                IRS_2016 = as.factor(IRS_2016_R1 + IRS_2016_R2),
                IRS_2017 = as.factor(IRS_2017_R1 + IRS_2017_R2)) %>%
  dplyr::mutate(across(starts_with("block_endm"), 
                              function(x) factor(x, 
                                          levels = c(0,1), 
                                          labels = c("Non-endemic","Endemic")))) %>%
  dplyr::select(vil_code, vl_affected_village, IRS_2015, IRS_2016, IRS_2017, block_endm_2015, block_endm_2016, block_endm_2017)

dat_winc %>%
  dplyr::inner_join(village_chars) -> dat_wchars
  
# dat_wchars %>%
#   st_drop_geometry() %>%
#   dplyr::select(patient_id, diag_date, days_fever, district, block, village, vil_code, vill_inc_2013:vill_inc_2020) -> for_graham
# 
# saveRDS(for_graham, here::here("data","delay_w_villinc.rds"))

# ---------------------------------------------------------------------------- #
# Add travel time values for each village

access.raster <- raster::raster(here::here("data","covariates","diag_facility_travel_time.tif"))
# malariaAtlas::autoplot_MAPraster(access.raster)

# Extract raster values at village points
dat_wchars$traveltime <- raster::extract(access.raster, dat_wchars)

# Define categorical 
dat_wchars <- dat_wchars %>%
  dplyr::mutate(travel_time_cat = cut(traveltime, 
                                      breaks = c(0, 10, 20, 30, 60, 124), 
                                      include.lowest = TRUE),
                travel_time_cat4 = cut(traveltime, 
                                       breaks = round(quantile(dat_wchars$traveltime, probs = c(0, 0.25, 0.5, 0.75, 1))), 
                                       include.lowest = TRUE,
                                       ordered_result = TRUE))

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
           vl_affected_village, IRS_2017, block_endm_2017, traveltime, vill_inc_2017, vill_inc_2018, vill_inc_2019,
           diag_year) %>%
  summarise(n_cases = n(),
            n_acd = sum(detection == "ACD"),
            p_acd = n_acd/n_cases,
            delay_med = median(days_fever, na.rm = T),
            delay_iqr = quantile(days_fever, p = 0.75, na.rm = T) - quantile(days_fever, p = 0.25, na.rm = T),
            gt90 = sum(gt90),
            gt30 = sum(gt30)) %>%
  ungroup() -> dat_village_yr

dat_village_yr %>%
  st_drop_geometry() %>%
  dplyr::select(district, block, village, vil_code, population, vill_inc_2017:vill_inc_2019, diag_year, delay_med) %>%
  pivot_wider(names_from = diag_year, values_from = delay_med, names_prefix = "med_delay_") -> vil_inc_delay

vil_inc_delay %>% 
  ggplot(aes(med_delay_2018, vill_inc_2019)) +
  geom_point(alpha = 0.5)

dat_wchars %>%
  group_by(district, block, village, vil_code, population,
           vl_affected_village, IRS_2017, block_endm_2017, traveltime, vill_inc_2017, vill_inc_2018, vill_inc_2019) %>%
  summarise(n_cases = n(),
            n_acd = sum(detection == "ACD"),
            p_acd = n_acd/n_cases,
            gt90 = sum(gt90),
            gt30 = sum(gt30)) -> dat_village_tot

# ---------------------------------------------------------------------------- #
# Save analysis datasets

saveRDS(dat_wchars, here::here("data","analysisdata_individual.rds"))

saveRDS(dat_village_yr, here::here("data","analysisdata_village_year.rds"))

saveRDS(dat_village_tot, here::here("data","analysisdata_village.rds"))

################################################################################
################################################################################
