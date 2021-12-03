################################################################################
# Description: Setup data set of CDF data from individuals included in Dubey et 
# al. matched to village GPS locations, for primary analysis
################################################################################
################################################################################

# Local data folder
datadir <- "C:/Users/phpuenig/Documents/VL/Data/KAMIS/Clean/linelist"

dat <- readRDS(file.path(datadir,"cdf_gps.rds")) %>%
  dplyr::filter(with_gps) %>% 
  dplyr::rename(marg_caste = caste4_r) %>%
  dplyr::mutate(age_child = factor(as.numeric(age_cdf < 16), levels = c(1, 0), labels = c("Under 16", "16+")),
                age_cat = factor(cut(age_cdf, breaks = c(0,15,35,100), include.lowest = TRUE)),
                occupation = recode(occ4_cat, 
                                    "0" = "Not working",
                                    "1" = "Unskilled",
                                    "2" = "Skilled",
                                    "3" = "Salaried.selfemployed"),
                prv_tx = factor((prv_tx_ka == "Yes" | prvtx_pkdl == "Yes"), levels = c(FALSE,TRUE), labels = c("No","Yes")),
                conslt_cat = factor(cut(num_conslt, breaks = c(0,2,5,8), include.lowest = TRUE)),
                num_conslt = as.factor(num_conslt),
                detection = recode(poss_acd, "0" = "PCD", "1" = "ACD"),
                diag_rainseason = factor(lubridate::month(diag_date) %in% 6:9, levels = c(FALSE,TRUE), labels = c("No","Yes")),
                block_endm_2017 = factor(block_endm_2017, levels = c(0,1), labels = c("Non-endemic","Endemic")),
                delay = days_fever - 14) %>% 
  dplyr::select(delay, days_fever, gt30, gt90, 
                age_cdf, age_child, age_cat, sex_cdf, hiv, marg_caste, occupation, 
                num_conslt, conslt_cat, prv_tx, detection, 
                diag_month, diag_year, diag_rainseason, blk_code, vil_code,
                data_entry_district, data_entry_block, patient_sc, patient_village, 
                population, latitude, longitude, vl_affect_1517, inc_2017, IRS_2017, block_endm_2017) 

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
boundary <- sf::st_union(blockmap)  
saveRDS(boundary, here::here("data","geography","boundary.rds"))

extent <- setNames(sf::st_bbox(blockmap), c("left","bottom","right","top"))
bh_lines <- ggmap::get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)

boundary %>%
  sf::st_transform(crs = sf::st_crs(7759)) -> boundary

# Transform to India projection so that buffer is calculated accurately
dat.sf <- sf::st_as_sf(dat, coords = c("longitude","latitude"), remove = FALSE) %>%
  sf::st_set_crs(4326) %>%
  sf::st_transform(7759)

buffer <- sf::st_buffer(boundary, 1e4)

ggplot() +
  geom_sf(data = boundary) +
  geom_sf(data = buffer, fill = NA) +
  geom_sf(data = dat.sf)

# Exclude likely erroneous points beyond a 10km buffer of the state border
dat_incl <- sf::st_intersection(dat.sf, buffer)

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


# Return final included locations to lat/long projection on KM scale
dat_clean <- dat_incl %>%
  sf::st_transform(4326) 

# ---------------------------------------------------------------------------- #
# Add annual village incidence

# Note: one patient in study data (IN/BI/ARR/BHA/2018/05/003) from BAIJUPATTI village in 
# DHANESHARI HSC which has missing population size - excluded.

# vil_inc <- readRDS(file.path(datadir,"kamis.rds")) %>%
#   dplyr::filter(!is.na(vil_code) & !is.na(population) & diag_year != 2021) %>%
#   dplyr::group_by(diag_year, vil_code, latitude, longitude, population) %>%
#   dplyr::count() %>%
#   dplyr::mutate(inc = n/population) %>%
#   dplyr::select(-n) %>%
#   dplyr::ungroup() %>%
#   tidyr::pivot_wider(names_from = "diag_year", values_from = inc, names_prefix = "vill_inc_") %>%
#   dplyr::mutate(across(starts_with("vill_inc"), function(x) replace_na(x, replace = 0)))

ggmap(bh_lines,
      base_layer = ggplot(data = filter(dat_clean, vl_affect_1517 > 0), aes(x = longitude, y = latitude))) +
  geom_point(cex = 0.3) +
  labs(x = "", y = "", title = "Bihar geo-tagged villages with non-zero VL incidence 2015-17")

ggsave(here::here("figures", "affected_villages.png"), height = 6, width = 8, units = "in")

# dat_clean %>%
#   dplyr::inner_join(vil_inc) -> dat_winc

# ---------------------------------------------------------------------------- #
# Add village IRS targeting and block endemicity

# village_chars <- readRDS(file.path(datadir,"village_irs_endemicity.rds")) %>% 
#   dplyr::mutate(vil_code = as.factor(vil_code),
#                 IRS_2015 = as.factor(IRS_2015_R1 + IRS_2015_R2),
#                 IRS_2016 = as.factor(IRS_2016_R1 + IRS_2016_R2),
#                 IRS_2017 = as.factor(IRS_2017_R1 + IRS_2017_R2)) %>%
#   dplyr::mutate(across(starts_with("block_endm"), 
#                               function(x) factor(x, 
#                                           levels = c(0,1), 
#                                           labels = c("Non-endemic","Endemic")))) %>%
#   dplyr::select(vil_code, vl_affected_village, IRS_2015, IRS_2016, IRS_2017, block_endm_2015, block_endm_2016, block_endm_2017)
# 
# dat_winc %>%
#   dplyr::inner_join(village_chars) -> dat_wchars
  
# dat_wchars %>%
#   st_drop_geometry() %>%
#   dplyr::select(patient_id, diag_date, days_fever, district, block, village, vil_code, vill_inc_2013:vill_inc_2020) -> for_graham
# 
# saveRDS(for_graham, here::here("data","delay_w_villinc.rds"))

# ---------------------------------------------------------------------------- #
# Add travel time values for each village

access.diag <- raster::raster(here::here("data","covariates","diag_facility_travel_time.tif"))
access.trt <- raster::raster(here::here("data","covariates","trt_facility_travel_time.tif"))

# malariaAtlas::autoplot_MAPraster(access.raster)

# Extract raster values at village points
dat_clean$traveltime <- raster::extract(access.diag, dat_clean)
dat_clean$traveltime_t <- raster::extract(access.trt, dat_clean)

# Define categorical 
dat_clean <- dat_clean %>%
  dplyr::mutate(travel_time_cat = cut(traveltime, 
                                      breaks = c(0, 10, 20, 30, 60, 150), 
                                      include.lowest = TRUE),
                travel_time_t_cat = cut(traveltime_t, 
                                      breaks = c(0, 10, 20, 30, 60, 150), 
                                      include.lowest = TRUE)
                # travel_time_cat4 = cut(traveltime, 
                #                        breaks = round(quantile(dat_clean$traveltime, probs = c(0, 0.25, 0.5, 0.75, 1))), 
                #                        include.lowest = TRUE,
                #                        ordered_result = TRUE)
                )

# waccess <- st_as_sf(dat_wchars, coords = c("longitude","latitude")) %>%
#   st_set_crs(4326)

# Check all points have an extracted raster value
ggplot() +
  geom_sf(data = boundary) +
  geom_sf(data = dat_clean, aes(col = is.na(traveltime)))

# Distribution of travel time
summary(dat_clean$traveltime)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    6.58   11.69   14.46   18.42  123.85 
summary(dat_clean$traveltime_t)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   10.63   17.94   20.91   27.16  130.55

ggplot(dat_clean, aes(traveltime)) + 
  geom_histogram(bins = 50) +
  labs(title = "Travel time to most accessible diagnosis facility",
       x = "Time (minutes)",
       y = "Count")
ggsave(here::here("figures", "covariates", "traveltime_hist.png"), height = 6, width = 8, units = "in")

ggplot(dat_clean, aes(traveltime_t)) + 
  geom_histogram(bins = 50) +
  labs(title = "Travel time to most accessible treatment facility",
       x = "Time (minutes)",
       y = "Count")
ggsave(here::here("figures", "covariates", "traveltime_trt_hist.png"), height = 6, width = 8, units = "in")

# ---------------------------------------------------------------------------- #
# Aggregate to village level

dat_clean <- dat_clean %>%
  dplyr::rename(district = data_entry_district,
         block = data_entry_block,
         hsc = patient_sc,
         village = patient_village)

dat_clean %>%
  dplyr::group_by(district, block, blk_code, hsc, village, vil_code, population,
           vl_affect_1517, IRS_2017, inc_2017, block_endm_2017, traveltime, traveltime_t,
           diag_year) %>%
  dplyr::summarise(n_cases = n(),
            n_acd = sum(detection == "ACD"),
            p_acd = n_acd/n_cases,
            delay_med = median(days_fever, na.rm = T),
            delay_iqr = quantile(days_fever, p = 0.75, na.rm = T) - quantile(days_fever, p = 0.25, na.rm = T),
            gt90 = sum(gt90),
            gt30 = sum(gt30)) %>%
  dplyr::ungroup() -> dat_village_yr

dat_village_yr %>%
  sf::st_drop_geometry() %>%
  dplyr::select(district, block, blk_code, village, vil_code, population, vl_affect_1517:traveltime_t, diag_year, delay_med) %>%
  tidyr::pivot_wider(names_from = diag_year, values_from = delay_med, names_prefix = "med_delay_") -> vil_inc_delay

# vil_inc_delay %>% 
#   ggplot(aes(med_delay_2018, inc_2017)) +
#   geom_point(alpha = 0.5)

dat_clean %>%
  dplyr::group_by(district, block, village, vil_code, population,
           vl_affect_1517, inc_2017, IRS_2017, block_endm_2017, traveltime, traveltime_t) %>%
  dplyr::summarise(n_cases = n(),
            n_acd = sum(detection == "ACD"),
            p_acd = n_acd/n_cases,
            gt90 = sum(gt90),
            gt30 = sum(gt30)) %>%
  dplyr::ungroup() -> dat_village_tot

# ---------------------------------------------------------------------------- #
# Save analysis datasets

saveRDS(dat_clean, file.path(datadir,"analysisdata_individual.rds"))

saveRDS(dat_village_yr, file.path(datadir,"analysisdata_village_year.rds"))

saveRDS(dat_village_tot, file.path(datadir,"analysisdata_village.rds"))

################################################################################
################################################################################
