################################################################################
# Description: Setup data set of CDF data from individuals included in Dubey et 
# al. matched to village GPS locations, for primary analysis
################################################################################
################################################################################

dat <- readRDS(here::here("data","kamis_cdf.rds")) %>%
  filter(with_gps) %>%
  dplyr::mutate(marg_caste_kamis = (caste_category %in% c("SC","ST"))) %>%
  dplyr::select(c(country:patient_id, ends_with("_kamis"), `treated_for_kala-azar_earlier`,
                caste_category, special_caste, basis_of_diagnosis,
                hiv_tests_done:tb,whether_other_place, ends_with("_cdf"),
                caste4, caste4_r, mahadalit4, prv_tx_ka, prvtx_pkdl, occ4_cat, adv_asha, 
                poss_acd, num_conslt)) 

# ---------------------------------------------------------------------------- #
# Add village IRS targeting and block endemicity

village_chars <- readRDS(here::here("data","village_irs_endemicity.rds")) %>%
  mutate(vil_code = as.factor(vil_code),
         IRS_2017 = factor(IRS_2017_R1 + IRS_2017_R2)) %>%
  dplyr::select(vil_code, IRS_2017, block_endm_2017)

dat %>%
  left_join(village_chars) -> dat
  
# ---------------------------------------------------------------------------- #
# Add travel time values for each village

access.raster <- raster::raster(here::here("data","covariates","diag_facility_travel_time.tif"))

# Make spdf
dat.spdf <- dat
coordinates(dat.spdf) <- ~ longitude + latitude

# Extract raster values at village points
dat$traveltime <- raster::extract(access.raster, dat.spdf)

# ---------------------------------------------------------------------------- #
# Save analysis dataset

saveRDS(dat, here::here("data","analysisdata.rds"))

################################################################################
################################################################################
