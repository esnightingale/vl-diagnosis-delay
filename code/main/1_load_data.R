################################################################################
# Description: Read geotagged case data, map context and access raster
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

# Read cleaned analysis data 
dat <- read_data()

# ---------------------------------------------------------------------------- #
# Exclude observations missing any covariate of interest

n_all <- nrow(dat)
dat_nona <- dat %>%
  dplyr::select(delay, gt30, dur_fev_r, age_s, sex, comorb1, caste4_r, occ4_cat, poss_acd, rain, prv_tx, 
                latitude, longitude, traveltime_s, traveltime_t_s, 
                inc_2017_gt0, IRS_2017, block_endm_2017, id, v, district, block, geometry) %>%
  drop_na() %>%
  st_as_sf()

n_nonmiss <- nrow(dat_nona)

# Compare delay between those with complete/incomplete covariate information
summary(dat$delay)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -10.00   11.00   16.00   31.04   44.00  496.00 
summary(dat_nona$delay)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -10.00   11.00   16.00   30.86   44.00  496.00

print(paste(n_all - n_nonmiss,"observations deleted due to missingness"))
# "84 observations deleted due to missingness"

dat <- dat_nona %>%
  dplyr::filter(delay >= 0)

print(paste(n_nonmiss - nrow(dat),"observations deleted due to negative delay"))

saveRDS(dat, here::here("data/analysis","dat_nona.rds"))

# ---------------------------------------------------------------------------- #
# Split fitting and validation data

# v.idx <- sample(1:nrow(dat), floor(nrow(dat)*val.size))
# 
# dat.fit <- dat[-v.idx,]
# dat.val <- dat[v.idx,]
# 
# dat.fit.df <- st_drop_geometry(dat.fit)
# dat.val.df <- st_drop_geometry(dat.val)
# 
# # Full dataset with outcome set to NA for validation points
# dat.fit.val <- dat
# dat.fit.val$days_fever[v.idx] <- NA
# 
# saveRDS(dat.fit, here::here("data/analysis","dat_fit.rds"))
# saveRDS(dat.val, here::here("data/analysis","dat_val.rds"))
# saveRDS(dat.fit.val, here::here("data/analysis","dat_fit_val.rds"))

################################################################################
################################################################################