################################################################################
# Description: Read geotagged case data, map context and access raster
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

# Local data folder
datadir <- "C:/Users/phpuenig/Documents/VL/Data/KAMIS/Clean/linelist"

# Proportion to withhold from fitting for final validation
val.size <- 0.25

# Read data and define final variables for modelling
dat <- readRDS(file.path(datadir,"analysisdata_individual.rds")) %>%
  dplyr::mutate(le30 = as.numeric(days_fever <= 30),
                id = row_number(),
                v = as.numeric(as.factor(vil_code)),
                rain = as.numeric(diag_rainseason),
                age_s = as.numeric(scale(age, center = T)),
                consult_gt0 = (num_conslt != "0"),
                traveltime_s = as.numeric(scale(traveltime, center = T)),
                traveltime_t_s = as.numeric(scale(traveltime_t, center = T)),
                IRS_2017 = factor((IRS_2017 != 0), labels = c("No","Yes")),
                inc_2017_t = inc_2017*1e3 + 1e-4,
                inc_2017_s = as.numeric(scale(inc_2017, center = T)),
                inc_2017_gt0 = factor((replace_na(inc_2017,0) > 0), labels = c("No","Yes")),
                vl_affect_1517 = factor(vl_affect_1517, levels = c(0,1), labels = c("No","Yes")))

# ---------------------------------------------------------------------------- #
# Exclude observations missing any covariate of interest

n_all <- nrow(dat)
dat <- dat %>%
  dplyr::select(delay, days_fever, age_s, sex, hiv, marg_caste, occupation, detection, prv_tx, 
                consult_gt0, latitude, longitude, traveltime, traveltime_s, traveltime_t_s, inc_2017_t,
                inc_2017_gt0, IRS_2017, block_endm_2017, id, v, district, block, rain, geometry) %>%
  drop_na() %>%
  st_as_sf()

print(paste(n_all - nrow(dat),"observations deleted due to missingness"))
# "84 observations deleted due to missingness"
saveRDS(dat, here::here("data/analysis","dat_nona.rds"))

# ---------------------------------------------------------------------------- #
# Split fitting and validation data

v.idx <- sample(1:nrow(dat), floor(nrow(dat)*val.size))

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