################################################################################
# Description: Setup data set of CDF data from individuals included in Dubey et 
# al. matched to village GPS locations, for primary analysis
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

library(tidyverse)
library(Rmisc)

# Local data folder
datadir <- "C:/Users/phpuenig/Documents/VL/Data"

# Raw linelist location
rawdir <- file.path(datadir,"KAMIS/Raw")

# Clean output location 
outdir <- file.path(datadir,"KAMIS/Clean")

# Start date
start <- lubridate::ymd("2013-01-01")

# ---------------------------------------------------------------------------- #
# Location IDs

loc_lookup <- readRDS(here::here(datadir, "KAMIS", "kamis_village_lookup.rds"))

# ---------------------------------------------------------------------------- #
# Patient Linelist

source(here::here("code/utils/clean_pat.R"))
clean_p <- clean_pat(here::here(rawdir,"state","pat.csv"),
                     start = start,
                     state_incl = c("BIHAR"), # c("BIHAR","JHARKHAND"),
                     log = here::here(outdir, "ll", "cleaning_log_pat.txt"))

# ---------------------------------------------------------------------------- #
# Match patients to village GPS

# Village incidence dataset with GPS for affected villages
village_gps <- read.csv(here::here(rawdir, "village", "village level data-Bi+Jh.csv"), header = T) %>%
  dplyr::filter(!is.na(population)) %>% 
  dplyr::mutate(inc_2017_gt0 = (X2017*1000/population > 0),
                IRS_2017 = ((insecticide_2017_R1 != "" | insecticide_2017_R2 != "")),
                block_endm_2017 = (block_endm_2017 > 0)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(vl_affect_1517 = (sum(c_across(X2015:X2017)) > 0)) %>%
  dplyr::ungroup() %>%
  dplyr::select(country:vil_code, vl_affect_1517, inc_2017_gt0, IRS_2017, block_endm_2017) %>%
  dplyr::mutate(across(where(is.character), as.factor)) %>%
  dplyr::filter(!(is.na(longitude) | is.na(latitude)))

# 12874 villages with GPS and known population size

village_gps %>%
  dplyr::right_join(clean_p, by = c("district" = "data_entry_district",
                                    "block" = "data_entry_block",
                                    "hsc" = "patient_sc",
                                    "village" ="patient_village")) -> pat_wgps

summary(!is.na(pat_wgps$longitude))
# Mode   FALSE    TRUE 
# logical    4129   38578 

# pat_wgps %>%
#  dplyr::filter(!is.na(longitude)) -> pat_wgps

# saveRDS(pat_wgps, file.path(outdir,"ll","vl_pat_wgps.rds"))

# ---------------------------------------------------------------------------- #
# Read CDF data from ACD study, from  SAS format
# Remove SAS formats and var labels, select vars of interest and rename/reformat

cdf <- haven::read_sas(here::here(datadir, "ACD evaluation study",
                                  "acddata_102620/acddata_102620.sas7bdat")) %>%
  haven::zap_formats() %>%
  haven::zap_label() %>%
  haven::zap_label() %>%
  dplyr::as_tibble() %>%
  dplyr::filter(included == 1) %>%
  dplyr::select(PID, AGE, SEX, HIV, Date_Diag,
                caste4_r, Prv_TX_KA, PrvTX_PKDL,
                occ4_cat, POSS_ACD, DUR_FEV_R) %>%
  dplyr::rename_with(tolower) %>%
  dplyr::rename(res_patient_code = pid,
                comorb = hiv) %>% #View()
  dplyr::mutate(# Identify missing value codes in factors
    across(c(sex, comorb, caste4_r:poss_acd),
           function(x) factor(x, exclude = c("","99","<NA>","NA",NA))),
    across(where(is.character), as.factor),
    across(c("sex","comorb","prv_tx_ka","prvtx_pkdl","caste4_r"), 
           function(x) relevel(x, ref = 2)),
    diag_year = as.factor(lubridate::year(date_diag)),
    diag_month = lubridate::floor_date(date_diag, "month"),
    rain = (lubridate::month(date_diag) %in% 6:9),
    age_s = as.numeric(scale(age, center = T)),
    age_cat = cut(age, quantile(age, probs = c(0,0.25,0.5,0.75,1)), include.lowest = TRUE),
    prv_tx = (prv_tx_ka == 1 | prvtx_pkdl == 1), #, labels = c("No","Yes")
    poss_acd = (poss_acd == 1), 
    delay = as.numeric(dur_fev_r) - 14,
    gt30 = (delay > 30),
    gt60 = (delay > 60),
    gt90 = (delay > 90),
    id = row_number()) %>%
  dplyr::select(-c(date_diag, prv_tx_ka, prvtx_pkdl))

# ---------------------------------------------------------------------------- #
# Join CDF data with patient village GPS

print("CDF patient not in KAMIS:")
cdf %>%
  dplyr::anti_join(pat_wgps) %>%
  nrow() # 2

pat_wgps %>%
  dplyr::select(country:patient_id) %>%
  inner_join(cdf) -> match

print("Matched patients, with GPS:")
print(summary(!is.na(match$longitude)))
# Mode   FALSE    TRUE 
# logical     649    4379 

match <- dplyr::filter(match, !is.na(longitude))

print("Matched patients, village population:")
print(summary(match$population[!is.na(match$longitude)]))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#      71    1609    2771    3729    5188   21273       1

# ---------------------------------------------------------------------------- #

match$v <- dplyr::group_by(match, vil_code) %>% group_indices()

dat <- match %>% 
  # dplyr::mutate(match_code = sample(1:n())) %>%
  dplyr::select(delay, gt30, gt60, dur_fev_r, 
                age, age_s, age_cat, sex, comorb, caste4_r, occ4_cat, 
                prv_tx, poss_acd, 
                diag_month, diag_year, rain,
                district, block,
                latitude, longitude, population, vl_affect_1517, 
                inc_2017_gt0, IRS_2017, block_endm_2017, res_patient_code,
                id, v) 

summary(dat)

# ---------------------------------------------------------------------------- #
# Exclude points beyond a small buffer of Bihar boundary

# Setup map context
blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) %>% 
  sf::st_set_crs(4326)
extent <- setNames(sf::st_bbox(blockmap), c("left","bottom","right","top"))
bh_lines <- ggmap::get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)

blockmap %>%
  sf::st_union() %>%
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
dat_clean <- sf::st_intersection(dat.sf, buffer)

nrow(dat.sf) - nrow(dat_clean)
# 1

# Plot with some geographic context
ggmap(bh_lines,
      base_layer = ggplot(data = filter(sf::st_transform(dat_clean, 4326) , vl_affect_1517 > 0), aes(x = longitude, y = latitude))) +
  geom_point(cex = 0.3) +
  # theme(axis.text = element_blank()) +
  labs(x = "", y = "", title = "Bihar geo-tagged villages with non-zero VL incidence 2015-17")

ggsave(here::here("figures", "affected_villages.png"), height = 6, width = 8, units = "in")

# ---------------------------------------------------------------------------- #
# Add travel time values for each village

access.diag <- raster::raster(here::here("data","covariates","diag_facility_travel_time.tif")) %>%
  raster::projectRaster(crs = 7759)
access.trt <- raster::raster(here::here("data","covariates","trt_facility_travel_time.tif")) %>%
  raster::projectRaster(crs = 7759)

# malariaAtlas::autoplot_MAPraster(access.diag)
# malariaAtlas::autoplot_MAPraster(access.trt)

# Extract raster values at village points
dat_clean$traveltime <- raster::extract(access.diag, dat_clean)
dat_clean$traveltime_t <- raster::extract(access.trt, dat_clean)

# Define categorical and scaled alternatives
dat_clean <- dat_clean %>%
  dplyr::mutate(travel_time_cat = cut(traveltime, 
                                      5,
                                      # breaks = c(0, 10, 20, 30, 60, 150), 
                                      include.lowest = TRUE),
                travel_time_t_cat = cut(traveltime_t, 
                                        5,
                                        # breaks = c(0, 10, 20, 30, 60, 150), 
                                        include.lowest = TRUE),
                traveltime_s = scale(traveltime),
                traveltime_t_s = scale(traveltime_t))

# Check all points have an extracted raster value
ggplot() +
  geom_sf(data = boundary) +
  geom_sf(data = dat_clean, aes(col = is.na(traveltime)))

# Distribution of travel time
summary(dat_clean$traveltime)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.2611   7.3979  12.3175  14.7635  18.2698 114.0283  
summary(dat_clean$traveltime_t)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.4199  11.4263  18.3357  21.1682  26.7491 119.8681 

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
# Exclude observations missing any covariate of interest

n_all <- nrow(dat_clean) #4378
n_all

# Count missingness in all vars
dat_clean %>% 
  sf::st_drop_geometry() %>% 
  dplyr::summarise(across(age:population, function(var) sum(is.na(var))*100/n_all)) %>% 
  t() 

# age        0.0000000
# sex        0.0000000
# comorb     1.2791229
# caste4_r   0.2512563
# occ4_cat   0.2284148
# prv_tx     0.1598904
# poss_acd   0.0000000
# rain       0.0000000
# latitude   0.0000000
# longitude  0.0000000
# population 0.0000000

dat_nona <- dat_clean %>%
  dplyr::select(delay, gt30, dur_fev_r, diag_month, age, age_s, age_cat, sex, comorb, caste4_r, occ4_cat, poss_acd, rain, prv_tx, 
                latitude, longitude, population, traveltime, traveltime_t, traveltime_s, traveltime_t_s, 
                inc_2017_gt0, IRS_2017, block_endm_2017, res_patient_code, id, v, district, block, geometry) %>%
  drop_na() %>%
  st_as_sf()

n_nonmiss <- nrow(dat_nona) #4294

print(paste(n_all - n_nonmiss,"observations deleted due to missingness"))
# "84 observations deleted due to missingness"

# Compare delay between those with complete/incomplete covariate information
summary(dat_clean$delay)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -10.00   11.00   16.00   31.04   44.00  496.00 
summary(dat_nona$delay)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -10.00   11.00   16.00   30.87   44.00  496.00 

summary(dat_clean$age)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   12.00   25.00   28.12   42.00   95.00 
summary(dat_nona$age)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   12.00   25.00   28.12   42.00   95.00

summary(dat_clean$sex)
# 2    1 
# 1874 2504 = 74.9%
summary(dat_nona$sex)
# 2    1 
# 1833 2461 = 74.5%

summary(dat_clean$poss_acd)
# FALSE    TRUE
# 2624    1754 = 66.8%
summary(dat_nona$poss_acd)
# FALSE    TRUE 
# 2574    1720 = 66.8%

summary(dat_clean$diag_month)
# Min.      1st Qu.       Median         Mean      3rd Qu.         Max.
# "2018-01-01" "2018-05-01" "2018-08-01" "2018-09-02" "2019-02-01" "2019-07-01"
summary(dat_nona$diag_month)
# Min.      1st Qu.       Median         Mean      3rd Qu.         Max. 
# "2018-01-01" "2018-05-01" "2018-08-01" "2018-09-03" "2019-02-01" "2019-07-01"

summary(dat_nona$dur_fev_r)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.00   25.00   30.00   44.87   58.00  510.00 

# ---------------------------------------------------------------------------- #
# Compare included with excluded patients

cdf %>% 
  dplyr::mutate(incl = (res_patient_code %in% dat_nona$res_patient_code)) %>% 
  dplyr::group_by(incl) %>% 
  dplyr::summarise(N = n(),
                   delay = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(delay),1)))),
                     # paste0(median(delay)," [",
                     #              paste(quantile(delay, 0.25), quantile(delay, 0.75), sep = "-"),
                     #              "]"),
                   age = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(age),1)))),
                     # paste0(median(age)," [",
                     #            paste(quantile(age, 0.25), quantile(age, 0.75), sep = "-"),
                     #            "]"),
                   p_female = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(sex == 2),2)))),
                   p_marg = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(na.omit(caste4_r) == 1),2)))),
                   p_marg_miss = round(100*mean(is.na(caste4_r)),1),
                   p_hiv = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(na.omit(comorb) == 1),2)))),
                   p_comorb_miss = round(100*mean(is.na(comorb)),1),
                   p_prv_tx = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(na.omit(prv_tx) == TRUE),2)))),
                   p_prvtx_miss = round(100*mean(is.na(prv_tx)),1),
                   p_unempl = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(na.omit(occ4_cat) == 0),2)))),
                   p_occ_miss = round(100*mean(is.na(occ4_cat)),1),
                   p_acd = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(poss_acd == TRUE),2))))
                   # p_villinc_gt0 = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(inc_2017_gt0 == TRUE),2)))),
                   # p_irs = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(IRS_2017 == TRUE),2)))),
                   # p_end = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(block_endm_2017 == TRUE),2)))),
                   # mean_travel = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(traveltime),2)))),
                   # mean_travel_t = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(traveltime_t),2))))
                   ) %>%
  # mutate(across(, function(x) round(x, 4))) %>% 
  tibble::column_to_rownames("incl") %>% 
  t() -> tab_excl

View(tab_excl)

write.csv(tab_excl, here::here("output","tables","compare_incl_excl.csv"), row.names = TRUE)

# ---------------------------------------------------------------------------- #
# Summarise patients with negative delay

dat_nona %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(delay = factor(delay < 0, 
                        levels = c(FALSE,TRUE), 
                        labels = c("> 14 days","<= 14 days"))) %>%
  dplyr::group_by(delay) %>%
  dplyr::summarise(N = n(),
                   mean_age = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(age),1)))),
                   p_child = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(age < 16),2)))),
                   p_female = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(sex == 2),2)))),
                   p_marg = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(caste4_r == 1),2)))),
                   p_comorb = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(comorb == 1),2)))),
                   p_prvtx = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(prv_tx == TRUE),2)))),
                   p_unempl = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(occ4_cat == 0),2)))),
                   p_acd = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(poss_acd == TRUE),2)))),
                   p_villinc_gt0 = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(inc_2017_gt0 == TRUE),2)))),
                   p_irs = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(IRS_2017 == TRUE),2)))),
                   p_end = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(block_endm_2017 == TRUE),2)))),
                   mean_travel = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(traveltime),2)))),
                   mean_travel_t = do.call(sprintf, c(fmt = "%s [%s, %s]", as.list(round(CI(traveltime_t),2))))) %>%
  tibble::column_to_rownames("delay") %>%
  t() -> tab_negative_delay

View(tab_negative_delay)
write.csv(tab_negative_delay, here::here("output","tables","pat_chars_negative_delay.csv"))

tab_negative_delay %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  tibble::rownames_to_column("Var") %>%
  dplyr::filter(Var != "N") %>%
  tidyr::pivot_longer(c(`> 14 days`, `<= 14 days`)) %>%
  tidyr::separate(value, into = c("low", "mean", "high"), sep = " - ", convert = TRUE) %>% 
  ggplot(aes(x = name, y = mean, ymin = low, ymax = high, col = name)) +
  geom_errorbar() +
  geom_point() +
  facet_wrap(~Var, scales = "free") +
  guides(col = "none")

ggsave(here::here("figures", "descriptive","covars_by_negdelay.png"), height = 10, width = 15)

# Exclude patients with negative delay
dat_final <- dat_nona %>%
  dplyr::filter(delay >= 0)

# Compare delay between those with complete/incomplete covariate information
summary(dat_nona$delay)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -10.00   11.00   16.00   30.87   44.00  496.00
summary(dat_final$delay)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   11.00   16.00   31.05   44.00  496.00 

print(paste(n_nonmiss - nrow(dat_final),"observations deleted due to negative delay"))

################################################################################
# Save analysis datasets

# Note: convert coordinates to km rather than m for more stable SPDE fitting later
st_geometry(dat_final) <- st_geometry(dat_final)/1000

# Separate for secure storage of village and individual chars 
decode <- dat_final %>% 
  sf::st_drop_geometry() %>% 
  dplyr::select(id, v)

pat <- dat_final %>%
  sf::st_drop_geometry() %>% 
  dplyr::select(delay:prv_tx, id)

vill <- dat_final %>% 
  # sf::st_drop_geometry() %>%
  dplyr::select(v, district, block, latitude:block_endm_2017) %>%
  distinct()

saveRDS(pat, file.path(datadir,"KAMIS/Clean/ll","analysisdata_pat.rds"))
saveRDS(vill, file.path(datadir,"KAMIS/Clean/ll","analysisdata_vill.rds"))
saveRDS(decode, "E:/vl-diagnosis-delay/data/decode.rds")

# ---------------------------------------------------------------------------- #

blockmap <- blockmap %>%
  sf::st_transform(7759)

# crs_km <- sp::CRS("+proj=utm +zone=18 +south +ellps=GRS80 +units=km +no_defs")

st_geometry(blockmap) <- st_geometry(blockmap)/1000

blockmap$kamis_master_block[blockmap$kamis_master_block == "UDA KISHAN GANJ"] <- "UDAKISHUNGANJ"
blockmap$kamis_master_block[blockmap$kamis_master_block == "KUSHESHWAR ASTHAN"] <- "KUSHESHWAR ASTHAN (SATIGHAT)"
blockmap$kamis_master_block[blockmap$kamis_master_dist == "MADHUBANI" & blockmap$kamis_master_block == "MADHUBANI"] <- "RAHIKA"
blockmap$kamis_master_block[blockmap$kamis_master_block == "PATNA RURAL"] <- "PATNA SADAR"

boundary <- blockmap %>%  
  sf::st_union() %>%
  st_set_crs(7759)

saveRDS(boundary, here::here("data","geography","boundary.rds"))
saveRDS(blockmap, here::here("data","geography","blockmap.rds"))

################################################################################
################################################################################
