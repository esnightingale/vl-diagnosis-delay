################################################################################
# Description: Clean raw KAMIS patient/diagnosis linelists and split by case type
################################################################################
################################################################################

library(tidyverse)

# Local data folder
datadir <- "C:/Users/phpuenig/Documents/VL/Data"

# Raw linelist location
rawdir <- file.path(datadir,"KAMIS/Raw")

# Clean output location 
outdir <- file.path(datadir,"KAMIS/Clean")

# Start date
start <- lubridate::ymd("2013-01-01")

redefine_fct <- function(x) {
  relevel(recode(x, "1" = "Yes", "2" = "No"), ref = "No")
}

# ---------------------------------------------------------------------------- #
# Location IDs

loc_lookup <- readRDS(here::here(datadir, "KAMIS", "kamis_village_lookup.rds"))

# ---------------------------------------------------------------------------- #
# Patient Linelist

source(here::here("code/utils/clean_pat.R"))
clean_p <- clean_pat(here::here(rawdir,"state","pat.csv"),
                      start = start,
                      state_incl = c("BIHAR","JHARKHAND"),
                      log = here::here(outdir, "linelist", "cleaning_log_pat.txt"))

saveRDS(clean_p, file.path(outdir, "linelist", "linelist_patients.rds"))

# ---------------------------------------------------------------------------- #
# Match patients to village GPS

# Village incidence dataset with GPS for affected villages
village_gps <- read.csv(here::here(rawdir, "village", "village level data-Bi+Jh.csv"), header = T) %>%
  dplyr::mutate(inc_2017 = replace_na(X2017*1000/population, 0),
                IRS_2017 = as.numeric((insecticide_2017_R1 != "" | insecticide_2017_R2 != ""))) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(vl_affect_1517 = as.numeric(sum(c_across(X2015:X2017)) > 0)) %>%
  dplyr::ungroup() %>%
  dplyr::select(country:vil_code, vl_affect_1517, inc_2017, IRS_2017, block_endm_2017) %>%
  dplyr::mutate(across(where(is.character), as.factor)) %>%
  dplyr::filter(!(is.na(longitude) | is.na(latitude)))

# 12894 villages with GPS
# vl_affected_village defined as any 

clean_p %>%
  dplyr::left_join(village_gps, by = c("data_entry_district" = "district",
                                       "data_entry_block" = "block",
                                       "patient_sc" = "hsc",
                                       "patient_village" = "village")) %>%
  dplyr::mutate(with_gps = !is.na(longitude)) -> pat_wgps

summary(pat_wgps$with_gps)
#    Mode   FALSE    TRUE 
# logical    4123   38584 

saveRDS(pat_wgps, file.path(outdir,"linelist","vl_patients_wgps.rds"))

# ---------------------------------------------------------------------------- #
# Read CDF data from ACD study, from  SAS format
# Remove SAS formats and var labels, select vars of interest and rename/reformat

cdf <- haven::read_sas(here::here(datadir, "ACD evaluation study",
                                  "acddata_102620/acddata_102620.sas7bdat")) %>%
  haven::zap_formats() %>%
  haven::zap_label() %>%
  haven::zap_label() %>%
  dplyr::as_tibble() %>% #View()
  dplyr::filter(included == 1) %>%
  dplyr::select(PID, DIST_CDF, BLOCK_CDF, VILL_CDF, AGE, SEX, HIV, Date_Diag,
                Trt_St_Dat, CASTE4, caste4_r, MAHADALIT4, Prv_TX_KA, PrvTX_PKDL,
                occ4_cat, ADV_ASHA, POSS_ACD, Refer, Num_conslt, DUR_FEV_R) %>%
  dplyr::rename_with(tolower) %>%
  dplyr::rename(patient_id = pid,
                diag_date = date_diag,
                days_fever = dur_fev_r,
                referred_by = refer,
                treatment_start_date = trt_st_dat) %>%
  dplyr::mutate(across(dist_cdf:vill_cdf, toupper),
                # Identify missing value codes in factors
                across(c(dist_cdf:vill_cdf,sex, hiv, caste4:referred_by),
                       function(x) factor(x, exclude = c("","99","<NA>","NA",NA))),
                across(where(is.character), as.factor),
                # Input data have indicators with 1=yes 2=no - recode to be more intuitive
                across(c("hiv","prv_tx_ka","prvtx_pkdl","caste4_r"), redefine_fct),
                sex = recode(sex, "1" = "Male", "2" = "Female"),
                age = as.numeric(age),
                diag_month = lubridate::floor_date(diag_date, "month"),
                diag_year = as.factor(lubridate::year(diag_date)),
                days_fever = as.numeric(days_fever),
                gt30 = (days_fever > 30),
                gt90 = (days_fever > 90))


saveRDS(cdf, here::here(datadir,"ACD evaluation study","cdf.rds"))

# ---------------------------------------------------------------------------- #
# Join CDF data with patient village GPS

print("CDF patient not in KAMIS:")
print(length(which(!cdf$patient_id %in% pat_wgps$res_patient_code))) # 2

cdf %>%
  inner_join(pat_wgps, by = c("patient_id" = "res_patient_code"), suffix = c("_kamis","_cdf")) -> match

print("Matched patients, with GPS:")
print(summary(match$with_gps))
#    Mode   FALSE    TRUE 
# logical     648    4380 

print("Matched patients, village population:")
print(summary(match$population[match$with_gps]))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#      71    1609    2771    3729    5188   21273       1

# ---------------------------------------------------------------------------- #
# Save merged CDF/GPS data

saveRDS(match, here::here(outdir,"linelist","cdf_gps.rds"))

################################################################################
# sink()
################################################################################