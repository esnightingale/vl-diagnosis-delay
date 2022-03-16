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
  dplyr::mutate(inc_2017_gt0 = (replace_na(X2017*1000/population, 0) > 0),
                IRS_2017 = ((insecticide_2017_R1 != "" | insecticide_2017_R2 != "")),
                block_endm_2017 = (block_endm_2017 > 0)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(vl_affect_1517 = (sum(c_across(X2015:X2017)) > 0)) %>%
  dplyr::ungroup() %>%
  dplyr::select(country:vil_code, vl_affect_1517, inc_2017_gt0, IRS_2017, block_endm_2017) %>%
  dplyr::mutate(across(where(is.character), as.factor)) %>%
  dplyr::filter(!(is.na(longitude) | is.na(latitude)))

# 12894 villages with GPS
# vl_affected_village defined as any 

village_gps %>%
  dplyr::right_join(clean_p, by = c("district" = "data_entry_district",
                                    "block" = "data_entry_block",
                                    "hsc" = "patient_sc",
                                    "village" ="patient_village")) -> pat_wgps

summary(!is.na(pat_wgps$longitude))
#    Mode   FALSE    TRUE 
# logical    4123   38584 

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
  dplyr::rename(patient_id = pid,
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
  dplyr::anti_join(pat_wgps, by = c("patient_id" ="res_patient_code")) %>%
  nrow() # 2

pat_wgps %>%
  dplyr::select(country:patient_id) %>%
  inner_join(cdf, 
             by = c("res_patient_code" = "patient_id")) -> match

print("Matched patients, with GPS:")
print(summary(!is.na(match$longitude)))
#    Mode   FALSE    TRUE 
# logical     648    4380 

match <- dplyr::filter(match, !is.na(longitude))

print("Matched patients, village population:")
print(summary(match$population[!is.na(match$longitude)]))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#      71    1609    2771    3729    5188   21273       1

# Define village index and drop unwanted columns
match <- match %>%
  dplyr::mutate(v = as.numeric(as.factor(as.character(vil_code))),
                match_code = sample(1:nrow(match))) %>%
  dplyr::select(-c(hsc, village, vil_code, patient_id)) 

# Separate for secure storage of village and individual chars 
code <- match %>% 
  dplyr::select(match_code, v)

pat <- match %>%
  dplyr::select(match_code, res_patient_code:id)

vill <- match %>%
  dplyr::select(v, country:block_endm_2017) %>%
  distinct()

# ---------------------------------------------------------------------------- #
# Save patient/village data and matching code

saveRDS(pat, here::here(outdir,"ll","pat.rds"))
saveRDS(vill, here::here(outdir,"ll","vill.rds"))
saveRDS(code, "E:/vl-diagnosis-delay/data/code.rds")

################################################################################
# sink()
################################################################################