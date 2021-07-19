################################################################################
# Description: Import Case Details Form data from Dubey et al. evaluation of
# ACD and compare with KAMIS
################################################################################
################################################################################

raw.path <- here::here("Data","KAMIS","Raw")
figdir <- here::here("Diagnosis delay","figures","CDF")
outdir <- here::here("Diagnosis delay","data")

# ---------------------------------------------------------------------------- #
# Read CDF data from ACD study, from  SAS format
# Remove SAS formats and var labels, select vars of interest and rename/reformat

cdf <- haven::read_sas(here::here("Data","ACD",
                                     "acddata_102620/acddata_102620.sas7bdat")) %>%
  haven::zap_formats() %>%
  haven::zap_label() %>%
  haven::zap_label() %>%
  dplyr::as_tibble() %>%
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
  dplyr::mutate(across(dist_cdf:vill_cdf, toupper)) %>%
  dplyr::mutate(across(c(dist_cdf:vill_cdf,sex, hiv, caste4:referred_by),
                       function(x) factor(x, exclude = c("","99")))) %>%
  dplyr::mutate(age = as.numeric(age),
                diag_month = lubridate::floor_date(diag_date, "month"),
                diag_year = as.factor(lubridate::year(diag_date)),
                days_fever = as.numeric(days_fever),
                gt30 = (days_fever > 30),
                gt90 = (days_fever > 90))

summary(cdf)

saveRDS(cdf, here::here(outdir,"cdf.rds"))

# ---------------------------------------------------------------------------- #
# Read linelist with matched GPS locations

kamis <- readRDS(here::here("Diagnosis delay","data",
                              "linelist_vl_wgps.rds")) %>%
  dplyr::rename(hiv = `hiv_positive?`,
                tb = `tb_positive?`,
                whether_other_place = `whether_other_place?...37`) %>%
  dplyr::mutate(age = as.numeric(age),
                marg_caste = (caste_category %in% c("SC","ST")),
                gt30 = (days_fever > 30),
                gt90 = (days_fever > 90))

summary(kamis)

saveRDS(kamis, here::here(outdir,"kamis.rds"))

# ---------------------------------------------------------------------------- #
# Check matching between full linelist and study data

length(which(!acddat$patient_id %in% ll_wgps$patient_id))
length(which(!ll_wgps$patient_id %in% acddat$patient_id))

# + Joy reports that two patients in Dubey data but not KAMIS were removed from the linelist. Reason is unknown.
# + 634 patients attributed to villages with no GPS data
# + One patient attributed to a village with also missing population estimate - Baijupatti village in Bhargama, Araria.

ll_wgps %>%
  anti_join(acddat, by = "patient_id") -> nonmatch1

acddat %>%
  anti_join(ll_wgps, by = "patient_id") -> nonmatch2

ll_wgps %>%
  inner_join(acddat, by = "patient_id", suffix = c("_kamis","_cdf")) -> match

summary(match$with_gps)
summary(match$population)

# ---------------------------------------------------------------------------- #
# Save merged KAMIS/GPS/CDF data for analysis

saveRDS(match, here::here(outdir,"kamis_cdf.rds"))

################################################################################
################################################################################
