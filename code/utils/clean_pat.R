################################################################################
################################################################################
#
# Import raw KAMIS patient linelist, clean and split by case type
#
################################################################################
################################################################################

clean_pat <- function(file.path,
                       start,
                       end = Sys.Date(),
                       unit_time = "month", #"year"
                       state_incl = c("BIHAR"),
                       log = "cleaning_log_pat.txt"){

  sink(file = log, type = "output")

  ll <- read_csv(file.path) %>%
    dplyr::mutate(across(c("data_entry_state","data_entry_district","data_entry_block",
                           "patient_sc", "patient_village"), toupper)) %>%
    filter(data_entry_state %in% state_incl) %>%
    dplyr::mutate(patient_state = toupper(data_entry_state),
                  patient_district = toupper(data_entry_district),
                  patient_block = toupper(data_entry_block),
                  block_id = patient_block_id,
                  hsc_id = patient_sc_id,
                  vid = patient_village_id)

  print("Duplicates:")
  print(
    ll %>%
      dplyr::summarise(res_patient_code = sum(duplicated(res_patient_code)),
                       patient_id = sum(duplicated(patient_id)))
  )

  delete_rows <- function(x) which(str_detect(x,"delete"))
  delete_rows_id <- unique(unlist(apply(ll,2,delete_rows)))
  print(paste("Entries marked 'DELETE': N =",length(delete_rows_id)))
  if (length(delete_rows_id) > 0) {ll <- ll[-delete_rows_id,]}

  # Tally of diagnoses per patient, per case type
  # tally <- ll %>%
  #   dplyr::group_by(res_patient_code, case_type) %>%
  #   dplyr::tally()

  n.b <- ll %>%
    dplyr::select(data_entry_district,data_entry_block) %>%
    unique() %>%
    nrow()

  n.v <- ll %>%
    dplyr::select(patient_sc, patient_village) %>%
    unique() %>%
    nrow()

  print(paste(nrow(ll),"Patients registered over", n.b,"unique blocks, residing in", n.v, "unique villages"))

  ll %>%
    dplyr::mutate(across(where(is.character), as.factor)) %>%
    dplyr::select(res_patient_code, patient_id, #sex, age, caste, special_caste,
                  # case_type, #treated_for_ka_earlier, referred_by,
                  patient_state, patient_district,
                  patient_block, patient_block_id,
                  # state, district,
                  patient_sc, #patient_sc_id,
                  patient_village, vid, #vil_code,
                  data_entry_state, sid,
                  data_entry_district, did,
                  data_entry_block, bid #, internal_date
                  ) -> final
  sink()

  return(final)

}
