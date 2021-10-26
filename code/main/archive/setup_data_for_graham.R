

kamis <- readRDS(here::here("data","kamis.rds")) %>%
  dplyr::select(patient_id, district, block, village, population, vil_code, diag_year, diag_date, days_fever)


villinc <- kamis %>%
  group_by(district, block, village, vil_code, population, diag_year) %>%
  count() %>%
  mutate(inc = n/population) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = diag_year, values_from = inc, names_prefix = "vil_inc_") 


cdf <- readRDS(here::here("data","cdf.rds")) %>%
  dplyr::select(patient_id, dist_cdf, block_cdf, vill_cdf, diag_year, diag_date, days_fever)

# names(cdf) <- gsub("_cdf", "",names(cdf))

saveRDS(kamis, here::here("data","kamis_patient_delay.rds"))
saveRDS(vil_inc, here::here("data","kamis_vill_inc.rds"))
saveRDS(cdf, here::here("data","dubey_patient_delay.rds"))


cdf %>% 
  inner_join(kamis, by = c("patient_id","diag_year"), suffix = c("_cdf","_kamis")) %>% #head()
  group_by(vil_code, diag_year) %>%
  summarise(delay_mean = mean(days_fever_cdf, na.rm = T),
            delay_sd = sd(days_fever_cdf, na.rm = T)) %>%
  ungroup() %>%
  left_join(villinc) -> vill_delay_inc

vill_delay_inc %>% 
  filter(diag_year == 2018) %>%
  ggplot(aes(delay_sd, vil_inc_2019)) +
  geom_point(alpha = 0.2)

# + 
#   scale_y_continuous(trans = "log2") +
#   scale_x_continuous(trans = "log2")

kamis %>%
  group_by(vil_code, diag_year) %>%
  summarise(delay_mean = median(days_fever, na.rm = T),
            delay_sd = sd(days_fever, na.rm = T)) %>%
  ungroup() %>%
  left_join(villinc) -> vill_delay_inc_kamis

vill_delay_inc_kamis %>% 
  filter(diag_year == 2018) %>%
  ggplot(aes(delay_mean, vil_inc_2019)) +
  geom_point(alpha = 0.2) +
  geom_smooth()
