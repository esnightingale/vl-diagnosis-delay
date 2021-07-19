################################################################################
# Description: Onset to diagnosis by patient characteristics
################################################################################
################################################################################

# Read linelist with matched GPS locations
ll_wgps <- readRDS(here::here("data","linelist_vl_wgps.rds")) %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  mutate(OD_missing = is.na(days_fever))

ll_filt <- ll_wgps %>%
  filter(diag_year %in% 2016:2019)

nrow(ll_wgps) #35432
nrow(ll_filt) #13714

# Village IRS targeting and block endemicity
village_chars <- readRDS(here::here("data","village_irs_endemicity.rds")) %>%
  mutate(vil_code = as.factor(vil_code))

################################################################################

block_endm <- village_chars %>%
  dplyr::select(vil_code, block_endm_2013:block_endm_2018) %>%
  pivot_longer(-vil_code, 
               values_to = "block_endemicity", 
               values_transform = list(block_endemicity = as.factor), 
               names_to = "diag_year", 
               names_pattern = "block_endm_(.*)",
               names_transform = list(diag_year = as.factor))

summary(block_endm)

vil_irs <- village_chars %>%
  dplyr::select(vil_code:IRS_2018_R1,total_sprays_1318) %>%
  pivot_longer(-c(vil_code, total_sprays_1318), 
               values_to = "irs_targeted", 
               names_to = c("diag_year","round"), 
               names_transform = list(diag_year = as.factor, 
                                      round = as.factor),
               names_pattern = "IRS_(.*)_(.*)")

# Note: 2018 R2 did not occur due to poor insecticide quality
summary(vil_irs)

vil_irs_yr <- vil_irs %>%
  filter(!(diag_year == 2018 & round == "R2")) %>%
  group_by(vil_code, total_sprays_1318, diag_year) %>%
  summarise(irs_targeted = as.factor(sum(irs_targeted, na.rm = T))) %>%
  ungroup()
# Note: 2018 R2 did not occur due to poor insecticide quality

summary(vil_irs_yr)


ll_wchars <- ll_wgps %>%
  # IRS/block endemicity included up to 2018
  filter(!diag_year %in% c(2019:2021)) %>%
  left_join(block_endm) %>%
  left_join(vil_irs_yr) 

# ---------------------------------------------------------------------------- #
# Plot distribution of delay by whether or not block considered endemic and 
# village targeted with IRS

ll_wchars %>%
  ggplot(aes(block_endemicity, days_fever)) +
  scale_y_continuous(trans = "log10") +
  geom_point(alpha = 0.2, position = position_jitter(width = 0.2, height = 0.1), col = "grey") +
  geom_boxplot(aes(col = block_endemicity), alpha = 0, lwd = 1)


ll_wchars %>%
  ggplot(aes(irs_targeted, days_fever)) +
  scale_y_continuous(trans = "log10") +
  geom_point(alpha = 0.2, position = position_jitter(width = 0.2, height = 0.1), col = "grey") +
  geom_boxplot(aes(col = irs_targeted), alpha = 0, lwd = 1)


