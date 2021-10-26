################################################################################
# Description: Onset to diagnosis by patient characteristics
################################################################################
################################################################################

# Read linelist with matched GPS locations
ll_wgps <- readRDS(here::here("KAMIS/Clean/linelist/linelist_vl_wgps.rds")) %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  mutate(OD_missing = is.na(days_fever),
         agegrp = cut(as.numeric(age), 5),
         excess_delay = (days_fever > 30),
         `hiv_positive?` = factor(`hiv_positive?`,levels = c("Yes","No")),
         `tb_positive?` = factor(`tb_positive?`,levels = c("Yes","No")),
         referred_by = as.factor(referred_by))

ll_filt <- ll_wgps %>%
  filter(diag_year %in% 2016:2019)

nrow(ll_wgps) #35432
nrow(ll_filt) #13714

summary(ll_filt$referred_by)
#              ASHA Medical Institute     Other Doctors         Other PHC              NA's
#              1080               341               283               217             11793

summary(ll_filt$referred_by)/nrow(ll_filt)
#              ASHA Medical Institute     Other Doctors         Other PHC              NA's
#        0.07875164        0.02486510        0.02063585        0.01582325        0.85992417

ggplot(ll_filt, aes(agegrp, days_fever)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log2")

ggplot(ll_filt, aes(caste_category, days_fever)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log2")

ggplot(ll_filt, aes(special_caste, days_fever)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log2")

ggplot(ll_filt, aes(`treated_for_kala-azar_earlier`, days_fever)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log2")

ggplot(ll_filt, aes(`hiv_positive?`, days_fever)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log2")
ggplot(ll_filt, aes(`tb_positive?`, days_fever)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log2")

ggplot(ll_filt, aes(`whether_other_place?`, days_fever)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log2")




dat <- ll_filt %>%
  dplyr::select(excess_delay, diag_year, sex, agegrp, caste_category, `treated_for_kala-azar_earlier`, `hiv_positive?`, `tb_positive?`, latitude, longitude) %>%
  na.omit()
m1 <- glm(excess_delay ~ diag_year + sex + agegrp + caste_category + `treated_for_kala-azar_earlier` + `hiv_positive?`,
    data = dat, family = "binomial")
summary(m1)

dat$res <- residuals(m1)

ggmap(bh_lines, base_layer = ggplot(dat, aes(x = longitude, y = latitude, col = res))) +
  geom_jitter(alpha = 0.5) +
  scale_colour_viridis_c(option = "viridis", na.value = "grey", direction = -1, end = 0.9) +
  labs(x = "", y = "", col = "Residual") -> map_res
map_res
