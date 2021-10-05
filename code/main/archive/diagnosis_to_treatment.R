################################################################################
# Description: Map onset to diagnosis by village
################################################################################
################################################################################

figdir <- "figures"

# Read linelist with matched GPS locations

ll_wgps <- readRDS(here::here("data","linelist_vl_wgps.rds")) %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  mutate(diag_trt = ,
         DT_missing = is.na(days_fever))

nrow(ll_wgps) #35432

# Setup map context

BH_blk <- readRDS(here::here("data","bihar_block.rds"))

centre <- c(lon = mean(ll_wgps$longitude), lat = mean(ll_wgps$latitude))
extent <- setNames(st_bbox(BH_blk), c("left","bottom","right","top"))
bh_lines <- get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)
# bh_road <- get_map(bbox = extent, source = "google", maptype = "roadmap", zoom = 8)
# bh_sat <- get_map(location = extent, source = "google", maptype = "satellite", zoom = 8)


################################################################################

# Missingness in OD delay

ll_wgps %>%
  group_by(diag_year) %>%
  summarise(N = n(),
            Yes = sum(!OD_missing)/N,
            No = sum(OD_missing)/N) %>%
  pivot_longer(c("Yes","No"), names_to = "With onset") %>%
ggplot(aes(x = diag_year, y = value, fill = `With onset`)) +
  geom_col() +
  labs(x = "Year of diagnosis", y = "Proportion") -> prop_wonset

ggsave(here::here(figdir,"vl_gps_withonset_byyear.png"),
       prop_wonset,
       height = 4, width = 6, units = "in")

# Locations of those with missing onset (2016-2018)

ll_filt <- ll_wgps %>%
  filter(diag_year %in% 2016:2019)

ll_filt %>%
  group_by(OD_missing) %>%
  tally()
# OD_missing     n
# FALSE      11622
# TRUE          98

ggmap(bh_lines, base_layer = ggplot(arrange(ll_filt,OD_missing), aes(x = longitude, y = latitude, col = OD_missing))) +
  geom_jitter(alpha = 0.8) +
  scale_colour_manual(values = c("grey","indianred")) +
  labs(x = "", y = "", col = "Missing onset", caption = "Diagnoses between 2016 and 2019") -> map_missOD

ggsave(here::here(figdir, "village_missOD.png"), map_missOD, height = 7, width = 10, units = "in")


# ---------------------------------------------------------------------------- #
# Observed delays by year

ll_wgps %>%
  group_by(diag_year) %>%
  summarise(mean = mean(days_fever, na.rm = T),
            median = median(days_fever, na.rm = T),
            q1 = quantile(days_fever, p = 0.05, na.rm = T),
            q3 = quantile(days_fever, p = 0.95, na.rm = T)) %>%
  ungroup() -> by_yr


by_yr %>%
  ggplot(aes(x = diag_year, group = 1, y = median, ymin = q1, ymax = q3)) +
  geom_errorbar(width = 0.5) +
  geom_point() +
  geom_line(aes(y = mean), col = "steelblue", lty = "dashed", lwd = 1.2) +
  scale_y_continuous(trans = "log2") +
  labs(x = "Year of diagnosis", y = "Onset to diagnosis (days)") -> OD_byyear

OD_byyear
ggsave(here::here(figdir,"OD_byyear.png"),
       OD_byyear,
       height = 4, width = 6, units = "in")

ll_wgps %>%
  filter(diag_year == 2019) %>%
  ggplot(aes(x = days_fever)) +
  geom_histogram(bins = 50) +
  xlim(c(0,365)) +
  labs(x = "Onset to diagnosis (days)", y = "Frequency")

pal <- c("grey",viridis::viridis(3, end = 0.8, direction = -1))
ll_wgps %>%
  group_by(diag_year) %>%
  summarise(N = n(),
            `> 90 days` = sum(days_fever > 90, na.rm = T)/N,
            `31-90 days` = sum(between(days_fever,31,90), na.rm = T)/N,
            `<= 30 days` = sum(days_fever <= 30, na.rm = T)/N,
            Missing = sum(OD_missing)/N) %>%
  pivot_longer(c(`> 90 days`,`31-90 days`,`<= 30 days`,Missing), names_to = "Delay") %>%
  mutate(Delay = factor(Delay, levels = c("Missing","<= 30 days", "31-90 days", "> 90 days"))) %>%
  ggplot(aes(x = diag_year, y = value, fill = Delay)) +
  geom_col() +
  scale_fill_manual(values = pal) +
  labs(x = "Year of diagnosis", y = "Proportion") -> prop_gt90days

prop_gt90days
ggsave(here::here(figdir,"prop_excessive_byyear.png"),
       prop_gt90days,
       height = 4, width = 6, units = "in")


# ---------------------------------------------------------------------------- #
# Map observed delays

ll_filt <- ll_filt %>%
  filter(!is.na(days_fever)) %>%
  mutate(ODcat5 = cut(days_fever, c(0,15,30,90,180,730), include.lowest = TRUE, ordered_result = TRUE),
         ODcat3 = cut(days_fever, c(0,15,90,730), include.lowest = TRUE, ordered_result = TRUE),
         ODexcess = (days_fever > 90)) %>%
  arrange(days_fever)

summary(ll_filt$days_fever)
summary(ll_filt$ODcat5)
summary(ll_filt$ODcat3)
summary(ll_filt$ODexcess)

# 38,891 diagnosed patients with village GPS

ggmap(bh_lines, base_layer = ggplot(ll_filt, aes(x = longitude, y = latitude, col = ODcat5))) +
  geom_jitter(alpha = 0.5) +
  scale_colour_viridis_d(option = "viridis", na.value = "grey", direction = -1, end = 0.9) +
  labs(x = "", y = "", col = "Delay") -> map_cat5

ggsave(here::here(figdir, "village_OD_cat5.png"), map_cat5, height = 7, width = 10, units = "in")

ggmap(bh_lines, base_layer = ggplot(ll_filt, aes(x = longitude, y = latitude, col = ODcat3))) +
  geom_jitter(alpha = 0.5) +
  scale_colour_viridis_d(option = "viridis", na.value = "grey", direction = -1, end = 0.9) +
  labs(x = "", y = "", col = "Delay") -> map_cat3

ggsave(here::here(figdir, "village_OD_cat3.png"), map_cat3, height = 7, width = 10, units = "in")

ggmap(bh_lines, base_layer = ggplot(arrange(ll_filt,ODexcess), aes(x = longitude, y = latitude, col = ODexcess))) +
  geom_jitter(alpha = 0.5) +
  scale_colour_manual(values = c("grey","indianred")) +
  labs(x = "", y = "", col = "Delay > 90 days", subtitle = "Diagnoses between 2016 and 2018") -> map_excess

ggsave(here::here(figdir, "village_OD_excess.png"), map_excess, height = 7, width = 10, units = "in")


# Median delay per village
ll_filt %>%
  group_by(vil_code, latitude, longitude) %>%
  summarise(N = n(),
            medianOD = median(days_fever, na.rm = T),
            propgt30 = mean((days_fever > 30), na.rm = T)*100,
            propgt90 = mean((days_fever > 90), na.rm = T)*100) %>%
  arrange(propgt90) -> by_village

ggmap(bh_lines, base_layer = ggplot(by_village, aes(x = longitude, y = latitude, col = medianOD))) +
  geom_point() +
  # scale_colour_gradient2(low = pal3[3], mid = pal3[2], high = pal3[1], midpoint = log(30,2), trans = "log2") +
  scale_colour_viridis_c(option = "viridis", na.value = "grey", direction = -1, trans = "log2", end = 0.9) +
  labs(x = "", y = "", col = "Median delay") -> map_medOD

ggsave(here::here(figdir, "village_medianOD.png"), map_medOD, height = 7, width = 10, units = "in")


ggmap(bh_lines, base_layer = ggplot(by_village, aes(x = longitude, y = latitude, col = propgt30))) +
  geom_point(alpha = 0.7) +
  scale_colour_viridis_c(option = "magma", na.value = "grey", direction = -1, end = 0.9) +
  labs(x = "", y = "", col = "% > 30 days") -> map_propgt30

ggsave(here::here(figdir, "village_propgt30.png"), map_propgt30, height = 7, width = 10, units = "in")


ggmap(bh_lines, base_layer = ggplot(by_village, aes(x = longitude, y = latitude, col = propgt90))) +
  geom_point(alpha = 0.7) +
  scale_colour_viridis_c(option = "magma", na.value = "grey", direction = -1, end = 0.9) +
  labs(x = "", y = "", col = "% > 90 days") -> map_propgt90

ggsave(here::here(figdir, "village_propgt90.png"), map_propgt90, height = 7, width = 10, units = "in")


# ---------------------------------------------------------------------------- #
# By block

matched_names <- read.csv(here::here("data","MatchedBlockNames.csv"), header = T) %>%
  mutate(across(everything(), toupper))
matched_names$kamis_master_block[matched_names$kamis_master_dist == "MADHUBANI" &
                                   matched_names$kamis_master_block == "MADHUBANI"] <- "RAHIKA"

ll_filt %>%
  group_by(district, block) %>%
  summarise(N = n(),
            medianOD = median(days_fever, na.rm = T),
            propgt90 = mean((days_fever > 90), na.rm = T)*100,
            propgt30 = mean((days_fever > 30), na.rm = T)*100,) %>%
  mutate(ODcat3 = cut(medianOD, c(0,15,90,730), include.lowest = TRUE, ordered_result = TRUE)) %>%
  ungroup() -> by_block

by_block <- BH_blk %>%
  left_join(matched_names) %>%
  left_join(by_block, by = c("kamis_master_dist" = "district", "kamis_master_block" = "block"))

pal3 <- viridis::viridis(3)
ggplot() +
  geom_sf(data = BH_blk, fill = "grey") +
  geom_sf(data = by_block, aes(geometry = geometry, fill = medianOD, alpha = N)) +
  # scale_fill_gradient2(low = pal3[3], mid = pal3[2], high = pal3[1], midpoint = log(30,2), trans = "log2") +
  scale_fill_viridis_c(option = "viridis", na.value = "grey", direction = -1, trans = "log2", end = 0.9) +
  labs(fill = "Median delay") -> blk_medOD

blk_medOD
ggsave(here::here(figures, "block_medianOD.png"), blk_medOD, height = 7, width = 10, units = "in")

# ggplot() +
#   geom_sf(data = BH_blk, fill = "grey") +
#   geom_sf(data = filter(by_block, !is.na(ODcat3)), aes(geometry = geometry, fill = ODcat3)) +
#   scale_fill_viridis_d(option = "viridis", na.value = "grey", direction = -1, end = 0.9) +
#   labs(fill = "Median delay") -> blk_ODcat3
#
# blk_ODcat3
# ggsave(here::here(figures, "block_OD_cat3.png"), blk_ODcat3, height = 7, width = 10, units = "in")

ggplot() +
  geom_sf(data = BH_blk) +
  geom_sf(data = by_block, aes(geometry = geometry, fill = propgt30, alpha = N)) +
  scale_fill_viridis_c(option = "magma", na.value = "grey", direction = -1, end = 0.9) +
  labs(fill = "% > 30 days") -> blk_probpgt30

blk_probpgt30
ggsave(here::here(figures, "block_propgt30.png"), blk_probpgt30, height = 7, width = 10, units = "in")


ggplot() +
  geom_sf(data = BH_blk) +
  geom_sf(data = by_block, aes(geometry = geometry, fill = propgt90)) +
  scale_fill_viridis_c(option = "magma", na.value = "grey", direction = -1, end = 0.9) +
  labs(fill = "% > 90 days") -> blk_probpgt90

blk_probpgt90
ggsave(here::here(figures, "block_propgt90.png"), blk_probpgt90, height = 7, width = 10, units = "in")

# Join two block maps
blk_probpgt30 <- blk_probpgt30 + theme(axis.text = element_blank())
blk_medOD <- blk_medOD + theme(axis.text = element_blank())
combined <- blk_medOD + blk_probpgt30 +
  plot_annotation(caption = paste0("Diagnoses between ",
                                   range(ll_filt$date_of_diagnosis)[1],
                                   " and ",
                                   range(ll_filt$date_of_diagnosis)[2]))
ggsave(here::here(figures, "OD_by_block.png"), combined, height = 7, width = 14, units = "in")

# ---------------------------------------------------------------------------- #
# By block and diagnosis year

ll_filt %>%
  group_by(district, block, diag_year) %>%
  summarise(N = n(),
            medianOD = median(days_fever, na.rm = T),
            propgt30 = mean((days_fever > 30), na.rm = T)*100,
            propgt90 = mean((days_fever > 90), na.rm = T)*100) %>%
  arrange(propgt90) %>%
  ungroup() -> by_block_yr

by_block_yr <- BH_blk %>%
  left_join(matched_names) %>%
  left_join(by_block_yr, by = c("kamis_master_dist" = "district", "kamis_master_block" = "block")) %>%
  filter(!is.na(diag_year))

ggplot() +
  geom_sf(data = BH_blk) +
  geom_sf(data = by_block_yr, aes(geometry = geometry, fill = medianOD)) +
  facet_wrap(~diag_year) +
  scale_fill_viridis_c(option = "viridis", na.value = "grey", direction = -1, trans = "log2", end = 0.9) +
  labs(fill = "Median delay") + theme(axis.text = element_blank()) -> blkyr_medOD

blkyr_medOD
ggsave(here::here(figures, "blockyr_medianOD.png"), blkyr_medOD, height = 7, width = 10, units = "in")

ggplot() +
  geom_sf(data = BH_blk) +
  geom_sf(data = by_block_yr, aes(geometry = geometry, fill = propgt30)) +
  facet_wrap(~diag_year) +
  scale_fill_viridis_c(option = "magma", na.value = "grey", direction = -1, end = 0.9) +
  labs(fill = "% > 30 days") + theme(axis.text = element_blank()) -> blkyr_probpgt30

blkyr_probpgt30
ggsave(here::here(figures, "blockyr_propgt30.png"), blkyr_probpgt30, height = 7, width = 10, units = "in")


################################################################################
################################################################################
