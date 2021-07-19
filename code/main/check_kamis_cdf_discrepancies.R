

match <- dplyr::filter(match, with_gps)

# Setup map context
blockmap <- readRDS(here::here("Diagnosis delay","data","bihar_block.rds"))

centre <- c(lon = mean(ll_wgps$longitude), lat = mean(ll_wgps$latitude))
extent <- setNames(st_bbox(blockmap), c("left","bottom","right","top"))
bh_lines <- get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)


# ---------------------------------------------------------------------------- #
# Check agreement between KAMIS and CDF on patient characteristics

match %>%
  dplyr::mutate(d_match = (district == dist_cdf)) %>%
  dplyr::pull(d_match) %>%
  summary()

match %>%
  dplyr::mutate(db_match = (district == dist_cdf &
                              block == block_cdf)) %>%
  dplyr::pull(db_match) %>%
  summary()

match %>%
  dplyr::mutate(dbv_match = (district == dist_cdf &
                               block == block_cdf &
                               village == vill_cdf)) %>%
  dplyr::pull(dbv_match) %>%
  summary()

match %>%
  group_by(sex_kamis, sex_cdf) %>%
  tally()

match %>%
  ggplot(aes(age_kamis, age_cdf)) +
  geom_jitter(alpha = 0.1) +
  geom_smooth() +
  labs(x = "KAMIS age", y = "CDF age") +
  theme_classic()

match %>%
  group_by(`treated_for_kala-azar_earlier`, prv_tx_ka) %>%
  tally()

match %>%
  group_by(hiv_kamis, hiv_cdf) %>%
  tally()

match %>%
  ggplot(aes(diag_date_kamis, diag_date_cdf)) +
  geom_jitter(alpha = 0.1) +
  geom_smooth() +
  labs(x = "KAMIS diag date", y = "CDF diag date") +
  theme_classic()

match %>%
  group_by(marg_caste, caste4_r) %>%
  tally()


match %>%
  group_by(referred_by_kamis, referred_by_cdf) %>%
  tally()

match %>%
  ggplot(aes(days_fever_kamis, days_fever_cdf)) +
  geom_jitter(alpha = 0.1) +
  geom_smooth() +
  labs(x = "KAMIS delay", y = "CDF delay") +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2") +
  theme_classic()

match %>%
  group_by(gt30_kamis, gt30_cdf) %>%
  tally()


match %>%
  group_by(gt90_kamis, gt90_cdf) %>%
  tally()

# ---------------------------------------------------------------------------- #
# Compare OD between KAMIS and CDF
# Continuous

match %>%
  mutate(diff = days_fever_kamis - days_fever_cdf) %>%
  dplyr::select(days_fever_kamis, days_fever_cdf, diff) %>%
  summary()

# days_fever_kamis days_fever_cdf        diff
# Min.   :  2.00   Min.   :  4.00   Min.   :-390.00
# 1st Qu.: 49.00   1st Qu.: 25.00   1st Qu.:   9.00
# Median : 68.00   Median : 30.00   Median :  29.00
# Mean   : 67.12   Mean   : 45.01   Mean   :  22.14
# 3rd Qu.: 71.00   3rd Qu.: 58.00   3rd Qu.:  41.00
# Max.   :159.00   Max.   :510.00   Max.   : 145.00
# NA's   :18                        NA's   :18

match %>%
  pivot_longer(c("days_fever_kamis","days_fever_cdf")) %>%
  ggplot(aes(value)) +
  geom_histogram() +
  geom_vline(aes(xintercept = mean(value, na.rm = T)), col = "red", lty = "dashed") +
  geom_vline(aes(xintercept = median(value, na.rm = T)), col = "dodgerblue", lty = "dashed") +
  facet_wrap(~name)

match %>%
  ggplot(aes(days_fever_kamis, days_fever_cdf)) +
  geom_jitter(alpha = 0.2) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "KAMIS value", y = "CDF value") -> compare_num
compare_num
ggsave(here::here(figdir, "days_fever_kamis_cdf.png"), compare_num, height = 6, width = 7, units = "in")


ggmap(bh_lines, base_layer = ggplot(arrange(match, days_fever_cdf), aes(x = longitude, y = latitude, col = days_fever_cdf))) +
  geom_jitter() +
  scale_colour_viridis_c(trans = "log10", direction = -1, option = "plasma") +
  labs(x = "", y = "", col = "Delay", subtitle = "Diagnoses between Jan 2018 - July 2019") -> map_dubey_delay
map_dubey_delay


# ---------------------------------------------------------------------------- #
# Categorical

gt30tab <- prop.table(table(match$gt30_kamis, match$gt30_cdf))
gt30tab
gt90tab <- prop.table(table(match$gt90_kamis, match$gt90_cdf))
gt90tab

match %>%
  # filter(!is.na(days_fever)) %>%
  group_by(gt30_kamis, gt30_cdf) %>%
  tally() %>%
  ungroup() %>%
  mutate(value = paste0(n, " (",round(n/nrow(match),2),")")) %>%
  pivot_wider(-n,values_from = value, names_from = gt30_cdf)
#           `<= 30`    `> 30`
# <= 30      2493 (0.5) 1169 (0.23)
# > 30       297 (0.06) 1049 (0.21)

match %>%
  # filter(!is.na(days_fever)) %>%
  group_by(gt90_kamis, gt90_cdf) %>%
  tally() %>%
  ungroup()  %>%
  mutate(value = paste0(n, " (", round(n/nrow(match),2),")")) %>%
  pivot_wider(-n,values_from = value, names_from = gt90_cdf)
#            `<= 90`    `> 90`
# <= 90      4539 (0.9) 248 (0.05)
# > 90       117 (0.02) 104 (0.02)

# ---------------------------------------------------------------------------- #
# Compare distribution of excessive (>90) delays

ggmap(bh_lines, base_layer = ggplot(arrange(match,gt90_cdf), aes(x = longitude, y = latitude, col = gt90_cdf))) +
  geom_jitter() +
  scale_colour_manual(values = c("grey", "indianred")) +
  guides(col = FALSE) +
  labs(x = "", y = "", col = "Delay > 90 days", title = "CDF") -> map_cdf_gt90
map_cdf_gt90
#
# # Filter KAMIS to same time period
# ll_filt <- ll_wgps %>%
#   filter(between(diag_date, min(match$Date_Diag), max(match$Date_Diag)) & !is.na(days_fever)) %>%
#   mutate(ODexcess = (days_fever > 90))

ggmap(bh_lines, base_layer = ggplot(arrange(match,gt90_kamis), aes(x = longitude, y = latitude, col = gt90_kamis))) +
  geom_jitter() +
  scale_colour_manual(values = c("grey","indianred")) +
  labs(x = "", y = "", col = "Delay > 90 days", title = "KAMIS", caption = "Diagnoses between Jan 2018-July 2019") -> map_kamis_gt90

compare_excess <- map_cdf_gt90 + map_kamis_gt90
compare_excess
ggsave(here::here(figdir, "compare_data_delaygt90.png"), compare_excess, height = 7, width = 16, units = "in")


# ---------------------------------------------------------------------------- #
# By block

match %>%
  group_by(district, block) %>%
  summarise(N = n(),
            medianOD = median(days_fever_cdf, na.rm = T),
            propgt90 = mean((days_fever_cdf > 90), na.rm = T)*100,
            propgt30 = mean((days_fever_cdf > 30), na.rm = T)*100) %>%
  mutate(ODcat3 = cut(medianOD, c(0,15,90,730), include.lowest = TRUE, ordered_result = TRUE)) %>%
  ungroup() -> by_block

by_block <- blockmap %>%
  left_join(by_block, by = c("kamis_master_dist" = "district", "kamis_master_block" = "block")) %>%
  rowwise() %>%
  mutate(pop = median(c_across(`2018`:`2019`)),
         inc = N*1e4/pop) %>%
  ungroup()

# pal3 <- viridis::viridis(3)
ggplot() +
  geom_sf(data = blockmap) +
  geom_sf(data = by_block, aes(geometry = geometry, fill = medianOD)) +
  scale_fill_viridis_c(option = "viridis", na.value = "white", direction = -1, trans = "log2", end = 0.9) +
  labs(fill = "Median delay") +
  theme(axis.text = element_blank()) -> blk_medOD

blk_medOD
ggsave(here::here(figdir, "block_medianOD.png"), blk_medOD, height = 7, width = 10, units = "in")

ggplot() +
  geom_sf(data = blockmap) +
  geom_sf(data = by_block, aes(geometry = geometry, fill = propgt30)) +
  scale_fill_viridis_c(option = "magma", na.value = "white", direction = -1, end = 0.8) +
  labs(fill = "% > 30 days") +
  theme(axis.text = element_blank()) -> blk_propgt30

blk_propgt30
ggsave(here::here(figdir, "block_propgt30.png"), blk_propgt30, height = 7, width = 10, units = "in")


ggplot() +
  geom_sf(data = blockmap) +
  geom_sf(data = by_block, aes(geometry = geometry, fill = propgt90)) +
  scale_fill_viridis_c(option = "magma", na.value = "white", direction = -1, end = 0.9) +
  labs(fill = "% > 90 days") +
  theme(axis.text = element_blank()) -> blk_propgt90

blk_propgt90
ggsave(here::here(figdir, "block_propgt90.png"), blk_propgt90, height = 7, width = 10, units = "in")

# Join two block maps
combined <- blk_medOD + blk_propgt30 +
  plot_annotation(caption = paste0("Diagnoses between ",
                                   range(match$Date_Diag)[1],
                                   " and ",
                                   range(match$Date_Diag)[2]))
ggsave(here::here(figdir, "OD_by_block.png"), combined, height = 7, width = 14, units = "in")

# ---------------------------------------------------------------------------- #
# Bivariate map of delay and incidence

by_block %>%
  filter(!is.na(propgt30)) %>%
  mutate(inc = replace_na(inc, 0)) -> plotdat
plotdat <- bi_class(plotdat, x = propgt30, y = inc, style = "quantile", dim = 3)

ggplot() +
  geom_sf(data = blockmap, fill = "white") +
  geom_sf(data = plotdat, aes(fill = bi_class), show.legend = FALSE) +
  bi_scale_fill(pal = "DkBlue", dim = 3, na.value = "white") +
  bi_theme() -> biv_map

legend <- bi_legend(pal = "DkBlue",
                    dim = 3,
                    xlab = "% > 30 days",
                    ylab = "Cases/10,000",
                    size = 10)

final <- cowplot::ggdraw() +
  cowplot::draw_plot(biv_map, 0, 0, 1, 1) +
  cowplot::draw_plot(legend, 0.6, 0.7, 0.3, 0.3)
final

ggsave(here::here(figdir, "block_propgt30_byinc_bv.png"), final, height = 7, width = 10, units = "in")


# ---------------------------------------------------------------------------- #
# Semi-variogram

dat <- filter(match, !is.na(longitude) & !is.na(days_fever_cdf))

# convert simple data frame into a spatial data frame object
coordinates(dat) <- ~ longitude + latitude

vg <- variogram(days_fever_cdf~1, data = dat)
plot(vg)

vgmod <- vgm(psill =  800, model = "Mat", nugget = 1000, range = 0.5)
plot(vg, model = vgmod)

vgfit <- fit.variogram(vg, model = vgmod)
plot(vg, model = vgfit)

vgfit
# model    psill    range kappa
# Nug 936.1536 0.000000   0.0
# Mat 959.9353 0.629266   0.5

png(here::here(figdir, "OD_semivariogram.png"), height = 500, width = 600)
plot(vg, model = vgfit)
dev.off()


