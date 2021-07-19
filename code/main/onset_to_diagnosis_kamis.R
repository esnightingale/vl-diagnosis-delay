################################################################################
# Description: Map onset to diagnosis by village
################################################################################
################################################################################

figdir <- "figures/KAMIS"

# Read linelist with matched GPS locations
ll_wgps <- readRDS(here::here("data","linelist_vl_wgps.rds")) %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  mutate(OD_missing = is.na(days_fever))

nrow(ll_wgps) #35432

# Filter to 2016-2019
ll_filt <- ll_wgps %>%
  filter(diag_year %in% 2016:2019)

# Setup map context
blockmap <- readRDS(here::here("data","bihar_block.rds"))

centre <- c(lon = mean(ll_wgps$longitude), lat = mean(ll_wgps$latitude))
extent <- setNames(st_bbox(blockmap), c("left","bottom","right","top"))
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


################################################################################
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


################################################################################
# Observed delays by location

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

# ---------------------------------------------------------------------------- #
# Individuals

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


# ---------------------------------------------------------------------------- #
# By village

ll_filt %>%
  group_by(vil_code, latitude, longitude) %>%
  summarise(N = n(),
            pop = unique(population),
            inc = replace_na(N*1e3/unique(pop),0),
            medianOD = median(days_fever, na.rm = T),
            propgt30 = mean((days_fever > 30), na.rm = T)*100,
            propgt90 = mean((days_fever > 90), na.rm = T)*100) %>%
  arrange(propgt90) -> by_village

ggmap(bh_lines, base_layer = ggplot(by_village, aes(x = longitude, y = latitude, col = medianOD))) +
  geom_point(pch = 19, alpha = 0.5) +
  # scale_colour_gradient2(low = pal3[3], mid = pal3[2], high = pal3[1], midpoint = log(30,2), trans = "log2") +
  scale_colour_viridis_c(option = "viridis", na.value = "grey", direction = -1, trans = "log2", end = 0.9) +
  scale_alpha_continuous(trans = "log10") +
  labs(x = "", y = "", col = "Median delay") -> map_medOD

map_medOD
ggsave(here::here(figdir, "village_medianOD.png"), map_medOD, height = 7, width = 10, units = "in")


ggmap(bh_lines, base_layer = ggplot(by_village, aes(x = longitude, y = latitude, col = propgt30))) +
  geom_point(alpha = 0.7) +
  scale_colour_viridis_c(option = "magma", na.value = "grey", direction = -1, end = 0.9) +
  scale_alpha_continuous(trans = "log10") +
  labs(x = "", y = "", col = "% > 30 days") -> map_propgt30

map_propgt30
ggsave(here::here(figdir, "village_propgt30.png"), map_propgt30, height = 7, width = 10, units = "in")


ggmap(bh_lines, base_layer = ggplot(by_village, aes(x = longitude, y = latitude, col = propgt90))) +
  geom_point(alpha = 0.7) +
  scale_colour_viridis_c(option = "magma", na.value = "grey", direction = -1, end = 0.9) +
  scale_alpha_continuous(trans = "log10") +
  labs(x = "", y = "", col = "% > 90 days") -> map_propgt90

map_propgt90
ggsave(here::here(figdir, "village_propgt90.png"), map_propgt90, height = 7, width = 10, units = "in")

# ---------------------------------------------------------------------------- #
# By block

ll_filt %>%
  group_by(district, block) %>%
  summarise(N = n(),
            medianOD = median(days_fever, na.rm = T),
            propgt90 = mean((days_fever > 90), na.rm = T)*100,
            propgt30 = mean((days_fever > 30), na.rm = T)*100) %>%
  mutate(ODcat3 = cut(medianOD, c(0,15,90,730), include.lowest = TRUE, ordered_result = TRUE)) %>%
  ungroup() -> by_block

by_block <- blockmap %>%
  left_join(by_block, by = c("kamis_master_dist" = "district", "kamis_master_block" = "block")) %>%
  rowwise() %>%
  mutate(pop = median(c_across(`2013`:`2021`)),
         inc = N*1e4/pop) %>%
  ungroup() 

# pal3 <- viridis::viridis(3)
ggplot() +
  geom_sf(data = blockmap, fill = "grey") +
  geom_sf(data = by_block, aes(geometry = geometry, fill = medianOD)) +
  # scale_fill_gradient2(low = pal3[3], mid = pal3[2], high = pal3[1], midpoint = log(30,2), trans = "log2") +
  scale_fill_viridis_c(option = "viridis", na.value = "white", direction = -1, trans = "log2", end = 0.9) +
  scale_alpha_continuous(trans = "log10") +
  labs(fill = "Median delay") -> blk_medOD

blk_medOD
ggsave(here::here(figdir, "block_medianOD.png"), blk_medOD, height = 7, width = 10, units = "in")

ggplot() +
  geom_sf(data = blockmap) +
  geom_sf(data = by_block, aes(geometry = geometry, fill = propgt30)) +
  scale_fill_viridis_c(option = "magma", na.value = "white", direction = -1, end = 0.8) +
  scale_alpha_continuous(trans = "log10") +
  labs(fill = "% > 30 days") -> blk_propgt30

blk_propgt30
ggsave(here::here(figdir, "block_propgt30.png"), blk_propgt30, height = 7, width = 10, units = "in")


ggplot() +
  geom_sf(data = blockmap) +
  geom_sf(data = by_block, aes(geometry = geometry, fill = propgt90)) +
  scale_fill_viridis_c(option = "magma", na.value = "white", direction = -1, end = 0.9) +
  scale_alpha_continuous(trans = "log10") +
  labs(fill = "% > 90 days") -> blk_propgt90

blk_propgt90
ggsave(here::here(figdir, "block_propgt90.png"), blk_propgt90, height = 7, width = 10, units = "in")

# Join two block maps
blk_propgt30 <- blk_propgt30 + theme(axis.text = element_blank())
blk_medOD <- blk_medOD + theme(axis.text = element_blank())
combined <- blk_medOD + blk_propgt30 +
  plot_annotation(caption = paste0("Diagnoses between ",
                                   range(ll_filt$date_of_diagnosis)[1],
                                   " and ",
                                   range(ll_filt$date_of_diagnosis)[2]))
ggsave(here::here(figdir, "OD_by_block.png"), combined, height = 7, width = 14, units = "in")

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

ggsave(here::here(figdir, "block_propgt30_byinc_bv.png"), final, height = 7, width = 14, units = "in")

# ---------------------------------------------------------------------------- #
# By block and diagnosis year

ll_wgps %>%
  filter(diag_year != 2021) %>%
  group_by(district, block, diag_year) %>%
  summarise(N = n(),
            medianOD = median(days_fever, na.rm = T),
            propgt30 = mean((days_fever > 30), na.rm = T)*100,
            propgt90 = mean((days_fever > 90), na.rm = T)*100) %>%
  arrange(propgt90) %>%
  ungroup() -> by_block_yr

by_block_yr <- blockmap %>%
  pivot_longer(`2013`:`2021`, names_to = "diag_year", values_to = "pop") %>%
  left_join(by_block_yr, by = c("kamis_master_dist" = "district", "kamis_master_block" = "block", "diag_year")) %>%
  mutate(inc = replace_na(N*1e4/pop, 0)) %>%
  filter(!is.na(diag_year))

ggplot() +
  geom_sf(data = blockmap) +
  geom_sf(data = filter(by_block_yr, diag_year %in% 2016:2019), aes(geometry = geometry, fill = medianOD)) +
  facet_wrap(~diag_year) +
  scale_fill_viridis_c(option = "viridis", na.value = "grey", direction = -1, trans = "log2", end = 0.9) +
  labs(fill = "Median delay") + theme(axis.text = element_blank()) -> blkyr_medOD

blkyr_medOD
ggsave(here::here(figdir, "blockyr_medianOD.png"), blkyr_medOD, height = 7, width = 10, units = "in")

ggplot() +
  geom_sf(data = blockmap) +
  geom_sf(data = filter(by_block_yr, diag_year %in% 2016:2019), aes(geometry = geometry, fill = propgt30)) +
  facet_wrap(~diag_year) +
  scale_fill_viridis_c(option = "magma", na.value = "grey", direction = -1, end = 0.9) +
  labs(fill = "% > 30 days") + theme(axis.text = element_blank()) -> blkyr_probpgt30

blkyr_probpgt30
ggsave(here::here(figdir, "blockyr_propgt30.png"), blkyr_probpgt30, height = 7, width = 10, units = "in")


# Bivariate map incidence and delay
by_block_yr %>%
  filter(!is.na(propgt30)) -> plotdat
plotdat <- bi_class(plotdat, x = propgt30, y = inc, style = "quantile", dim = 3) 

for (y in 2013:2020) {
  
ggplot() +
  geom_sf(data = blockmap, aes(geometry = geometry), fill = "white") +
  geom_sf(data = filter(plotdat, diag_year == y), aes(geometry = geometry, fill = bi_class), show.legend = FALSE) +
  labs(title = y) +
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

ggsave(here::here(figdir, paste0("block_propgt30_byinc_bv_",y,".png")), final, height = 7, width = 14, units = "in")

}

# ---------------------------------------------------------------------------- #
# Semi-variogram

dat <- ll_filt

# convert simple data frame into a spatial data frame object
coordinates(dat) <- ~ longitude + latitude

vg <- variogram(days_fever~1, data = dat)
plot(vg)

vgmod <- vgm(psill = 600, model = "Mat", nugget = 1200, range = 0.2)
plot(vg, model = vgmod)

vgfit <- fit.variogram(vg, model = vgmod)    
plot(vg, model = vgfit)

vgfit

png(here::here(figdir, "OD_semivariogram.png"), height = 500, width = 600)
plot(vg, model = vgfit)
dev.off()

################################################################################
################################################################################
