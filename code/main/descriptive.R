################################################################################
# Description: Import Case Details Form data from Dubey et al. and compare with 
# KAMIS
################################################################################
################################################################################

figdir <- "figures/descriptive"

# Load analysis data
dat <- readRDS(here::here("data","analysisdata_individual.rds")) %>%
  dplyr::mutate(delay_cat5 = cut(days_fever, c(0,15,30,90,180,730), include.lowest = TRUE, ordered_result = TRUE),
                delay_cat3 = cut(days_fever, c(0,30,90,730), include.lowest = TRUE, ordered_result = TRUE),
                detection = factor(poss_acd, labels = c("PCD","ACD"))) %>%
  arrange(days_fever)

# Setup map context
blockmap <- readRDS(here::here("data","geography","bihar_block.rds"))

# centre <- c(lon = mean(match$longitude), lat = mean(match$latitude))
extent <- setNames(st_bbox(blockmap), c("left","bottom","right","top"))
bh_lines <- get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)

################################################################################
# Overall summary

summary(dat$days_fever)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.00   25.00   30.00   45.01   58.00  510.00 
summary(dat$delay_cat5)
# [0,15]   (15,30]   (30,90]  (90,180] (180,730] 
# 198      2284      1607       251        54 
summary(dat$delay_cat3)
# [0,30]  (30,90] (90,730] 
# 2482     1607      305 
summary(dat$gt90)
# FALSE    TRUE 
# 4089     305

################################################################################
# Observed delays over time

# ---------------------------------------------------------------------------- #
# By year

dat %>%
  ggplot(aes(x = days_fever)) +
  geom_histogram(bins = 50) +
  xlim(c(0,365)) +
  facet_wrap(~diag_year) +
  labs(x = "Onset to diagnosis (days)", y = "Frequency") -> delay_hist
delay_hist

ggsave(here::here(figdir,"delay_hist.png"),
       delay_hist,
       height = 4, width = 8, units = "in")

# ---------------------------------------------------------------------------- #
# By month

dat %>%
  st_drop_geometry() %>%
  group_by(diag_month) %>%
  summarise(N = n(),
            mean = mean(days_fever, na.rm = T),
            median = median(days_fever, na.rm = T),
            q1 = quantile(days_fever, p = 0.05, na.rm = T),
            q3 = quantile(days_fever, p = 0.95, na.rm = T),
            `> 90 days` = sum(days_fever > 90, na.rm = T)/N,
            `31-90 days` = sum(between(days_fever,31,90), na.rm = T)/N,
            `<= 30 days` = sum(days_fever <= 30, na.rm = T)/N) %>%
  ungroup() -> by_mth

by_mth %>%
  ggplot(aes(x = diag_month, group = 1, y = median, ymin = q1, ymax = q3)) +
  geom_errorbar(width = 0.5) +
  geom_point() +
  geom_line(aes(y = mean), col = "steelblue", lty = "dashed", lwd = 1.2) +
  scale_y_continuous(trans = "log2") +
  labs(x = "Month of diagnosis", y = "Onset to diagnosis (days)") -> delay_bymth
delay_bymth

ggsave(here::here(figdir,"delay_bymth.png"),
       delay_bymth,
       height = 4, width = 6, units = "in")


pal <- viridis::viridis(3, end = 0.8, direction = -1)
by_mth %>%
  pivot_longer(c(`> 90 days`,`31-90 days`,`<= 30 days`), names_to = "Delay") %>%
  mutate(Delay = factor(Delay, levels = c("<= 30 days", "31-90 days", "> 90 days"))) %>%
  ggplot(aes(x = diag_month, y = value, fill = Delay)) +
  geom_col() +
  scale_fill_manual(values = pal) +
  labs(x = "Month of diagnosis", y = "Proportion") -> prop_excessive
prop_excessive

ggsave(here::here(figdir,"prop_excessive_bymth.png"),
       prop_excessive,
       height = 4, width = 6, units = "in")

################################################################################
# Observed delays by location

# ---------------------------------------------------------------------------- #
# Individuals

## Continuous
ggmap(bh_lines, 
      base_layer = ggplot(arrange(st_jitter(dat), 
                                  days_fever), 
                          aes(col = days_fever))) +
  geom_sf(alpha = 0.5, cex = 0.8) +
  scale_colour_viridis_c(trans = "log10", direction = -1, option = "plasma") +
  labs(x = "", y = "", col = "Delay", subtitle = "Diagnoses between Jan 2018 - July 2019") -> map_village
map_village

ggsave(here::here(figdir,"delay_byvil.png"),
       map_village,
       height = 7, width = 10, units = "in")


## 5 categories
ggmap(bh_lines, 
      base_layer = ggplot(match, 
                          aes(x = longitude, y = latitude, col = ODcat5))) +
  geom_jitter(alpha = 0.8) +
  scale_colour_viridis_d(option = "viridis", na.value = "grey", direction = -1, end = 0.9) +
  labs(x = "", y = "", col = "Delay") -> map_cat5
map_cat5

ggsave(here::here(figdir, "village_OD_cat5.png"), 
       map_cat5, 
       height = 7, width = 10, units = "in")

## 3 categories
ggmap(bh_lines, 
      base_layer = ggplot(match, 
                          aes(x = longitude, y = latitude, col = ODcat3))) +
  geom_jitter(alpha = 0.8) +
  scale_colour_viridis_d(option = "viridis", na.value = "grey", direction = -1, end = 0.9) +
  labs(x = "", y = "", col = "Delay") -> map_cat3
map_cat3

ggsave(here::here(figdir, "village_OD_cat3.png"), 
       map_cat3, 
       height = 7, width = 10, units = "in")


## Binary > 90 days
ggmap(bh_lines, 
      base_layer = ggplot(arrange(match,gt90_cdf), 
                          aes(x = longitude, y = latitude, col = gt90_cdf))) +
  geom_jitter(alpha = 0.8) +
  scale_colour_manual(values = c("grey","indianred")) +
  labs(x = "", y = "", col = "Delay > 90 days", subtitle = "Diagnoses between Jan 2018 - July 2019") -> map_excess
map_excess

ggsave(here::here(figdir, "village_OD_excess.png"), 
       map_excess, 
       height = 7, width = 10, units = "in")

# ---------------------------------------------------------------------------- #
# By village

match %>%
  group_by(vil_code, latitude, longitude) %>%
  summarise(N = n(),
            pop = unique(population),
            inc = replace_na(N*1e3/unique(pop),0),
            medianOD = median(days_fever_cdf, na.rm = T),
            propgt30 = mean((days_fever_cdf > 30), na.rm = T)*100,
            propgt90 = mean((days_fever_cdf > 90), na.rm = T)*100) %>%
  arrange(propgt90) -> by_village

## Median
ggmap(bh_lines, 
      base_layer = ggplot(by_village, 
                          aes(x = longitude, y = latitude, col = medianOD))) +
  geom_point(pch = 19, alpha = 0.5) +
  # scale_colour_gradient2(low = pal3[3], mid = pal3[2], high = pal3[1], midpoint = log(30,2), trans = "log2") +
  scale_colour_viridis_c(option = "viridis", na.value = "grey", direction = -1, trans = "log2", end = 0.9) +
  scale_alpha_continuous(trans = "log10") +
  labs(x = "", y = "", col = "Median delay") -> map_medOD
map_medOD

ggsave(here::here(figdir, "village_medianOD.png"), 
       map_medOD, 
       height = 7, width = 10, units = "in")


## % greater than 30 days
ggmap(bh_lines, 
      base_layer = ggplot(by_village, 
                          aes(x = longitude, y = latitude, col = propgt30))) +
  geom_point(alpha = 0.7) +
  scale_colour_viridis_c(option = "magma", na.value = "grey", direction = -1, end = 0.9) +
  scale_alpha_continuous(trans = "log10") +
  labs(x = "", y = "", col = "% > 30 days") -> map_propgt30
map_propgt30

ggsave(here::here(figdir, "village_propgt30.png"), 
       map_propgt30, 
       height = 7, width = 10, units = "in")


## % greater than 90 days
ggmap(bh_lines, 
      base_layer = ggplot(by_village, 
                          aes(x = longitude, y = latitude, col = propgt90))) +
  geom_point(alpha = 0.7) +
  scale_colour_viridis_c(option = "magma", na.value = "grey", direction = -1, end = 0.9) +
  scale_alpha_continuous(trans = "log10") +
  labs(x = "", y = "", col = "% > 90 days") -> map_propgt90
map_propgt90

ggsave(here::here(figdir, "village_propgt90.png"), 
       map_propgt90, 
       height = 7, width = 10, units = "in")


# ---------------------------------------------------------------------------- #
# By block

match %>%
  group_by(district, block) %>%
  summarise(N = n(),
            medianOD = median(days_fever_cdf, na.rm = T),
            pgt90 = mean(gt90_cdf, na.rm = T)*100,
            pgt30 = mean(gt30_cdf, na.rm = T)*100) %>%
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
  geom_sf(data = by_block, aes(geometry = geometry, fill = pgt30)) +
  scale_fill_viridis_c(option = "magma", na.value = "white", direction = -1, end = 0.8) +
  labs(fill = "% > 30 days") + 
  theme(axis.text = element_blank()) -> blk_pgt30
blk_pgt30

ggsave(here::here(figdir, "block_propgt30.png"), blk_pgt30, height = 7, width = 10, units = "in")


ggplot() +
  geom_sf(data = blockmap) +
  geom_sf(data = by_block, aes(geometry = geometry, fill = pgt90)) +
  scale_fill_viridis_c(option = "magma", na.value = "white", direction = -1, end = 0.9) +
  labs(fill = "% > 90 days") +
  theme(axis.text = element_blank()) -> blk_pgt90
blk_pgt90

ggsave(here::here(figdir, "block_propgt90.png"), blk_pgt90, height = 7, width = 10, units = "in")

# Join two block maps
combined <- blk_pgt30 + blk_pgt90 +
  plot_annotation(title = "Proportion of cases diagnosed with excess delays",
                  caption = paste0("Diagnoses between ",
                                   range(match$diag_date_cdf)[1],
                                   " and ",
                                   range(match$diag_date_cdf)[2]))
combined

ggsave(here::here(figdir, "excess_delay_byblock.png"), combined, height = 7, width = 14, units = "in")

# ---------------------------------------------------------------------------- #
# Bivariate map of delay and incidence

by_block %>%
  filter(!is.na(pgt30)) %>%
  mutate(inc = replace_na(inc, 0)) -> plotdat

bv <- create_bivmap(blockmap, plotdat,
                   xvar = "pgt30", yvar = "inc",
                   xlab = "% > 30 days", ylab = "Cases/10,000",
                   pal = "DkBlue")
bv

ggsave(here::here(figdir, "block_propgt30_byinc_bv.png"), bv, height = 7, width = 10, units = "in")

# ---------------------------------------------------------------------------- #
# Semi-variogram

dat <- filter(match, !is.na(longitude) & !is.na(days_fever_cdf))

# convert simple data frame into a spatial data frame object
coordinates(dat) <- ~ longitude + latitude
crs(dat) <- "+proj=longlat +datum=WGS84 +units=km +no_defs"

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

