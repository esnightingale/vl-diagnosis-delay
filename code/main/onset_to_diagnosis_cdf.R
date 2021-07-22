################################################################################
# Description: Import Case Details Form data from Dubey et al. and compare with 
# KAMIS
################################################################################
################################################################################

figdir <- "figures/CDF"

# Load merged KAMIS/GPS/CDF data 
match <- readRDS(here::here("data","kamis_cdf.rds"))

# Setup map context
blockmap <- readRDS(here::here("data","bihar_block.rds"))

centre <- c(lon = mean(ll_wgps$longitude), lat = mean(ll_wgps$latitude))
extent <- setNames(st_bbox(blockmap), c("left","bottom","right","top"))
bh_lines <- get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)

################################################################################
# Observed delays over time

# ---------------------------------------------------------------------------- #
# By year

match %>%
  # filter(diag_year_cdf == 2019) %>%
  ggplot(aes(x = days_fever_cdf)) +
  geom_histogram(bins = 50) +
  xlim(c(0,365)) +
  facet_wrap(~diag_year_cdf) +
  labs(x = "Onset to diagnosis (days)", y = "Frequency") -> OD_hist
OD_hist

ggsave(here::here(figdir,"OD_hist.png"),
       OD_hist,
       height = 4, width = 8, units = "in")

# ---------------------------------------------------------------------------- #
# By month

match %>%
  group_by(diag_month_cdf) %>%
  summarise(mean = mean(days_fever_cdf, na.rm = T),
            median = median(days_fever_cdf, na.rm = T),
            q1 = quantile(days_fever_cdf, p = 0.05, na.rm = T),
            q3 = quantile(days_fever_cdf, p = 0.95, na.rm = T)) %>%
  ungroup() -> by_mth

by_mth %>%
  ggplot(aes(x = diag_month_cdf, group = 1, y = median, ymin = q1, ymax = q3)) +
  geom_errorbar(width = 0.5) +
  geom_point() +
  geom_line(aes(y = mean), col = "steelblue", lty = "dashed", lwd = 1.2) +
  scale_y_continuous(trans = "log2") +
  labs(x = "Month of diagnosis", y = "Onset to diagnosis (days)") -> OD_bymth
OD_bymth

ggsave(here::here(figdir,"OD_bymth.png"),
       OD_bymth,
       height = 4, width = 6, units = "in")


pal <- c("grey",viridis::viridis(3, end = 0.8, direction = -1))
match %>%
  group_by(diag_month_cdf) %>%
  summarise(N = n(),
            `> 90 days` = sum(days_fever_cdf > 90, na.rm = T)/N,
            `31-90 days` = sum(between(days_fever_cdf,31,90), na.rm = T)/N,
            `<= 30 days` = sum(days_fever_cdf <= 30, na.rm = T)/N,
            Missing = sum(is.na(days_fever_cdf))/N) %>%
  pivot_longer(c(`> 90 days`,`31-90 days`,`<= 30 days`,Missing), names_to = "Delay") %>%
  mutate(Delay = factor(Delay, levels = c("Missing","<= 30 days", "31-90 days", "> 90 days"))) %>%
  ggplot(aes(x = diag_month_cdf, y = value, fill = Delay)) +
  geom_col() +
  scale_fill_manual(values = pal) +
  labs(x = "Month of diagnosis", y = "Proportion") -> prop_excessive
prop_excessive

ggsave(here::here(figdir,"prop_excessive_bymth.png"),
       prop_excessive,
       height = 4, width = 6, units = "in")


################################################################################
# Observed delays by location

match <- match %>%
  filter(!is.na(days_fever_cdf)) %>%
  mutate(ODcat5 = cut(days_fever_cdf, c(0,15,30,90,180,730), include.lowest = TRUE, ordered_result = TRUE),
         ODcat3 = cut(days_fever_cdf, c(0,30,90,730), include.lowest = TRUE, ordered_result = TRUE),
         ODexcess = (days_fever_cdf > 90)) %>%
  arrange(days_fever_cdf)

summary(match$days_fever_cdf)
summary(match$ODcat5)
summary(match$ODcat3)
summary(match$ODexcess)

# ---------------------------------------------------------------------------- #
# Individuals

## Continuous
ggmap(bh_lines, 
      base_layer = ggplot(arrange(match, days_fever_cdf), 
                          aes(x = longitude, y = latitude, col = days_fever_cdf))) +
  geom_jitter(alpha = 0.8) +
  scale_colour_viridis_c(trans = "log10", direction = -1, option = "plasma") +
  labs(x = "", y = "", col = "Delay", subtitle = "Diagnoses between Jan 2018 - July 2019") -> map_village
map_village

ggsave(here::here(figdir,"OD_byvil.png"),
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

ggsave(here::here(figdir, "village_OD_cat3.png"), 
       map_cat3, 
       height = 7, width = 10, units = "in")


## Binary > 90 days
ggmap(bh_lines, 
      base_layer = ggplot(arrange(match,gt90_cdf), 
                          aes(x = longitude, y = latitude, col = gt90_cdf))) +
  geom_jitter(alpha = 0.8) +
  scale_colour_manual(values = c("grey","indianred")) +
  labs(x = "", y = "", col = "Delay > 90 days", subtitle = "Diagnoses between 2016 and 2018") -> map_excess

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
            propgt90 = mean(gt90_cdf, na.rm = T)*100,
            propgt30 = mean(gt30_cdf, na.rm = T)*100) %>%
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
plotdat <- bi_class(plotdat, x = propgt30, y = inc, style = "fisher", dim = 3) 

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

dat <- filter(match, !is.na(longitude) & !is.na(DUR_FEV_R))

# convert simple data frame into a spatial data frame object
coordinates(dat) <- ~ longitude + latitude

vg <- variogram(DUR_FEV_R~1, data = dat)
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



# ---------------------------------------------------------------------------- #
# Compare OD between KAMIS and CDF
# Continuous

match %>% 
  mutate(diff = days_fever - days_fever_cdf) %>%
  dplyr::select(days_fever, days_fever_cdf, diff) %>% 
  summary()

# days_fever     days_fever_cdf       diff         
# Min.   :  1.00   Min.   :  4.0   Min.   :-662.000  
# 1st Qu.: 20.00   1st Qu.: 25.0   1st Qu.: -15.000  
# Median : 30.00   Median : 30.0   Median :   0.000  
# Mean   : 37.74   Mean   : 45.5   Mean   :  -7.738  
# 3rd Qu.: 35.00   3rd Qu.: 60.0   3rd Qu.:   0.000  
# Max.   :546.00   Max.   :690.0   Max.   : 526.000  
# NA's   :20                       NA's   :20  

match %>%
  pivot_longer(c("days_fever","days_fever_cdf")) %>%
  ggplot(aes(value)) +
  geom_histogram() +
  geom_vline(aes(xintercept = mean(value, na.rm = T)), col = "red", lty = "dashed") +
  geom_vline(aes(xintercept = median(value, na.rm = T)), col = "dodgerblue", lty = "dashed") +
  facet_wrap(~name)

match %>%
  ggplot(aes(days_fever, days_fever_cdf)) +
  geom_jitter(alpha = 0.2) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "KAMIS value", y = "CDF value") -> compare_num
ggsave(here::here("figures", "days_fever_kamis_cdf.png"), compare_num, height = 6, width = 7, units = "in")

# ---------------------------------------------------------------------------- #
# Categorical

gt30tab <- prop.table(table(match$kamis_gt30, match$cdf_gt30))
gt30tab
gt90tab <- prop.table(table(match$kamis_gt90, match$cdf_gt90))
gt90tab

match %>%
  filter(!is.na(days_fever)) %>%
  group_by(kamis_gt30, cdf_gt30) %>%
  tally() %>%
  ungroup() %>%
  mutate(value = paste0(n, " (",round(n/nrow(match),2),")")) %>% 
  pivot_wider(-n,values_from = value, names_from = cdf_gt30)
#           `<= 30`    `> 30`     
# <= 30      2493 (0.5) 1169 (0.23)
# > 30       297 (0.06) 1049 (0.21)

match %>%
  filter(!is.na(days_fever)) %>%
  group_by(kamis_gt90, cdf_gt90) %>%
  tally() %>%
  ungroup()  %>%
  mutate(value = paste0(n, " (", round(n/nrow(match),2),")")) %>% 
  pivot_wider(-n,values_from = value, names_from = cdf_gt90)
#            `<= 90`    `> 90`    
# <= 90      4539 (0.9) 248 (0.05)
# > 90       117 (0.02) 104 (0.02)


pal <- viridis::viridis(3, end = 0.9, direction = -1)
ggmap(bh_lines, base_layer = ggplot(arrange(match,Delay), aes(x = longitude, y = latitude, col = Delay))) +
  geom_jitter(alpha = 0.5) +
  scale_colour_manual(values = pal) +
  labs(x = "", y = "", col = "Delay > 90 days", subtitle = "Diagnoses between Jan 2018-July 2019") -> map_delay_cat
map_delay_cat

# ---------------------------------------------------------------------------- #
# Compare distribution of excessive (>90) delays

ggmap(bh_lines, base_layer = ggplot(arrange(match,cdf_gt90), aes(x = longitude, y = latitude, col = cdf_gt90))) +
  geom_jitter() +
  scale_colour_manual(values = c("grey", "indianred")) +
  guides(col = FALSE) +
  labs(x = "", y = "", col = "Delay > 90 days", title = "CDF") -> map_cdf_gt90
map_cdf_gt90

# Filter KAMIS to same time period
ll_filt <- ll_wgps %>%
  filter(between(diag_date, min(match$Date_Diag), max(match$Date_Diag)) & !is.na(days_fever)) %>%
  mutate(ODexcess = (days_fever > 90))

ggmap(bh_lines, base_layer = ggplot(arrange(ll_filt,ODexcess), aes(x = longitude, y = latitude, col = ODexcess))) +
  geom_jitter() +
  scale_colour_manual(values = c("grey","indianred")) +
  labs(x = "", y = "", col = "Delay > 90 days", title = "KAMIS", caption = "Diagnoses between Jan 2018-July 2019") -> map_kamis_gt90

compare_excess <- map_cdf_gt90 + map_kamis_gt90

ggsave(here::here(figdir, "compare_data_delaygt90.png"), compare_excess, height = 7, width = 16, units = "in")



ggmap(bh_lines, base_layer = ggplot(arrange(match,kamis_gt90), aes(x = longitude, y = latitude, col = kamis_gt90))) +
  geom_jitter() +
  scale_colour_manual(values = c("grey","indianred")) +
  labs(x = "", y = "", col = "Delay > 90 days", title = "KAMIS", caption = "Diagnoses between Jan 2018-July 2019") -> map_kamis_gt90

compare_excess <- map_cdf_gt90 + map_kamis_gt90
compare_excess
ggsave(here::here(figdir, "compare_data_delaygt90.png"), compare_excess, height = 7, width = 16, units = "in")


# ---------------------------------------------------------------------------- #
# Compare KAMIS and CDF

match %>%
  group_by(district, block) %>%
  summarise(N = n(),
            KAMIS = mean((days_fever > 30), na.rm = T)*100,
            Dubey = mean((DUR_FEV_R > 30), na.rm = T)*100) %>%
  ungroup() %>%
  pivot_longer(c("KAMIS","Dubey")) -> by_block_both

by_block_both <- blockmap %>%
  left_join(by_block_both, by = c("kamis_master_dist" = "district", "kamis_master_block" = "block")) %>%
  rowwise() %>%
  mutate(pop = median(c_across(`2018`:`2019`)),
         inc = N*1e4/pop) %>%
  ungroup() 

by_block_both %>%
  filter(!is.na(value)) %>%
  mutate(inc = replace_na(inc, 0)) -> plotdat2
plotdat2 <- bi_class(plotdat2, x = value, y = inc, style = "quantile", dim = 3) 

ggplot() +
  geom_sf(data = blockmap, fill = "white") +
  geom_sf(data = plotdat2, aes(fill = bi_class), show.legend = FALSE) +
  facet_wrap(~name) +
  bi_scale_fill(pal = "DkBlue", dim = 3, na.value = "white") -> biv_map_both

final <- cowplot::ggdraw() +
  cowplot::draw_plot(biv_map_both, 0, 0, 1, 1) +
  cowplot::draw_plot(legend, 0.78, 0.65, 0.25, 0.25)
final

ggsave(here::here(figdir, "compare_kamis_dubey_bv.png"), final, height = 7, width = 15, units = "in")


