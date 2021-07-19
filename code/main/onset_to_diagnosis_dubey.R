################################################################################
# Description: Import Case Details Form data from Dubey et al. and compare with 
# KAMIS
################################################################################
################################################################################

figdir <- "figures/CDF"

library(haven)
acddat <- read_sas("~/VL/Data/ACD/acddata_102620/acddata_102620.sas7bdat", NULL) %>%
  as_tibble() %>%
  filter(included == 1) %>%
  mutate(days_fever_cdf = as.numeric(DUR_FEV_R),
         cdf_gt30 = (days_fever_cdf > 30),
         cdf_gt90 = (days_fever_cdf > 90))

summary(acddat)

# Read linelist with matched GPS locations
ll_wgps <- readRDS("~/VL/Diagnosis delay/data/linelist_vl_wgps.rds") %>% #here::here("data","linelist_vl_wgps.rds")
  mutate(OD_missing = is.na(days_fever),
         kamis_gt30 = (days_fever > 30),
         kamis_gt90 = (days_fever > 90))

# Setup map context
blockmap <- readRDS(here::here("data","bihar_block.rds"))

centre <- c(lon = mean(ll_wgps$longitude), lat = mean(ll_wgps$latitude))
extent <- setNames(st_bbox(blockmap), c("left","bottom","right","top"))
bh_lines <- get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)

# ---------------------------------------------------------------------------- #
# Check matching between full linelist and study data

length(which(!acddat$PID %in% ll_wgps$patient_id))
length(which(!ll_wgps$patient_id %in% acddat$PID))

# ll_wgps %>%
#   anti_join(acddat, by = c("patient_id" = "PID")) -> nonmatch1
# 
# acddat %>%
#   anti_join(ll_wgps, by = c("PID" = "patient_id")) -> nonmatch2
# 
# acddat %>%
#   anti_join(linelist_diag_vl, by = c("PID" = "res_patient_code")) -> nonmatch3
# 
# acddat %>%
#   anti_join(linelist_patients, by = c("PID" = "res_patient_code")) -> nonmatch4

ll_wgps %>%
  right_join(acddat, by = c("patient_id" = "PID")) %>%
  mutate(Delay = factor(
    case_when(DUR_FEV_R <= 30 ~ "<= 30 days",
              DUR_FEV_R > 30 & DUR_FEV_R <= 90 ~ "30 - 90 days",
              DUR_FEV_R > 90 ~ "> 90 days"),
    levels = c("<= 30 days","30 - 90 days","> 90 days"))) -> match

summary(is.na(match$longitude))
summary(is.na(match$population))

# Save merged KAMIS/GPS/CDF data for analysis
saveRDS(match, here::here("data","kamis_cdf.rds"))

# ---------------------------------------------------------------------------- #
# Plot locations in Dubey et al. data

ggmap(bh_lines, base_layer = ggplot(ll_wgps, aes(x = longitude, y = latitude))) +
  geom_jitter(alpha = 0.1, cex = 0.8) -> ll_loc

ggmap(bh_lines, base_layer = ggplot(match, aes(x = longitude, y = latitude))) +
  geom_jitter(alpha = 0.1, cex = 0.8) -> acd_loc

compare <- ll_loc + acd_loc
ggsave(here::here("figures", "compare_locations.png"), compare, height = 7, width = 16, units = "in")

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


ggmap(bh_lines, base_layer = ggplot(arrange(match, DUR_FEV_R), aes(x = longitude, y = latitude, col = DUR_FEV_R))) +
  geom_jitter() +
  scale_colour_viridis_c(trans = "log10", direction = -1, option = "plasma") +
  labs(x = "", y = "", col = "Delay", subtitle = "Diagnoses between Jan 2018 - July 2019") -> map_dubey_delay
map_dubey_delay


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
# By block

match %>%
  group_by(district, block) %>%
  summarise(N = n(),
            medianOD = median(DUR_FEV_R, na.rm = T),
            propgt90 = mean((DUR_FEV_R > 90), na.rm = T)*100,
            propgt30 = mean((DUR_FEV_R > 30), na.rm = T)*100) %>%
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


