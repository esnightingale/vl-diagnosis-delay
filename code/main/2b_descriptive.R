################################################################################
# Description: Import Case Details Form data from Dubey et al. and compare with 
# KAMIS
################################################################################
################################################################################

figdir <- "figures/descriptive"
datadir <- "~/VL/Data/Analysis"

# Load analysis data
dat <- readRDS(here::here("data/analysis","dat_nona.rds")) %>%
  dplyr::filter(delay >= 0) %>%
  dplyr::mutate(delay_cat5 = cut(delay, c(0,15,30,90,180,730), include.lowest = TRUE, ordered_result = TRUE),
                delay_cat3 = cut(delay, c(0,30,90,730), include.lowest = TRUE, ordered_result = TRUE),
                delay_gt30 = (delay > 30),
                delay_gt90 = (delay > 90)) %>%
  arrange(delay)

# Setup map context
blockmap <- readRDS(here::here("data","geography","bihar_block.rds"))

# centre <- c(lon = mean(match$longitude), lat = mean(match$latitude))
extent <- setNames(st_bbox(blockmap), c("left","bottom","right","top"))
bh_lines <- get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)

################################################################################
# Overall summary

summary(dat$delay)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   11.00   16.00   31.26   46.00  496.00 
summary(dat$delay_cat5)
# [0,15]   (15,30]   (30,90]  (90,180] (180,730] 
# 1475      1406      1208       221        43 
summary(dat$delay_cat3)
# [0,30]  (30,90] (90,730] 
# 2881     1208      264 
summary(dat$delay_gt90)
#    Mode   FALSE    TRUE 
# logical    4089     264

################################################################################
# Observed delays over time

# ---------------------------------------------------------------------------- #
# By year

dat %>%
  ggplot(aes(x = delay, after_stat(density))) +
  geom_histogram(bins = 50) +
  xlim(c(0,365)) +
  labs(x = "Onset to diagnosis (days)", y = "Density") -> delay_hist
delay_hist

ggsave(here::here(figdir,"delay_hist.png"),
       delay_hist,
       height = 4, width = 8, units = "in")

# ---------------------------------------------------------------------------- #
# By month

dat %>%
  st_drop_geometry() %>%
  group_by(diag_month) %>%
  dplyr::summarise(N = n(),
            mean = mean(delay, na.rm = T),
            median = median(delay, na.rm = T),
            q1 = Rmisc::CI(delay)[3],
            q3 = Rmisc::CI(delay)[1],
            `> 90 days` = sum(delay > 90, na.rm = T)/N,
            `31-90 days` = sum(between(delay,31,90), na.rm = T)/N,
            `<= 30 days` = sum(delay <= 30, na.rm = T)/N) %>%
  ungroup() -> by_mth

by_mth %>%
  ggplot(aes(x = diag_month, group = 1, y = mean, ymin = q1, ymax = q3)) +
  geom_ribbon(alpha = 0.1) +
  geom_line(lty = "solid") +
  geom_line(aes(y = median), col = "steelblue", lty = "solid") +
  labs(x = "Month of diagnosis", y = "Delay beyond 14 days fever") -> delay_bymth
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
  labs(x = "Month of diagnosis", y = "Proportion", fill = "Excess delay") -> prop_excessive
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
                                  delay), 
                          aes(col = delay))) +
  geom_sf(alpha = 0.5, cex = 0.8) +
  scale_colour_viridis_c(trans = "log10", direction = -1, option = "plasma") +
  labs(x = "", y = "", col = "Delay", subtitle = "Diagnoses between Jan 2018 - July 2019") -> map_village
map_village

ggsave(here::here(figdir,"delay_byvil.png"),
       map_village,
       height = 7, width = 10, units = "in")


## 5 categories
ggmap(bh_lines, 
      base_layer = ggplot(dat, 
                          aes(x = longitude, y = latitude, col = delay_cat5))) +
  geom_jitter(alpha = 0.8) +
  scale_colour_viridis_d(option = "viridis", na.value = "grey", direction = -1, end = 0.9) +
  labs(x = "", y = "", col = "Delay") -> map_cat5
map_cat5

ggsave(here::here(figdir, "village_delay_cat5.png"), 
       map_cat5, 
       height = 7, width = 10, units = "in")

## 3 categories
ggmap(bh_lines, 
      base_layer = ggplot(dat, 
                          aes(x = longitude, y = latitude, col = delay_cat3))) +
  geom_jitter(alpha = 0.8) +
  scale_colour_viridis_d(option = "viridis", na.value = "grey", direction = -1, end = 0.9) +
  labs(x = "", y = "", col = "Delay") -> map_cat3
map_cat3

ggsave(here::here(figdir, "village_delay_cat3.png"), 
       map_cat3, 
       height = 7, width = 10, units = "in")


## Binary > 90 days
ggmap(bh_lines, 
      base_layer = ggplot(arrange(dat, delay_gt90), 
                          aes(x = longitude, y = latitude, col = delay_gt90))) +
  geom_jitter(alpha = 0.8) +
  scale_colour_manual(values = c("grey","indianred")) +
  labs(x = "", y = "", col = "Delay > 90 days", subtitle = "Diagnoses between Jan 2018 - July 2019") -> map_excess
map_excess

ggsave(here::here(figdir, "village_delay_excess.png"), 
       map_excess, 
       height = 7, width = 10, units = "in")

# ---------------------------------------------------------------------------- #
# By village

dat %>%
  group_by(vil_code, latitude, longitude) %>%
  dplyr::summarise(N = n(),
            pop = unique(population),
            inc = replace_na(N*1e3/unique(pop),0),
            med_delay = median(delay, na.rm = T),
            propgt30 = mean((delay > 30), na.rm = T)*100,
            propgt90 = mean((delay > 90), na.rm = T)*100) %>%
  arrange(propgt90) -> by_village

## Median
ggmap(bh_lines, 
      base_layer = ggplot(by_village, 
                          aes(x = longitude, y = latitude, col = med_delay))) +
  geom_point(pch = 19, alpha = 0.5) +
  scale_colour_viridis_c(option = "viridis", na.value = "grey", direction = -1, trans = "log2", end = 0.9) +
  # scale_alpha_continuous(trans = "log10") +
  labs(x = "", y = "", col = "Median delay") -> map_med_delay
map_med_delay

ggsave(here::here(figdir, "village_med_delay.png"), 
       map_med_delay, 
       height = 7, width = 10, units = "in")


## % greater than 30 days
ggmap(bh_lines, 
      base_layer = ggplot(by_village, 
                          aes(x = longitude, y = latitude, col = propgt30))) +
  geom_point(alpha = 0.7) +
  scale_colour_viridis_c(option = "magma", na.value = "grey", direction = -1, end = 0.9) +
  # scale_alpha_continuous(trans = "log10") +
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
  # scale_alpha_continuous(trans = "log10") +
  labs(x = "", y = "", col = "% > 90 days") -> map_propgt90
map_propgt90

ggsave(here::here(figdir, "village_propgt90.png"), 
       map_propgt90, 
       height = 7, width = 10, units = "in")


# ---------------------------------------------------------------------------- #
# By block

dat %>%
  st_drop_geometry() %>%
  group_by(district, block) %>%
  dplyr::summarise(N = n(),
            med_delay = median(delay, na.rm = T),
            pgt90 = mean(delay_gt90, na.rm = T)*100,
            pgt30 = mean(delay_gt30, na.rm = T)*100) %>%
  mutate(delay_cat3 = cut(med_delay, c(0,15,90,730), include.lowest = TRUE, ordered_result = TRUE)) %>%
  ungroup() -> by_block

by_block <- blockmap %>%
  left_join(by_block, by = c("kamis_master_dist" = "district", "kamis_master_block" = "block")) %>%
  rowwise() %>%
  dplyr::mutate(pop = median(c_across(`2018`:`2019`)),
         inc = N*1e4/pop) %>%
  ungroup() 

# pal3 <- viridis::viridis(3)
ggplot() +
  geom_sf(data = blockmap) +
  geom_sf(data = by_block, aes(geometry = geometry, fill = med_delay)) +
  scale_fill_viridis_c(option = "viridis", na.value = "white", direction = -1, trans = "sqrt", end = 0.9) +
  labs(fill = "Median delay") + 
  theme(axis.text = element_blank()) -> blk_med_delay
blk_med_delay

ggsave(here::here(figdir, "block_med_delay.png"), blk_med_delay, height = 7, width = 10, units = "in")

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
                                   range(dat$diag_month)[1],
                                   " and ",
                                   range(dat$diag_month)[2]))
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

dat.sp <- dat %>%
   st_drop_geometry() %>%
   filter(!is.na(longitude) & !is.na(delay))

# convert simple data frame into a spatial data frame object
coordinates(dat.sp) <- ~ longitude + latitude

vg <- variogram(delay~1, data = dat)
plot(vg)

vgmod <- vgm(psill =  800, model = "Mat", nugget = 1000, range = 100)
plot(vg, model = vgmod)

vgfit <- fit.variogram(vg, model = vgmod)    
plot(vg, model = vgfit)

vgfit
# model    psill    range kappa
# Nug 927.8130  0.00000   0.0
# Mat 935.2576 60.21329   0.5

png(here::here(figdir, "delay_semivariogram.png"), height = 500, width = 600)
plot(vg, model = vgfit)
dev.off()

# ---------------------------------------------------------------------------- #
# Tabulate

make_tab <- function(dat, varname){
  
  dat %>%
    dplyr::mutate(Value = !!sym(varname)) %>%
    dplyr::group_by(Value) %>%
    dplyr::summarise(N = n(),
                     mean = round(mean(delay),1),
                     ci = paste(round(mean + sd(delay)/sqrt(N) * qnorm(p = c(0.025, 0.975)),1), collapse = ","),
                     med = median(delay),
                     iqr = paste(quantile(delay, p = c(0.25, 0.75)), collapse = ","),
                     n_gt30 = sum(delay_gt30),
                     p_gt30 = round(n_gt30/N,2)) %>%
    dplyr::mutate(`Delay, mean [95% CI]` = paste0(mean, " [",ci,"]"),
                  `Delay, median [IQR]` = paste0(med, " [",iqr,"]"),
                  `> 30 days, N (%)` = paste0(n_gt30, " (",p_gt30,")"),
                  Variable = varname) %>%
    dplyr::select(Variable, Value, N, `Delay, mean [95% CI]`, `> 30 days, N (%)`) %>%
    dplyr::arrange(Value) -> tab
  
  return(tab)
  
}

make_plot <- function(t) {
  
  t %>%
    filter(Value != "NA") %>%
    separate(`Delay, mean [95% CI]`, into = c("mean","ll","ul"), sep = "[\\[\\,\\]]", convert = TRUE) %>%
    ggplot(aes(Value, mean, ymin = ll, ymax = ul)) +
    geom_linerange() +
    geom_point() +
    geom_hline(yintercept = mean(dat.df$delay), col = "grey") +
    labs(subtitle = unique(t$Variable), x = "", y =  "") %>%
    return()
  
}

dat.df %>%
  dplyr::rename(Sex = sex,
                Age = age_cat,
                `Scheduled caste or tribe` = marg_caste,
                Occupation = occupation,
                `HIV status` = hiv,
                `Previous VL/PKDL treatment` = prv_tx,
                `Number of prior consultations` = conslt_cat,
                Detection = detection,
                `Block endemic in 2017` = block_endm_2017,
                `Village IRS targeted in 2017` = IRS_2017_1,
                `Village incidence > 0 in 2017` = inc_2017_gt0,
                `Travel time to nearest diagnosis facility` = travel_time_cat,
                `Travel time to nearest treatment facility` = travel_time_t_cat) -> dat.tab

# varlist <- list("sex","age_child", "marg_caste","hiv", "prv_tx_ka", "num_conslt_4","detection","block_endm_2017" ,"IRS_2017_1", "inc_2017_gt0", "travel_time_cat")
varlist <- list("Sex","Age","Scheduled caste or tribe","Occupation","HIV status", 
                "Previous VL/PKDL treatment","Number of prior consultations",
                "Detection","Block endemic in 2017",
                "Village IRS targeted in 2017","Village incidence > 0 in 2017", 
                "Travel time to nearest diagnosis facility",
                "Travel time to nearest treatment facility")

tabs.list <- lapply(varlist, make_tab, dat = dat.tab)

plots <- lapply(tabs.list, make_plot)

png(here::here(figdir, "covariate_mean_ci_plot.png"), height = 10, width = 15, units = "in", res = 300)
do.call("grid.arrange", plots)
dev.off()

tab <- bind_rows(lapply(tabs.list, function(t) mutate(t, Value = as.character(Value))))
View(tab)

write.csv(tab, here::here("output","table_descriptive.csv"), row.names = FALSE)

################################################################################
################################################################################