################################################################################
# Description: Import Case Details Form data from Dubey et al. and compare with 
# KAMIS
################################################################################
################################################################################

figdir <- "figures/descriptive"

# Load analysis data
dat <- read_data() %>%
  dplyr::mutate(delay_cat5 = gtools::quantcut(delay, 5),
                delay_cat3 = gtools::quantcut(delay, 3),
                gt90 = (delay > 90),
                gt60 = (delay > 60)) %>%
  arrange(delay)

# Setup map context
blockmap <- readRDS(here::here("data","geography","blockmap.rds")) %>% 
  sf::st_set_crs(7759)
boundary <- readRDS(here::here("data","geography","boundary.rds")) %>% 
  sf::st_set_crs(7759)
boundary.spdf <- as_Spatial(boundary)

################################################################################
# Overall summary

summary(dat$delay)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   11.00   16.00   31.05   44.00  496.00
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   11.00   16.00   31.26   46.00  496.00 
summary(dat$delay_cat5)
# [0,8]   (8,16]  (16,21]  (21,46] (46,496] 
# 909     1498      228      960      676 
summary(dat$delay_cat3)
# [0,15]  (15,31] (31,496] 
# 1449     1590     1232 
summary(dat$gt90)
# Mode   FALSE    TRUE 
# logical    4014     257 

dat %>%
  ggplot(aes(x = dur_fev_r)) +
  geom_histogram(bins = 50) +
  xlim(c(0,250)) +
  geom_vline(xintercept = 14, lty = "dashed") +
  annotate("text", x = 200, y = 300, size = 2, 
           label = paste0("N = 24", #Axis cut at 365 days\n
                          # sum(dat$days_fever < 14),
                          " cases with < 14 days fever before diagnosis\nN = ",
                          sum(dat$dur_fev_r > 250),
                          " cases with > 250 days"
                          )) +
  # facet_wrap(~diag_year) +
  labs(x = "Onset to diagnosis (days)", y = "Frequency") -> days_fever_hist
days_fever_hist

ggsave(here::here(figdir,"days_fever_hist.png"),
       days_fever_hist,
       height = 4, width = 5, units = "in")

plot_vgm(dat$delay, dat, title = "Observed delay")
#   model    psill    range kappa
# 1   Nug 331.7657  0.00000   0.0
# 2   Mat 348.9317 49.96614   0.5

################################################################################
# Observed delays over time

dat %>%
  st_drop_geometry() %>%
  group_by(diag_month) %>%
  dplyr::summarise(N = n(),
            mean = mean(delay, na.rm = T),
            median = median(delay, na.rm = T),
            q1 = Rmisc::CI(delay)[3],
            q3 = Rmisc::CI(delay)[1],
            `> 90 days` = sum(dur_fev_r > 90, na.rm = T)/N,
            `31-90 days` = sum(between(dur_fev_r,31,90), na.rm = T)/N,
            `<= 30 days` = sum(dur_fev_r <= 30, na.rm = T)/N) %>%
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
# Visualise the spatial range of the data and pattern of delays

ggplot() +
  geom_sf(data = blockmap, fill = NA) +
  geom_sf(data = st_jitter(dat), aes(col = delay_cat5), alpha = 0.9, cex = 1.0) +
  scale_colour_viridis_d(direction = -1, end = 0.9) + #trans = "log2"
  labs(col = "Delay (days)") + 
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9,0.85)) 

ggsave(here::here("figures/descriptive/fig1_block.png"), height = 6, width = 8, units = "in")

ggplot() +
  geom_sf(data = boundary, fill = NA) +
  geom_sf(data = st_jitter(dat), aes(col = delay_cat5), alpha = 0.9, cex = 1.0) +
  scale_colour_viridis_d(direction = -1, end = 0.9) + #trans = "log2"
  labs(col = "Delay (days)") + 
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9,0.85)) 
# +
#   annotation_scale(location = "br") +
#   annotation_north_arrow(location = "tr", pad_x = unit(2, "cm"), pad_y = unit(1, "cm"))

ggsave(here::here("figures/descriptive/fig1_state.png"), height = 6, width = 8, units = "in")

# ---------------------------------------------------------------------------- #
# By village

dat %>%
  group_by(v, latitude, longitude) %>%
  dplyr::summarise(N = n(),
            med_delay = median(delay, na.rm = T),
            propgt30 = mean(gt30, na.rm = T)*100,
            propgt60 = mean(gt60, na.rm = T)*100,
            propgt90 = mean(gt90, na.rm = T)*100) %>%
  arrange(propgt90) -> by_village

## Median
ggplot() +
  geom_sf(data = boundary, fill = NA) +
  geom_sf(data = by_village, aes(col = med_delay), alpha = 0.9, cex = 1.0) +
  geom_point(pch = 19, alpha = 0.9) +
  scale_colour_viridis_c(option = "viridis", na.value = "grey", direction = -1, trans = "identity", end = 0.9) +
  # scale_alpha_continuous(trans = "log10") +
  labs(x = "", y = "", col = "Median delay") + 
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9,0.85)) -> map_med_delay
map_med_delay

ggsave(here::here(figdir, "village_med_delay.png"),
       map_med_delay,
       height = 7, width = 10, units = "in")


## % greater than 30 days
ggplot() +
  geom_sf(data = boundary, fill = NA) +
  geom_sf(data = by_village, aes(col = propgt30), alpha = 0.9, cex = 1.0) +
  geom_point(alpha = 0.7) +
  scale_colour_viridis_c(option = "viridis", na.value = "grey", direction = -1, end = 0.9) +
  # scale_alpha_continuous(trans = "log10") +
  labs(x = "", y = "", col = "% > 30 days") + 
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9,0.85)) -> map_propgt30
map_propgt30

ggsave(here::here(figdir, "village_propgt30.png"),
       map_propgt30,
       height = 7, width = 10, units = "in")


## % greater than 90 days
ggplot() +
  geom_sf(data = boundary, fill = NA) +
  geom_sf(data = by_village, aes(col = propgt90), alpha = 0.9, cex = 1.0) +
  geom_point(alpha = 0.7) +
  scale_colour_viridis_c(option = "viridis", na.value = "grey", direction = -1, end = 0.9) +
  # scale_alpha_continuous(trans = "log10") +
  labs(x = "", y = "", col = "% > 90 days") + 
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9,0.85)) -> map_propgt90
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
            pgt90 = mean(gt90, na.rm = T)*100,
            pgt60 = mean(gt60, na.rm = T)*100,
            pgt30 = mean(gt30, na.rm = T)*100,
            pACD = mean(poss_acd*100)) %>%
  mutate(delay_cat3 = cut(med_delay, c(0,15,90,730), include.lowest = TRUE, ordered_result = TRUE)) %>%
  ungroup() -> by_block

by_block <- blockmap %>%
  left_join(by_block, by = c("kamis_master_dist" = "district", "kamis_master_block" = "block")) %>%
  rowwise() %>%
  dplyr::mutate(pop = median(c_across(`2018`:`2019`)),
         inc = N*1e4/pop) %>%
  ungroup() 

ggplot() +
  geom_sf(data = by_block, aes(geometry = geometry, fill = inc)) +
  geom_sf(data = blockmap, fill = NA) +
  scale_fill_viridis_c(option = "plasma", na.value = "white", direction = 1, trans = "log10", end = 0.9) +
  labs(fill = "Incidence\nper 10,000") + 
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9,0.85)) -> blk_inc
blk_inc

saveRDS(by_block, here::here("data","geography","by_block.rds"))
ggsave(here::here(figdir, "block_incidence.png"), blk_inc, height = 6, width = 8, units = "in")

# pal3 <- viridis::viridis(3)
ggplot() +
  geom_sf(data = blockmap) +
  geom_sf(data = by_block, aes(geometry = geometry, fill = med_delay)) +
  scale_fill_viridis_c(option = "viridis", na.value = "white", direction = -1, trans = "sqrt", end = 0.9) +
  labs(fill = "Median delay") + 
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9,0.85)) -> blk_med_delay
blk_med_delay

ggsave(here::here(figdir, "block_med_delay.png"), blk_med_delay, height = 6, width = 8, units = "in")

ggplot() +
  geom_sf(data = blockmap) +
  geom_sf(data = by_block, aes(geometry = geometry, fill = pgt30)) +
  scale_fill_viridis_c(option = "viridis", na.value = "white", direction = -1, end = 0.9) +
  labs(fill = "% > 30 days") + 
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9,0.85)) -> blk_pgt30
blk_pgt30

ggsave(here::here(figdir, "block_propgt30.png"), blk_pgt30, height = 6, width = 8, units = "in")


ggplot() +
  geom_sf(data = blockmap) +
  geom_sf(data = by_block, aes(geometry = geometry, fill = pgt60)) +
  scale_fill_viridis_c(option = "viridis", na.value = "white", direction = -1, end = 0.9) +
  labs(fill = "% > 60 days") +
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9,0.85))  -> blk_pgt60
blk_pgt60

ggsave(here::here(figdir, "block_propgt60.png"), blk_pgt60, height = 6, width = 8, units = "in")

ggplot() +
  geom_sf(data = blockmap) +
  geom_sf(data = by_block, aes(geometry = geometry, fill = pgt90)) +
  scale_fill_viridis_c(option = "viridis", na.value = "white", direction = -1, end = 0.9) +
  labs(fill = "% > 90 days") +
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9,0.85)) -> blk_pgt90
blk_pgt90

ggsave(here::here(figdir, "block_propgt90.png"), blk_pgt90, height = 6, width = 8, units = "in")

# Join two block maps
combined <- blk_pgt30 + blk_pgt60 +
  plot_annotation(title = "Proportion of cases diagnosed with excess delays",
                  caption = paste0("Diagnoses between ",
                                   range(dat$diag_month)[1],
                                   " and ",
                                   range(dat$diag_month)[2]))
combined

ggsave(here::here(figdir, "excess_delay_byblock.png"), combined, height = 7, width = 14, units = "in")

# ---------------------------------------------------------------------------- #
# % ACD by block

ggplot() +
  geom_sf(data = blockmap) +
  geom_sf(data = by_block, aes(geometry = geometry, fill = pACD)) +
  scale_fill_viridis_c(option = "magma", na.value = "white", direction = -1, end = 0.9) +
  labs(fill = "% via ACD") +
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9,0.85)) -> blk_pACD
blk_pACD

ggsave(here::here(figdir, "block_propACD.png"), blk_pACD, height = 6, width = 8, units = "in")

# ---------------------------------------------------------------------------- #
# Bivariate map of delay and incidence

# by_block %>%
#   filter(!is.na(pgt30)) %>%
#   mutate(inc = replace_na(inc, 0)) -> plotdat
# 
# bv <- create_bivmap(blockmap, plotdat,
#                    xvar = "pgt30", yvar = "inc",
#                    xlab = "% > 30 days", ylab = "Cases/10,000",
#                    pal = "DkBlue")
# bv
# 
# ggsave(here::here(figdir, "block_propgt30_byinc_bv.png"), bv, height = 7, width = 10, units = "in")

# ---------------------------------------------------------------------------- #
# Bivariate map of delay and % ACD
# 
# by_block %>%
#   filter(!is.na(pgt30)) -> plotdat
# 
# bv <- create_bivmap(blockmap, plotdat,
#                     xvar = "pgt30", yvar = "pACD",
#                     xlab = "% > 30 days", ylab = "% via ACD",
#                     pal = "DkBlue")
# bv
# 
# ggsave(here::here(figdir, "block_propgt30_byACD_bv.png"), bv, height = 7, width = 10, units = "in")

# ---------------------------------------------------------------------------- #
# Semi-variogram

dat.sp <- dat %>%
   st_drop_geometry() %>%
   filter(!is.na(longitude) & !is.na(delay))

# convert simple data frame into a spatial data frame object
coordinates(dat.sp) <- ~ longitude + latitude

vg <- variogram(delay~1, data = dat, cressie = TRUE)
plot(vg)

vgmod <- vgm(psill =  300, model = "Mat", nugget = 400, range = 75)
plot(vg, model = vgmod)

vgfit <- fit.variogram(vg, model = vgmod, fit.kappa = FALSE)    
plot(vg, model = vgfit)

vgfit
# model    psill    range kappa
# Nug 331.6653  0.00000   0.0
# Mat 348.8659 49.93885   0.5

png(here::here(figdir, "delay_semivariogram.png"), height = 500, width = 600)
plot(vg, model = vgfit, xlab = "Distance (km)")
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

dat %>%
  dplyr::rename(Sex = sex,
                Age = age_cat,
                `Scheduled caste or tribe` = caste4_r,
                Occupation = occ4_cat,
                `HIV status` = comorb,
                `Previous VL/PKDL treatment` = prv_tx,
                Detection = poss_acd,
                `Block endemic in 2017` = block_endm_2017,
                `Village IRS targeted in 2017` = IRS_2017,
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