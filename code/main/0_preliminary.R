################################################################################
# Description: Basic plotting of diagnosis delay against potential covariates
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/descriptive/covariate effects"

dat <- readRDS(here::here("data","analysisdata_individual.rds"))
dat.df <- st_drop_geometry(dat) %>%
  dplyr::mutate(vill_inc_2017_gt0 = factor(as.numeric(replace_na(vill_inc_2017,0) > 0), 
                                             levels = c(0,1), labels = c("No","Yes")),
                vill_inc_2017_t = vill_inc_2017*1e3 + 1e-4,
                IRS_2017_1 = factor(as.numeric(IRS_2017 != 0), 
                                    levels = c(0,1), labels = c("No","Yes")),
                diag_quarter = as.factor(lubridate::quarter(diag_month)),
                num_conslt = as.factor(num_conslt),
                # Add 0.5 to all times to avoid zeros in log transformation
                traveltime_adj = traveltime + 0.5)

village <- readRDS(here::here("data","analysisdata_village.rds")) 

# Setup map context
blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
  st_transform(4326)
boundary <- sf::st_union(blockmap)  
boundary.spdf <- as_Spatial(boundary)

extent <- setNames(st_bbox(blockmap), c("left","bottom","right","top"))
bh_lines <- get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)

# Travel time raster
access <- raster::raster(here::here("data","covariates","diag_facility_travel_time.tif"))

# ---------------------------------------------------------------------------- #
# Delay versus patient characteristics
# + Age
# + Sex
# + Marginalised caste (SC/ST)
# + HIV status
# + Number of consultations prior to formal diagnosis (<=2, 3-5, 6-8)
# + Occupation (None, unskilled, skilled, selfempl/salaried)
# + Season of diagnosis (month or rain/dry)
# + Active/passive detection

dat.df %>%
  ggplot(aes(age, days_fever, col = detection)) +
  geom_point(alpha = 0.3, cex = 0.8) +
  geom_smooth() +
  scale_y_continuous(trans = "log2") +
  labs(x = "Age (years)", y = "Days fever", col = "") -> plot_age

create_raincloud(dat.df, "age_child", "days_fever", 
                 xlab = "", ylab = "Age", # ylab = "Days fever"
                 y_trans = "log2", col_by = "age_child") -> plot_child

create_raincloud(dat.df, "age_cat", "days_fever", 
                 xlab = "", ylab = "Age", # ylab = "Days fever"
                 y_trans = "log2", col_by = "age_cat") -> plot_agecat

create_raincloud(dat.df, "sex", "days_fever", 
                 xlab = "Sex", ylab = "Days fever", 
                 y_trans = "log2", col_by = "sex") -> plot_sex

create_raincloud(dat.df, "marg_caste", "days_fever", 
                 xlab = "Marginalised caste", ylab = "Days fever", 
                 y_trans = "log2", col_by = "marg_caste", 
                 drop_na = TRUE) -> plot_caste

create_raincloud(dat.df, "occ4_cat", "days_fever", 
                 xlab = "Occupation", ylab = "Days fever", 
                 y_trans = "log2", col_by = "occ4_cat", 
                 drop_na = TRUE) -> plot_occ

create_raincloud(dat.df, "hiv", "days_fever", 
                 xlab = "HIV positive", ylab = "Days fever", 
                 y_trans = "log2", col_by = "hiv", 
                 drop_na = TRUE) -> plot_hiv

create_raincloud(dat.df, "prv_tx", "days_fever", 
                 xlab = "Previous treatment for VL or PKDL", ylab = "Days fever", 
                 y_trans = "log2", col_by = "prv_tx", 
                 drop_na = TRUE) -> plot_ptvl

create_raincloud(dat.df, "num_conslt", "days_fever", 
                 xlab = "No. consultations prior to diagnosis", 
                 ylab = "Days fever", 
                 y_trans = "log2", col_by = "num_conslt", 
                 drop_na = TRUE) -> plot_conslt

create_raincloud(dat.df, "conslt_cat", "days_fever", 
                 xlab = "No. consultations prior to diagnosis", 
                 ylab = "Days fever", 
                 y_trans = "log2", col_by = "conslt_cat", 
                 drop_na = TRUE) -> plot_consltcat

create_raincloud(dat.df, "diag_quarter", "days_fever", 
                 xlab = "Quarter of diagnosis", 
                 ylab = "Days fever", 
                 y_trans = "log2", col_by = "diag_quarter", 
                 drop_na = TRUE) -> plot_quarter

create_raincloud(dat.df, "diag_rainseason", "days_fever", 
                 xlab = "Diagnosis during rainy season", 
                 ylab = "Days fever", 
                 y_trans = "log2", col_by = "diag_rainseason", 
                 drop_na = TRUE) -> plot_season

create_raincloud(dat.df, "detection", "days_fever", 
                 xlab = "Detection route", ylab = "Days fever", 
                 y_trans = "log2", col_by = "detection", 
                 drop_na = TRUE) -> plot_acd

png(here::here(figdir, "delay_vs_patchar.png"), 
    height = 10, width = 16, units = "in", res = 300)
gridExtra::grid.arrange(plot_age,
                        plot_child,
                        plot_agecat,
                        plot_sex,
                        plot_caste,
                        plot_hiv,
                        plot_occ,
                        plot_ptvl,
                        plot_consltcat,
                        plot_quarter,
                        plot_season,
                        plot_acd,
                        nrow = 3)
dev.off()

pdf(here::here(figdir, "delay_vs_patchar.pdf"), height = 5, width = 6)

plot_age
plot_child
plot_agecat
plot_sex
plot_caste
plot_hiv
plot_occ
plot_ptvl
plot_conslt
plot_consltcat
plot_quarter
plot_season
plot_acd

dev.off()

# ---------------------------------------------------------------------------- #
# Delay versus village characteristics
# + Recent incidence
# + IRS
# + Block endemicity

dat.df %>%
  ggplot(aes(vill_inc_2017*1000, days_fever, col = detection)) +
  geom_point(alpha = 0.3, cex = 0.8) +
  geom_smooth() +
  scale_y_continuous(trans = "log2") +
  scale_x_continuous(trans = "log10") +
  labs(x = "Village incidence per 1,000 (2017)", 
       y = "Days fever", 
       col = "") -> plot_vilinc

create_raincloud(dat.df, "vill_inc_2017_gt1e3", "days_fever", 
                 xlab = "Greater than 1/1,000 incidence in 2017", ylab = "Days fever", 
                 y_trans = "log2", col_by = "vill_inc_2017_gt1e3") -> plot_vilinc1

create_raincloud(dat.df, "vill_inc_2017_gt0", "days_fever", 
                 xlab = "Non-zero incidence in 2017", ylab = "Days fever", 
                 y_trans = "log2", col_by = "vill_inc_2017_gt0") -> plot_vilinc1

create_raincloud(dat.df, "block_endm_2017", "days_fever", 
                 xlab = "Block status in 2017", ylab = "Days fever", 
                 y_trans = "log2", col_by = "block_endm_2017") -> plot_blockendm

# Number of rounds
create_raincloud(dat.df, "IRS_2017", "days_fever", 
                 xlab = "Rounds of IRS in 2017", ylab = "Days fever", 
                 y_trans = "log2", col_by = "IRS_2017") -> plot_irs
# Any rounds
create_raincloud(dat.df, "IRS_2017_1", "days_fever", 
                 xlab = "Any IRS in 2017", ylab = "Days fever", 
                 y_trans = "log2", col_by = "IRS_2017_1") -> plot_irs1

png(here::here(figdir, "delay_vs_vilchar.png"), 
    height = 8, width = 10, units = "in", res = 300)
gridExtra::grid.arrange(plot_vilinc,
                        plot_vilinc1,
                        plot_blockendm,
                        plot_irs1,
                        nrow = 2)
dev.off()

pdf(here::here(figdir, "delay_vs_vilchar.pdf"), height = 5, width = 6)
plot_vilinc
plot_vilinc1
plot_blockendm
plot_irs
plot_irs1
dev.off()

# ---------------------------------------------------------------------------- #
# Delays versus travel time

summary(dat.df$travel_time_cat)
# [0,10]  (10,20]  (20,30]  (30,60] (60,124] 
# 1889     1576      545      332       51 

summary(dat.df$travel_time_cat4)
# [0,7]   (7,12]  (12,18] (18,124] 
# 1196     1071      981     1145 

pdf(here::here(figdir, "delay_vs_traveltime.pdf"), height = 5, width = 6)

ggplot(dat.df, aes(x = traveltime_adj)) +
  geom_histogram(binwidth = 5) +
  # scale_x_continuous(trans = "log2") +
  labs(title = "Travel time to most accessible facility",
       subtitle = paste0("Median [IQR]: ", round(median(dat.df$traveltime)), " [", 
                         round(quantile(dat.df$traveltime, 0.25)),
                         ", ",
                         round(quantile(dat.df$traveltime, 0.75)),
                         "] minutes."),
       x = "Travel time (minutes)", y = "Count")

create_raincloud(dat.df, "travel_time_cat", "days_fever", xlab = "Travel time (minutes)", ylab = "Days fever", 
                 y_trans = "log2", col_by = "travel_time_cat", drop_na = TRUE)

create_raincloud(dat.df, "travel_time_cat4", "days_fever", xlab = "Travel time (minutes)", ylab = "Days fever", 
                 y_trans = "log2", col_by = "travel_time_cat4", drop_na = TRUE)

dat.df %>%
  ggplot(aes(traveltime_adj, days_fever)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  scale_y_continuous(trans = "log2") +
  scale_x_continuous(trans = "log2") +
  labs(y = "Delay (days)", x =  "Travel time (minutes)")

dat.df %>%
  ggplot(aes(traveltime_adj, days_fever, col = age_child)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  scale_y_continuous(trans = "log2") +
  scale_x_continuous(trans = "log2") +
  labs(y = "Delay (days)", x =  "Travel time (minutes)", col = "")

dat.df %>%
  ggplot(aes(traveltime_adj, days_fever, col = sex)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  scale_y_continuous(trans = "log2") +
  scale_x_continuous(trans = "log2") +
  labs(y = "Delay (days)", x =  "Travel time (minutes)", col = "Sex")

dat.df %>%
  filter(!is.na(marg_caste)) %>%
  ggplot(aes(traveltime_adj, days_fever, col = marg_caste)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  scale_y_continuous(trans = "log2") +
  scale_x_continuous(trans = "log2") +
  labs(y = "Delay (days)", x =  "Travel time (minutes)", col = "Marginalised\ncaste")

dat.df %>%
  filter(!is.na(occ4_cat)) %>%
  ggplot(aes(traveltime_adj, days_fever, col = occ4_cat)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  scale_y_continuous(trans = "log2") +
  scale_x_continuous(trans = "log2") +
  labs(y = "Delay (days)", x =  "Travel time (minutes)", col = "Occupation")

dat.df %>%
  filter(!is.na(hiv)) %>%
  ggplot(aes(traveltime_adj, days_fever, col = hiv)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  scale_y_continuous(trans = "log2") +
  scale_x_continuous(trans = "log2") +
  labs(y = "Delay (days)", x =  "Travel time (minutes)", col = "HIV status")

dat.df %>%
  ggplot(aes(traveltime_adj, days_fever, col = diag_rainseason)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  scale_y_continuous(trans = "log2")  +
  scale_x_continuous(trans = "log2") +
  labs(y = "Delay (days)", x =  "Travel time (minutes)", col = "Rainy season")

dat.df %>%
  ggplot(aes(traveltime_adj, days_fever, col = detection)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  scale_y_continuous(trans = "log2")  +
  scale_x_continuous(trans = "log2") +
  labs(y = "Delay (days)", x =  "Travel time (minutes)", col = "Detection")

dat.df %>%
  ggplot(aes(gt30, traveltime_adj)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log2")  +
  labs(x = "Delay > 30 days", y =  "Travel time (minutes)")

dat.df %>%
  ggplot(aes(gt90, traveltime_adj)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log2")  +
  labs(x = "Delay > 90 days", y =  "Travel time (minutes)")

dev.off()

access.df <- as.data.frame(raster::rasterToPoints(access)) %>%
  # Add 0.01 to all times to avoid zeros in log transformation
  mutate(traveltime_adj = diag_facility_travel_time + 1)

ggplot() +
  geom_tile(data = access.df, aes(x = x, y = y, fill = traveltime_adj)) +
  geom_sf(data = boundary, fill = NA, lty = "dashed") +
  geom_sf(data = dat, col = "white", cex = 0.5) +
  scale_fill_viridis_c(trans = "log2") +
  labs(fill = "Time (minutes)") + 
  coord_sf(crs = st_crs(4326))

ggsave(here::here(figdir, "villages_vs_traveltime.png"), height = 8, width = 10, units = "in")


# ---------------------------------------------------------------------------- #
# Distribution of delays by detection route

dat %>%
  st_drop_geometry() %>%
  ggplot(aes(x = detection, y = days_fever)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log2") +
  labs(x = "Detection route",
       y = "Days fever",
       title = "Days fever prior to diagnosis, by detection route")

ggsave(here::here(figdir, "delay_vs_acd.png"), height = 5, width = 6, units = "in")

# ---------------------------------------------------------------------------- #
# Map delays by detection route

ggmap(bh_lines, 
      base_layer = ggplot(arrange(st_jitter(dat), 
                                  detection, days_fever), 
                          aes(col = days_fever))) +
  geom_sf(alpha = 0.5, cex = 0.8) +
  scale_colour_viridis_c(trans = "log2", direction = -1, option = "plasma") +
  facet_wrap(~detection) +
  labs(x = "", y = "", col = "Delay", 
       title = "Onset to diagnosis delay, by method of detection",
       caption = "Diagnoses between Jan 2018 - July 2019") -> map_byacd
map_byacd

ggsave(here::here(figdir,"map_delay_acd_pcd.png"), height = 6, width = 8, units = "in")

# ---------------------------------------------------------------------------- #
# Proportion of cases detected via ACD - by month of diagnosis

dat %>%
  st_drop_geometry() %>%
  group_by(diag_month) %>%
  summarise(N = n(),
            propacd = mean(detection == "ACD", na.rm = T),
            se = sqrt((propacd*(1 - propacd))/N),
            ll = propacd + se*qnorm(p = 0.025),
            ul = propacd + se*qnorm(p = 0.975),
            ACD = sum(detection == "ACD", na.rm = T)/N,
            PCD = sum(detection != "ACD", na.rm = T)/N) %>%
  ungroup() -> by_mth

by_mth %>%
  ggplot(aes(x = diag_month, group = 1, y = propacd, ymin = ll, ymax = ul)) +
  geom_errorbar(width = 0.5) +
  geom_point() +
  labs(x = "Month of diagnosis", y = "Proportion diagnosed via ACD") -> acd_bymth
acd_bymth

ggsave(here::here(figdir,"propacd_bymth.png"),
       acd_bymth,
       height = 4, width = 6, units = "in")


pal <- viridis::viridis(2, end = 0.8, direction = -1)
by_mth %>%
  pivot_longer(c("ACD","PCD"), names_to = "Detection") %>%
  ggplot(aes(x = diag_month, y = value, fill = Detection)) +
  geom_col() +
  scale_fill_manual(values = pal) +
  labs(x = "Month of diagnosis", y = "Proportion", fill = "Detected via") -> prop_acd
prop_acd

ggsave(here::here(figdir,"prop_acd_bymth.png"),
       prop_acd,
       height = 4, width = 6, units = "in")

# ---------------------------------------------------------------------------- #
# Proportion of cases detected via ACD - by village

ggmap(bh_lines, 
      base_layer = ggplot(village,
                          aes(col = p_acd*100))) +
  geom_sf(alpha = 0.5) +
  scale_colour_viridis_c(direction = -1, option = "viridis", trans = "sqrt") +
  labs(x = "", y = "", col = "% ACD", 
       title = "% village cases detected by ACD",
       caption = "Diagnoses between Jan 2018 - July 2019") -> map_p_acd
map_p_acd

ggsave(here::here(figdir,"map_perc_acd.png"), height = 6, width = 8, units = "in")

# ---------------------------------------------------------------------------- #
# Proportion of cases detected via ACD - by block

dat %>%
  st_drop_geometry() %>%
  group_by(district, block) %>%
  summarise(N = n(),
            n_acd = sum(detection == "ACD"),
            p_acd = n_acd/N) %>%
  ungroup() -> by_block

by_block <- blockmap %>%
  left_join(by_block, by = c("kamis_master_dist" = "district", "kamis_master_block" = "block")) %>%
  rowwise() %>%
  mutate(pop = median(c_across(`2018`:`2019`)),
         inc = replace_na(N*1e4/pop, 0)) %>%
  ungroup() 


p <- create_bivmap(blockmap, filter(by_block, !is.na(p_acd)),
                   yvar = "p_acd", xvar = "inc",
                   ylab = "% via ACD", xlab = "Cases/10,000",
                   pal = "DkViolet")
p

ggsave(here::here(figdir, "block_propACD_byinc_bv.png"), height = 7, width = 10, units = "in")

# ---------------------------------------------------------------------------- #
# Tabulate

make_tab <- function(dat, varname){

  dat %>%
    dplyr::mutate(Value = !!sym(varname)) %>%
    dplyr::group_by(Value) %>%
    dplyr::summarise(N = n(),
                     mean = round(mean(days_fever),1),
                     ci = paste(round(mean + sd(days_fever)/sqrt(N) * qnorm(p = c(0.025, 0.975)),1), collapse = ","),
                     med = median(days_fever),
                     iqr = paste(quantile(days_fever, p = c(0.25, 0.75)), collapse = ","),
                     n_gt30 = sum(gt30),
                     p_gt30 = round(n_gt30/N,2)) %>%
    dplyr::mutate(`Days fever, mean [95% CI]` = paste0(mean, " [",ci,"]"),
                  `Days fever, median [IQR]` = paste0(med, " [",iqr,"]"),
                  `> 30 days, N (%)` = paste0(n_gt30, " (",p_gt30,")"),
                  Variable = varname) %>%
    dplyr::select(Variable, Value, N, `Days fever, mean [95% CI]`, `> 30 days, N (%)`) %>%
    dplyr::arrange(Value) -> tab

  return(tab)
  
}

make_plot <- function(t) {
  
  t %>%
    filter(Value != "NA") %>%
    separate(`Days fever, mean [95% CI]`, into = c("mean","ll","ul"), sep = "[\\[\\,\\]]", convert = TRUE) %>%
    ggplot(aes(Value, mean, ymin = ll, ymax = ul)) +
    geom_linerange() +
    geom_point() +
    geom_hline(yintercept = 45, col = "grey") +
    labs(subtitle = unique(t$Variable), x = "", y =  "") %>%
    return()
  
}

dat.df %>%
  dplyr::rename(Sex = sex,
         Age = age_cat,
         `Scheduled caste or tribe` = marg_caste,
         Occupation = occ4_cat,
         `HIV status` = hiv,
         `Previous VL/PKDL treatment` = prv_tx,
         `Number of prior consultations` = conslt_cat,
         Detection = detection,
         `Block endemic in 2017` = block_endm_2017,
         `Village IRS targeted in 2017` = IRS_2017_1,
         `Village incidence > 0 in 2017` = vill_inc_2017_gt0,
         `Travel time to nearest facility` = travel_time_cat) -> dat.tab

# varlist <- list("sex","age_child", "marg_caste","hiv", "prv_tx_ka", "num_conslt_4","detection","block_endm_2017" ,"IRS_2017_1", "vill_inc_2017_gt0", "travel_time_cat")
varlist <- list("Sex","Age","Scheduled caste or tribe","Occupation","HIV status", 
                "Previous VL/PKDL treatment","Number of prior consultations",
                "Detection","Block endemic in 2017",
                "Village IRS targeted in 2017","Village incidence > 0 in 2017", 
                "Travel time to nearest facility")

tabs.list <- lapply(varlist, make_tab, dat = dat.tab)

plots <- lapply(tabs.list, make_plot)

png(here::here(figdir, "covariate_mean_ci_plot.png"), height = 7, width = 10, units = "in", res = 300)
do.call("grid.arrange", plots)
dev.off()

tab <- bind_rows(lapply(tabs.list, function(t) mutate(t, Value = as.character(Value))))
View(tab)

write.csv(tab, here::here("output","table_descriptive.csv"), row.names = FALSE)

################################################################################
################################################################################

