################################################################################
# Description: Basic plotting of diagnosis delay against potential covariates
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/descriptive"
datadir <- "~/VL/Data/KAMIS/Clean/ll"

dat <- read_data() %>%
  dplyr::mutate(inc_2017_gt0 = factor(as.numeric(replace_na(inc_2017,0) > 0), 
                                             levels = c(0,1), labels = c("No","Yes")),
                inc_2017_gt1e3 = factor(as.numeric(replace_na(inc_2017,0) > 1), 
                                      levels = c(0,1), labels = c("No","Yes")),
                inc_2017_t = inc_2017*1e3 + 1e-4,
                IRS_2017 = factor(as.numeric(IRS_2017 != 0), 
                                    levels = c(0,1), labels = c("No","Yes")),
                # Add 0.5 to all times to avoid zeros in log transformation
                traveltime_adj = traveltime + 0.5,
                delay_gt30 = (delay > 30),
                delay_gt90 = (delay > 90))

dat.df <- sf::st_drop_geometry(dat) 

# Setup map context
blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
  st_transform(4326)
boundary <- sf::st_union(blockmap)
boundary.spdf <- as_Spatial(boundary)

extent <- setNames(st_bbox(blockmap), c("left","bottom","right","top"))
bh_lines <- get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)

# Travel time rasters
access_d <- raster::raster(here::here("data","covariates","diag_facility_travel_time.tif"))
access_t <- raster::raster(here::here("data","covariates","trt_facility_travel_time.tif"))

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

dat.df <- filter(dat.df, delay >= 0)

dat.df %>%
  ggplot(aes(age, delay, col = detection)) +
  geom_point(alpha = 0.3, cex = 0.8) +
  geom_smooth() +
  scale_y_continuous(trans = "log2") +
  labs(x = "Age (years)", y = "Delay", col = "") -> plot_age

create_raincloud(dat.df, "age_child", "delay", 
                 xlab = "", ylab = "Delay", 
                 y_trans = "log2", col_by = "age_child") -> plot_child

create_raincloud(dat.df, "age_cat", "delay", 
                 xlab = "", ylab = "Delay",
                 y_trans = "log2", col_by = "age_cat") -> plot_agecat

create_raincloud(dat.df, "sex", "delay", 
                 xlab = "Sex", ylab = "Delay", 
                 y_trans = "log2", col_by = "sex") -> plot_sex

create_raincloud(dat.df, "marg_caste", "delay", 
                 xlab = "Marginalised caste", ylab = "Delay", 
                 y_trans = "log2", col_by = "marg_caste", 
                 drop_na = TRUE) -> plot_caste

create_raincloud(dat.df, "occupation", "delay", 
                 xlab = "Occupation", ylab = "Delay", 
                 y_trans = "log2", col_by = "occupation", 
                 drop_na = TRUE) -> plot_occ

create_raincloud(dat.df, "hiv", "delay", 
                 xlab = "HIV positive", ylab = "Delay", 
                 y_trans = "log2", col_by = "hiv", 
                 drop_na = TRUE) -> plot_hiv

create_raincloud(dat.df, "prv_tx", "delay", 
                 xlab = "Previous treatment for VL or PKDL", ylab = "Delay", 
                 y_trans = "log2", col_by = "prv_tx", 
                 drop_na = TRUE) -> plot_ptvl

create_raincloud(dat.df, "num_conslt", "delay", 
                 xlab = "No. consultations prior to diagnosis", 
                 ylab = "Delay", 
                 y_trans = "log2", col_by = "num_conslt", 
                 drop_na = TRUE) -> plot_conslt

create_raincloud(dat.df, "conslt_cat", "delay", 
                 xlab = "No. consultations prior to diagnosis", 
                 ylab = "Delay", 
                 y_trans = "log2", col_by = "conslt_cat", 
                 drop_na = TRUE) -> plot_consltcat

create_raincloud(dat.df, "diag_rainseason", "delay", 
                 xlab = "Diagnosis during rainy season", 
                 ylab = "Delay", 
                 y_trans = "log2", col_by = "diag_rainseason", 
                 drop_na = TRUE) -> plot_season

create_raincloud(dat.df, "detection", "delay", 
                 xlab = "Detection route", ylab = "Delay", 
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
plot_season
plot_acd

dev.off()

# ---------------------------------------------------------------------------- #
# Delay versus village characteristics
# + Recent incidence
# + IRS
# + Block endemicity

dat.df %>%
  ggplot(aes(inc_2017*1000, delay, col = detection)) +
  geom_point(alpha = 0.3, cex = 0.8) +
  geom_smooth() +
  scale_y_continuous(trans = "log2") +
  scale_x_continuous(trans = "log10") +
  labs(x = "Village incidence per 1,000 (2017)", 
       y = "Delay", 
       col = "") -> plot_vilinc

create_raincloud(dat.df, "inc_2017_gt1e3", "delay", 
                 xlab = "Greater than 1/1,000 incidence in 2017", ylab = "Delay", 
                 y_trans = "log2", col_by = "inc_2017_gt1e3") -> plot_vilinc1

create_raincloud(dat.df, "inc_2017_gt0", "delay", 
                 xlab = "Non-zero incidence in 2017", ylab = "Delay", 
                 y_trans = "log2", col_by = "inc_2017_gt0") -> plot_vilinc1

create_raincloud(dat.df, "block_endm_2017", "delay", 
                 xlab = "Block status in 2017", ylab = "Delay", 
                 y_trans = "log2", col_by = "block_endm_2017") -> plot_blockendm

# Any rounds of IRS
create_raincloud(dat.df, "IRS_2017", "delay", 
                 xlab = "Any IRS in 2017", ylab = "Delay", 
                 y_trans = "log2", col_by = "IRS_2017") -> plot_irs

png(here::here(figdir, "delay_vs_vilchar.png"), 
    height = 8, width = 10, units = "in", res = 300)
gridExtra::grid.arrange(plot_vilinc,
                        plot_vilinc1,
                        plot_blockendm,
                        plot_irs,
                        nrow = 2)
dev.off()

pdf(here::here(figdir, "delay_vs_vilchar.pdf"), height = 5, width = 6)
plot_vilinc
plot_vilinc1
plot_blockendm
plot_irs
dev.off()

# ---------------------------------------------------------------------------- #
# Delays versus travel time

summary(dat.df$travel_time_cat)
# [0,10]  (10,20]  (20,30]  (30,60] (60,124] 
# 1889     1576      545      332       51 

summary(dat.df$travel_time_t_cat)
# [0,10]  (10,20]  (20,30]  (30,60] (60,150] 
# 999     1485     1022      772      101 

pdf(here::here(figdir, "delay_vs_traveltime.pdf"), height = 5, width = 6)

ggplot(dat.df, aes(x = traveltime)) +
  geom_histogram(binwidth = 5) +
  # scale_x_continuous(trans = "log2") +
  labs(title = "Travel time to most accessible diagnosis facility",
       subtitle = paste0("Median [IQR]: ", round(median(dat.df$traveltime)), " [", 
                         round(quantile(dat.df$traveltime, 0.25)),
                         ", ",
                         round(quantile(dat.df$traveltime, 0.75)),
                         "] minutes."),
       x = "Travel time (minutes)", y = "Count")

ggplot(dat.df, aes(x = traveltime_t)) +
  geom_histogram(binwidth = 5) +
  # scale_x_continuous(trans = "log2") +
  labs(title = "Travel time to most accessible treatment facility",
       subtitle = paste0("Median [IQR]: ", round(median(dat.df$traveltime_t)), " [", 
                         round(quantile(dat.df$traveltime, 0.25)),
                         ", ",
                         round(quantile(dat.df$traveltime, 0.75)),
                         "] minutes."),
       x = "Travel time (minutes)", y = "Count")

create_raincloud(dat.df, "travel_time_cat", "delay", xlab = "Travel time (minutes)", ylab = "Delay", 
                 y_trans = "log2", col_by = "travel_time_cat", drop_na = TRUE)

create_raincloud(dat.df, "travel_time_t_cat", "delay", xlab = "Travel time (minutes)", ylab = "Delay", 
                 y_trans = "log2", col_by = "travel_time_t_cat", drop_na = TRUE)

dat.df %>%
  ggplot(aes(traveltime_t, delay)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  labs(y = "Delay (days)", x =  "Travel time (minutes)")

dat.df %>%
  ggplot(aes(traveltime_t, delay, col = age_child)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  labs(y = "Delay (days)", x =  "Travel time (minutes)", col = "")

dat.df %>%
  ggplot(aes(traveltime_t, delay, col = sex)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  labs(y = "Delay (days)", x =  "Travel time (minutes)", col = "Sex")

dat.df %>%
  filter(!is.na(marg_caste)) %>%
  ggplot(aes(traveltime_t, delay, col = marg_caste)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  labs(y = "Delay (days)", x =  "Travel time (minutes)", col = "Marginalised\ncaste")

dat.df %>%
  filter(!is.na(occupation)) %>%
  ggplot(aes(traveltime_t, delay, col = occupation)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  labs(y = "Delay (days)", x =  "Travel time (minutes)", col = "Occupation")

dat.df %>%
  filter(!is.na(hiv)) %>%
  ggplot(aes(traveltime_t, delay, col = hiv)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  labs(y = "Delay (days)", x =  "Travel time (minutes)", col = "HIV status")

dat.df %>%
  ggplot(aes(traveltime_t, delay, col = diag_rainseason)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  labs(y = "Delay (days)", x =  "Travel time (minutes)", col = "Rainy season")

dat.df %>%
  ggplot(aes(traveltime_t, delay, col = detection)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  labs(y = "Delay (days)", x =  "Travel time (minutes)", col = "Detection")

dat.df %>%
  ggplot(aes(delay_gt30, traveltime_t)) +
  geom_boxplot() +
  labs(x = "Delay > 30 days", y =  "Travel time (minutes)")

dat.df %>%
  ggplot(aes(delay_gt90, traveltime_t)) +
  geom_boxplot() +
  labs(x = "Delay > 90 days", y =  "Travel time (minutes)")

dev.off()

access.df <- as.data.frame(raster::rasterToPoints(access_d)) %>%
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

dat.df %>%
  ggplot(aes(x = detection, y = delay)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log2") +
  labs(x = "Detection route",
       y = "Delay",
       title = "Delay to diagnosis, by detection route")

ggsave(here::here(figdir, "delay_vs_acd.png"), height = 5, width = 6, units = "in")

# ---------------------------------------------------------------------------- #
# Map delays by detection route

ggmap(bh_lines, 
      base_layer = ggplot(arrange(st_jitter(dat), 
                                  detection, delay), 
                          aes(col = delay))) +
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
  dplyr::summarise(N = n(),
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
  dplyr::summarise(N = n(),
            n_acd = sum(detection == "ACD"),
            p_acd = n_acd/N) %>%
  ungroup() -> by_block

by_block <- blockmap %>%
  left_join(by_block, by = c("kamis_master_dist" = "district", "kamis_master_block" = "block")) %>%
  rowwise() %>%
  dplyr::mutate(pop = median(c_across(`2018`:`2019`)),
         inc = replace_na(N*1e4/pop, 0)) %>%
  ungroup() 


p <- create_bivmap(blockmap, filter(by_block, !is.na(p_acd)),
                   yvar = "p_acd", xvar = "inc",
                   ylab = "% via ACD", xlab = "Cases/10,000",
                   pal = "DkViolet")
p

ggsave(here::here(figdir, "block_propACD_byinc_bv.png"), height = 7, width = 10, units = "in")

################################################################################
################################################################################

