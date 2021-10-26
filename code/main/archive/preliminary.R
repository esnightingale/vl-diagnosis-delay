################################################################################
# Description: Descriptive analyses of diagnosis delays given village location,
# recent incidence and active or passive detection.
################################################################################
################################################################################

library(tidyverse)
library(INLA)
library(inlabru)
library(raster)
library(sf)
library(ggmap)
library(mapr)

theme_set(theme_minimal())
figdir <- "figures/descriptive"

dat <- readRDS(here::here("data","analysisdata_individual.rds")) %>%
  dplyr::mutate(detection = factor(poss_acd, labels = c("PCD","ACD")))

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
# Compare delays by detection route

dat %>%
  ggplot(aes(x = poss_acd, y = days_fever_cdf)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log2") +
  labs(x = "Detection through active surveillance",
       y = "Days fever",
       title = "Days fever prior to diagnosis, by detection route")

ggsave(here::here(figdir, "delay_vs_acd.png"), height = 8, width= 10, units = "in")

# ---------------------------------------------------------------------------- #
# Map delays by detection route

ggmap(bh_lines, 
      base_layer = ggplot(arrange(st_jitter(dat), 
                                  poss_acd, days_fever_cdf), 
                          aes(col = days_fever_cdf))) +
  geom_sf(alpha = 0.5, cex = 0.8) +
  scale_colour_viridis_c(trans = "log2", direction = -1, option = "plasma") +
  facet_wrap(~detection) +
  labs(x = "", y = "", col = "Delay", 
       title = "Onset to diagnosis delay, by method of detection",
       caption = "Diagnoses between Jan 2018 - July 2019") -> map_byacd

map_byacd

ggsave(here::here(figdir,"map_delay_acd_pcd.png"), height = 7, width = 10, units = "in")


# ---------------------------------------------------------------------------- #
# Proportion of cases detected via ACD - by village

ggmap(bh_lines, 
      base_layer = ggplot(village,
                          aes(col = p_acd*100))) +
  geom_sf(alpha = 0.5) +
  scale_colour_viridis_c(direction = -1, option = "viridis") +
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
            n_acd = sum(poss_acd),
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

ggsave(here::here(figdir, "block_propACD_byinc_bv.png"), final, height = 7, width = 10, units = "in")

# ---------------------------------------------------------------------------- #

access.df <- as.data.frame(rasterToPoints(access)) %>%
  mutate(traveltime_adj = diag_facility_travel_time + 0.01)

# Add 0.01 to all times to avoid zeros in log transformation
ggplot() +
  geom_tile(data = access.df, aes(x = x, y = y, fill = traveltime_adj)) +
  geom_sf(data = boundary, fill = NA, lty = "dashed") +
  geom_sf(data = dat, col = "white", cex = 0.5) +
  scale_fill_viridis_c(trans = "log2") +
  labs(fill = "Time (minutes)") + 
  coord_sf(crs = st_crs(4326))

ggsave(here::here(figdir, "villages_vs_traveltime.png"), height = 8, width= 10, units = "in")
