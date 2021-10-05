################################################################################
# Description: Load travel time raster from Malaria Atlas Project
################################################################################
################################################################################

# Required packages
library(malariaAtlas)
library(gdistance)
library(abind)

# ---------------------------------------------------------------------------- #

# Get shapefile for Bihar in India - transform projection from lat/long

ind.shp <- malariaAtlas::getShp(ISO = "IND", admin_level = "admin1")
bh.shp <- ind.shp[ind.shp@data$name_1 == "Bihar",] 

# Retrieve friction surface from MAP
friction <- malariaAtlas::getRaster(
  surface = "A global friction surface enumerating land-based travel speed for a nominal year 2015",
  extent = bbox(bh.shp))

malariaAtlas::autoplot_MAPraster(friction)

# ---------------------------------------------------------------------------- #

# Convert friction surface to transition matrix
Tmat <- gdistance::transition(friction, function(x) 1/mean(x), 8) 

# Correct for 3-dimensional surface of the earth
Tmat.GC <- gdistance::geoCorrection(Tmat) 

# ---------------------------------------------------------------------------- #

# Read in health facility locations
facility <- read.csv(here::here("data","covariates",
                                    "Health facilties coordinates-revised-7-mar-2018.csv"), 
                         header = TRUE) %>%
  dplyr::filter(State == "Bihar") %>%
  tidyr::separate(Coordinates, into = c("y", "x"), sep = ",", convert = TRUE) %>%
  dplyr::mutate(diag = toupper(Diagnosis.Facility.Available),
                trt = toupper(LAmB.Treatment.Available),
                Capacity = as.factor(
                  case_when(diag == "YES" & trt == "NO"~ "Diagnosis",
                            trt == "YES" & diag == "NO" ~ "LAmB Treatment",
                            diag == "YES" & trt == "YES" ~ "Diagnosis + treatment",
                            TRUE ~ "None")))

sp::coordinates(facility) <- ~ x + y
sp::proj4string(facility) <- sp::proj4string(bh.shp)

# Keep only point coordinates within the shapefile bounds (allow a small buffer)
overlap <- sp::over(facility, raster::buffer(bh.shp, 1))
facility <- facility[!is.na(overlap),]

# sp::plot(bh.shp, main = "Health facilities with VL diagnosis/treatment capacity in Bihar")
# sp::plot(facility, add = T)

ggplot() +
  geom_sf(data = st_as_sf(bh.shp), aes(geometry = geometry), fill = "white") +
  geom_sf(data = st_as_sf(facility), aes(geometry = geometry, col = Capacity)) +
  labs(title = "Health facilities with VL capacity in endemic districts of Bihar") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()) -> map_facilities

ggsave(here::here("figures","covariates","facilities.png"), map_facilities, height = 5, width = 7, units = "in")

# Extract only those with diagnosis capacity
diag.facility <- facility[facility@data$diag == "YES",]

# Convert to matrix of lat/long
points <- as.matrix(diag.facility@coords)

# ---------------------------------------------------------------------------- #

# Calculate travel times via the "accumulated cost surface" algorithm
access.raster <- gdistance::accCost(Tmat.GC, diag.facility)

# Save raster of travel time to facility
raster::writeRaster(access.raster, here::here("data","covariates","diag_facility_travel_time.tif"), overwrite = TRUE)

# ---------------------------------------------------------------------------- #

# Redefine infinite values as NA to avoid showing in plot
access.raster@data@values[is.infinite(access.raster@data@values)] <- NA

# Plot constructed surface
p <- malariaAtlas::autoplot_MAPraster(access.raster, 
                                      shp_df = bh.shp, printed = F)

full_plot <- p[[1]] +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(fill = NA, colour = "white"),
        panel.background = element_rect(fill = "white", colour = "white")) +
  scale_fill_viridis_c(trans = "sqrt") +
  labs(title = "Travel time to most accessible health facility",
       fill = "Minutes \nof travel")

print(full_plot)
ggsave(here::here("figures","covariates","traveltime.png"), full_plot, height = 5, width = 7, units = "in")

# Overlay facility locations
full_plot <- full_plot +
  geom_point(data = data.frame(diag.facility@coords), aes(x = x, y = y),
           col = "white", cex = 0.8) 
ggsave(here::here("figures","covariates","traveltime_wfacilities.png"), full_plot, height = 5, width = 7, units = "in")

# ---------------------------------------------------------------------------- #
# 
# villages <- filter(match, !is.na(longitude))
# coordinates(villages) <- ~ longitude + latitude
# 
# vil.access <- data.frame(coordinates(villages),
#                          vil_code = villages$vil_code,
#                          sex = factor(villages$SEX),
#                          agecat = factor(villages$agecat),
#                          caste = factor(villages$caste4_r),
#                          hiv = factor(villages$HIV, levels = c(2,1), labels = c("Neg","Pos")),
#                          prev_trt_ka = factor(villages$Prv_TX_KA),
#                          acd = factor(villages$POSS_ACD),
#                          days_fever = villages$DUR_FEV_R,
#                          fever_gt30 = (villages$DUR_FEV_R > 30),
#                          fever_gt90 = (villages$DUR_FEV_R > 90),
#                          traveltime = raster::extract(access.raster, villages))
# 
# ggplot(vil.access, aes(traveltime, days_fever)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth() +
#   scale_x_continuous(trans = "log2") +
#   scale_y_continuous(trans = "log2")
# 
# ggplot(vil.access, aes(fever_gt90, traveltime)) +
#   geom_boxplot() +
#   scale_y_continuous(trans = "log2")
# 
# ggplot(vil.access, aes(sex, days_fever)) +
#   geom_boxplot() +
#   scale_y_continuous(trans = "log2")
# 
# ggplot(vil.access, aes(agecat, days_fever)) +
#   geom_boxplot() +
#   scale_y_continuous(trans = "log2")
# 
# ggplot(vil.access, aes(caste, days_fever)) +
#   geom_boxplot() +
#   scale_y_continuous(trans = "log2")
# 
# ggplot(vil.access, aes(hiv, days_fever)) +
#   geom_boxplot() +
#   scale_y_continuous(trans = "log2")
# 
# ggplot(vil.access, aes(prev_trt_ka, days_fever)) +
#   geom_boxplot() +
#   scale_y_continuous(trans = "log2")
# 
# ggplot(vil.access, aes(acd, days_fever)) +
#   geom_boxplot() +
#   scale_y_continuous(trans = "log2")
# 
# 
# # ---------------------------------------------------------------------------- #
# # Exploratory model
# 
# dat <- vil.access %>%
#   dplyr::select(longitude, latitude, days_fever, fever_gt30, fever_gt30, sex, agecat, caste, hiv, prev_trt_ka, acd, traveltime) %>%
#   na.omit()
# 
# coordinates(dat) <- ~ longitude + latitude
# 
# mod <- glm(days_fever ~ sex + agecat + caste + hiv + prev_trt_ka + acd + traveltime, 
#            data = dat,
#            family = "poisson")
# 
# summary(mod)
# plot(mod)
# 
# # ---------------------------------------------------------------------------- #
# # Semi-variogram of model residuals
# 
# dat$resid <- residuals(mod, "pearson")
# 
# vg <- variogram(resid~1, data = dat)
# plot(vg)
# 
# vgmod <- vgm(psill =  15, model = "Mat", nugget = 20, range = 0.5)
# plot(vg, model = vgmod)
# 
# vgfit <- fit.variogram(vg, model = vgmod)    
# plot(vg, model = vgfit)
# 
# vgfit
# # model    psill     range kappa
# # Nug 19.04896 0.0000000   0.0
# # Mat 21.70701 0.8020592   0.5
# 
# png(here::here("figures","fit", "resid_variogram.png"), height = 600, width = 700, res = 150)
# plot(vg, model = vgfit)
# dev.off()
# 
################################################################################
################################################################################