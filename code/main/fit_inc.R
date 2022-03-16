################################################################################
# Description: Univariate regression on covariates and assessment of empirical 
# variogram for spatial dependence
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit/exploratory/"
outdir <- "output/exploratory"

dat <- read_data()
mesh <- readRDS(here::here("data/analysis","mesh.rds"))

# Setup map context
boundary <- readRDS(here::here("data","geography","boundary.rds"))

boundary.spdf <- as_Spatial(boundary)

blockmap <- readRDS(here::here("data","geography","blockmap.rds")) %>% 
  sf::st_set_crs(7759)

by_block <- readRDS(here::here("data","geography","by_block.rds"))

#------------------------------------------------------------------------------#
# Define formulae

spde <- inla.spde2.pcmatern(mesh = mesh, 
                            prior.range = c(10, 0.01), # P(range < U) = a
                            prior.sigma = c(2, 0.01), # P(sigma > U) = a
                            constr = TRUE)
   
indexs <- inla.spde.make.index("s", spde$n.spde)

#------------------------------------------------------------------------------#
# Aggregate detection data by village

dat %>% 
  group_by(district, block, v, latitude, longitude, traveltime, traveltime_t, 
           inc_2017_gt0, IRS_2017, block_endm_2017, population) %>% 
  dplyr::summarise(n_detected = n(),
                   n_acd = sum(poss_acd == TRUE),
                   population = unique(population)) %>% 
  ungroup() -> dat_agg

summary(dat_agg)

#------------------------------------------------------------------------------#
# Generate a projection matrix A

# coo <- st_coordinates(dat)
coo <- st_coordinates(dat_agg)

A <- inla.spde.make.A(mesh = mesh, loc = coo)

dim(A)
# 2302 3574

# ---------------------------------------------------------------------------- #
# Make stack

stk <- inla.stack(
  data = list(y = dat_agg$n_detected,
              E = dat_agg$population),
  # data = list(y = as.numeric(dat$poss_acd)),
  A = list(A, 1),
  effects = list(s = indexs,  # the spatial index,
                 data.frame(
                   Intercept = rep(1, nrow(dat_agg)))
  )
)

# ---------------------------------------------------------------------------- #
# Fit model

f <- y ~ -1 + Intercept + f(s, model = spde)

fit.inc <- inla(f,
                 family = "poisson",
                 E = E,
                 data = inla.stack.data(stk, spde = spde),
                 control.predictor = list(
                         compute = TRUE, A = inla.stack.A(stk)),
                 control.compute = list(dic = FALSE, 
                                              waic = TRUE, 
                                              config = TRUE),
                 verbose = TRUE)

saveRDS(fit.inc, here::here(outdir, "fit_inc.rds"))
fit.inc$summary.hyperpar[,c(1,3,5)]
#                  mean 0.025quant 0.975quant
# Range for s 7.6117561  5.6146997 10.0431934
# Stdev for s 0.8426712  0.7329724  0.9731476

# ---------------------------------------------------------------------------- #
# Project fitted SPDE

p <- plot_spde_mean(fit.inc, title = "", trans = TRUE, palopt = "plasma")

# no_inc <- which(is.na(by_block$inc))
# blkfill <- rep(NA, nrow(blockmap))
# blkfill[no_inc] <- "white"

p +
  geom_sf(data = blockmap) + #, fill = blkfill
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = c(0.9, 0.85)) +
  scale_fill_viridis_c(option = "plasma",direction = 1) +
  labs(#title = "Predicted mean of days fever prior to diagnosis", 
    x = "", y = "", fill = "Incidence rate") -> spde_inc

ggsave(here::here(figdir,"spde_inc.png"), spde_acd, height = 6, width = 8, units = "in")

################################################################################
################################################################################