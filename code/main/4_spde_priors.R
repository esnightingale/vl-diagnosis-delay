################################################################################
# Description: Univariate regression on covariates and assessment of empirical 
# variogram for spatial dependence
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit/exploratory"
outdir <- "output/exploratory"

dat <- readRDS(here::here("data/analysis","dat_nona.rds"))  %>%
  st_transform(7759)

mesh <- readRDS(here::here("data/analysis","mesh.rds"))

# Setup map context
blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
  sf::st_transform(7759)

boundary <- blockmap %>%
  sf::st_union()

boundary.spdf <- as_Spatial(boundary)

# extent <- setNames(st_bbox(blockmap), c("left","bottom","right","top"))
# bh_lines <- get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)

# ---------------------------------------------------------------------------- #
# Define SPDE

# Shorter range - 5km
spde1 <- inla.spde2.pcmatern(mesh = mesh, 
                            prior.range = c(5000, 0.01), # P(range < U) = a
                            prior.sigma = c(1.5, 0.1), # P(sigma > U) = a
                            constr = TRUE)

# Longer range - 50km
spde2 <- inla.spde2.pcmatern(mesh = mesh, 
                             prior.range = c(5e5, 0.01), # P(range < U) = a
                             prior.sigma = c(1.5, 0.1), # P(sigma > U) = a
                             constr = TRUE)

# Smaller sigma
spde3 <- inla.spde2.pcmatern(mesh = mesh, 
                             prior.range = c(5000, 0.01), # P(range < U) = a
                             prior.sigma = c(1.1, 0.1), # P(sigma > U) = a
                             constr = TRUE)

# Greater sigma
spde4 <- inla.spde2.pcmatern(mesh = mesh, 
                             prior.range = c(5000, 0.01), # P(range < U) = a
                             prior.sigma = c(2, 0.1), # P(sigma > U) = a
                             constr = TRUE)

spde.list <- list(`R5e3-Sig1.5` = spde1, 
                  `R5e5-Sig1.5` = spde2, 
                  `R5e3-Sig1.1` = spde3, 
                  `R5e3-Sig2` = spde4)

#------------------------------------------------------------------------------#
# Generate the index set for this SPDE

indexs <- inla.spde.make.index("s", spde1$n.spde)

#------------------------------------------------------------------------------#
# Generate a projection matrix A - projects the spatially continuous GRF from
# the observations to the mesh nodes

# For training points
coo <- st_coordinates(dat)
A <- inla.spde.make.A(mesh = mesh, loc = coo)

dim(A)
nrow(dat)

# ---------------------------------------------------------------------------- #
# Define model matrix based on all covariates of interest, removing automatic 
# intercept

# Covariates of interest
covs <- c("age_s","sex","hiv","prv_tx",
          "marg_caste","occupation",
          "block_endm_2017", "IRS_2017_1","vill_inc_2017_gt0", 
          "traveltime_s", "rain",
          "detection")

X <- model.matrix(as.formula(paste("~ ",paste(covs, collapse = " + "))), 
                  data = dat)[,-1] 

# ---------------------------------------------------------------------------- #
# Make stack

stk <- inla.stack(
  data = list(y = dat$delay),
  A = list(A, 1, 1),
  effects = list(s = indexs,  # the spatial index,
                 v = dat$v,
                 data.frame(
                   Intercept = 1,
                   X)
  )
)

# saveRDS(stk, here::here("data/analysis","stack.rds"))

# ---------------------------------------------------------------------------- #
# Fit baseline IID model 

fit.iid <- init_inla(f = y ~ -1 + Intercept + f(v, model = 'iid',
                                                 prior = 'pc.prec', 
                                                 param = c(1, 0.01)),
                     data.stack = stk,
                     family = "nbinomial")

# Check spatial correlation in fitted IID effects
plot_vgm(values = fit.iid$fit$summary.random$v$mean,
         loc = distinct(dplyr::select(dat,v, geometry)))

# ---------------------------------------------------------------------------- #
# Fit baseline models with SPDE and IID to explain spatial structure
# Compare priors for SPDE

# Define formulae

f.list <- list(
  f1 = y ~ -1 + Intercept + f(v, model = 'iid',
                                 prior = 'pc.prec', 
                                 param = c(1, 0.01)) +
                                 f(s, model = spde1),
  f2 = y ~ -1 + Intercept + f(v, model = 'iid',
                                 prior = 'pc.prec', 
                                 param = c(1, 0.01)) +
                                 f(s, model = spde2),
  f3 = y ~ -1 + Intercept + f(v, model = 'iid',
                                  prior = 'pc.prec', 
                                  param = c(1, 0.01)) +
                                f(s, model = spde3),
  f4 = y ~ -1 + Intercept + f(v, model = 'iid',
                                  prior = 'pc.prec', 
                                  param = c(1, 0.01)) +
                                f(s, model = spde4)
)

# Fit models
fits <- lapply(f.list, init_inla, data.stack = stk, family = "nbinomial")

plyr::llply(fits, function(x) x$fit$dic$dic)
# $f1
# [1] 38552.04
# 
# $f2
# [1] 38551.69
# 
# $f3
# [1] 38548.86
# 
# $f4
# [1] 38549.47

plyr::llply(fits, function(x) summary(x$fit))

plyr::llply(fits, function(x) x$fit$summary.hyperpar[,c(1,3,5)])

# Project fitted SPDEs
pdf(here::here(figdir, "compare_spde_priors.pdf"), height = 7, width = 10) 
purrr::imap(fits, function(x, nm) plot_spde(x$fit, nm, limit1 = c(-1,1), limit2 = c(0,0.5))) 
dev.off()

names(fits) <- c("R5e3-Sig1.5", 
                 "R5e5-Sig1.5", 
                 "R5e3-Sig1.1", 
                 "R5e3-Sig2")
p.list <- purrr::imap(fits, function(x, nm) plot_spde(x$fit, nm, limit1 = c(-1,1), limit2 = c(0,0.5))) 

png(here::here(figdir, "compare_spde_priors.png"), height = 4000, width = 8000, res = 350) 
gridExtra::grid.arrange(grobs = p.list) #plyr::llply(fits, function(x) plot_spde(x$fit))
dev.off()

saveRDS(fits, here::here(outdir, "fits_compare_priors.rds"))

saveRDS(spde1, here::here("data/analysis","spde.rds"))

################################################################################
################################################################################