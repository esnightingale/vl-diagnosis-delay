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

#------------------------------------------------------------------------------#
# Define formulae

spde <- inla.spde2.pcmatern(mesh = mesh, 
                            prior.range = c(10, 0.01), # P(range < U) = a
                            prior.sigma = c(1, 0.01), # P(sigma > U) = a
                            constr = TRUE)

f <- "y ~ -1 + Intercept + f(s, model = spde)"
    
indexs <- inla.spde.make.index("s", spde$n.spde)

#------------------------------------------------------------------------------#
# Generate a projection matrix A

coo <- st_coordinates(dat)
A <- inla.spde.make.A(mesh = mesh, loc = coo)

dim(A)
# 4271 3574
nrow(dat)

# ---------------------------------------------------------------------------- #
# Define model matrix based on all covariates of interest, removing automatic 
# intercept

# Covariates of interest
covs <- c("age_s","sex","comorb","prv_tx",
          "caste4_r","occ4_cat","poss_acd",
          "block_endm_2017", "IRS_2017","inc_2017_gt0",
          "traveltime_s", "traveltime_t_s", "rain")

X <- model.matrix(as.formula(paste("~ ",paste(covs, collapse = " + "))),
                  data = dat)[,-1]

# ---------------------------------------------------------------------------- #
# Make stack

stk <- inla.stack(
  data = list(y = as.numeric(dat$poss_acd)),
  A = list(A, 1),
  effects = list(s = indexs,  # the spatial index,
                 # id = dat$id, # observation level index
                 data.frame(
                   Intercept = rep(1, nrow(dat)))
  )
)

# ---------------------------------------------------------------------------- #

fit.acd <- inla(f,
                 family = "binomial",
                 Ntrials = 1,
                 data = inla.stack.data(stk),
                 control.predictor = list(
                         compute = TRUE, link = 1,
                         A = inla.stack.A(stk)),
                 control.compute = list(dic = FALSE, 
                                              waic = TRUE, 
                                              config = TRUE),
                 verbose = TRUE)

fit.acd$summary.hyperpar[,c(1,3,5)]

# Project fitted SPDE
png(here::here(figdir, "spde_acd.png"), height = 6, width = 9, units = "in") 
plot_spde_mean(fit.acd, nm, limits = c(-1,1))
dev.off()

################################################################################
################################################################################