################################################################################
# Description: Setup SPDE prior and data stack
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
# SPDE prior specification and index set

# + 10km is the minimum edge-length of the mesh so shorter range unlikely
# + Assuming most delays between 0-40 days, SD unlikely to be greater than 10 
#   + Looser upper bound on SD preferably to overly stringent

spde <- inla.spde2.pcmatern(mesh = mesh, 
                            prior.range = c(10, 0.01), # P(range < U) = a
                            prior.sigma = c(2, 0.1), # P(sigma > U) = a; on the log scale
                            constr = TRUE)

indexs <- inla.spde.make.index("s", spde$n.spde)

#------------------------------------------------------------------------------#
# Generate a projection matrix A - projects the spatially continuous GRF from
# the observations to the mesh nodes

# For training points
coo <- st_coordinates(dat)
A <- inla.spde.make.A(mesh = mesh, loc = coo)

dim(A)
# 4271 3674
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
  data = list(y = dat$delay),
  A = list(A, 1, 1),
  effects = list(s = indexs,  # the spatial index,
                 id = dat$id, # observation level index
                 data.frame(
                   Intercept = 1,
                   X) # covariate model matrix
  )
)

saveRDS(stk, here::here("data/analysis","stack.rds"))

# Also make stack for binomial output for later sensitivity analysis
stk_bin <- inla.stack(
  data = list(y = as.numeric(dat$gt30)),
  A = list(A, 1, 1),
  effects = list(s = indexs,  # the spatial index,
                 id = dat$id, # observation level index
                 data.frame(
                   Intercept = 1,
                   X) # covariate model matrix
  )
)

saveRDS(stk_bin, here::here("data/analysis","stack_bin.rds"))

################################################################################
################################################################################