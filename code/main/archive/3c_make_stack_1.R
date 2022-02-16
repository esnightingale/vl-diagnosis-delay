################################################################################
# Define data stacks
################################################################################

dat <- readRDS(here::here("data/analysis","dat_nona.rds")) %>%
  st_transform(7759) %>%
  filter(delay >= 0) %>%
  mutate(delay_gt30 = as.numeric(delay > 30))

spde <- readRDS(here::here("data/analysis","spde.rds"))
mesh <- readRDS(here::here("data/analysis","mesh.rds"))
coop <- readRDS(here::here("data/analysis","coop.rds"))

# Covariates of interest
covs <- c("age_s","sex","hiv","prv_tx",
          "marg_caste","occupation",
          "block_endm_2017", "IRS_2017","inc_2017_gt0", 
          "traveltime_s", "traveltime_t_s", "rain",
          "detection")

#------------------------------------------------------------------------------#
# Generate the index set for this SPDE

indexs <- inla.spde.make.index("s", spde$n.spde)
lengths(indexs)

#------------------------------------------------------------------------------#
# Generate a projection matrix A - projects the spatially continuous GRF from
# the observations to the mesh nodes

# For training points
coo <- st_coordinates(dat)
A <- inla.spde.make.A(mesh = mesh, loc = coo)

dim(A)
nrow(dat)

# Define model matrix based on all covariates of interest, removing automatic 
# intercept
X <- model.matrix(as.formula(paste("~ ",paste(covs, collapse = " + "))), 
                  data = dat)[,-1] 

# Training stack
stk.train <- inla.stack(
  tag = "train",
  data = list(y = dat$delay), delay_gt30
  A = list(A, 1, 1, 1),
  effects = list(s = indexs,  # the spatial index,
                 v = dat$v, # village index
                 id = dat$id, # individual index
                 data.frame(  # covariates
                   Intercept = 1, 
                   X)
  )
)

# Prediction points
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)

# Stack for smooth prediction from intercept and fitted spatial field (no covariates)
stk.pred <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(Ap, 1),
  effects = list(s = indexs,
                 data.frame(
                   Intercept = rep(1, nrow(coop)))
  )
)
 
stk.full <- inla.stack(stk.train, stk.pred)

saveRDS(stk.full, here::here("data/analysis","stack.rds"))

################################################################################
################################################################################