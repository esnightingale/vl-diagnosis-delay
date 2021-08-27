################################################################################
# Define data stacks
################################################################################

dat <- readRDS(here::here("data/analysis","dat_nona.rds")) 
spde <- readRDS(here::here("data/analysis","spde.rds"))
mesh <- readRDS(here::here("data/analysis","mesh.rds"))
coop <- readRDS(here::here("data/analysis","coop.rds"))

# Covariates of interest
covs <- c("age_s","sex","hiv",
          "marg_caste","occ4_cat",
          "block_endm_2017", "IRS_2017_1","vill_inc_2017_gt0", "prv_tx",
          "traveltime_s", "rain",
          "detection","consult_gt1")

#------------------------------------------------------------------------------#
# Generate the index set for this SPDE

indexs <- inla.spde.make.index("s", spde$n.spde)
lengths(indexs)

#------------------------------------------------------------------------------#
# Split training and testing data

# Random partition
# test.idx <- sample(1:nrow(dat), floor(nrow(dat)*test.size))
# 
# dat.train <- dat[-test.idx,]
# dat.test <- dat[test.idx,]

# Spatial partition
make_partition <- function(data, r = 0.1) {
  
  # Randomly sample one observation from data
  samp <- sample_n(data, 1)
  
  # Define buffer of radius r around observation
  buff <- st_buffer(samp, dist = r)
  
  # Intersect full dataset with buffer
  index <- st_intersects(data, buff, sparse = FALSE)
  
  train <- data[-index,]
  test <- data[index,]
  
  return(list(train = train, test = test))
  
}

M = 10
partitions <- lapply(1:M, make_partition, data = dat)
saveRDS(partitions, here::here("data/analysis","partitions.rds"))

plot_partition <- function(part) {
  ggplot() +
    gg(as_Spatial(boundary)) +
    gg(as_Spatial(part$train)) +
    gg(as_Spatial(part$test), col = "red")
}
lapply(partitions, plot_partition)

# Define function to make stacks from given train/test partition
make_stacks <- function(dat.partition) {

  dat.train <- dat.partition$train
  dat.test <- dat.partition$test
  
#------------------------------------------------------------------------------#
# Generate a projection matrices A - projects the spatially continuous GRF from
# the observations to the mesh nodes

  # For training points
  coo <- st_coordinates(dat.train)
  A <- inla.spde.make.A(mesh = mesh, loc = coo)
  
  dim(A)
  nrow(dat.train)
 
  # Testing
  coot <- st_coordinates(dat.test)
  At <- inla.spde.make.A(mesh = mesh, loc = coot)
  
  dim(At)
  nrow(dat.test)

  # Prediction
  Ap <- inla.spde.make.A(mesh = mesh, loc = coop)

#------------------------------------------------------------------------------#
# Organise the data, projection matrices and fixed/random effects into stacks

  # Define model matrix based on all covariates of interest, removing automatic 
  # intercept
  X1 <- model.matrix(as.formula(paste("~ ",paste(covs, collapse = " + "))), 
                    data = dat.train)[,-1] 
  X2 <- model.matrix(as.formula(paste("~ ",paste(covs, collapse = " + "))), 
                    data = dat.test)[,-1] 
  
  # Training stack
  stk.train <- inla.stack(
    tag = "train",
    data = list(y = dat.train$days_fever),
    A = list(A, 1),
    effects = list(v = indexv,  # the spatial index,
                   data.frame(  # covariates
                     Intercept = 1, 
                     X1)
                   )
  )
  
  # Testing stack
  stk.test <- inla.stack(
    tag = "test",
    data = list(y = dat.test$days_fever),
    A = list(A, 1),
    effects = list(v = indexv,  # the spatial index,
                   data.frame(  # covariates
                     Intercept = 1, 
                     X2)
    )
  )
  
  # Stack for smooth prediction from intercept and fitted spatial field (no covariates)
  stk.pred <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(Ap, 1),
    effects = list(v = indexv,
                   data.frame(
                     Intercept = rep(1, nrow(coop)))
                   )
  )
  
  # Combine stacks
  stk.full <- inla.stack(stk.train, stk.test, stk.pred)

  return(stk.full)
}


################################################################################
################################################################################