################################################################################
# Define data stacks
################################################################################

dat <- readRDS(here::here("data/analysis","dat_nona.rds")) %>%
  sf::st_transform(7759)
spde <- readRDS(here::here("data/analysis","spde.rds"))
mesh <- readRDS(here::here("data/analysis","mesh.rds"))
coop <- readRDS(here::here("data/analysis","coop.rds"))

# Covariates of interest
covs <- c("age_s","sex","hiv",
          "marg_caste","occupation",
          "block_endm_2017", "IRS_2017_1","vill_inc_2017_gt0", "prv_tx",
          "traveltime_s", "rain",
          "detection")

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
make_spatial_partition <- function(data, r = 1e5) {
  
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

make_village_partition <- function(data, p = 0.1) {
  
  # Randomly sample (100*p)% of villages from data
  n_samp <- round(n_distinct(data$v)*p)
  print(n_samp)
  samp <- sample(unique(data$v), n_samp)
  
  # Exclude all observations in those villages
  train <- data[!data$v %in% samp,]
  test <- data[data$v %in% samp,]
  
  return(list(train = train, test = test))
  
}

M = 10
partitions <- lapply(1:M, make_village_partition, data = dat, p = 0.1)
saveRDS(partitions, here::here("data/analysis",paste0("village_partition_",M,".rds")))

plot_partition <- function(part) {
  ggplot() +
    gg(as_Spatial(boundary)) +
    gg(as_Spatial(part$train)) +
    gg(as_Spatial(part$test), col = "red") -> p
  
  print(p)
}

pdf(here::here("figures","fit","spatial_partition.pdf"))
lapply(partitions, plot_partition)
dev.off()

# ---------------------------------------------------------------------------- #

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
  
  # Define model matrix based on all covariates of interest, removing automatic 
  # intercept
  X1 <- model.matrix(as.formula(paste("~ ",paste(covs, collapse = " + "))), 
                     data = dat.train)[,-1] 
 
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
  
  # Prediction points
  Ap <- inla.spde.make.A(mesh = mesh, loc = coop)
  
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
  
  if (nrow(dat.test) > 0){
    # Testing
    coot <- st_coordinates(dat.test)
    At <- inla.spde.make.A(mesh = mesh, loc = coot)
    
    dim(At)
    nrow(dat.test)
    
    X2 <- model.matrix(as.formula(paste("~ ",paste(covs, collapse = " + "))), 
                       data = dat.test)[,-1] 
    
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

  # Combine stacks
  stk.full <- inla.stack(stk.train, stk.test, stk.pred)

  }else{  
    stk.full <- inla.stack(stk.train, stk.pred)
    
  }
  return(stk.full)
}

dat.full <- list(train = dat, test = NULL)
stk.full <- make_stacks(dat.full)

################################################################################
################################################################################