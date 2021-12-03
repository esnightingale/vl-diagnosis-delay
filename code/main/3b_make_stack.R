################################################################################
# Define data stacks
################################################################################

dat <- readRDS(here::here("data/analysis","dat_nona.rds")) %>%
  dplyr::filter(delay >= 0) %>%
  sf::st_transform(7759)
spde <- readRDS(here::here("data/analysis","spde.rds"))
mesh <- readRDS(here::here("data/analysis","mesh.rds"))
coop <- readRDS(here::here("data/analysis","coop.rds"))

# Setup map context
blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
  sf::st_transform(7759)

boundary <- blockmap %>%
  sf::st_union()

# Covariates of interest
covs <- c("age_s","sex","hiv",
          "marg_caste","occupation",
          "block_endm_2017", "IRS_2017","inc_2017_gt0", "prv_tx",
          "traveltime_s", "traveltime_t_s", "rain",
          "detection")

#------------------------------------------------------------------------------#
# Generate the index set for this SPDE

indexs <- inla.spde.make.index("s", spde$n.spde)
lengths(indexs)

#------------------------------------------------------------------------------#
# Split training and testing data

# # Spatial partition
# make_spatial_partition <- function(data, r = 1e4) { # sample 1 village, default 10km radius
#   
#   # villages <- distinct(dplyr::select(data, v, geometry))
#   
#   # Randomly sample n.start observations
#   samp <- sample_n(data, 1)
#   # print(paste0("Sampled village ID: ",samp$v))
#   
#   # Define buffer of radius r around sampled observation(s)
#   buff <- st_buffer(samp, dist = r)
#   
#   # Intersect full dataset with buffer
#   p.index <- st_intersects(data, buff, sparse = FALSE)
#   print(summary(p.index))
#   
#   return(p.index)
#   
#   # train <- data[!index,]
#   # test <- data[index,]
#   # 
#   # return(list(train = train, test = test))
#   
# }
# 
# 
# # make_village_partition <- function(data, p = NULL, n = NULL) {
# # 
# #   # Randomly sample fixed number or proportion (100*p)% of villages from data
# #   if (is.null(n)){
# #     n_samp <- round(n_distinct(data$v)*p)
# #     print(n_samp)
# #   }elseif (!is.null(p)){ n_samp <- n 
# #   }else {print("Provide n or p")}
# # 
# #   samp <- sample(unique(data$v), n_samp)
# #   # print(length(samp))
# # 
# #   # Exclude all observations in those villages
# #   train <- dat[!dat$v %in% samp,]
# #   test <- dat[dat$v %in% samp,]
# # 
# #   return(list(train = train, test = test))
# # 
# # }
# 
# M = 10
# partitions <- lapply(1:M, function(x) make_spatial_partition(data = dat, r = 5e3))
# # partitions <- lapply(1:M, make_village_partition, data = dat, n = 10) 
# saveRDS(partitions, here::here("data/analysis",paste0("village_radius_partition_",M,".rds")))
# 
# plot_partition <- function(data, p.index) {
#   
#   train <- data[!p.index,]
#   test <- data[p.index,]
#   
#   ggplot() +
#     gg(as_Spatial(boundary)) +
#     gg(as_Spatial(train), cex = 1, pch = 1) +
#     gg(as_Spatial(test), col = "red", cex = 1) -> p
#   
#   print(p)
# }
# 
# pdf(here::here("figures","fit","cross-validation","spatial_partition.pdf"),
#     height = 9, width = 12)
# lapply(partitions, function(p.index) plot_partition(dat, p.index))
# dev.off()

# ---------------------------------------------------------------------------- #

# Define function to make stacks from given train/test partition index
make_stacks <- function(data, p.index) {

  dat.train <- data[!p.index,]
  dat.test <- data[p.index,]
  
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
    data = list(y = dat.train$delay),
                A = list(A, 1, 1),
                effects = list(s = indexs,  # the spatial index,
                               v = dat.train$v,
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
                effects = list(s = indexs,
                               data.frame(
                                 Intercept = rep(1, nrow(coop)))
                )
  )
  
  if (nrow(dat.test) > 0){
    print(nrow(dat.test))
    
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
      data = list(y = NA),
                  A = list(At, 1, 1),
                  effects = list(s = indexs,  # the spatial index,
                                 v = dat.test$v,
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
  return(list(stack = stk.full, p.index = p.index))
}

# stk.partitions <- lapply(partitions, function(p.index) make_stacks(data = dat, p.index = p.index))

# saveRDS(stk.partitions, 
#         here::here("data/analysis","stack_partitions.rds"))

# Setup full stack with no testing data
stk.full <- make_stacks(dat, p.index = rep(FALSE, nrow(dat)))
saveRDS(stk.full, here::here("data/analysis","stack.rds"))

################################################################################
################################################################################