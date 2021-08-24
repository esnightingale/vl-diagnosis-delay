################################################################################
# Define data stacks
################################################################################

# Organise the data, projection matrices and fixed/random effects into a "stack".
# Need separate ones for estimation at observed locations and prediction across
# defined grid

# List A contains 1 to indicate fixed effects are mapped one to one, and a projection
# matrix for the random effects

covs_pat <- c("age_s","sex","hiv","marg_caste","occ4_cat","prv_tx")
covs_choice <- c("detection","conslt_gt1", "traveltime_s")
covs_vil <- c("block_endm_2017", "IRS_2017_1","vill_inc_2017_gt0")

covs_pat_f <- paste(covs_pat, collapse = " + ")
covs_vil_f <- paste(covs_vil, collapse = " + ")

# Rename levels of num_conslt to avoid special character errors in later formulae
# dat.fit$num_conslt_cat <- factor(dat.fit$num_conslt_cat, labels = c("0_1","1_3","3_5","5_8"))

# Define model matrix based on selected covariates from previous baseline/IID modelling, removing automatic intercept
X <- model.matrix(as.formula(paste("~ ",covs_pat_f," + ", covs_vil_f)), data = dat.train)[,-1] 
Xt <- model.matrix(as.formula(paste("~ ",covs_pat_f," + ", covs_vil_f)), data = dat.test)[,-1] 
# Xp <- model.matrix(~ traveltime_s, data = dat.pred)[,-1] 


# Make the stack
stk.train <- inla.stack(
  tag = "train",
  data = list(y = dat.train$days_fever),
  A = list(A, 1, 1),
  effects = list(s = indexs,  #the spatial index,
                 id = dat.train$id, # IID index
                 data.frame(Intercept = 1, # Covariates
                            X))
)

# stack for validation at witheld test points stk.v
stk.test <- inla.stack(
  tag = "test",
  data = list(y = NA),
  A = list(At, 1, 1),
  effects = list(s = indexs, 
                 id = dat.test$id, 
                 data.frame(Intercept = 1, 
                            Xt))
)

# stack for prediction stk.p - only predicting from intercept and fitted spatial field
# stk.p <- inla.stack(
#   tag = "pred",
#   data = list(y = NA),
#   A = list(Ap, 1),
#   effects = list(s = indexs,
#                  data.frame(Intercept = rep(1, nrow(coop)),
#                             X = Xp))
# )

# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.train, stk.test)


# # Stack for spatio-temporal 
# 
# # Specifying a new set of SPDE components ####
# 
# Groups = "q"
# 
# NGroups <- length(unique(dat.fit[,Groups])) 
# 
# A2 <- inla.spde.make.A(mesh, # Leave
#                        loc = coo, # Leave
#                        group = as.numeric(as.factor(dat.fit.df[,Groups])), # this must be a numeric value counting from 1. If the groups variable is a factor, this will happen by default.
#                        n.group = NGroups) 
# 
# w.Host2 <- inla.spde.make.index(
#   name    = 'w', 
#   n.spde  = Hosts.spde$n.spde,
#   n.group = NGroups)  
# 
# StackHost2 <- inla.stack( 
#   data = list(y = TestHosts[,resp]), # Leave
#   
#   A = list(1, 1, 1, HostA2), # Change the A matrix to the new one
#   
#   effects = list(
#     Intercept = rep(1, N), # Leave
#     X = X, # Leave
#     ID = TestHosts$ID, # Leave
#     
#     w = w.Host2)) # CHANGE
# 
# # Stack for estimation/validation/prediction
#  
# # Prediction stack has the response set as NA
# 
# # stack for estimation stk.e
# stk.e <- inla.stack(
#   tag = "est", 
#   data = list(y = dat.fit$days_fever),
#   A = list(A, 1, 1),
#   effects = list(s = indexs, i = dat.fit.df$i, 
#                  list(data.frame(b0 = 1, 
#                                  # patient chars
#                                  age = dat.fit.df$age,
#                                  sex = dat.fit.df$sex,
#                                  hiv = dat.fit.df$hiv,
#                                  marg_caste = dat.fit.df$marg_caste,
#                                  detection = dat.fit.df$detection,
#                                  prv_tx_ka = dat.fit.df$prv_tx_ka,
#                                  # village chars
#                                  num_conslt_cat = dat.fit.df$num_conslt_cat,
#                                  vill_inc_2017_gt0 = dat.fit.df$vill_inc_2017_gt0,
#                                  IRS_2017_1 = dat.fit.df$IRS_2017_1,
#                                  block_endm_2017 = dat.fit.df$block_endm_2017,
#                                  traveltime = dat.fit.df$traveltime)))
# )
# 
# 
# # stack for validation at witheld test points stk.v
# stk.v <- inla.stack(
#   tag = "val",
#   data = list(y = dat.val.df$days_fever),
#   A = list(Av, 1, 1),
#   effects = list(s = indexs, i = dat.val.df$i, 
#                  list(data.frame(b0 = 1, 
#                                  # patient chars
#                                  age = dat.val.df$age,
#                                  sex = dat.val.df$sex,
#                                  hiv = dat.val.df$hiv,
#                                  marg_caste = dat.val.df$marg_caste,
#                                  detection = dat.val.df$detection,
#                                  prv_tx_ka = dat.val.df$prv_tx_ka,
#                                  # village chars
#                                  num_conslt_cat = dat.val.df$num_conslt_cat,
#                                  vill_inc_2017_gt0 = dat.val.df$vill_inc_2017_gt0,
#                                  IRS_2017_1 = dat.val.df$IRS_2017_1,
#                                  block_endm_2017 = dat.val.df$block_endm_2017,
#                                  traveltime = dat.val.df$traveltime)))
# )
# 
# # stack for prediction stk.p
# stk.p <- inla.stack(
#   tag = "pred",
#   data = list(y = NA),
#   A = list(1, Ap),
#   effects = list(data.frame(b0 = rep(1, nrow(coop))), s = indexs)
# )
# 
# # stk.full has stk.e and stk.p
# stk.full <- inla.stack(stk.e, stk.v, stk.p)
