################################################################################
# Description: 
# 
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit"
outdir <- "output"

covs_pat <- c("age_s","comorb1", "poss_acdTRUE") #"sexFemale","prv_txYes","marg_casteYes","occupationUnskilled","occupationSkilled", "`occupationSalaried.selfemployed`",
covs_vil_aware <- c("block_endm_2017TRUE", "inc_2017_gt0TRUE") # "IRS_2017_1Yes",
covs_vil_access <- c("traveltime_t_s")

covs.list <- list(None = NULL, 
                  `Patient only` = covs_pat,
                  `Awareness only` = covs_vil_aware,
                  `Access only` = covs_vil_access,
                  `Patient + village awareness` = c(covs_pat,covs_vil_aware),
                  `Patient + village access` = c(covs_pat,covs_vil_access),
                  `Village awareness + access` = c(covs_vil_aware,covs_vil_access),
                  All = c(covs_pat, covs_vil_aware, covs_vil_access),
                  Interaction = c(covs_pat, covs_vil_aware, covs_vil_access, "poss_acdTRUE:block_endm_2017TRUE")) 

dat <- read_data()
mesh <- readRDS(here::here("data/analysis","mesh.rds"))

boundary <- readRDS(here::here("data","geography","boundary.rds"))
boundary.spdf <- as_Spatial(boundary)

#------------------------------------------------------------------------------#
# SPDE prior specification and index set

# + 10km is the minimum edge-length of the mesh so shorter range unlikely
# + Assuming most delays between 0-40 days, SD unlikely to be greater than 10 
#   + However SD on internal scale
#   + Looser upper bound on SD preferably to overly stringent

spde <- inla.spde2.pcmatern(mesh = mesh, 
                            prior.range = c(10, 0.01), # P(range < U) = a
                            prior.sigma = c(2, 0.05), # P(sigma > U) = a; on the log scale
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
covs <- c("age_s","comorb","poss_acd", "block_endm_2017","inc_2017_gt0","traveltime_t_s")

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

#------------------------------------------------------------------------------#
# Initialise each model with all data

fit_covs <- function(covs.list) {
  
# Define formula
  f <- as.formula(paste0("y ~ -1 + Intercept +", 
                         paste0(covs.list, collapse = " + "), 
                         "+ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(1, 0.01)) +
                            f(s, model = spde)"))

  print(f)
  
# Fit full models 
  fit <- inla(f,
              family = "poisson",
              data = inla.stack.data(stk),
              control.predictor = list(
                compute = TRUE, link = 1,
                A = inla.stack.A(stk)),
              control.compute = list(waic = TRUE, 
                                     config = TRUE,
                                     cpo = FALSE,
                                     return.marginals.predictor = TRUE),
              control.fixed = list(mean = 0, 
                                   prec = 0.1, 
                                   mean.intercept = 0, 
                                   prec.intercept = 0.1),
              verbose = TRUE)
  
  res <- list(f = f, fit = fit)
  return(res)
  
}

fits.main <- plyr::llply(covs.list, fit_covs)
# saveRDS(fits.main, here::here(outdir, "fits_covs_main.rds"))

fit_int <- fit_covs(covs.list$Interaction)

#------------------------------------------------------------------------------#
# Refit full model without OLRE

f.spde <- as.formula(paste0("y ~ -1 + Intercept +", 
                       paste0(covs.list$All, collapse = " + "), 
                       "+ f(s, model = spde)"))

refit.spde <- inla(f.spde,
                  family = "poisson",
                  data = inla.stack.data(stk),
                  control.predictor = list(
                    compute = TRUE, link = 1,
                    A = inla.stack.A(stk)),
                  control.compute = list(waic = TRUE,
                                         config = TRUE,
                                         cpo = FALSE,
                                         return.marginals.predictor = TRUE),
                  control.mode = list(result = fits.main$All, restart = TRUE),
                  control.fixed = list(mean = 0,
                                       prec = 0.1,
                                       mean.intercept = 0,
                                       prec.intercept = 0.1),
                  verbose = TRUE)

refit.spde <- list(f = f.spde, fit = refit.spde)

# saveRDS(refit.spde, here::here(outdir, "fit_covs_spde.rds"))

#------------------------------------------------------------------------------#
# Refit full model without SPDE

f.iid <- as.formula(paste0("y ~ -1 + Intercept +", 
                           paste0(covs.list$All, collapse = " + "), 
                           "+ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(1, 0.01))"))

refit.iid<- inla(f.iid,
                 family = "poisson",
                 data = inla.stack.data(stk),
                 control.predictor = list(
                   compute = TRUE, link = 1,
                   A = inla.stack.A(stk)),
                 control.compute = list(waic = TRUE,
                                        config = TRUE,
                                        cpo = FALSE,
                                        return.marginals.predictor = TRUE),
                 control.mode = list(result = fits.main$All, restart = TRUE),
                 control.fixed = list(mean = 0,
                                      prec = 0.1,
                                      mean.intercept = 0,
                                      prec.intercept = 0.1),
                 verbose = TRUE)

refit.iid <- list(f = f.iid, fit = refit.iid)

# saveRDS(refit.iid, here::here(outdir, "fit_covs_iid.rds"))

#------------------------------------------------------------------------------#
# Refit full model without any random effects

f.fixed <- as.formula(paste0("y ~ -1 + Intercept +", 
                           paste0(covs.list$All, collapse = " + ")))

refit.fixed <- inla(f.fixed,
                 family = "poisson",
                 data = inla.stack.data(stk),
                 control.predictor = list(
                   compute = TRUE, link = 1,
                   A = inla.stack.A(stk)),
                 control.compute = list(waic = TRUE,
                                        config = TRUE,
                                        cpo = FALSE,
                                        return.marginals.predictor = TRUE),
                 # control.mode = list(result = fits.main$All, restart = TRUE),
                 control.fixed = list(mean = 0,
                                      prec = 0.1,
                                      mean.intercept = 0,
                                      prec.intercept = 0.1),
                 verbose = TRUE)

refit.fixed <- list(f = f.fixed, fit = refit.fixed)

# saveRDS(refit.iid, here::here(outdir, "fit_covs_iid.rds"))

#------------------------------------------------------------------------------#
# Refit full model as binomial for delay > 30 days

f.all <- as.formula(paste0("y ~ -1 + Intercept +", 
                           paste0(covs.list$All, collapse = " + "), 
                           "+ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(1, 0.01)) +
                            f(s, model = spde)"))

refit.bin <- inla(f.all,
                 family = "binomial",
                 data = inla.stack.data(stk_bin),
                 control.predictor = list(
                   compute = TRUE, link = 1,
                   A = inla.stack.A(stk_bin)),
                 control.compute = list(waic = TRUE,
                                        config = TRUE,
                                        cpo = FALSE,
                                        return.marginals.predictor = TRUE),
                 control.fixed = list(mean = 0,
                                      prec = 0.1,
                                      mean.intercept = 0,
                                      prec.intercept = 0.1),
                 verbose = TRUE)
refit.bin <- list(f = f.all, fit = refit.bin)

# saveRDS(refit.iid, here::here(outdir, "fit_full_bin.rds"))

#------------------------------------------------------------------------------#
# Refit with SPDE replicated by ACD

# indexs <- inla.spde.make.index(name = "s",
#                                   n.spde = spde$n.spde,
#                                   n.repl = 2)
# 
# m <- nrow(dat)
# A = inla.spde.make.A(mesh,
#                      loc = st_coordinates(dat), 
#                      index = rep(1:m, times = 2),
#                      repl = rep(1:2, each = m))
# 
# # Covariates of interest
# covs <- c("age_s","comorb","poss_acd",
#           "block_endm_2017", "inc_2017_gt0", 
#           "traveltime_t_s")
# 
# X <- model.matrix(as.formula(paste("~ ",paste(covs, collapse = " + "))), 
#                   data = dat)[,-1] 
# 
# stk_rep <- inla.stack(
#   data = list(y = dat$delay),
#   A = list(A, 1, 1),
#   effects = list(s = indexs,  # the spatial index,
#                  id = dat$id, # observation level index
#                  data.frame(
#                    Intercept = 1,
#                    X) # covariate model matrix
#   )
# )
# 
# saveRDS(stk, here::here("data/analysis","stack_rep.rds"))
# 
# f.rep <- as.formula(paste0("y ~ -1 + Intercept +", 
#                            paste0(covs.list$All, collapse = " + "), 
#                            "+ f(id, model = 'iid',
#                              prior = 'pc.prec', 
#                              param = c(10, 0.01)) +
#                             f(s, model = spde, replicate = s.repl)"))
# 
# refit.rep <- inla(f.rep,
#                   family = "poisson",
#                   data = inla.stack.data(stk),
#                   control.predictor = list(
#                     compute = TRUE, link = 1,
#                     A = inla.stack.A(stk)),
#                   control.compute = list(waic = TRUE,
#                                          config = TRUE,
#                                          cpo = FALSE,
#                                          return.marginals.predictor = FALSE),
#                   control.fixed = list(mean = 0,
#                                        prec = 0.1,
#                                        mean.intercept = 0,
#                                        prec.intercept = 0.1),
#                   verbose = TRUE)
# refit.rep <- list(f = f.rep, fit = refit.rep)
# 
# saveRDS(refit.rep, here::here(outdir, "fit_full_repACD.rds"))

#------------------------------------------------------------------------------#
# Final model list

fits.all <- append(fits.main,
                   list(`Fixed effects only` = refit.fixed,
                        `All (IID only)` = refit.iid,
                        `All (SPDE only)` = refit.spde,
                        `All (Binomial)` = refit.bin))
saveRDS(fits.all, here::here(outdir, "fits_final.rds"))

################################################################################
################################################################################

