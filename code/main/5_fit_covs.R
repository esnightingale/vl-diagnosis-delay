################################################################################
# Description: 
# 
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit"
outdir <- "output"

covs_pat <- c("age_s","sexFemale","hivYes","prv_txYes","marg_casteYes","occupationUnskilled",
              "occupationSkilled", "`occupationSalaried.selfemployed`", "detectionACD","rain")
covs_vil_aware <- c("block_endm_2017Endemic", "IRS_2017_1Yes","vill_inc_2017_gt0Yes")
covs_vil_access <- c("traveltime_s*rain")

covs.list <- list(None = NULL, 
                  `Patient only` = covs_pat,
                  `Awareness only` = covs_vil_aware,
                  `Access only` = covs_vil_access,
                  `Patient + village awareness` = c(covs_pat,covs_vil_aware),
                  `Patient + village access` = c(covs_pat[-10],covs_vil_access),
                  `Village awareness + access` = c(covs_vil_aware,covs_vil_access),
                  All = c(covs_pat[-10],covs_vil_aware, covs_vil_access)) 

dat <- readRDS(here::here("data/analysis","dat_nona.rds")) %>%
  st_transform(crs = st_crs(7759))
mesh <- readRDS(here::here("data/analysis","mesh.rds"))
spde <- readRDS(here::here("data/analysis","spde.rds"))
stk <- readRDS(here::here("data/analysis","stack_incl_lt14.rds"))[[1]]

blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
  st_transform(7759)

boundary <- sf::st_union(blockmap)
boundary.spdf <- as_Spatial(boundary)

#------------------------------------------------------------------------------#
# Initialise each model with all data

fit_covs <- function(covs.list) {
  
# Define formula
  f <- as.formula(paste0("y ~ -1 + Intercept +", 
                         paste0(covs.list, collapse = " + "), 
                         "+ f(v, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(1, 0.01)) +
                            f(s, model = spde)"))

  print(f)
  
# Fit full models 
  fit <- inla(f,
              family = "nbinomial",
              data = inla.stack.data(stk),
              control.predictor = list(
                compute = TRUE, link = 1,
                A = inla.stack.A(stk)),
              control.compute = list(dic = TRUE, 
                                     waic = TRUE, 
                                     config = TRUE,
                                     cpo = TRUE),
              control.fixed = list(mean = 0, 
                                   prec = 0.1, 
                                   mean.intercept = 0, 
                                   prec.intercept = 0.1),
              verbose = TRUE)
  
  res <- list(f = f, fit = fit)
  return(res)
  
}

fits.init <- plyr::llply(covs.list, fit_covs)
# saveRDS(fits.init, here::here(outdir, "fits_covs_init.rds"))

#------------------------------------------------------------------------------#
# Refit full model without IID

f.spde <- as.formula(paste0("y ~ -1 + Intercept +", 
                       paste0(covs.list$All, collapse = " + "), 
                       "+ f(s, model = spde)"))

refit.spde <- inla(f.spde,
                  family = "nbinomial",
                  data = inla.stack.data(stk),
                  control.predictor = list(
                    compute = TRUE, link = 1,
                    A = inla.stack.A(stk)),
                  control.compute = list(dic = TRUE,
                                         waic = TRUE,
                                         config = TRUE,
                                         cpo = TRUE),
                  control.mode = list(result = fits.init$All, restart = TRUE),
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
                           "+ f(v, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(1, 0.01))"))

refit.iid<- inla(f.iid,
                 family = "nbinomial",
                 data = inla.stack.data(stk),
                 control.predictor = list(
                   compute = TRUE, link = 1,
                   A = inla.stack.A(stk)),
                 control.compute = list(dic = TRUE,
                                        waic = TRUE,
                                        config = TRUE,
                                        cpo = TRUE),
                 control.mode = list(result = fits.init$All, restart = TRUE),
                 control.fixed = list(mean = 0,
                                      prec = 0.1,
                                      mean.intercept = 0,
                                      prec.intercept = 0.1),
                 verbose = TRUE)

refit.iid <- list(f = f.iid, fit = refit.iid)

# saveRDS(refit.iid, here::here(outdir, "fit_covs_iid.rds"))

#------------------------------------------------------------------------------#
# Final model list

fits.all <- append(fits.init, 
                   list(`All (IID only)` = refit.iid,
                        `All (SPDE only)` = refit.spde))
saveRDS(fits.all, here::here(outdir, "fits_covs_all.rds"))

################################################################################
################################################################################

