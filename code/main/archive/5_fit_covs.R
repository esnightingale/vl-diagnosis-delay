################################################################################
# Description: 
# 
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit/poisson"
outdir <- "output/poisson"

covs_pat <- c("age_s","hivYes", "detectionACD") #"sexFemale","prv_txYes","marg_casteYes","occupationUnskilled","occupationSkilled", "`occupationSalaried.selfemployed`",
covs_vil_aware <- c("block_endm_2017Endemic", "inc_2017_gt0Yes") # "IRS_2017_1Yes",
covs_vil_access <- c("traveltime_t_s")

covs.list <- list(None = NULL, 
                  `Patient only` = covs_pat,
                  `Awareness only` = covs_vil_aware,
                  `Access only` = covs_vil_access,
                  `Patient + village awareness` = c(covs_pat,covs_vil_aware),
                  `Patient + village access` = c(covs_pat,covs_vil_access),
                  `Village awareness + access` = c(covs_vil_aware,covs_vil_access),
                  All = c(covs_pat, covs_vil_aware, covs_vil_access)) 

# dat <- readRDS(here::here("data/analysis","dat_nona.rds")) %>%
#   st_transform(crs = st_crs(7759))
mesh <- readRDS(here::here("data/analysis","mesh.rds"))
spde <- readRDS(here::here("data/analysis","spde.rds"))
stk <- readRDS(here::here("data/analysis","stack.rds"))
stk_bin <- readRDS(here::here("data/analysis","stack_bin.rds"))

boundary <- readRDS(here::here("data","geography","boundary.rds"))
boundary.spdf <- as_Spatial(boundary)

#------------------------------------------------------------------------------#
# Initialise each model with all data

fit_covs <- function(covs.list) {
  
# Define formula
  f <- as.formula(paste0("y ~ -1 + Intercept +", 
                         paste0(covs.list, collapse = " + "), 
                         "+ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(10, 0.01)) +
                            f(s, model = spde)"))

  print(f)
  
# Fit full models 
  fit <- inla(f,
              family = "poisson",
              data = inla.stack.data(stk),
              control.predictor = list(
                compute = TRUE, link = 1,
                A = inla.stack.A(stk)),
              control.compute = list(dic = TRUE, 
                                     waic = TRUE, 
                                     config = TRUE,
                                     cpo = TRUE,
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
saveRDS(fits.main, here::here(outdir, "fits_covs_main.rds"))

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
                  control.compute = list(dic = TRUE,
                                         waic = TRUE,
                                         config = TRUE,
                                         cpo = TRUE,
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
                             param = c(10, 0.01))"))

refit.iid<- inla(f.iid,
                 family = "poisson",
                 data = inla.stack.data(stk),
                 control.predictor = list(
                   compute = TRUE, link = 1,
                   A = inla.stack.A(stk)),
                 control.compute = list(dic = TRUE,
                                        waic = TRUE,
                                        config = TRUE,
                                        cpo = TRUE,
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
# Refit full model as binomial for delay > 30 days

f.all <- as.formula(paste0("y ~ -1 + Intercept +", 
                           paste0(covs.list$All, collapse = " + "), 
                           "+ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(10, 0.01)) +
                            f(s, model = spde)"))

refit.bin <- inla(f.all,
                 family = "binomial",
                 data = inla.stack.data(stk_bin),
                 control.predictor = list(
                   compute = TRUE, link = 1,
                   A = inla.stack.A(stk_bin)),
                 control.compute = list(dic = TRUE,
                                        waic = TRUE,
                                        config = TRUE,
                                        cpo = TRUE,
                                        return.marginals.predictor = TRUE),
                 control.fixed = list(mean = 0,
                                      prec = 0.1,
                                      mean.intercept = 0,
                                      prec.intercept = 0.1),
                 verbose = TRUE)
refit.bin <- list(f = f.all, fit = refit.bin)

# saveRDS(refit.iid, here::here(outdir, "fit_full_bin.rds"))

#------------------------------------------------------------------------------#
# Final model list

fits.all <- append(fits.main, 
                   list(`All (IID only)` = refit.iid,
                        `All (SPDE only)` = refit.spde,
                        `All (Binomial)` = refit.bin))
saveRDS(fits.all, here::here(outdir, "fits_final.rds"))

################################################################################
################################################################################

