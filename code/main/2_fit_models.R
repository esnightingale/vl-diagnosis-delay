################################################################################
# Description: 
# 
# Initialise IID and SPDE models using all fitting data
# 
# Refit using subsampled training data and assess prediction for remaining test
# data.
# 
# Want models with all covariates and with covariates in each of four domains
# 
# Compare domains, select strongest predictors from each
# - Do I want this? Want to interpret each covariate or just domain contribution?
# 
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit"
outdir <- "output"

covs_pat <- c("age_s","sex","hiv")
covs_ses <- c("marg_caste","occ4_cat")
covs_aware <- c("block_endm_2017", "IRS_2017_1","vill_inc_2017_gt0","prv_tx")
covs_phys <- c("traveltime_s","rain")
covs_ctl <- "detection"
covs_choice <- "consult_gt1"

covs.list <- list(None = NULL, 
                  Health = covs_pat, 
                  SES = covs_ses, 
                  Awareness = covs_aware,
                  Choice = covs_choice,
                  Physical = covs_phys,
                  Control = covs_ctl,
                  All = c(covs_pat, covs_ses, covs_aware, covs_choice, covs_phys, covs_ctl))

dat.fit <- readRDS(here::here("data/analysis","dat_fit.rds"))
mesh <- readRDS(here::here("data/analysis","mesh.rds"))
spde <- readRDS(here::here("data/analysis","spde.rds"))

#------------------------------------------------------------------------------#
# Initialise each model with all training data

fit_covs <- function(covs.list) {
  
# Define formula
  f <- as.formula(paste0("days_fever ~ ", paste0(covs.list, collapse = " + "), "+ f(v, model = spde)"))

# Fit full models 
  fit <- init_inla(f, 
                   family = "poisson",
                   data = dat.fit)
  
  res <- list(f = f, fit = fit)
  return(res)
  
}

fits.init <- plyr::llply(covs.list, fit_covs)
saveRDS(fits.init, here::here(outdir, "/base/fits_init.rds"))

# Compare fitted regression estimates
ggregplot::Efxplot(plyr::llply(fits.init, function(x) x$fit$fit),
                   ModelNames = names(covs.list),
                   Intercept = FALSE)
ggsave(here::here(figdir, "fits_init_efx.png"), height = 6, width = 8, units = "in", dpi = 320)


#----------------#
# Spatial LOO CV #
#----------------#

out.field <- INLA::inla.spde2.result(fits.init$None$fit$fit,'v', spde, do.transf = TRUE)
range.out <- INLA::inla.emarginal(function(x) x, out.field$marginals.range.nominal[[1]])
# 0.6886794
  
run_sloo <- function(res.init) {
  
  sloo <- INLAutils::inlasloo(st_drop_geometry(dat.fit),
                              long = "longitude",
                              lat = "latitude",
                              y = "days_fever",
                              family = "poisson",
                              ss = 10,
                              rad = 0.1,
                              modform = res.init$f,
                              mesh = mesh,
                              print = TRUE,
                              plot = TRUE,
                              control.mode = list(result = res.init$fit, restart = TRUE))
  
  return(sloo)
  
}

fits.sloo <- plyr::llply(covs.list, run_sloo)

saveRDS(here::here(outdir, "cross-validation/fits_sloo.rds"))

for (m in seq_along(fits.init)) {
  
  sloo <- run_sloo(fits.init[[m]])
  saveRDS(here::here(outdir, paste0("cross-validation/sloo_",names(fits.init)[m],".rds")))
  
}

#--------------------------------------#
# Refit with sub-sampled training data #
#--------------------------------------#

## Randomly resample for each iteration/model

# Starting from initial fits using all data, randomly sub-sample a test set and 
# refit with these outcome values set to NA
run_crossval <- function(res.init, M = 10){
  res <- lapply(1:M, 
                function(x) {
                  print(x)
                  train_inla(formula = res.init$f, 
                             fit.init = res.init$fit, 
                             test.percent = 0.2)
                  }
                )
  return(res)
}
  
  # In list, tends to crash (parallellise?)
  # fits.xval <- plyr::llply(fits.init, run_crossval, M = M)
  # saveRDS(fits.xval, here::here(outdir, "fits_xval.rds"))
  
  # In loop, save as go along
  for (model in c("IID",
                  # "SPDE",
                  "IID_SPDE")) {
    
    saveRDS(run_crossval(fits.init[[model]], M = M), 
            here::here(outdir, 
                       paste0("fits_xval_",model,".rds")))
  }

## Pre-define M resampled training/test sets and refit all models with these

test.percent = 0.2
lapply(1:M, function(x) sample(1:nrow(dat.fit), floor(nrow(dat.fit)*test.percent)))
test.idx <- sample(1:nrow(dat.fit), floor(nrow(dat.fit)*test.percent))

dat.train <- dat.fit
dat.train$days_fever[test.idx] <- NA

#------------------------------------------------------------------------------#
# Spatial LOOCV

spatial.field <- INLA::inla.spde2.result(fit.base.spde$fit,'v', spde, do.transf = TRUE)
range.spatial <- INLA::inla.emarginal(function(x) x, spatial.field$marginals.range.nominal[[1]])

f1 <- days_fever ~ 1 + f(v, model = 'iid') 
f2 <- days_fever ~ 1 + f(v, model = spde)

INLAutils::inlasloo(st_drop_geometry(dat.fit),
                    long = "longitude",
                    lat = "latitude",
                    y = "days_fever",
                    family = "poisson",
                    ss = 10,
                    rad = range.spatial*0.15,
                    modform = list(f1,f2),
                    mesh = mesh,
                    mae = TRUE,
                    ds = TRUE,
                    print = TRUE,
                    plot = TRUE)


################################################################################
################################################################################

