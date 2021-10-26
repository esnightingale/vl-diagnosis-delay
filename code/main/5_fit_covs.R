################################################################################
# Description: 
# 
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit"
outdir <- "output"

covs_pat <- c("age_s","sex","hiv","prv_tx")
covs_ses <- c("marg_caste","occupation")
# covs_pat_aware <- c("prv_tx") #,"consult_gt0"
covs_vil_aware <- c("block_endm_2017", "IRS_2017_1","vill_inc_2017_gt0")
covs_access <- c("traveltime_s","rain")
covs_ctl <- "detection"
# covs_choice <- "consult_gt0"

covs.list <- list(None = NULL, 
                  Health = covs_pat, 
                  SES = covs_ses, 
                  Awareness = covs_pat_aware,
                  # Choice = covs_choice,
                  Access = covs_access,
                  Control = covs_ctl,
                  Patient = c(covs_pat, covs_ses, "detection"),
                  Village = covs_vil_aware,
                  All = c(covs_pat, covs_ses, "detection", covs_vil_aware, covs_access)) #covs_choice, covs_phys

dat <- readRDS(here::here("data/analysis","dat_fit_val.rds"))
mesh <- readRDS(here::here("data/analysis","mesh.rds"))
spde <- readRDS(here::here("data/analysis","spde.rds"))

blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
  st_transform(7759)

boundary <- sf::st_union(blockmap)
boundary.spdf <- as_Spatial(boundary)

#------------------------------------------------------------------------------#
# Initialise each model with all training data

fit_covs <- function(covs.list) {
  
# Define formula
  f <- as.formula(paste0("days_fever ~ ", 
                         paste0(covs.list, collapse = " + "), 
                         "+ f(v, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(1, 0.01)) 
                          + f(v2, model = spde)"))

# Fit full models 
  fit <- inla(f,
              family = "nbinomial",
              data = dat,
              control.predictor = list(
                compute = TRUE, link = 1),
              control.compute = list(dic = TRUE, 
                                     waic = TRUE, 
                                     config = TRUE),
              control.fixed = list(mean = 0, 
                                   prec = 0.1, 
                                   mean.intercept = 0, 
                                   prec.intercept = 0.1),
              verbose = TRUE)
  
  res <- list(f = f, fit = fit)
  return(res)
  
}

covs.list <- covs.list[c("None","Patient","Village","Access","All")]
fits.init <- plyr::llply(covs.list, fit_covs)
saveRDS(fits.init, here::here(outdir, "fits_covs_init.rds"))

plyr::llply(fits.init, function(x) summary(x$fit))

# Compare fitted regression estimates
ggregplot::Efxplot(plyr::llply(fits.init, function(x) x$fit),
                   ModelNames = names(covs.list),
                   Intercept = FALSE)
ggsave(here::here(figdir, "fits_init_efx.png"), height = 6, width = 8, units = "in", dpi = 320)

#------------------------------------------------------------------------------#
# Plot the fitted spatial fields

pdf(here::here(figdir, "map_fitted_spdes.pdf"), height = 4, width = 15)
plyr::llply(fits.init, function(x) print(plot_spde(x)))
dev.off()

#------------------------------------------------------------------------------#

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

