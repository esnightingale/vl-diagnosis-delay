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

figdir <- "figures/fit/cross-validation/no_covs"
outdir <- "output/cross-validation/no_covs"

covs_pat <- c("age_s","sex","hiv")
covs_ses <- c("marg_caste","occ4_cat")
covs_aware <- c("detection","consult_gt1","prv_tx")
covs_phys <- c("traveltime_s","rain")
covs_vilctl <- c("block_endm_2017", "IRS_2017_1","vill_inc_2017_gt0")

covs_all <- c(covs_pat, covs_ses, covs_aware, covs_phys, covs_vilctl)

dat.fit <- readRDS(here::here("data/analysis","dat_fit.rds"))
spde <- readRDS(here::here("data/analysis","spde.rds"))

#------------------------------------------------------------------------------#
# Initialise IID and SPDE models with all data

covs.list <- NULL # covs_all
M <- 3

# run_crossval <- function(covs.list, outdir, M = 10) {
  
# Define formulae
  f1 <- as.formula(paste0("days_fever ~ ", paste0(covs.list, collapse = " + "), " +  f(id, model = 'iid')"))
  f2 <- as.formula(paste0("days_fever ~ ", paste0(covs.list, collapse = " + "), " + f(v, model = spde)"))
  f3 <- as.formula(paste0("days_fever ~ ", paste0(covs.list, collapse = " + "), " + f(id, model = 'iid') + f(v, model = spde)"))

# Fit full models 
  fits.init <- plyr::llply(list(IID = f1, 
                                # SPDE = f2, 
                                IID_SPDE = f3), 
                           init_inla, 
                           family = "poisson",
                           data = dat.fit)
  saveRDS(fits.init, here::here(outdir, "fits_init.rds"))

  # lapply(fits.init, plyr::llply(fits.init, function(x) summary(x$res)))
  plyr::llply(fits.init, function(x) x$fit$dic$dic)
  # $IID
  # [1] 23643.69
  # 
  # $IID_SPDE
  # [1] 23643.97
  
# Compare fitted regression estimates

  ggregplot::Efxplot(plyr::llply(fits.init, function(x) x$fit),
                     ModelNames = c("IID",
                                    # "SPDE",
                                    "Both"),
                     Intercept = FALSE)
  ggsave(here::here(figdir, "fits_init_efx.png"), height = 6, width = 8, units = "in", dpi = 320)

#----------------#
# Spatial LOO CV #
#----------------#
  
INLAutils::inlasloo(st_drop_geometry(dat.fit),
                    long = "longitude",
                    lat = "latitude",
                    y = "days_fever",
                    family = "poisson",
                    ss = 4,
                    rad = 0.1,
                    modform = f3,
                    mesh = mesh,
                    print = TRUE,
                    plot = TRUE)

#-------------------------------------#
# Refit with subsampled training data #
#-------------------------------------#

  run_crossval <- function(fit.init, M){
    res <- lapply(1:M, 
                  function(x) {
                    print(x)
                    train_inla(formula = fit.init$f, fit.init = fit.init$fit, test.percent = 0.2)
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

################################################################################
################################################################################

