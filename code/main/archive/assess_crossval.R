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

figdir <- "figures/fit/cross-validation/covs_all"
outdir <- "output/cross-validation/covs_all"

covs_pat <- c("age_s","sex","hiv")
covs_ses <- c("marg_caste","occ4_cat")
covs_choice <- c("detection","consult_gt1","prv_tx")
covs_phys <- c("traveltime_s","rain")
covs_vilctl <- c("block_endm_2017", "IRS_2017_1","vill_inc_2017_gt0")

covs_all <- c(covs_pat, covs_ses, covs_choice, covs_phys, covs_vilctl)

dat.fit <- readRDS(here::here("data/analysis","dat_fit.rds"))

# Read cross-validation outputs for three model structures (10 reps each)
fits.xval <- lapply(c("IID","SPDE","IID_SPDE"), 
                    function(model) readRDS(here::here(outdir, paste0("fits_xval_",model,".rds"))))

#------------------------------------------------------------------------------#
# Compare regression coefficients between reps/models

fits.xval.all <- unlist(fits.xval, recursive = FALSE)

gen_Efxplot <- function(fits){
  ggregplot::Efxplot(
    plyr::llply(fits, 
                function(res) return(res$fit)),
    Intercept = FALSE) +
    guides(colour = "none") -> p
  return(p)
}

p.list <- plyr::llply(fits.xval, gen_Efxplot)
gridExtra::grid.arrange(grobs = p.list)

ggsave(here::here(figdir, "fits_xval_efx.png"), height = 7, width = 9, units = "in")

#---------------------------------------------#
# Correlation of fitted regression parameters #
#---------------------------------------------#

get_coefs <- function(res){
  coefs <- data.frame(Covariate = res$fit$names.fixed, Estimate = res$fit$summary.fixed$mean)
  return(coefs)
}

coef_reps <- bind_rows(lapply(fits.xval.both, get_coefs), .id = "Replicate") %>%
  pivot_wider(names_from = "Covariate", values_from = "Estimate") %>%
  column_to_rownames("Replicate")

coef_corr <- cor(as.matrix(coef_reps))

corrplot::corrplot(coef_corr, method = "color", type = "lower") #order = "hclust"

#------------------------------------------------------------------------------#
# Check prediction at test observations

# Extract fitted marginals for sampled test observations defined by test.index
# Calculate probability of exceeding 30 days from these marginal distributions
# Compare to observed delay from full dataset, dat.fit

calc_excprob <- function(res, exc = 30){
  
  idx <- res$test.idx
  test.marg <- res$fit$marginals.fitted.values[idx]
  test.fit <- res$fit$summary.fitted.values$mean[idx]
  
  dat.fit %>% 
    slice(idx) %>%
    mutate(exceed = (days_fever > exc),
           excprob = sapply(test.marg,
                            FUN = function(marg){1 - inla.pmarginal(q = exc, marginal = marg)}),
           sqerr = (days_fever - test.fit)^2) %>%
    return()
}

out.xval <- bind_rows(
  plyr::llply(fits.xval,
              function(res){
                bind_rows(
                  plyr::llply(res, 
                              calc_excprob, 
                              exc = 30),
                  .id = "Model")
              }),
  .id = "Replicate")

RMSE <- 
  RMSE

png(here::here(figdir, "xval_exc30_val.png"), height = 500, width = 800)
create_raincloud(exceed30, xvar = "exceed", yvar = "excprob", 
                 xlab = paste("Observed delay > 30 days"),
                 ylab = paste("Fitted probability of delay > 30 days"), 
                 col_by = "exceed", 
                 drop_na = TRUE) 
# + facet_wrap(~Replicate) 
dev.off()

# ---------------------------------------------------------------------------- #
# Testing for residual spatial correlation after accounting for patient and 
# village covariates

# dat.fit$Residuals <- inlatools::residuals(final.iid)[-v.idx]
# dat.fit$IID <- final.iid$summary.random$i$mean[-v.idx]
# # dat.fit <- dat.fit %>% mutate(Total = Residuals + IID)
# 
# vg <- variogram(Residuals~1, data = dat.fit)
# print(plot(vg))
# 
# vgmod <- vgm(psill =  0.01, model = "Mat", nugget = 0.06, range = 100)
# plot(vg, model = vgmod)
# 
# vgfit <- fit.variogram(vg, model = vgmod)
# plot(vg, model = vgfit, main = "Residuals") -> p_res
# 
# print(vgfit)
# #   model      psill    range kappa
# # 1   Nug 0.06119009  0.00000   0.0
# # 2   Mat 0.01280572 20.14779   0.5
# p_res
# 
# vg <- variogram(IID~1, data = dat.fit)
# print(plot(vg))
# 
# vgmod <- vgm(psill =  0.01, model = "Mat", nugget = 0.06, range = 100)
# plot(vg, model = vgmod)
# 
# vgfit <- fit.variogram(vg, model = vgmod)
# plot(vg, model = vgfit, main = "IID") -> p_iid
# 
# print(vgfit)
# #   model      psill    range kappa
# # 1   Nug 0.20746644  0.00000   0.0
# # 2   Mat 0.07844214 38.15515   0.5
# p_iid
# 
# png(here::here(figdir, "/fit/iid_varselect_semivariogram.png"), height = 500, width = 1000)
# gridExtra::grid.arrange(p_res, p_iid, nrow = 1)
# dev.off()

################################################################################
################################################################################