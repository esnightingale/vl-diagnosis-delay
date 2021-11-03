################################################################################
# Description: Preliminary variable selection on patient characteristics, with
# non-spatial IID effects. 
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit/exploratory"
outdir <- "output/exploratory"

# covs_pat <- list("age_s","sex","hiv","prv_tx","marg_caste","occupation", "detection")
covs_pat <- list("age_s","sexFemale","hivYes","prv_txYes","marg_casteYes",c("occupationUnskilled",
              "occupationSkilled", "`occupationSalaried.self.employed`"), "detectionACD", "rain")
covs_pat2 <- list("age_s","sexFemale","hivYes","prv_txYes","marg_casteYes","occupationUnskilled",
                  "occupationSkilled", "`occupationSalaried.selfemployed`", "detectionACD", "rain")

dat <- readRDS(here::here("data/analysis","dat_nona.rds"))
stk <- readRDS(here::here("data/analysis","stack.rds"))

# ---------------------------------------------------------------------------- #
# Univariate fits

fit_covs_uni <- function(cov) {
  
  # Define formula
  f <- as.formula(paste0("y ~ ", 
                         paste0(cov, collapse = " + "),
                         "+ f(v, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(1, 0.01)) +
                            f(s, model = spde)"))
  
  # Fit model
  fit <- init_inla(f, data.stack = stk, family = "nbinomial", cpo = F)
  
  return(fit)
  
}

fits.uni <- plyr::llply(covs_pat, fit_covs_uni)
saveRDS(fits.uni, here::here(outdir, "fits_pat_univar.rds"))

# ---------------------------------------------------------------------------- #
# Multivariate fit

f <- as.formula(paste0("y ~ ", 
                       paste0(covs_pat2, collapse = " + "), 
                       "+ f(v, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(1, 0.01)) +
                            f(s, model = spde)"))

fit.multi <- init_inla(f, data.stack = stk, family = "nbinomial", cpo = F)
saveRDS(fit.multi, here::here(outdir, "fit_pat_multivar.rds"))

# ---------------------------------------------------------------------------- #
# Compare fitted regression estimates

fits.all <- rlist::list.append(fits.uni, fit.multi)

ggregplot::Efxplot(plyr::llply(fits.all, function(x) x$fit), 
                   Intercept = FALSE) +
  scale_colour_manual(values = c(rep("indianred",length(fits.uni)),
                                 "steelblue")) +
  theme(legend.position = "none")

ggsave(here::here(figdir, "pat_uni_multi_efx.png"), height = 6, width = 8, units = "in", dpi = 300)

# ---------------------------------------------------------------------------- #
# Patient covariate selection with IID 

pat.select <- ggregplot::INLAModelSel("days_fever", 
                                      c("age_s","sex","hiv","prv_tx","marg_caste","occupation", "detection","rain"), 
                                      "v", "iid", 
                                      "nbinomial", 
                                      dat)
saveRDS(pat.select, file = here::here(outdir, "pat_varselect_iid.rds"))

summary(pat.select$FinalModel)
# Fixed effects:
#   mean    sd 0.025quant 0.5quant 0.975quant   mode kld
# (Intercept)                      3.810 0.016      3.778    3.810      3.842  3.810   0
# age_s                            0.075 0.010      0.056    0.075      0.094  0.075   0
# hivYes                           0.233 0.051      0.133    0.233      0.334  0.233   0
# occupationUnskilled             -0.036 0.022     -0.079   -0.036      0.008 -0.036   0
# occupationSkilled                0.040 0.039     -0.036    0.040      0.116  0.040   0
# occupationSalaried.selfemployed -0.015 0.038     -0.090   -0.015      0.060 -0.015   0
# detectionACD                    -0.180 0.020     -0.218   -0.180     -0.141 -0.180   0
# 
# Random effects:
#   Name	  Model
# v IID model
# 
# Model hyperparameters:
#   mean    sd 0.025quant 0.5quant 0.975quant mode
# size for the nbinomial observations (1/overdispersion) 4.36 0.131       4.10     4.35       4.62 4.35
# Precision for v                                        7.18 0.468       6.32     7.16       8.16 7.11
# 
# Expected number of effective parameters(stdev): 1006.77(42.47)
# Number of equivalent replicates : 4.26 
# 
# Deviance Information Criterion (DIC) ...............: 38516.63
# Deviance Information Criterion (DIC, saturated) ....: 5377.06
# Effective number of parameters .....................: 1018.10
# 
# Marginal log-Likelihood:  -19522.03 

# Compare fitted regression estimates
ggregplot::Efxplot(pat.select$FinalModel,
                   Intercept = FALSE)
ggsave(here::here(figdir, "pat_varselect_efx.png"), height = 6, width = 8, units = "in", dpi = 320)


# Assess spatial correlation in fitted IID effects
plot_vgm(values = pat.select$FinalModel$summary.random$v$mean,
         loc = distinct(dplyr::select(dat,v, geometry)))

covs.final <- pat.select$Removed[[length(pat.select$Removed)]]
# [1] "age_s"      "hiv"        "occupation" "detection" 

# f.select <- as.formula(paste0("y", " ~ ", 
#                               paste(covs.final, collapse = " + "), 
#                               "+ f(v, model = 'iid',
#                              prior = 'pc.prec', 
#                              param = c(1, 0.01))")) 
# 
# fit.select <- inla(f.select,
#                    family = "nbinomial",
#                    data = dat,
#                    control.predictor = list(
#                      compute = TRUE, link = 1),
#                    control.compute = list(dic = TRUE, 
#                                           waic = TRUE, 
#                                           config = TRUE),
#                    control.fixed = list(mean = 0, 
#                                         prec = 0.1, 
#                                         mean.intercept = 0, 
#                                         prec.intercept = 0.1),
#                    verbose = TRUE) 
# 
# summary(fit.select)

################################################################################
################################################################################
