################################################################################
# Description: Preliminary variable selection on patient characteristics, with
# non-spatial IID effects. 
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit/univariate/poisson"
outdir <- "output/univariate/poisson"

covs_all <- list("age_s","sexMale","hivYes","prv_txYes","marg_casteYes",c("occupationUnskilled",
              "occupationSkilled", "`occupationSalaried.selfemployed`"), "detectionACD",
             "block_endm_2017Endemic", "IRS_2017Yes","inc_2017_gt0Yes",
             "traveltime_s","traveltime_t_s", "rain")

covs2 <- list("age_s","sexMale","hivYes","prv_txYes","marg_casteYes","occupationUnskilled",
              "occupationSkilled", "`occupationSalaried.selfemployed`") 
covs3 <- list("age_s","sexMale","hivYes","prv_txYes","marg_casteYes","occupationUnskilled",
              "occupationSkilled", "`occupationSalaried.selfemployed`", "detectionACD")
covs_full <- list("age_s","sexMale","hivYes","prv_txYes","marg_casteYes","occupationUnskilled",
              "occupationSkilled", "`occupationSalaried.selfemployed`", "detectionACD",
             "block_endm_2017Endemic", "IRS_2017Yes","inc_2017_gt0Yes",
             "traveltime_s","rain") #"traveltime_t_s", 

covs_select <- list("age_s","hivYes", "detectionACD",
              "block_endm_2017Endemic","inc_2017_gt0Yes","traveltime_s","rain") 

mesh <- readRDS(here::here("data/analysis","mesh.rds")) 
spde <- readRDS(here::here("data/analysis","spde.rds"))
stk <- readRDS(here::here("data/analysis","stack.rds"))
idx <- stk$data$index$train

# Original sf object for plotting vgm
dat <- readRDS(here::here("data/analysis","dat_nona.rds")) %>%
  filter(delay >= 0) %>%
  st_transform(crs = st_crs(7759))

# ---------------------------------------------------------------------------- #
# Univariate fits

family <- "poisson"

fit_covs_uni <- function(cov) {
  
  # Define formula
  f <- as.formula(paste0("y ~ -1 + Intercept +", 
                         paste0(cov, collapse = " + "), 
                         "+ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(10, 0.01))"))
  
  print(f)
  
  # Fit model
  fit <- init_inla(f, data.stack = stk, family = family, cpo = F)$fit
  
  res <- list(f = f, fit = fit)
  return(res)
  
}

fits.uni <- plyr::llply(covs, fit_covs_uni)
saveRDS(fits.uni, here::here(outdir, "fits_covs_univar_nonspatial.rds"))

# ---------------------------------------------------------------------------- #
# Patient-only

f <- as.formula(paste0("y ~ ", 
                       paste0(covs2, collapse = " + "), 
                       "+ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(10, 0.01))"))

fit.pat <- init_inla(f, data.stack = stk, family = "poisson", cpo = F)

# saveRDS(fit.pat, here::here(outdir, "fit_covs_pat_nonspatial.rds"))

p_vgm_pat <- plot_vgm(fit.pat$fit$summary.random$id$mean, dat, title = "Fitted IID effects: Patient-only model")
#   model     psill    range kappa
# 1   Nug 0.6696714     0.00   0.0
# 2   Mat 0.2294127 28905.71   0.5

# ---------------------------------------------------------------------------- #
# Patient + detection
f <- as.formula(paste0("y ~ ", 
                       paste0(covs3, collapse = " + "), 
                       "+ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(10, 0.01))"))

fit.pat2 <- init_inla(f, data.stack = stk, family = "poisson", cpo = F)

# saveRDS(fit.pat2, here::here(outdir, "fit_covs_pat_nonspatial.rds"))

p_vgm_pat2 <- plot_vgm(fit.pat2$fit$summary.random$id$mean, dat, title = "Fitted IID effects: Patient + detection model")
#   model    psill    range kappa
# 1   Nug 0.672013     0.00   0.0
# 2   Mat 0.211144 29801.58   0.5

# ---------------------------------------------------------------------------- #
# Full multivariate fit

f <- as.formula(paste0("y ~ ", 
                       paste0(covs_full, collapse = " + "), 
                       "+ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(10, 0.01))"))

fit.multi <- init_inla(f, data.stack = stk, family = "poisson", cpo = T)

# saveRDS(fit.multi, here::here(outdir, "fit_covs_multivar_nonspatial.rds"))

p_vgm_full <- plot_vgm(fit.multi$fit$summary.random$id$mean, dat, title = "Fitted IID effects: Full model")
# model     psill    range kappa
# 1   Nug 0.6771872     0.00   0.0
# 2   Mat 0.1932114 27455.23   0.5

# ---------------------------------------------------------------------------- #
# Compare univariate versus multivariate coeffcients

fits.all <- rlist::list.append(fits.uni, fit.pat, fit.pat2, fit.multi)

c("Intercept",
  "Age (std)",
  "Sex: Male",
  "HIV positive",
  "Prv. VL/PKDL",
  "Marginalised caste",
  "Occupation: Unskilled",
  "Skilled",
  "Salaried/self-employed",
  "Detection: ACD",
  "Block endemic (2017)",
  "Village IRS targeting (2017)",
  "Village incidence > 0 (2017)",
  "Travel time - diagnosis facility (std)",
  "Travel time - treatment facility (std)",
  "Diagnosis season: rain",
  "Intercept"
) -> axis.labs

Efxplot(plyr::llply(fits.all, function(x) x$fit), 
        VarNames = rev(axis.labs),
        Intercept = FALSE, Size = 1.5) +
  scale_colour_manual(values = c(rep("indianred",length(fits.uni)),
                                 "goldenrod",
                                 "forestgreen",
                                 "steelblue")) +
  theme(legend.position = "none") + 
  scale_y_continuous(limits = c(-0.4, 0.6)) +
  labs(title = "Estimated covariate effects: non-spatial models")

# Look at spatial correlation after including patient characteristics and detection
# png(here::here(figdir, "vgm_iid_multi_select_nonspatial.png"), height = 500, width = 1000)
gridExtra::grid.arrange(p_vgm_pat, p_vgm_pat2, p_vgm_full, nrow = 1)
# dev.off()

# ---------------------------------------------------------------------------- #
# Multivariate fit with selection

f <- as.formula(paste0("y ~ ", 
                       paste0(covs_select, collapse = " + "), 
                       "+ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(10, 0.01))"))

fit.multi2 <- init_inla(f, data.stack = stk, family = "poisson", cpo = F)

saveRDS(fit.multi2, here::here(outdir, "fit_covs_multivar_select_nonspatial.rds"))

p_vgm2 <- plot_vgm(fit.multi2$fit$summary.random$id$mean, dat, title = "Fitted IID effects: Full model with selection")

png(here::here(figdir, "vgm_iid_multi_select_nonspatial.png"), height = 500, width = 1000)
gridExtra::grid.arrange(p_vgm, p_vgm2, nrow = 1)
dev.off()

# ---------------------------------------------------------------------------- #
# Compare fitted regression estimates

fits.all <- rlist::list.append(fits.uni, fit.multi, fit.multi2)
saveRDS(fits.all, here::here(outdir, "fits_covs_all_nonspatial.rds"))

c("Intercept",
  "Age (std)",
  "Sex: Male",
  "HIV positive",
  "Prv. VL/PKDL",
  "Marginalised caste",
  "Occupation: Unskilled",
  "Skilled",
  "Salaried/self-employed",
  "Detection: ACD",
  "Block endemic (2017)",
  "Village IRS targeting (2017)",
  "Village incidence > 0 (2017)",
  "Travel time - diagnosis facility (std)",
  "Travel time - treatment facility (std)",
  "Diagnosis season: rain",
  "Intercept"
  ) -> axis.labs

Efxplot(plyr::llply(fits.all, function(x) x$fit), 
                   VarNames = rev(axis.labs),
                   Intercept = FALSE, Size = 1.5) +
  scale_y_continuous(limits = c(-0.4, 0.6)) +
  scale_colour_manual(values = c(rep("indianred",length(fits.uni)),
                                 "steelblue",
                                 "forestgreen")) +
  theme(legend.position = "none") + 
  labs(title = "Estimated covariate effects: non-spatial models",
       caption = "Univariate estimate in red, multivariate in blue and multivariate with selected covariates only in green.")

ggsave(here::here(figdir, "covs_uni_multi_efx_nonspatial.png"), height = 6, width = 8, units = "in", dpi = 300)

################################################################################
################################################################################
