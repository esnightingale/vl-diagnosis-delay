################################################################################
# Description: Preliminary variable selection on patient characteristics, with
# non-spatial IID effects. 
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit/univariate"
outdir <- "output/univariate"

covs <- list("age_s","sexMale","hivYes","prv_txYes","marg_casteYes",c("occupationUnskilled",
              "occupationSkilled", "`occupationSalaried.selfemployed`"), "detectionACD",
             "block_endm_2017Endemic", "IRS_2017Yes","inc_2017_gt0Yes",
             "traveltime_s","traveltime_t_s", "rain")
covs2 <- list("age_s","sexMale","hivYes","prv_txYes","marg_casteYes","occupationUnskilled",
              "occupationSkilled", "`occupationSalaried.selfemployed`", "detectionACD",
             "block_endm_2017Endemic", "IRS_2017Yes","inc_2017_gt0Yes",
             "traveltime_s","rain") #"traveltime_t_s", 
covs3 <- list("age_s","hivYes", "detectionACD",
              "block_endm_2017Endemic","inc_2017_gt0Yes") 

spde <- readRDS(here::here("data/analysis","spde.rds"))
stk <- readRDS(here::here("data/analysis","stack.rds"))

# ---------------------------------------------------------------------------- #
# Univariate fits

fit_covs_uni <- function(cov) {
  
  # Define formula
  f <- as.formula(paste0("y ~ -1 + Intercept +", 
                         paste0(cov, collapse = " + "), 
                         "+ f(v, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(1, 0.01)) +
                            f(s, model = spde)"))
  
  print(f)
  
  # Fit model
  fit <- init_inla(f, data.stack = stk, family = "nbinomial", cpo = F)$fit
  
  res <- list(f = f, fit = fit)
  return(res)
  
}

fits.uni <- plyr::llply(covs, fit_covs_uni)
saveRDS(fits.uni, here::here(outdir, "fits_covs_univar.rds"))

ggregplot::Efxplot(plyr::llply(fits.uni, function(x) x$fit), 
                   Intercept = FALSE) +
  theme(legend.position = "none")

# ---------------------------------------------------------------------------- #
# Full multivariate fit

f <- as.formula(paste0("y ~ ", 
                       paste0(covs2, collapse = " + "), 
                       "+ f(v, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(1, 0.01)) +
                            f(s, model = spde)"))

fit.multi <- init_inla(f, data.stack = stk, family = "nbinomial", cpo = T)

saveRDS(fit.multi, here::here(outdir, "fit_covs_multivar.rds"))

plot_obsvpred(fit.multi$fit, stk)

# p_le30 <- pnbinom(q = 30, prob = fit.nb2$fit$summary.fitted.values$mean[idx], size = 1)
# 
# ggplot(data.frame(le30 = (stk$data$data$y[idx] <= 30),
#                   p = p_le30),
#        aes(le30, p)) +
#   geom_boxplot()
# 
# predy <- sapply(fit.nb2$fit$summary.fitted.values$mean[idx], function(p) rnbinom(n = 100, prob = p, size = 1))
# 
# 
# preds <- data.frame(obs = stk$data$data$y[idx],
#                     p = apply(predy,2,mean),
#                     lo = apply(predy, 2, quantile, 0.25),
#                     hi = apply(predy, 2, quantile, 0.75))
# ggplot(preds,
#        aes(obs, p, ymin = lo, ymax = hi)) +
#   geom_errorbar(alpha = 0.5, col = "grey") +
#   geom_point(alpha = 0.1) + 
#   geom_abline(slope = 1) +
#   scale_x_continuous(trans = "sqrt") +
#   scale_y_continuous(trans = "sqrt") 

ggregplot::Efxplot(fit.nb2$fit, 
                   Intercept = FALSE) +
  theme(legend.position = "none")

# ---------------------------------------------------------------------------- #
# Multivariate fit with selection

f <- as.formula(paste0("y ~ ", 
                       paste0(covs3, collapse = " + "), 
                       "+ f(v, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(1, 0.01)) +
                            f(s, model = spde)"))

fit.multi2 <- init_inla(f, data.stack = stk, family = "nbinomial", cpo = F)

saveRDS(fit.multi2, here::here(outdir, "fit_covs_multivar_select.rds"))

# ---------------------------------------------------------------------------- #
# Compare fitted regression estimates

fits.all <- rlist::list.append(fits.uni, fit.multi, fit.multi2)

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

lapply(fits.all, function(x) summary(x$fit))

ggregplot::Efxplot(plyr::llply(fits.all, function(x) x$fit), 
                   VarNames = rev(axis.labs),
                   Intercept = FALSE
                   ) +
  scale_colour_manual(values = c(rep("indianred",length(fits.uni)),
                                 "steelblue",
                                 "forestgreen")) +
  theme(legend.position = "none") + 
  labs(caption = "Univariate estimate in red, multivariate in blue and multivariate with selected covariates only in green.")

ggsave(here::here(figdir, "covs_uni_multi_efx.png"), height = 6, width = 8, units = "in", dpi = 300)

################################################################################
################################################################################
