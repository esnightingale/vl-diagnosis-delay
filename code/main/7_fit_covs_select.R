################################################################################
# Description: Preliminary variable selection on patient characteristics, with
# non-spatial IID effects. 
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit/univariate"
outdir <- "output/univariate"

covs_all <- list("age_s","sex1","comorb1","prv_txTRUE","caste4_r1",c("occ4_cat1",
              "occ4_cat2", "occ4_cat3"), "poss_acdTRUE",
             "block_endm_2017TRUE", "IRS_2017TRUE","inc_2017_gt0TRUE",
             "traveltime_s","traveltime_t_s", "rainTRUE")
#c("age_cat.12.25.","age_cat.25.42.","age_cat.42.95.") 
covs_pat1 <- list("age_s","sex1","comorb1","prv_txTRUE","caste4_r1","occ4_cat1",
                                                                   "occ4_cat2", "occ4_cat3") 
covs_pat <- c(covs_pat1, list("poss_acdTRUE"))
covs_aware <- list("block_endm_2017TRUE", "IRS_2017TRUE","inc_2017_gt0TRUE")
covs_accessD <- list("traveltime_s", "rainTRUE")
covs_accessT <- list("traveltime_t_s", "rainTRUE") 

covs.multi.list <- list(Patient = covs_pat1, 
                        `Patient + detection` = covs_pat, 
                        Awareness = covs_aware, 
                        `Access: diag` = covs_accessD, 
                        `Access: trt` = covs_accessT)

mesh <- readRDS(here::here("data/analysis","mesh.rds")) 
stk <- readRDS(here::here("data/analysis","stack.rds"))

# Original sf object for plotting vgm
dat <- read_data() 

# Specify likelihood family = NB or Poisson
family <- "poisson"

# ---------------------------------------------------------------------------- #
# Null fit

f <- as.formula(paste0("y ~ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(10, 0.01))"))

fit.null <- inla(f,
                  family = family,
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
fit.null <- list(f = f, fit = fit.null)

p_vgm_null <- plot_vgm(fit.null$fit$summary.random$id$mean, dat, title = "Fitted IID effects: Null model",ylim = c(0.5,1))
#   model     psill    range kappa
# 1   Nug 0.5942502  0.00000   0.0
# 2   Mat 0.2564780 30.99597   0.5
# model     psill    range kappa
# 1   Nug 0.6274121  0.00000   0.0
# 2   Mat 0.1936387 18.62789   0.5

# Estimated variogram not sensitive to choice of prior precision

saveRDS(fit.null, here::here(outdir, "fit_null.rds"))

# ---------------------------------------------------------------------------- #
# Univariate fits

fit_covs <- function(cov) {
  
  # Define formula
  f <- as.formula(paste0("y ~", 
                         paste0(cov, collapse = " + "), 
                         "+ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(10, 0.01))"))
  
  print(f)
  
  # Fit model
  fit <- init_inla(f, data.stack = stk, family = family)$fit
  
  res <- list(f = f, fit = fit)
  return(res)
  
}

fits.uni <- plyr::llply(covs_all, fit_covs)
# saveRDS(fits.uni, here::here(outdir, "fits_covs_univar_nonspatial.rds"))

# ---------------------------------------------------------------------------- #
# Domain fits

fits.domain <- plyr::llply(covs.multi.list, fit_covs)

p_vgms_domain <- purrr::imap(fits.domain, function(x, nm) plot_vgm(x$fit$summary.random$id$mean, dat, 
                                                                  title = paste("Fitted IID effects:", nm),
                                                                  ylim = c(0.5,1)))

# model     psill    range kappa
# 1   Nug 0.7179050  0.00000   0.0
# 2   Mat 0.1600005 21.48299   0.5
# model     psill    range kappa
# 1   Nug 0.7156381  0.00000   0.0
# 2   Mat 0.1482520 20.83197   0.5
# model     psill    range kappa
# 1   Nug 0.6844554  0.00000   0.0
# 2   Mat 0.1854343 20.17741   0.5
# model     psill    range kappa
# 1   Nug 0.6564771  0.00000   0.0
# 2   Mat 0.1908631 18.72602   0.5
# model     psill   range kappa
# 1   Nug 0.6666980  0.0000   0.0
# 2   Mat 0.1901168 18.9717   0.5

# Compare travel time to diagnosis vs treatment
plyr::llply(fits.domain, function(x) summary(x$fit))
plyr::llply(fits.domain, function(x) x$fit$waic$waic) 
# $Patient
# [1] 28578.2
# 
# $`Patient + detection`
# [1] 28578.51
# 
# $Awareness
# [1] 28586.48
# 
# $`Access: diag`
# [1] 28585.52
# 
# $`Access: trt`
# [1] 28584.96

# Treatment facility travel time has a somewhat larger coefficient, but neither 
# significant on 95%CrI
  
# ---------------------------------------------------------------------------- #
# Full multivariate fit

covs_full <- c(covs_pat, covs_aware, covs_accessT) 

f <- as.formula(paste0("y ~ ", 
                       paste0(covs_full, collapse = " + "), 
                       "+ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(10, 0.01))"))

fit.full <- init_inla(f, data.stack = stk, family = family)

# saveRDS(fit.full, here::here(outdir, "fit_covs_multivar_nonspatial.rds"))

p_vgm_full <- plot_vgm(fit.full$fit$summary.random$id$mean, dat, title = "Fitted IID effects: Full model", ylim = c(0.5,1))
# model     psill    range kappa
# 1   Nug 0.7238787  0.00000   0.0
# 2   Mat 0.1316704 20.82969   0.5

# ---------------------------------------------------------------------------- #
# Compare univariate versus multivariate coefficients

fits.all <- rlist::list.append(c(fits.uni,
                                 setNames(fits.domain[-1], NULL)),
                               fit.full)
rlist::list.save(fits.all, here::here(outdir, "fits_covs_full_nonspatial.rds"))

c("Intercept",
  "Age (std)",
  "Sex: Male",
  "HIV positive",
  "Prv. VL/PKDL",
  "Marginalised caste (SC/ST)",
  "Occupation: Unskilled",
  "Skilled",
  "Salaried/self-employed",
  "Detection: ACD",
  "Block endemic (2017)",
  "Village IRS targeting (2017)",
  "Village incidence > 0 (2017)",
  "Travel time - diagnosis facility (std)",
  "Travel time - treatment facility (std)",
  "Diagnosis season: rain"
) -> axis.labs

Efxplot(plyr::llply(fits.all, function(x) x$fit),
        ModelGroups = c(rep("Univariate",length(fits.uni)),
                        rep("Domain multivariate",length(fits.domain)-1),
                        "Full multivariate"),
        VarNames = rev(axis.labs),
        Intercept = FALSE,
        Size = 1.5,
        # exp = FALSE
        ) +
  geom_vline(xintercept = c(3.5, 6.5), col = "lightgrey", lwd = 0.5) +
  scale_colour_viridis_d(option = "plasma", direction = -1, end = 0.8) +
  theme(
    # legend.position = "none",
    axis.text.y = element_text(face = rev(c("plain","bold","plain","bold",rep("plain",5),"bold","bold","plain","bold","plain","bold","plain"))),
    # legend.text = element_text(face = c("plain","bold","plain"))
  ) + 
  labs(title = "Estimated covariate effects: non-spatial models",
       colour = "")

ggsave(here::here(figdir, "covs_domain_efx_nonspatial.png"), height = 6, width = 8, units = "in", dpi = 300)

# Look at spatial correlation after including patient characteristics and detection
png(here::here(figdir, "vgms_domain_nonspatial.png"), height = 400, width = 2500)
gridExtra::grid.arrange(grobs = list.append(list.prepend(p_vgms_domain,p_vgm_null), 
                                            p_vgm_full), nrow = 1)
dev.off()

# ---------------------------------------------------------------------------- #
# Multivariate fit with selection

covs_select <- c("age_s","comorb1", "poss_acdTRUE",
                 "block_endm_2017TRUE","inc_2017_gt0TRUE","traveltime_t_s") 

f <- as.formula(paste0("y ~ ", 
                       paste0(covs_select, collapse = " + "), 
                       "+ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(10, 0.01))"))

fit.select <- init_inla(f, data.stack = stk, family = family)

# saveRDS(fit.select, here::here(outdir, "fit_covs_multivar_select_nonspatial.rds"))

p_vgm_select <- plot_vgm(fit.select$fit$summary.random$id$mean, dat, title = "Fitted IID effects: Full model with selection", ylim = c(0.5,1))
# model     psill   range kappa
# 1   Nug 0.7231392  0.0000   0.0
# 2   Mat 0.1328259 20.0026   0.5

png(here::here(figdir, "vgms_select_nonspatial.png"), height = 1000, width = 500)
gridExtra::grid.arrange(p_vgm_null, p_vgm_select, nrow = 2)
dev.off()

# ---------------------------------------------------------------------------- #
# Compare univariate, multivariate and multivariate with selection

fits.all <- rlist::list.append(c(fits.uni, 
                                 setNames(fits.domain[-1], NULL)),
                               fit.select)
rlist::list.save(fits.all, here::here(outdir, "fits_covs_select_all_nonspatial.rds"))

c("Intercept",
  "Age (std)",
  "Sex: Male",
  "HIV positive",
  "Prv. VL/PKDL treatment",
  "Marginalised caste (SC/ST)",
  "Occupation: Unskilled",
  "Skilled",
  "Salaried/self-employed",
  "Detection: ACD",
  "Block endemic (2017)",
  "Village IRS targeting (2017)",
  "Village incidence > 0 (2017)",
  "Travel time - diagnosis facility (std)",
  "Travel time - treatment facility (std)",
  "Diagnosis season: rain"
) -> axis.labs


png(here::here(figdir, "covs_select_efx_nonspatial.png"), height = 8, width = 8, unit = "in", res = 320)
Efxplot(plyr::llply(fits.all, function(x) x$fit),
                   ModelGroups = c(rep("Univariate",length(fits.uni)), 
                                   rep("Domain multivariate",length(fits.domain)-1), 
                                   "Selected multivariate"),
                   VarNames = rev(axis.labs),
                   Intercept = FALSE, Size = 1.5) +
  geom_vline(xintercept = c(3.5, 6.5), col = "lightgrey", lwd = 1) +
  scale_colour_viridis_d(option = "plasma", end= 0.8, direction = -1) +
  scale_y_continuous(breaks = seq(0.6,1.8, by = 0.2)) +
  theme(legend.position = c(0.75,0.2),
        legend.title = element_blank(),
        legend.box.background = element_rect(color="white", size=2),
        axis.text.y = element_text(face = rev(c("plain","bold","plain","bold",rep("plain",5),"bold","bold","plain","bold","plain","bold","plain"))),
        # legend.text = element_text(face = c("plain","bold","plain"))
        ) + 
  labs(#title = "Estimated covariate effects: non-spatial models",
       # colour = "",
       caption = "Selected covariates from each domain are highlighted in bold.")
dev.off()

################################################################################
################################################################################
