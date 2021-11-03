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
stk <- readRDS(here::here("data/analysis","stack.rds"))

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

#------------------------------------------------------------------------------#
# Simple fit summaries

mod.compare <- data.frame(DIC = sapply(fits.all, function(x) x$fit$dic$dic),
                          MAE = sapply(fits.all, function(x) get_mae(summ_pred(x$fit, data.stack = stk))),
                          RMSE = sapply(fits.all, function(x) get_rmse(summ_pred(x$fit, data.stack = stk))),
                          Logscore = sapply(fits.all, function(x) -mean(log(x$fit$cpo$cpo), na.rm = T)),
                          cpoFail = sapply(fits.all, function(x) sum(x$fit$cpo$failure > 0, na.rm = T)))

write.csv(mod.compare, here::here(outdir, "mod_compare.csv"), row.names = FALSE)

pdf(here::here(figdir, "pit_histograms.pdf"), height = 7, width = 10) 
purrr::imap(fits.all, function(x, nm) pit_hist(x$fit, title = nm))
dev.off()

#------------------------------------------------------------------------------#
# Compare fitted regression estimates

c("Intercept",
  "Age (std)",
  "Sex: Female",
  "HIV positive",
  "Prv. VL/PKDL",
  "Marginalised caste",
  "Occupation: Unskilled",
  "Skilled",
  "Salaried/self-employed",
  "Detection: ACD",
  "Diagnosis season: rain",
  "Block endemic (2017)",
  "Village IRS targeting (2017)",
  "Village incidence > 0 (2017)",
  "Travel time (std)",
  "Travel time*rainy season") -> axis.labs

png(here::here(figdir, "fitted_effects.png"), height = 2000, width = 2400, res = 300)
ggregplot::Efxplot(lapply(fits.all, function(x) x$fit),
                   ModelNames = names(fits.all),
                   Intercept = FALSE,
                   VarNames = rev(axis.labs))  
dev.off()

#------------------------------------------------------------------------------#
# Compare fitted SPDEs

# Exclude IID-only model here
spde.fits <- fits.all[-7]

## Compare range and SD posteriors:

spde.range <- bind_rows(plyr::llply(spde.fits, function(x) x$fit$summary.hyperpar["Range for s",c(1,3,5)])) %>%
  mutate(Model = 1:length(spde.fits))

png(here::here(figdir, "mod_compare_spde_range.png"), height = 1500, width = 3000, res = 300)
ggplot(spde.range, aes(Model, mean, ymin = `0.025quant`, ymax = `0.975quant`)) +
  geom_errorbar(width = 0.5) +
  geom_point() +
  labs(y = "Range") +
  scale_x_discrete(limits = 1:9, labels = names(spde.fits)) +
  scale_y_continuous(trans = "log10",labels = scales::number_format(accuracy = 1000))
dev.off()

spde.sd <- bind_rows(plyr::llply(spde.fits, function(x) x$fit$summary.hyperpar["Stdev for s",c(1,3,5)])) %>%
  mutate(Model = 1:length(spde.fits)) 

png(here::here(figdir, "mod_compare_spde_stdev.png"), height = 1500, width = 3000, res = 300)
ggplot(spde.sd, aes(Model, mean, ymin = `0.025quant`, ymax = `0.975quant`)) +
  geom_errorbar(width = 0.5) +
  geom_point() +
  labs(y = "Stdev") +
  scale_x_discrete(limits = 1:9, labels = names(spde.fits)) +
  scale_y_continuous(trans = "log10")
dev.off()


## Plot projection of the fitted SPDEs:

pdf(here::here(figdir, "compare_spde_pat_all_3cat.pdf"), height = 7, width = 10) 
purrr::imap(fits.list, function(x, nm) plot_spde(x$fit, nm, limit.mean = c(-0.6,0.6), limit.sd = c(0,0.5)))
dev.off()


## Summary measures of variation in SPDE/IID effects (mean absolute values:

tab_random_mav <- bind_rows(purrr::imap(fits.all, function(x, nm) random_mav(x$fit, name = nm)))

# Calculate difference from baseline model
base_mav_iid <- tab_random_mav$mean_abs_iid[tab_random_mav$Model == "None"]
base_mav_spde <- tab_random_mav$mean_abs_spde[tab_random_mav$Model == "None"]
tab_random_mav$pdiff_base_iid <- (tab_random_mav$mav_iid - base_mav_iid)*100/base_mav_iid
tab_random_mav$pdiff_base_spde <- (tab_random_mav$mav_spde - base_mav_spde)*100/base_mav_spde

write.csv(tab_random_mav, here::here(outdir, "tab_random_mav.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------- #
# Percent difference between each SPDE and baseline, over space

# Define a projection over the region
rang <- apply(mesh$loc[, c(1, 2)], 2, range)
proj <- inla.mesh.projector(mesh, 
                            xlim = rang[, 1], 
                            ylim = rang[, 2], 
                            dims = c(300, 300))

# Evaluate fitted SPDE from null model over this projection
mean_base <- inla.mesh.project(proj, fits.all$None$fit$summary.random$s$mean)
sd_base <- inla.mesh.project(proj, fits.all$None$fit$summary.random$s$sd)

# Set up a data frame for plotting
df <- expand.grid(x = proj$x, y = proj$y)
df$mean_base <- as.vector(mean_base)
df$sd_base <- as.vector(sd_base)

pdf(here::here(figdir, "baseline_diff_spde.pdf"), height = 7, width = 10) 
purrr::imap(fits.all[-c(1,7)], function(x, nm) plot_diff_spde(x$fit, name = nm))
dev.off()

################################################################################
################################################################################

