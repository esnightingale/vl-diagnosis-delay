################################################################################
# Description: 
# 
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit/excl_lt14"
outdir <- "output/excl_lt14"

dat <- readRDS(here::here("data/analysis","dat_nona.rds")) %>%
  st_transform(crs = st_crs(7759))
mesh <- readRDS(here::here("data/analysis","mesh.rds"))
spde <- readRDS(here::here("data/analysis","spde.rds"))
stk <- readRDS(here::here("data/analysis","stack_excl_lt14.rds"))[[1]]
coop <- readRDS(here::here("data/analysis","coop.rds"))

blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
  st_transform(7759)

boundary <- sf::st_union(blockmap)
boundary.spdf <- as_Spatial(boundary)

fits.all <- readRDS(here::here(outdir,"fits_covs_all.rds"))

#------------------------------------------------------------------------------#
# Simple fit summaries

mod.compare <- data.frame(DIC = sapply(fits.all, function(x) x$fit$dic$dic),
                          WAIC = sapply(fits.all, function(x) x$fit$waic$waic),
                          MAE = sapply(fits.all, function(x) get_mae(summ_pred(x$fit, data.stack = stk))),
                          RMSE = sapply(fits.all, function(x) get_rmse(summ_pred(x$fit, data.stack = stk))),
                          Logscore = sapply(fits.all, function(x) -mean(log(x$fit$cpo$cpo), na.rm = T)),
                          cpoFail = sapply(fits.all, function(x) sum(x$fit$cpo$failure > 0, na.rm = T))) %>%
  rownames_to_column("Model")

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

png(here::here(figdir, "fitted_effects_main.png"), height = 2000, width = 2400, res = 300)
ggregplot::Efxplot(lapply(fits.all[-c(9,10)], function(x) x$fit),
                   ModelNames = names(fits.all)[-c(9,10)],
                   Intercept = FALSE,
                   VarNames = rev(axis.labs))  
dev.off()

order <- c("age_s","sexFemale","hivYes","prv_txYes","marg_casteYes","occupationUnskilled",
              "occupationSkilled", "`occupationSalaried.selfemployed`", "detectionACD","rain",
              "block_endm_2017Endemic", "IRS_2017_1Yes","vill_inc_2017_gt0Yes","traveltime_s","traveltime_s:rain")

png(here::here(figdir, "fitted_effects_supp.png"), height = 2000, width = 2400, res = 300)
ggregplot::Efxplot(lapply(fits.all[8:10], function(x) x$fit),
                   ModelNames = c("Full model","IID only","SPDE only"),
                   Intercept = FALSE,
                   VarOrder= rev(order),
                   VarNames = rev(axis.labs))  
dev.off()


#------------------------------------------------------------------------------#
# Compare fitted SPDEs

# Exclude IID-only model here
spde.fits <- fits.all[-which(names(fits.all) == "All (IID only)")] 

# Compare range and SD posteriors:
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

pdf(here::here(figdir, "compare_spde_quants.pdf"), height = 7, width = 10) 
# purrr::imap(spde.fits, function(x, nm) plot_spde(x$fit, nm, stat = "mean", limit1 = c(-1,1), limit2 = c(0,0.5)))
purrr::imap(spde.fits, function(x, nm) plot_spde(x$fit, nm, stat = "quantile", limit1 = c(-1,1), limit2 = c(-1,1)))
dev.off()

plot_spde(spde.fits$All$fit, "All", stat = "quantile", limit1 = NULL, limit2 = NULL) #limit1 = c(-1,1), limit2 = c(-1,1))

## Summary measures of variation in SPDE/IID effects (mean absolute values:

tab_random_mav <- bind_rows(purrr::imap(fits.all, function(x, nm) random_mav(x$fit, name = nm)))

# Calculate difference from baseline model
base_mav_iid <- tab_random_mav$mav_iid[tab_random_mav$Model == "None"]
base_mav_spde <- tab_random_mav$mav_spde[tab_random_mav$Model == "None"]
tab_random_mav$pdiff_mav_iid <- (tab_random_mav$mav_iid - base_mav_iid)*100/base_mav_iid
tab_random_mav$pdiff_mav_spde <- (tab_random_mav$mav_spde - base_mav_spde)*100/base_mav_spde

base_msv_iid <- tab_random_mav$msv_iid[tab_random_mav$Model == "None"]
base_msv_spde <- tab_random_mav$msv_spde[tab_random_mav$Model == "None"]
tab_random_mav$pdiff_msv_iid <- (tab_random_mav$msv_iid - base_msv_iid)*100/base_msv_iid
tab_random_mav$pdiff_msv_spde <- (tab_random_mav$msv_spde - base_msv_spde)*100/base_msv_spde

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

pdf(here::here(figdir, "baseline_diff_spde_cat.pdf"), height = 7, width = 10) 
plot_spde(spde.fits[[1]]$fit, "Null", stat = "mean", limit1 = c(-1,1), limit2 = c(0,0.5))
purrr::imap(spde.fits[-c(1,9)], function(x, nm) plot_diff_spde_cat(x$fit, name = nm))
dev.off()

pdf(here::here(figdir, "baseline_diff_spde.pdf"), height = 7, width = 10) 
purrr::imap(spde.fits[-c(1,9)], function(x, nm) plot_diff_spde(x$fit, name = nm))
dev.off()

pdf(here::here(figdir, "baseline_diff_spde_abs.pdf"), height = 7, width = 10) 
purrr::imap(spde.fits[-c(1,9)], function(x, nm) plot_diff_spde(x$fit, name = nm, absolute = TRUE, limit.mean = c(0,0.1), limit.sd = c(0,0.01)))
dev.off()

pdf(here::here(figdir, "baseline_diff_spde_only.pdf"), height = 7, width = 10) 
plot_diff_spde(spde.fits[[9]]$fit, name = "All (SPDE only) ", limit.mean = NULL, limit.sd = NULL)
dev.off()

# ---------------------------------------------------------------------------- #
# Summarise fitted versus observed

idx <- stk$data$index$train
C <- 30 

compare_fitted <- function(fit){
  
  fitted <- data.frame(ll = fit$summary.fitted.values[idx, "0.025quant"],
                       med = fit$summary.fitted.values[idx, "0.5quant"],
                       pred = fit$summary.fitted.values[idx, "mean"],
                       ul = fit$summary.fitted.values[idx, "0.975quant"],
                       exc.prob = sapply(fit$marginals.fitted.values[idx],
                                         FUN = function(marg){1-inla.pmarginal(q = C, marginal = marg)}),
                       # Match to observed values from the full stack at the specified indices
                       obs = stk$data$data$y[idx]) %>%
    dplyr::mutate(exc.obs = (obs > C),
                  abs.err = abs(pred - obs),
                  sq.err = (pred - obs)^2)

  return(fitted)

}

obs_vs_fitted <- lapply(fits.all, function(x) compare_fitted(x$fit))

pdf(here::here(figdir,"obs_vs_fitted.pdf"), height = 7, width = 8)
purrr::imap(obs_vs_fitted, function(x, nm) plot_preds(x, name = nm))
dev.off()

pdf(here::here(figdir,"obs_vs_fitted_exc40.pdf"), height = 7, width = 8)
purrr::imap(obs_vs_fitted, function(x, nm) plot_exc(x, name = nm))
dev.off()

pdf(here::here(figdir,"obs_vs_fitted.pdf"), height = 7, width = 8)
purrr::imap(obs_vs_fitted, function(x, nm) plot_resids(x, name = nm))
dev.off()

# ---------------------------------------------------------------------------- #
# Compare predictions from covariate-adjusted SPDEs to null SPDE 

# Extract fitted values at prediction points

pidx <- stk$data$index$pred

# Evaluate fitted SPDE from null model over this projection
mean_base <- spde.fits$None$fit$summary.fitted.values[pidx, "mean"]
sd_base <- spde.fits$None$fit$summary.fitted.values[pidx, "sd"]

# Set up a data frame for plotting
df <- data.frame(x = coop[,1], y = coop[,2])
df$mean_base <- as.vector(mean_base)
df$sd_base <- as.vector(sd_base)

# Plot null prediction
gnull <- ggplot() + 
  geom_raster(data = df, aes(x = x, y = y, fill = mean_base)) +
  gg(boundary.spdf, fill = "transparent") +
  scale_fill_viridis_c(option = "magma", direction = 1, na.value = "transparent", limits = NULL) +
  labs(title = "Null",
       subtitle = "Mean",
       fill = "",
       x = "", y = "") +
  coord_fixed(ratio = 1) + 
  theme_bw()

pdf(here::here(figdir, "baseline_diff_pred_cat.pdf"), height = 7, width = 10) 
gnull
purrr::imap(spde.fits[-c(1,9)], function(x, nm) plot_diff_pred(x$fit, name = nm))
dev.off()

################################################################################
################################################################################
