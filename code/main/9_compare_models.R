################################################################################
# Description: 
# 
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit"
outdir <- "output"

dat <- read_data()
mesh <- readRDS(here::here("data/analysis","mesh.rds"))
spde <- readRDS(here::here("data/analysis","spde.rds"))
stk <- readRDS(here::here("data/analysis","stack.rds"))
# idx <- stk$data$index$train

boundary <- readRDS(here::here("data","geography","boundary.rds"))
boundary.spdf <- as_Spatial(boundary)

# Spatial models
fits <- readRDS(here::here(outdir,"fits_final.rds"))

# null model
fit.null <- readRDS(here::here(outdir,"univariate", "fit_null.rds"))

fits <- rlist::list.prepend(fits[-9], fit.null) # exclude interaction fit
names(fits)[1] <- "Null"

# Cut off for exceedances
C = 30 

#------------------------------------------------------------------------------#
# Simple fit summaries

fit_summs <- lapply(fits, function(x) summ_fit(x$fit, data.stack = stk, C = C))

mod.compare <- data.frame(WAIC = sapply(fits, function(x) x$fit$waic$waic),
                          mLL = sapply(fits, function(x) x$fit$mlik["log marginal-likelihood (Gaussian)",]),
                          MSE = sapply(fit_summs, get_mse),
                          brier = sapply(fit_summs, function(x) mean((x$exc.prob - as.numeric(x$exc.obs))^2))) %>%
  rownames_to_column("Model")
write.csv(mod.compare, here::here(outdir, "mod_compare_all.csv"), row.names = FALSE)

mod.compare.main <- filter(mod.compare, Model %in% c("Null","All (IID only)", "None", "All", "All (SPDE only)"))
write.csv(mod.compare.main, here::here(outdir, "mod_compare.csv"), row.names = FALSE)

# pdf(here::here(figdir, "pit_histograms.pdf"), height = 7, width = 10) 
# purrr::imap(fits, function(x, nm) pit_hist(x$fit, title = nm))
# dev.off()

# ---------------------------------------------------------------------------- #
# Summarise fitted versus observed

pdf(here::here(figdir,"obs_vs_fitted.pdf"), height = 7, width = 8)
purrr::imap(fit_summs, function(x, nm) plot_preds(x, name = nm))
dev.off()

pdf(here::here(figdir,"obs_vs_fitted_exc30.pdf"), height = 7, width = 8)
purrr::imap(fit_summs, function(x, nm) plot_exc(x, name = nm))
dev.off()

# pdf(here::here(figdir,"obs_vs_fitted.pdf"), height = 7, width = 8)
# purrr::imap(obs_vs_fitted, function(x, nm) plot_resids(x, name = nm))
# dev.off()

#------------------------------------------------------------------------------#
# Precision of IID effects and spatial range/SD

fits$Null$fit$summary.hyperpar
#                      mean         sd 0.025quant 0.5quant 0.975quant     mode
# Precision for id 1.020131 0.02418332  0.9724507 1.019738   1.068401 1.017056

fits$`All (IID only)`$fit$summary.hyperpar
#                      mean         sd 0.025quant 0.5quant 0.975quant     mode
# Precision for id 1.077161 0.02572244   1.026203 1.076932   1.128319 1.075659

fits$All$fit$summary.hyperpar
#                        mean          sd 0.025quant   0.5quant 0.975quant       mode
# Precision for id  1.1652845  0.03121173  1.1014187  1.1664947  1.2237411  1.1714110
# Range for s      47.2586851 14.85399471 26.1310378 44.5197495 83.5768644 39.5595050
# Stdev for s       0.3239117  0.04817799  0.2297483  0.3248052  0.4163401  0.3308667

#------------------------------------------------------------------------------#
# Spatial correlation in fitted IID effects

p_vgm_null <- plot_vgm(fits$Null$fit$summary.random$id$mean, dat, title = "Null model",ylim = c(0.5,1))

p_vgm_select <- plot_vgm(fits$`All (IID only)`$fit$summary.random$id$mean, dat, title = "Selected covariates, non-spatial model",ylim = c(0.5,1))

p_vgm_final <- plot_vgm(fits$All$fit$summary.random$id$mean, dat, title = "Selected covariates, spatial model",ylim = c(0.5,1))

png(here::here(figdir, "vgms_null_select_final.png"), height = 500, width = 1500)
gridExtra::grid.arrange(p_vgm_null, p_vgm_select, p_vgm_final, nrow = 1)
dev.off()

coo <- sf::st_coordinates(dat)
png(here::here(figdir, "MI_null_select_final.png"), height = 400, width = 1200)
par(mfrow = c(1,3))
lapply(list(Null = fits$Null$fit$summary.random$id$mean,
            IID = fits$`All (IID only)`$fit$summary.random$id$mean,
            Final = fits$All$fit$summary.random$id$mean),
       function(x) calc_MI(x,coo, range = 30))
dev.off()

#------------------------------------------------------------------------------#
# Summary measures of variation in SPDE/IID effects (mean absolute values:

tab_random_mav <- bind_rows(purrr::imap(fits[-which(names(fits) == "All (Binomial)")], function(x, nm) random_mav(x$fit, name = nm)))

# Calculate difference from baseline model
base_mav_iid <- tab_random_mav$mav_iid[tab_random_mav$Model == "None"]
base_mav_spde <- tab_random_mav$mav_spde[tab_random_mav$Model == "None"]
tab_random_mav$pdiff_mav_iid <- (tab_random_mav$mav_iid - base_mav_iid)*100/base_mav_iid
tab_random_mav$pdiff_mav_spde <- (tab_random_mav$mav_spde - base_mav_spde)*100/base_mav_spde

base_msv_iid <- tab_random_mav$msv_iid[tab_random_mav$Model == "None"]
base_msv_spde <- tab_random_mav$msv_spde[tab_random_mav$Model == "None"]
tab_random_mav$pdiff_msv_iid <- (tab_random_mav$msv_iid - base_msv_iid)*100/base_msv_iid
tab_random_mav$pdiff_msv_spde <- (tab_random_mav$msv_spde - base_msv_spde)*100/base_msv_spde

write.csv(tab_random_mav, here::here(outdir, "tab_random_mav_refnull.csv"), row.names = FALSE)

#------------------------------------------------------------------------------#
# Compare fitted regression estimates

c("Intercept",
  "Age (std)",
  # "Sex: Female",
  "HIV positive",
  # "Prv. VL/PKDL",
  # "Marginalised caste",
  # "Occupation: Unskilled",
  # "Skilled",
  # "Salaried/self-employed",
  "Detection: ACD",
  # "Diagnosis season: rain",
  "Block endemic (2017)",
  # "Village IRS targeting (2017)",
  "Village incidence > 0 (2017)",
  "Travel time - treatment facility (std)"
  # "Diagnosis season: rain",
  # "Travel time*rainy season"
  ) -> axis.labs

png(here::here(figdir, "fitted_effects_main.png"), height = 6, width = 8, unit = "in", res = 320)
Efxplot(lapply(fits[c(10, 11, 9)], function(x) x$fit),
        ModelNames = c("Fixed effects only","Non-spatial random effects","Spatial random effects"),
        Intercept = FALSE,
        exp = TRUE,
        # VarOrder= rev(order),
        VarNames = rev(axis.labs)) +
  scale_colour_viridis_d(option = "plasma", begin = 0.2, end = 0.8, direction = 1) +
  theme(legend.position = c(0.75,0.2),
        legend.title = element_blank(),
        legend.box.background = element_rect(color="white", size=2))
# theme(legend.position = c(0.8,0.2))
dev.off()


png(here::here(figdir, "fitted_effects_supp1.png"), height = 6, width = 8, unit = "in", res = 320)
Efxplot(lapply(fits[3:9], function(x) x$fit),
                   ModelNames = names(fits)[3:9],
                   Intercept = FALSE,
                   exp = TRUE,
                   VarNames = rev(axis.labs)) +
  scale_colour_viridis_d(option = "plasma", end = 0.8, direction = -1) +
  theme(legend.position = c(0.8,0.25),
        legend.title = element_blank(),
        legend.box.background = element_rect(color="white", size=2))
dev.off()

png(here::here(figdir, "fitted_effects_supp2.png"), height = 6, width = 8, unit = "in", res = 320)
Efxplot(lapply(fits[c(11,9,13)], function(x) x$fit),
                   ModelNames = c("Non-spatial (IID)", "Spatial (IID+SPDE)","Spatial (Binomial likelihood)"),
                   Intercept = FALSE,
                   exp = TRUE,
                   # VarOrder= rev(order),
                   VarNames = rev(axis.labs)) +
  scale_colour_viridis_d(option = "plasma", begin = 0.2, end = 0.8, direction = 1) +
  theme(legend.position = c(0.75,0.2),
        legend.title = element_blank(),
        legend.box.background = element_rect(color="white", size=2))
  # theme(legend.position = c(0.8,0.2))
dev.off()

#------------------------------------------------------------------------------#
# Compare fitted SPDEs

# Exclude IID-only model here
spde.fits <- fits[-which(names(fits) %in% c("Null","All (IID only)", "Fixed effects only"))] 

names <- gsub("\\+", "+\n", names(spde.fits))

# Compare range and SD posteriors:
spde.range <- bind_rows(plyr::llply(spde.fits, function(x) x$fit$summary.hyperpar["Range for s",c(1,3,5)])) %>% 
  dplyr::mutate(Model = 1:length(spde.fits))

png(here::here(figdir, "mod_compare_spde_range.png"), height = 1500, width = 3500, res = 300)
ggplot(spde.range, aes(Model, mean, ymin = `0.025quant`, ymax = `0.975quant`)) +
  geom_errorbar(width = 0.5) +
  geom_point() +
  labs(y = "Range", x = "") +
  scale_x_discrete(limits = factor(1:10), 
                   labels = names) #+
  # scale_y_continuous(trans = "log10",labels = scales::number_format(accuracy = 1000))
dev.off()

spde.sd <- bind_rows(plyr::llply(spde.fits, function(x) x$fit$summary.hyperpar["Stdev for s",c(1,3,5)])) %>%
  mutate(Model = 1:length(spde.fits)) 

png(here::here(figdir, "mod_compare_spde_stdev.png"), height = 1500, width = 3500, res = 300)
ggplot(spde.sd, aes(Model, mean, ymin = `0.025quant`, ymax = `0.975quant`)) +
  geom_errorbar(width = 0.5) +
  geom_point() +
  labs(y = "Stdev", x = "") +
  scale_x_discrete(limits = factor(1:10), labels = names) +
  scale_y_continuous(trans = "log10")
dev.off()

#------------------------------------------------------------------------------#
## Plot projection of the fitted SPDEs:

# Project fitted SPDEs
lims <- range(spde.fits$None$fit$summary.random$s$mean)
p.list <- purrr::imap(spde.fits[1:8], function(x, nm) plot_spde_mean(x$fit, nm, limits = lims)) 
png(here::here(figdir, "compare_spde_fits.png"), height = 1000, width = 3000, res = 350) 
gridExtra::grid.arrange(grobs = p.list[c(1,8)], nrow = 1) 
dev.off()

png(here::here(figdir, "compare_spde_fits_all.png"), height = 4000, width = 9000, res = 350) 
gridExtra::grid.arrange(grobs = p.list, nrow = 2) 
dev.off()

lims2 <- range(spde.fits$`All (Binomial)`$fit$summary.random$s$mean)
p.list.supp <- purrr::imap(spde.fits[c(8,10)], function(x, nm) plot_spde_mean(x$fit, nm, limits = lims2))
png(here::here(figdir, "compare_spde_fits_supp.png"), height = 1000, width = 3000, res = 350) 
gridExtra::grid.arrange(grobs = p.list, nrow = 1) 
dev.off()

png(here::here(figdir, "spde_only_fit_supp.png"), height = 1000, width = 1300, res = 350) 
plot_spde_mean(spde.fits$`All (SPDE only)`$fit, "All (SPDE only)")
dev.off()

# pdf(here::here(figdir, "compare_spde_quants.pdf"), height = 7, width = 10) 
# # purrr::imap(spde.fits, function(x, nm) plot_spde(x$fit, nm, stat = "mean", limit1 = c(-1,1), limit2 = c(0,0.5)))
# purrr::imap(spde.fits, function(x, nm) plot_spde(x$fit, nm, stat = "quantile", limit1 = c(-1,1), limit2 = c(-1,1)))
# dev.off()
# 
# plot_spde(spde.fits$All$fit, "All", stat = "quantile", limit1 = c(-1.5,1.5), limit2 = c(-1.5,1.5)) #limit1 = c(-1,1), limit2 = c(-1,1))

# ---------------------------------------------------------------------------- #
# Percent difference between each SPDE and baseline, over space
# 
# # Define a projection over the region
# rang <- apply(mesh$loc[, c(1, 2)], 2, range)
# proj <- inla.mesh.projector(mesh, 
#                             xlim = rang[, 1], 
#                             ylim = rang[, 2], 
#                             dims = c(300, 300))
# 
# # Evaluate fitted SPDE from null model over this projection
# mean_base <- inla.mesh.project(proj, fits$None$fit$summary.random$s$mean)
# sd_base <- inla.mesh.project(proj, fits$None$fit$summary.random$s$sd)
# 
# # Set up a data frame for plotting
# df <- expand.grid(x = proj$x, y = proj$y)
# df$mean_base <- as.vector(mean_base)
# df$sd_base <- as.vector(sd_base)
# 
# pdf(here::here(figdir, "baseline_diff_spde_cat.pdf"), height = 7, width = 10) 
# plot_spde(spde.fits[[1]]$fit, "Null", stat = "mean", limit1 = c(-1,1), limit2 = c(0,0.5))
# purrr::imap(spde.fits[-1], function(x, nm) plot_diff_spde_cat(x$fit, name = nm))
# dev.off()
# 
# pdf(here::here(figdir, "baseline_diff_spde.pdf"), height = 7, width = 10) 
# purrr::imap(spde.fits[-1], function(x, nm) plot_diff_spde(x$fit, name = nm))
# dev.off()

################################################################################
################################################################################
