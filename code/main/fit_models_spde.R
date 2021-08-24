################################################################################
# Description: Fit spde model with chosen covariates and compare to base/IID
################################################################################
################################################################################

figdir <- "figures/fit"
outdir <- "output"

#------------------------------------------------------------------------------#
# Fit three models: fixed effect only, with IID random effects and with spde on village


f1 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + ")))
f2 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + "), " +  f(id, model = 'iid')"))
f3 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + "), " +  f(id, model = 'iid') + f(s, model = spde)"))

fit.fixed <- inla(f1,
                family = "poisson",
                data = inla.stack.data(stk.full),
                control.predictor = list(
                  compute = TRUE, link = 1,
                  A = inla.stack.A(stk.full)),
                control.compute = list(dic = TRUE, 
                                       waic = TRUE, 
                                       # cpo = TRUE,
                                       config = TRUE),
                control.fixed = list(mean = 0, prec = 0.1, 
                                     mean.intercept = 0, prec.intercept = 0.001),
                verbose = TRUE
)

fit.iid <- inla(f2,
                 family = "poisson",
                 data = inla.stack.data(stk.full),
                 control.predictor = list(
                   compute = TRUE, link = 1,
                   A = inla.stack.A(stk.full)),
                 control.compute = list(dic = TRUE, 
                                        waic = TRUE, 
                                        # cpo = TRUE,
                                        config = TRUE),
                 control.fixed = list(mean = 0, prec = 0.1, 
                                      mean.intercept = 0, prec.intercept = 0.001),
                 verbose = TRUE
)

fit.spde <- inla(f3,
                 family = "poisson",
                 data = inla.stack.data(stk.full),
                 control.predictor = list(
                   compute = TRUE, link = 1,
                   A = inla.stack.A(stk.full)),
                 control.compute = list(dic = TRUE, 
                                        waic = TRUE, 
                                        # cpo = TRUE,
                                        config = TRUE),
                 control.fixed = list(mean = 0, prec = 0.1, 
                                      mean.intercept = 0, prec.intercept = 0.001),
                 verbose = TRUE
)

saveRDS(fit.fixed, here::here(outdir,"fit_final_fixed.rds"))
saveRDS(fit.iid, here::here(outdir,"fit_final_iid.rds"))
saveRDS(fit.spde, here::here(outdir,"fit_final_spde.rds"))

#------------------------------------------------------------------------------#
# Compare models

mod.list <- list(Fixed = fit.fixed, IID = fit.iid, SPDE = fit.spde)

lapply(mod.list, summary)
# lapply(mod.list, pit_hist)

# all.mod.list <- c(mod.list, vil.mod.list)
# names(all.mod.list) <- c("Base", "IID", "SPDE", "Inc17", "IRS17", "BlkEndm", "Travel", "All")

sapply(mod.list, function(f) f$dic$dic)
#     Base      IID     SPDE
# 76288.29 23664.60 23651.42 23720.88 23720.73 23719.25 23720.43 23719.40
# No added village covariates reduce DIC

png(here::here(figdir, "final_models_dic.png"), height = 5, width = 8, units = "in", res = 300)
ggregplot::INLADICFig(mod.list, ModelNames = c("Fixed", "IID", "SPDE"))
dev.off()

png(here::here(figdir, "final_covariate_effects.png"), height = 5, width = 8, units = "in", res = 300)
ggregplot::Efxplot(mod.list, ModelNames = c("Fixed", "IID", "SPDE"), Intercept = FALSE)
dev.off()

# ---------------------------------------------------------------------------- #
# Assess spatial independence of residuals from SPDE model

# Identify indices which correspond to validation points
index.e <- inla.stack.index(stack = stk.full, tag = "est")$data

dat.fit$resids.base <- (dat.fit$days_fever - fit.base$summary.fitted.values$mean[index.e])^2
vg <- variogram(resids.base~1, 
                data = dat.fit)
plot(vg)

dat.fit$resids.iid <- (dat.fit$days_fever - fit.iid$summary.fitted.values$mean[index.e])^2
vg <- variogram(resids.iid~1, 
                data = dat.fit)
plot(vg)

vg <- variogram(fit.iid$summary.random$id$mean[index.e]~1,
                data = dat.fit)
plot(vg)

vgmod <- vgm(psill = 3e7, model = "Mat", nugget = 1e7, range = 100)
plot(vg, model = vgmod)

vgfit <- fit.variogram(vg, model = vgmod)    
vgfit
# model    psill    range kappa
# 1   Nug 10583277  0.00000   0.0
# 2   Mat 29683113 66.16962   0.5

png(here::here(figdir, "allcovs_fixed_semivariogram.png"), height = 500, width = 600)
plot(vg, model = vgfit, main = "Semivariogram of fixed effect model residuals")
dev.off()

#------------------------------------------------------------------------------#
# Plot the fitted spatial field

rang <- apply(mesh$loc[, c(1, 2)], 2, range)

proj <- inla.mesh.projector(mesh, 
                            xlim = rang[, 1], 
                            ylim = rang[, 2], 
                            dims = c(300, 300))

mean_i <- inla.mesh.project(proj, fit.spde$summary.random$s$mean)
sd_i <- inla.mesh.project(proj, fit.spde$summary.random$s$sd)

df <- expand.grid(x = proj$x, y = proj$y)
df$mean_i <- as.vector(mean_i)
df$sd_i <- as.vector(sd_i)

pal <- viridis(2)
gmean <- ggplot() + 
  geom_raster(data = df, aes(x = x, y = y, fill = mean_i)) +
  gg(boundary.spdf, fill = NA) +
  scale_fill_gradient2(na.value = "transparent", low = pal[1], mid = "white", high = pal[2]) +
  coord_fixed(ratio = 1) + 
  theme_bw()

gsd <- ggplot() + 
  geom_raster(data = df, aes(x = x, y = y, fill = sd_i)) +
  gg(boundary.spdf, fill = NA) +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + 
  theme_bw()

cowplot::plot_grid(gmean, gsd)

# ggregplot::INLARange(list(fit.spde), MaxRange = 1, MeshList = list(mesh))

#------------------------------------------------------------------------------#
# Check prediction at validation points

# Identify indices which correspond to validation points
index.v <- inla.stack.index(stack = stk.full, tag = "val")$data

# Extract fitted marginals
# Calculate probability of exceeding 30 days from this marginal distribution

calc_excprob <- function(res, dat, idx, exc = 30){
  dat %>% 
    mutate(exceed = (days_fever > exc),
           excprob = sapply(res$marginals.fitted.values[idx],
                            FUN = function(marg){1 - inla.pmarginal(q = exc, marginal = marg)})) %>%
    return()
}

exceed30 <- bind_rows(lapply(mod.list, calc_excprob, dat = dat.val.df, idx = index.v, exc = 30),
                      .id = "model")

exc = 30
png(here::here(figdir, "compare_val_exc30_allcovs.png"), height = 500, width = 1000)
create_raincloud(exceed30, xvar = "exceed", yvar = "excprob", 
                 xlab = paste("Observed delay >", exc, "days"),
                 ylab = paste("Fitted probability of delay >", exc, "days"), 
                 col_by = "exceed", 
                 drop_na = TRUE) +
  facet_wrap(~model) 
dev.off()

#------------------------------------------------------------------------------#
# Map across prediction points

# Identify indices which correspond to validation points
index.p <- inla.stack.index(stack = stk.full, tag = "pred")$data

# Extract summary stats of fitted values at these indices
mean <- fit.spde$summary.fitted.values[index.p, "mean"]
ll <- fit.spde$summary.fitted.values[index.p, "0.025quant"]
ul <- fit.spde$summary.fitted.values[index.p, "0.975quant"]

pred <- data.frame(x = coop[,1], y = coop[,2], Mean = mean, Lower = ll, Upper = ul, iqr = ul - ll) 

ggmap(bh_lines, 
      base_layer = ggplot(data = pred)) + # , alpha = 1/iqr , col = Mean
  geom_tile(aes(x = x, y = y, fill = Mean), alpha = 0.8) +
  scale_fill_viridis_c(option = "viridis", direction = -1) +
  gg(boundary.spdf, aes(fill = NULL)) +
  # geom_sf(data = dat.fit, pch = "+") +
  # scale_color_viridis_c(option = "viridis", direction = -1, trans = "log2") +
  labs(title = "Predicted mean of days fever prior to diagnosis", 
       x = "", y = "", fill = "Delay (days)") 

ggsave(here::here(figdir,"pred_allcovs_spde.png"), height = 6, width = 9, units = "in")

pred.long <- pred %>%
  tidyr::pivot_longer(-c("x","y","iqr"))

ggmap(bh_lines, 
      base_layer = ggplot(data = pred.long)) +
  geom_tile(aes(x = x, y = y, fill = value)) +
  gg(boundary.spdf, aes(fill = NULL)) +
  # geom_sf(data = dat.fit, pch = "+") +
  scale_fill_viridis_c(option = "viridis", direction = -1) +
  labs(title = "Fitted mean and 2.5-97.5 quantiles for days fever prior to diagnosis", 
       x = "", y = "", fill = "Delay (days)") +
  facet_wrap(~name) +
  coord_fixed(ratio = 1)

ggsave(here::here(figdir,"pred_quants_allcovs_spde.png"), height = 6, width = 9, units = "in")

#------------------------------------------------------------------------------#
# Exceedance probabilities

# Extract fitted marginals
# Calculate probability of exceeding 30 days from this marginal distribution
pred$excprob <- sapply(fit.spde$marginals.fitted.values[index.p],
                       FUN = function(marg){1 - inla.pmarginal(q = 30, marginal = marg)})

ggmap(bh_lines, 
      base_layer = ggplot(data = pred)) +
  geom_tile(aes(x = x, y = y, fill = excprob)) +
  gg(boundary.spdf, aes(fill = NULL)) +
  # geom_sf(data = dat.fit, pch = "+") +
  scale_fill_viridis_c(direction = 1, option = "plasma") + #trans = "log10", 
  labs(title = "Predicted probability of delay exceeding 30 days",
       x = "", y = "", fill = "P(delay > 30)") +
  coord_fixed(ratio = 1) -> map_exc
map_exc

ggsave(here::here(figdir,"exc30_allcovs_spde.png"), map_exc, height = 6, width = 9, units = "in")

# ---------------------------------------------------------------------------- #
# Backward covariate selection based on DIC increase of at least Delta

select <- ggregplot::INLAModelSel("days_fever", c(covs_pat, covs_vil), 
                                  Random = "id", RandomModel = "iid", 
                                  Family = "poisson", 
                                  Data = dat.fit.df, 
                                  Delta = 2)

covs.select <- select$Removed[[length(select$Removed)]]
# "age_s"          "detection"      "num_conslt_cat"

# Same covariates selected when including village effects
