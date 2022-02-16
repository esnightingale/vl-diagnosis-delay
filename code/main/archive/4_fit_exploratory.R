################################################################################
# Description: Univariate regression on covariates and assessment of empirical 
# variogram for spatial dependence
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit/exploratory/poisson"
outdir <- "output/exploratory/poisson"

dat <- readRDS(here::here("data/analysis","dat_nona.rds"))  %>%
  mutate(delay_cat = gtools::quantcut(delay, 5)) 

mesh <- readRDS(here::here("data/analysis","mesh.rds")) 
stk <- readRDS(here::here("data/analysis","stack.rds"))

blockmap <- readRDS(here::here("data","geography","blockmap.rds"))
boundary <- readRDS(here::here("data","geography","boundary.rds"))
boundary.spdf <- as_Spatial(boundary)

# ---------------------------------------------------------------------------- #
# Investigate raw spatial dependence 

idx <- stk$data$index$train

vgm <- variogram(delay ~ 1, #log(days_fever)
                 dat,
                     cressie = TRUE)

fit.vgm <- fit.variogram(vgm, vgm("Mat"), fit.kappa = FALSE)
fit.vgm

# model    psill   range kappa
# 1   Nug 331.6658  0.0000   0.0
# 2   Mat 348.8682 49.9397   0.5


plot(vgm, model = fit.vgm,
     xlab = "Distance (km)", main = "Unadjusted delay beyond 14 days") -> plot.vgm.raw 
plot.vgm.raw

ggsave(here::here("figures/descriptive/vgm_delay.png"), height = 6, width = 8, units = "in")

# ---------------------------------------------------------------------------- #
# Fit baseline IID models with different priors

family <- "poisson"

make_f <- function(u){
  f <- as.formula(sprintf('y ~ -1 + Intercept + 
                                f(id, model = "iid", prior = "pc.prec",
                                  param = c(%s, 0.01))', u))
  return(f)
}

# PC prior on precision of IID effect: P[1/sqrt(prec) > u] = alpha
# Upper limit u is on sd scale

u <- list(`0.001` = 0.001, `0.01` = 0.001,`0.1` = 0.1, `1` = 1, `2` = 2, `5` = 5, `10` = 10)
f.list <- lapply(u, make_f)

# Fit models
fits.base <- lapply(f.list, init_inla, data.stack = stk, family = family, cpo = FALSE)

plyr::llply(fits.base, function(x) x$fit$waic$waic)
# $`0.001`
# [1] 29511.74
# 
# $`0.01`
# [1] 29511.8
# 
# $`0.1`
# [1] 28595.03
# 
# $`1`
# [1] 28586.54
# 
# $`2`
# [1] 28586.03
# 
# $`5`
# [1] 28586.04
#
# $`10`
# [1] 28585
# 
# $`100`
# [1] 28585.34

plyr::llply(fits.base, function(x) x$fit$summary.hyperpar)
# $`0.001`
# mean         sd 0.025quant 0.5quant 0.975quant     mode
# Precision for id 2.113357 0.03534019   2.040762 2.113707   2.182083 2.111306
# 
# $`0.01`
# mean        sd 0.025quant 0.5quant 0.975quant   mode
# Precision for id 2.113436 0.0352481    2.04146 2.113746   2.182108 2.1113
# 
# $`0.1`
# mean         sd 0.025quant 0.5quant 0.975quant     mode
# Precision for id 1.03341 0.02434065  0.9827425 1.033709   1.080984 1.031974
# 
# $`1`
# mean         sd 0.025quant 0.5quant 0.975quant     mode
# Precision for id 1.021218 0.02449294  0.9704758 1.022103    1.06868 1.021758
#
# $`2`
# mean         sd 0.025quant 0.5quant 0.975quant     mode
# Precision for id 1.020791 0.02492144   0.969656 1.021552   1.068504 1.020752
# 
# $`5`
# mean         sd 0.025quant 0.5quant 0.975quant     mode
# Precision for id 1.020278 0.02466197   0.969232 1.020721   1.068006 1.019074
#
# $`10`
# mean         sd 0.025quant 0.5quant 0.975quant     mode
# Precision for id 1.019686 0.02594099   0.969133 1.021609   1.068966 1.026123
# 
# $`100`
# mean         sd 0.025quant 0.5quant 0.975quant     mode
# Precision for id 1.019771 0.02518603  0.9689675 1.021095   1.067851 1.021942

# pc.prec with u >= 1 very similar, better on WAIC than smaller limits.

# ---------------------------------------------------------------------------- #
# Fit baseline models with SPDE and IID to explain spatial structure

# Define formulae
f1 <- y ~ -1 + Intercept + f(id, model = 'iid',
                             prior = "pc.prec",
                             param = c(10, 0.01)) 
f2 <- y ~ -1 + Intercept + f(s, model = spde) 
f3 <- y ~ -1 + Intercept + f(id, model = 'iid',
                               prior = "pc.prec",
                               param = c(10, 0.01)) + f(s, model = spde) 
f4 <- y ~ -1 + Intercept + f(v, model = 'iid',
                             prior = "pc.prec",
                             param = c(10, 0.01))
f5 <- y ~ -1 + Intercept + f(v, model = 'iid',
                             prior = "pc.prec",
                             param = c(10, 0.01)) + f(s, model = spde) 

f.list <- list(IID = f1, SPDE = f2, IID_SPDE = f3, vIID = f4, vIID_SPDE = f5)

# Fit models
fits.base <- lapply(f.list, init_inla, data.stack = stk, family = family, cpo = TRUE)

plyr::llply(fits.base, function(x) x$fit$dic$dic)
 # $IID
 # [1] 28906.32
 # 
 # $SPDE
 # [1] 57468.55
 # 
 # $IID_SPDE
 # [1] 28902.2
 # 
 # $vIID
 # [1] 57313.92
 # 
 # $vIID_SPDE
 # [1] 57313.92

plyr::llply(fits.base, function(x) summary(x$fit))

par(mfrow = c(2,3))
plyr::llply(fits.base, function(x) hist(x$fit$cpo$pit, breaks = 50, ylim = c(0,3), prob = T))

# Plot obs vs fitted
idx <- stk$data$index$train
png(here::here(figdir, "obsvpred_basemodels.png"), height = 1200, width = 2400, res = 150)
inlabru::multiplot(plotlist = plyr::llply(fits.base, function(f) plot_obsvpred(f$fit, stk, idx = idx)), cols = 3)
dev.off()

saveRDS(fits.base, here::here(outdir, "fits_base.rds"))

# ---------------------------------------------------------------------------- #
# Variogram of fitted IID effects

## From IID-only model
hist(fits.base[["IID"]]$fit$summary.random$id$mean, breaks = 30, prob = TRUE)

vgm.iid <- variogram(fits.base[["IID"]]$fit$summary.random$id$mean ~ 1,
                     dat,
                     cressie = TRUE)

fit.vgm.iid <- fit.variogram(vgm.iid, vgm("Mat"), fit.kappa = FALSE)
fit.vgm.iid
# model     psill    range kappa
# 1   Nug 0.1786862     0.00   0.0
# 2   Mat 0.1019022 45498.06   0.5

plot(vgm.iid,
     fit.vgm.iid,
     # ylim = c(0,1),
     xlab = "Distance (km)", main = "Fitted IID effects: IID model") -> plot.vgm.iid
plot.vgm.iid

## From IID+SPDE model
hist(fits.base[["IID_SPDE"]]$fit$summary.random$id$mean, breaks = 30, prob = TRUE)

vgm.iid.spde <- variogram(fits.base[["IID_SPDE"]]$fit$summary.random$id$mean ~ 1,
                          filter(dat, delay >= 0),
                     cressie = TRUE)
fit.vgm.iid.spde <- fit.variogram(vgm.iid.spde, vgm("Mat"), fit.kappa = FALSE)
fit.vgm.iid.spde
# model      psill range kappa
# 1   Nug 0.18353093     0   0.0
# 2   Mat 0.09734597 57303   0.3

plot(vgm.iid.spde,
     fit.vgm.iid.spde,
     # ylim = c(0,1),
     xlab = "Distance (km)", main = "Fitted IID effects: IID + SPDE model"
) -> plot.vgm.iid.spde

## From vIDD model
hist(fits.base[["vIID"]]$fit$summary.random$v$mean, breaks = 30, prob = TRUE)

vgm.viid <- variogram(fits.base[["vIID"]]$fit$summary.random$v$mean ~ 1,
                          unique(dplyr::select(dat, v, geometry)),
                          cressie = TRUE)
fit.vgm.viid <- fit.variogram(vgm.viid, vgm("Mat"), fit.kappa = FALSE)
fit.vgm.viid
# model     psill    range kappa
# 1   Nug 0.5328215  0.00000   0.0
# 2   Mat 0.1060925 13.27141   0.5

plot(vgm.viid,
     fit.vgm.viid,
     # ylim = c(0,1),
     xlab = "Distance (km)", main = "Fitted IID effects: vIID model"
) -> plot.vgm.viid

## From vIID+SPDE model
hist(fits.base[["vIID_SPDE"]]$fit$summary.random$v$mean, breaks = 30, prob = TRUE)

vgm.viid.spde <- variogram(fits.base[["vIID_SPDE"]]$fit$summary.random$v$mean ~ 1,
                           unique(dplyr::select(dat, v, geometry)),
                          cressie = TRUE)
fit.vgm.viid.spde <- fit.variogram(vgm.viid.spde, vgm("Mat"), fit.kappa = FALSE)
fit.vgm.viid.spde
# model     psill    range kappa
# 1   Nug 0.5328215  0.00000   0.0
# 2   Mat 0.1060925 13.27141   0.5

plot(vgm.viid.spde,
     fit.vgm.viid.spde,
     # ylim = c(0,1),
     xlab = "Distance (km)", main = "Fitted IID effects: vIID + SPDE model"
) -> plot.vgm.viid.spde

png(here::here(figdir, "fits_iid_variograms.png"), height = 800, width = 1200)
gridExtra::grid.arrange(plot.vgm.iid, plot.vgm.iid.spde, plot.vgm.viid, plot.vgm.viid.spde, nrow = 2)
dev.off()


## Residuals from SPDE model
resids <- dat$delay - fitted(fits.base$SPDE$fit)[idx]
hist(resids, breaks = 30, prob = TRUE)

vgm.spde <- variogram(resids ~ 1,
                      dat,
                      cressie = TRUE)

fit.vgm.spde <- fit.variogram(vgm.spde, vgm("Mat"), fit.kappa = FALSE)
fit.vgm.spde
# model      psill range kappa
# 1   Nug 0.18353093     0   0.0
# 2   Mat 0.09734597 57303   0.3

plot(vgm.spde,
     fit.vgm.spde,
     # ylim = c(0,1),
     xlab = "Distance (km)", main = "Residuals: SPDE model"
) -> plot.vgm.spde

## Residuals from vIID + SPDE model
resids <- dat$delay - fitted(fits.base$vIID_SPDE$fit)[idx]
hist(resids, breaks = 30, prob = TRUE)

vgm.viid.spde <- variogram(resids ~ 1,
                           dat,
                           cressie = TRUE)

fit.vgm.viid.spde <- fit.variogram(vgm.viid.spde, vgm("Mat"), fit.kappa = FALSE)
fit.vgm.viid.spde
# model      psill range kappa
# 1   Nug 0.18353093     0   0.0
# 2   Mat 0.09734597 57303   0.3

plot(vgm.viid.spde,
     fit.vgm.viid.spde,
     # ylim = c(0,1),
     xlab = "Distance (km)", main = "Residuals: vIID + SPDE model"
) -> plot.vgm.viid.spde


# ---------------------------------------------------------------------------- #
# Variogram of fitted residuals

get_inla_resid <- function(res, observed){
  df <- data.frame(predicted = res$fit$summary.fitted.values$mean[1:length(observed)], 
                   lower = res$fit$summary.fitted.values$`0.025quant`[1:length(observed)], 
                   upper = res$fit$summary.fitted.values$`0.975quant`[1:length(observed)], 
                   observed = observed) %>%
    dplyr::mutate(residual = predicted - observed,
                  standardResidual = residual/stats::sd(residual),
                  standardUpper = (upper - observed)/stats::sd(residual),
                  standardLower = (lower - observed)/stats::sd(residual))
  return(df)
}

resids <- plyr::llply(fits.base, get_inla_resid, observed = dat.prj$delay)
hist(resids[["IID"]]$standardResidual, breaks = 30, prob = TRUE)
hist(resids[["IID_SPDE"]]$standardResidual, breaks = 30, prob = TRUE)

vgm.resid1 <- variogram(resids[["IID"]]$standardResidual ~ 1, 
                        dat.prj, 
                        cressie = TRUE)
fit.vgm.resid1 <- fit.variogram(vgm.resid1, vgm("Mat"))
fit.vgm.resid1
# model     psill    range kappa
# 1   Nug 0.3429789     0.00   0.0
# 2   Mat 0.2206443 38150.47   0.5

plot(vgm.resid1, model = fit.vgm.resid1, 
     ylim = c(0,0.6),
     xlab = "Distance (m)",
     main = "Residuals: IID model") -> plot.vgm.resid1
plot.vgm.resid1

vgm.resid2 <- variogram(resids[["IID_SPDE"]]$standardResidual ~ 1, 
                        dat.prj, 
                        cressie = TRUE)
fit.vgm.resid2 <- fit.variogram(vgm.resid2, vgm("Mat"))
fit.vgm.resid2
# model     psill    range kappa
# 1   Nug 0.3359520     0.00   0.0
# 2   Mat 0.2542666 40706.37   0.5

plot(vgm.resid2, 
     model = fit.vgm.resid2,
     ylim = c(0,0.6),
     xlab = "Distance (m)",
     main = "Residuals: IID + SPDE model") -> plot.vgm.resid2
plot.vgm.resid2

vgm <- variogram(delay ~ 1, 
                 filter(dat.prj, delay >= 0),
                 cressie = TRUE)

fit.vgm <- fit.variogram(vgm, vgm("Mat"), fit.kappa = FALSE)
plot(vgm, model = fit.vgm, 
     # ylim = c(0,0.5),
     xlab = "Distance (m)", main = "Unadjusted delay beyond 14 days") -> plot.vgm.raw
plot.vgm.raw


png(here::here(figdir, "raw_vs_fitresid_variograms.png"), height = 500, width = 1500)
gridExtra::grid.arrange(plot.vgm.raw, plot.vgm.resid1, plot.vgm.resid2, nrow = 1)
dev.off()

INLAutils::ggplot_inla_residuals(fits.base[["IID"]]$fit, observed = dat.prj$delay)
INLAutils::ggplot_inla_residuals(fits.base[["IID_SPDE"]]$fit, observed = dat.prj$delay)

#------------------------------------------------------------------------------#
# Plot the fitted spatial field

pdf(here::here(figdir, "map_fitted_spdes.pdf"), height = 8, width = 12)
plot_spde(fits.base[["SPDE"]]$fit)
plot_spde(fits.base[["IID_SPDE"]]$fit)
plot_spde(fits.base[["vIID_SPDE"]]$fit)
dev.off()

#------------------------------------------------------------------------------#
# Plot the fitted values

# Identify indices of summary.fitted.values which correspond to the predictions
# via inla.stack.index tagged with "pred".
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

coop_proj <- coop %>%
  st_multipoint() %>%
  st_sfc(crs = 7759) %>%
  st_transform(4326) %>%
  st_coordinates()

plot_fitted <- function(res){
  
  fit <- res$fit
  
  # Extract summary stats of fitted values at these indices
  med <- fit$summary.fitted.values[index, "0.5quant"]
  ll <- fit$summary.fitted.values[index, "0.025quant"]
  ul <- fit$summary.fitted.values[index, "0.975quant"]
  
  pred <- data.frame(x = coop_proj[,1], y = coop_proj[,2], 
                     med = med, 
                     ll = ll, ul = ul) %>%
    pivot_longer(-x:-y)
  
  # Map fitted value + quantiles
  # ggmap(bh_lines, 
  #       base_layer = ggplot(data = pred,
  #                           aes(x = x, y = y, col = value))) +
  ggplot(data = pred,
         aes(x = x, y = y, col = value)) +
    # gg(boundary.spdf, fill = "transparent") +
    geom_point(alpha = 0.5) +
    scale_colour_viridis_c(direction = -1, option = "plasma") +
    labs(x = "", y = "", col = "Median") +
    coord_fixed(ratio = 1) +
    facet_wrap(~name) -> p
  
  return(p)
}

plot_fit_space <- plyr::llply(fits.base, plot_fitted)

plot_fit_space[[1]]

ggsave(here::here(figdir,"base_fit.png"), height = 3, width = 10, units = "in")

# ---------------------------------------------------------------------------- #
# Plot expected days prior to diagnosis (NB2 model)

t.index <- inla.stack.index(stack = stk, tag = "train")$data

margs <- fits.base$SPDE$fit$marginals.fitted.values[t.index]

tmargs <- lapply(margs,
                    function(marg) inla.tmarginal(function(p) (1-p)/p, marg))

exp_delay <- lapply(tmargs, 
                    function(marg) inla.emarginal(function(m) m, marg))

hpd_delay <- lapply(tmargs, 
                    function(marg) inla.hpdmarginal(0.95, marg))

data.frame(obs = stk$data$data$y[t.index],
           E = unlist(exp_delay),
           low = hpd.df$low,
           high = hpd.df$high) %>%
  ggplot(aes(obs, E, ymin = low, ymax = high)) +
  geom_point(alpha = 0.3) +
  geom_errorbar(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, col = "grey") +
  geom_smooth() +
  scale_x_continuous(trans = "sqrt") +
  scale_y_continuous(trans = "sqrt")

# ---------------------------------------------------------------------------- #
# Correlation between all covariates

vars <- dat %>%
  st_drop_geometry() %>%
  dplyr::select(-longitude, -latitude, -v, -id, -vill_inc_2017_t, -traveltime) %>%
  dplyr::mutate(across(everything(), as.numeric))

M <- cor(vars)
corrplot::corrplot(M, method = "color", type = "lower", order = "hclust")

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# matrix of the p-value of the correlation
p.mat <- cor.mtest(vars)

png(here::here("figures","descriptive", "covariate effects", "corrplot.png"), 
    height = 8, width = 8, units = "in", res = 300)
corrplot::corrplot(M, method = "color", type = "lower", order = "hclust", 
                   p.mat = p.mat, sig.level = 0.01)
dev.off()

################################################################################