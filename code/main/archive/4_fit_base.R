################################################################################
# Description: Univariate regression on covariates and assessment of empirical 
# variogram for spatial dependence
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit/exploratory"
outdir <- "output/exploratory"

dat <- readRDS(here::here("data/analysis","dat_nona.rds")) 
mesh <- readRDS(here::here("data/analysis","mesh.rds"))
spde <- readRDS(here::here("data/analysis","spde.rds"))

# Setup map context
blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
  sf::st_transform(7759)

boundary <- blockmap %>%
  sf::st_union()

boundary.spdf <- as_Spatial(boundary)

################################################################################
# Define data stack
################################################################################

#------------------------------------------------------------------------------#
# Generate a projection matrix A - projects the spatially continuous GRF from
# the observations to the mesh nodes

# For training points
coo <- st_coordinates(dat)
A <- inla.spde.make.A(mesh = mesh, loc = coo)

dim(A)
nrow(dat)

# Define model matrix based on all covariates of interest, removing automatic 
# intercept
X <- model.matrix(as.formula(paste("~ ",paste(covs, collapse = " + "))), 
                  data = dat)[,-1] 

# Training stack
stk.train <- inla.stack(
  tag = "train",
  data = list(y = dat$days_fever),
  A = list(A, 1, 1),
  effects = list(s = indexs,  # the spatial index,
                 v = dat$v,
                 data.frame(  # covariates
                   Intercept = 1, 
                   X)
  )
)

# Prediction points
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)

# Stack for smooth prediction from intercept and fitted spatial field (no covariates)
stk.pred <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(Ap, 1),
  effects = list(s = indexs,
                 data.frame(
                   Intercept = rep(1, nrow(coop)))
  )
)

stk.full <- inla.stack(stk.train, stk.pred)


################################################################################
################################################################################
# ---------------------------------------------------------------------------- #
# Fit baseline models with SPDE and IID to explain spatial structure

# Define formulae
f1 <- y ~ -1 + Intercept + f(v, model = 'iid',
                             prior = "pc.prec", 
                             param = c(1, 0.01)) 
f2 <- y ~ -1 + Intercept + f(s, model = spde)
f3 <- y ~ -1 + Intercept + f(v, model = 'iid',
                             prior = "pc.prec", 
                             param = c(1, 0.01)) +  f(s, model = spde)
f.list <- list(IID = f1, SPDE = f2, BOTH = f3)

# Fit models
fits.base <- lapply(f.list, init_inla, data.stack = stk.full, family = "nbinomial")

plyr::llply(fits.base, function(x) x$fit$dic$dic)
# $IID
# [1] 38651.08
# 
# $SPDE
# [1] 38762.94

plyr::llply(fits.base, function(x) summary(x$fit))

saveRDS(fits.base, here::here(outdir, "fits_base.rds"))

# ---------------------------------------------------------------------------- #

f4 <- y ~ -1 + Intercept + age_s + sexFemale + hivYes + prv_txYes + detectionACD +
  f(v, model = 'iid', prior = "pc.prec", param = c(1, 0.01)) +  
  f(s, model = spde)

fit_patient <- inla(f4,
                    family = "nbinomial",
                    data = inla.stack.data(stk.full),
                    control.predictor = list(
                      compute = TRUE, link = 1,
                      A = inla.stack.A(stk.full)),
                    control.compute = list(dic = TRUE, 
                                           waic = TRUE, 
                                           config = TRUE),
                    control.fixed = list(mean = 0, 
                                         prec = 0.1, 
                                         mean.intercept = 0, 
                                         prec.intercept = 0.1),
                    verbose = TRUE)

summary(fit_patient)

# ---------------------------------------------------------------------------- #
# Variogram of fitted IID effects

## From IID-only model
# hist(fits.base[["IID"]]$fit$summary.random$v$mean, breaks = 30, prob = TRUE)
# 
# vgm.iid <- variogram(fits.base[["IID"]]$fit$summary.random$v$mean ~ 1,
#                      unique(dplyr::select(dat, v, geometry)),
#                      cressie = TRUE)
# 
# fit.vgm.iid <- fit.variogram(vgm.iid, vgm("Mat"), fit.kappa = TRUE)
# fit.vgm.iid
# # model     psill    range kappa
# # 1   Nug 0.1786862     0.00   0.0
# # 2   Mat 0.1019022 45498.06   0.5
# 
# plot(vgm.iid,
#      fit.vgm.iid,
#      ylim = c(0,0.35),
#      xlab = "Distance (m)", main = "Fitted IID effects: IID model") -> plot.vgm.iid
# 
# 
# ## From IID+SPDE model
# hist(fits.base[["SPDE"]]$fit$summary.random$id$mean, breaks = 30, prob = TRUE)
# 
# vgm.spde <- variogram(fits.base[["SPDE"]]$fit$summary.random$id$mean ~ 1,
#                      dat.fit,
#                      cressie = TRUE)
# fit.vgm.spde <- fit.variogram(vgm.spde, vgm("Mat"), fit.kappa = TRUE)
# fit.vgm.spde
# # model      psill range kappa
# # 1   Nug 0.18353093     0   0.0
# # 2   Mat 0.09734597 57303   0.3
# 
# plot(vgm.spde,
#      fit.vgm.spde,
#      ylim = c(0,0.35),
#      xlab = "Distance (km)", main = "Fitted IID effects: IID + SPDE model"
#      ) -> plot.vgm.spde
# 
# 
# png(here::here(figdir, "raw_vs_fitIID_variograms.png"), height = 500, width = 1500)
# gridExtra::grid.arrange(plot.vgm.raw, plot.vgm.iid, plot.vgm.spde, nrow = 1)
# dev.off()

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

resids <- plyr::llply(fits.base, get_inla_resid, observed = dat$days_fever)
hist(resids[["IID"]]$standardResidual, breaks = 30, prob = TRUE)
hist(resids[["SPDE"]]$standardResidual, breaks = 30, prob = TRUE)

vgm.resid1 <- variogram(resids[["IID"]]$standardResidual ~ 1, 
                        dat, 
                        cressie = TRUE)
fit.vgm.resid1 <- fit.variogram(vgm.resid1, vgm("Mat"))
fit.vgm.resid1
# model     psill    range kappa
# 1   Nug 0.3669700     0.00   0.0
# 2   Mat 0.1671981 53555.08   0.5

plot(vgm.resid1, #model = fit.vgm.resid1, 
     # ylim = c(0,0.5),
     xlab = "Distance (km)",
     main = "Residuals: IID model") -> plot.vgm.resid1
plot.vgm.resid1

vgm.resid2 <- variogram(resids[["SPDE"]]$standardResidual ~ 1, 
                        dat, 
                        cressie = TRUE)
# fit.vgm.resid2 <- fit.variogram(vgm.resid2, vgm("Exp"))
# fit.vgm.resid2
#   model     psill    range kappa
# 1   Nug 0.6812076     0.00   0.0
# 2   Mat 0.2800832 36687.68   0.5

plot(vgm.resid2, 
     # model = fit.vgm.resid2, 
     # ylim = c(0,0.5),
     xlab = "Distance (km)",
     main = "Residuals: SPDE model") -> plot.vgm.resid2
plot.vgm.resid2

png(here::here(figdir, "raw_vs_fitresid_variograms.png"), height = 500, width = 1500)
gridExtra::grid.arrange(plot.vgm.raw, plot.vgm.resid1, plot.vgm.resid2, nrow = 1)
dev.off()

INLAutils::ggplot_inla_residuals(fits.base[["IID"]]$fit, observed = dat$days_fever)
INLAutils::ggplot_inla_residuals(fits.base[["SPDE"]]$fit, observed = dat$days_fever)

#------------------------------------------------------------------------------#
# Plot the fitted spatial field

png(here::here(figdir, "map_fitted_spde.png"), height = 500, width = 1500)
plot_spde(fits.base[["SPDE"]])
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