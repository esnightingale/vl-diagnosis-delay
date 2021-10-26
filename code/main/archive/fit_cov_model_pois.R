################################################################################
# Description: Setup and fit geostatistical model to diagnosis delay data from
# case details forms.
# Define spatial and temporal effect structure only, without covariate effects.
################################################################################
################################################################################

# source(here::here("code","setup_env.R"))

figdir <- "figures/fit"
outdir <- "output"

#------------------------------------------------------------------------------#

# Model formula
f_pat <- y ~ 0 + b0 + 
  age_s + sex + hiv + marg_caste + detection + prv_tx_ka + num_conslt_cat + 
  f(s, model = spde) + 
  f(i, model = "iid") 

# f_vil <- y ~ 0 + b0 + 
#   log(vill_inc_2017_t,10) + IRS_2017_1 + block_endm_2017 + traveltime_s +
#   # f(s, model = spde) + 
#   f(i, model = "iid")
# 
# f_vil <- y ~ 0 + b0 + 
#   vill_inc_2017_gt0 + IRS_2017_1 + block_endm_2017 + traveltime_s +
#   # f(s, model = spde) + 
#   f(i, model = "iid") 
# 
# f_all <- y ~ 0 + b0 + 
#   age_s + sex + hiv + marg_caste + detection + prv_tx_ka + num_conslt + 
#   log(vill_inc_2017_t,10) + IRS_2017_1 + block_endm_2017 + 
#   traveltime_s +
#   # f(s, model = spde) + 
#   f(i, model = "iid") 

#------------------------------------------------------------------------------#

# Fit fixed effects (~1min)
res.fixed <- inla(y ~ 0 + b0 + 
                    age_s + sex + hiv + marg_caste + detection + prv_tx_ka + num_conslt_cat,
            family = "poisson",
            control.family = list(link = "log"),
            data = inla.stack.data(stk.full),
            control.predictor = list(
              compute = TRUE, link = 1,
              A = inla.stack.A(stk.full)
            ),
            control.compute = list(dic = TRUE, 
                                   waic = TRUE, 
                                   cpo = TRUE,
                                   config = TRUE),
            control.fixed = list(mean = 0, prec = 0.1, 
                                 mean.intercept = 0, prec.intercept = 0.001),
            verbose = TRUE
)

summary(res.fixed)

# Fixed effects:
#   mean    sd 0.025quant 0.5quant 0.975quant   mode kld
# b0                   2.400 1.826     -1.185    2.400      5.981  2.400   0
# age_s                0.075 0.002      0.071    0.075      0.080  0.075   0
# sexMale              1.192 1.826     -2.393    1.192      4.774  1.192   0
# sexFemale            1.208 1.826     -2.377    1.208      4.789  1.208   0
# hivYes               0.296 0.010      0.276    0.296      0.316  0.296   0
# marg_casteYes        0.048 0.005      0.038    0.048      0.058  0.048   0
# detectionACD        -0.199 0.005     -0.208   -0.199     -0.189 -0.199   0
# prv_tx_kaYes        -0.002 0.008     -0.018   -0.002      0.014 -0.002   0
# num_conslt_cat(1,3]  0.130 0.008      0.115    0.130      0.145  0.130   0
# num_conslt_cat(3,5]  0.425 0.008      0.409    0.425      0.441  0.425   0
# num_conslt_cat(5,8]  0.817 0.010      0.797    0.817      0.837  0.817   0
# 
# Expected number of effective parameters(stdev): 13.50(0.00)
# Number of equivalent replicates : 320.21 
# 
# Deviance Information Criterion (DIC) ...............: 105278.20
# Deviance Information Criterion (DIC, saturated) ....: 81783.43
# Effective number of parameters .....................: 13.50
# 
# Watanabe-Akaike information criterion (WAIC) ...: 105677.16
# Effective number of parameters .................: 436.45
# 
# Marginal log-Likelihood:  -52745.97 


# Fitting (~4mins)
res.spde <- inla(f_pat,
                 family = "poisson",
                 control.family = list(link = "log"),
                 data = inla.stack.data(stk.full),
                 control.predictor = list(
                  compute = TRUE, link = 1,
                  A = inla.stack.A(stk.full)
                 ),
                 control.compute = list(dic = TRUE, 
                                        waic = TRUE, 
                                        cpo = TRUE,
                                        config = TRUE),
                 control.fixed = list(mean = 0, prec = 0.1, 
                                      mean.intercept = 0, prec.intercept = 0.001),
                 verbose = TRUE
)

summary(res.spde)

# autoplot(res)
# INLAutils::ggplot_inla_residuals(res, dat.fit$days_fever)

summary(res.spde$cpo$failure)
pit_hist(res.spde)

saveRDS(res.fixed,here::here(outdir,"fit_pat_fixed.rds"))
saveRDS(res.spde,here::here(outdir,"fit_pat_spde.rds"))

# ---------------------------------------------------------------------------- #
# Assess spatial independence of fitted IID effects, after accounting for patient
# covariates

# Identify indices which correspond to fit
index.e <- inla.stack.index(stack = stk.full, tag = "est")$data

# Calculate residuals from fixed effect model
dat.fit$resid.fixed <- res.fixed$summary.fitted.values$mean[index.e] - dat.fit$days_fever

vg <- variogram(resid.fixed~1, data = dat.fit)
plot(vg)

vgmod <- vgm(psill = 700, model = "Sph", nugget = 800, range = 100)
plot(vg, model = vgmod)

vgfit <- fit.variogram(vg, model = vgmod)    
plot(vg, model = vgfit)

vgfit
# model    psill    range
# 1   Nug 816.8027  0.00000
# 2   Sph 668.4870 87.72689

png(here::here(figdir, "inla_fixed_semivariogram.png"), height = 500, width = 600)
plot(vg, model = vgfit)
dev.off()

#------------------------------------------------------------------------------#
# Plot the fitted spatial field

rang <- apply(mesh$loc[, c(1, 2)], 2, range)

proj <- inla.mesh.projector(mesh, 
                            xlim = rang[, 1], 
                            ylim = rang[, 2], 
                            dims = c(300, 300))

mean_i <- inla.mesh.project(proj, res$summary.random$i$mean[index.e])
sd_i <- inla.mesh.project(proj, res$summary.random$i$sd[index.e])

df <- expand.grid(x = proj$x, y = proj$y)
df$mean_i <- as.vector(mean_i)
df$sd_i <- as.vector(sd_i)

gmean <- ggplot(df, aes(x = x, y = y, fill = mean_i)) + 
  geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + 
  theme_bw()

gsd <- ggplot(df, aes(x = x, y = y, fill = sd_i)) + 
  geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + 
  theme_bw()

plot_grid(gmean, gsd)

#------------------------------------------------------------------------------#

# Plot the results

# Identify indices of summary.fitted.values which correspond to the predictions
# via inla.stack.index tagged with "pred".
index.p <- inla.stack.index(stack = stk.full, tag = "pred")$data

# Extract summary stats of fitted values at these indices
mean <- res$summary.fitted.values[index.p, "mean"]
ll <- res$summary.fitted.values[index.p, "0.025quant"]
ul <- res$summary.fitted.values[index.p, "0.975quant"]

pred <- data.frame(x = coop[,1], y = coop[,2], Mean = mean, Lower = ll, Upper = ul, iqr = ul - ll) 

ggmap(bh_lines, 
      base_layer = ggplot(data = pred,
                          aes(x = x, y = y, fill = Mean))) + # , alpha = 1/iqr , col = Mean
  geom_tile(alpha = 0.6) +
  scale_fill_viridis_c(option = "viridis", direction = -1) +
  # scale_color_viridis_c(option = "viridis", direction = -1, trans = "log2") +
  labs(title = "Predicted mean of days fever prior to diagnosis", 
       x = "", y = "", fill = "Delay (days)") +
  coord_fixed(ratio = 1)

ggsave(here::here(figdir,"cov_fit_pois.png"), height = 6, width = 9, units = "in")


pred.long <- pred %>%
  tidyr::pivot_longer(-c("x","y","iqr"))

ggmap(bh_lines, 
      base_layer = ggplot(data = pred.long,
                          aes(x = x, y = y, fill = value))) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = -1) +
  labs(title = "Fitted mean and 2.5-97.5 quantiles for days fever prior to diagnosis", 
       x = "", y = "", fill = "Delay (days)") +
  facet_wrap(~name) +
  coord_fixed(ratio = 1)

ggsave(here::here(figdir,"cov_fit_pois_quants.png"), height = 6, width = 9, units = "in")

#------------------------------------------------------------------------------#
# Exceedance probabilities

# Again at prediction locations
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

# Extract fitted marginals
# Calculate probability of exceeding 30 days from this marginal distribution
pred$excprob <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){1 - inla.pmarginal(q = 30, marginal = marg)})

ggmap(bh_lines, 
      base_layer = ggplot(data = pred,
                          aes(x = x, y = y, fill = excprob))) +
  geom_tile() +
  scale_fill_viridis_c(direction = 1, option = "plasma") + #trans = "log10", 
  labs(title = "Predicted probability of delay exceeding 30 days",
       x = "", y = "", fill = "P(delay > 30)") +
  coord_fixed(ratio = 1) -> map_exc
map_exc

ggsave(here::here(figdir,"cov_pois_exc30.png"), map_exc, height = 6, width = 9, units = "in")


################################################################################
################################################################################
