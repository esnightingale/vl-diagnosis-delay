################################################################################
# Description: Univariate regression on covariates and assessment of empirical 
# variogram for spatial dependence
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit/exploratory"
outdir <- "output/exploratory"

dat <- read_data()

mesh <- readRDS(here::here("data/analysis","mesh.rds")) 
spde <- readRDS(here::here("data/analysis","spde.rds"))
stk <- readRDS(here::here("data/analysis","stack.rds"))

# ---------------------------------------------------------------------------- #
# Simplest model: individual level IID random effects
# Does this capture variation in the observed data?

# Formulae with different prior standard deviations
f1 <- y ~ -1 + Intercept + f(id, model = 'iid',
                             prior = "pc.prec",
                             param = c(0.1, 0.01)) 
f2 <- y ~ -1 + Intercept + f(id, model = 'iid',
                             prior = "pc.prec",
                             param = c(1, 0.01)) 
f3 <- y ~ -1 + Intercept + f(id, model = 'iid',
                             prior = "pc.prec",
                             param = c(10, 0.01)) 

f.list <- list(low = f1, mid = f2, high = f3)

# Fit models
fits.iid <- lapply(f.list, init_inla, data.stack = stk, family = "poisson", cpo = TRUE)


plyr::llply(fits.iid, function(x) x$fit$waic$waic)
# $low
# [1] 28586.93
# 
# $mid
# [1] 28585.23
# 
# $high
# [1] 28585.25

plyr::llply(fits.iid, function(x) summary(x$fit))
# Fitted precision not very sensitive to choice of PC prior

par(mfrow = c(3,1))
plyr::llply(fits.iid, function(x) hist(x$fit$cpo$pit, breaks = 30, ylim = c(0,10), prob = T))
# Of course indiv IID only model will be bad on PIT as this is a LOO measure and each 
# IID effect doesn't inform on any other observation (so only predict avg via intercept)

# Plot obs vs fitted
idx <- stk$data$index$train
png(here::here(figdir, "compare_indivIID_priorprec.png"), height = 900, width = 400)
inlabru::multiplot(plotlist = plyr::llply(fits.iid, function(f) plot_obsvpred(f$fit, stk = stk, idx = idx)))
dev.off()

# ---------------------------------------------------------------------------- #
# Compare likelihoods

# f_covs <- "delay ~ 1 + age_s + sex + hiv + prv_tx + marg_caste + occupation + 
#   detection + block_endm_2017 + IRS_2017 + inc_2017_gt0 + traveltime_s + rain"
# formula(f_covs),

f.iid <- delay ~ 1 + f(id, model = 'iid',
                             prior = "pc.prec",
                             param = c(10, 0.01)) 

fam.list <- list(Po = "poisson", NB = "nbinomial", NB1 = "nbinomial2")
# fits.fam <- lapply(fam.list, function(fam) init_inla(f = f2, data.stack = stk, family = fam, cpo = TRUE))

fits.fam <- lapply(fam.list, function(fam) inla(formula = f.iid,
                                                data = dat,
                                                family = fam,
                                                control.predictor = list(compute = TRUE),
                                                control.compute = list(dic = TRUE,
                                                                       waic = TRUE,
                                                                       config = TRUE,
                                                                       cpo = TRUE,
                                                                       return.marginals.predictor = TRUE)))


plyr::llply(fits.fam, function(x) x$waic$waic)
# $Po
# [1] 28585.97
# 
# $NB
# [1] 36772.69
# 
# $NB1
# [1] 37929.18

plyr::llply(fits.fam, function(x) summary(x))

par(mfrow = c(2,2))
plyr::llply(fits.fam, function(x) hist(x$cpo$pit, breaks = 30))

disp.check <- plyr::llply(fits.fam[c(1,2)], function(x) inlatools::dispersion_check(x))
lapply(disp.check, plot)

dist.check <- inlatools::fast_distribution_check(fits.fam[c(1,2)])
plot(dist.check)

# Plot obs vs fitted

png(here::here(figdir, "compare_indivIID_family.png"), height = 900, width = 400)
inlabru::multiplot(plotlist = plyr::llply(fits.fam,
                                          function(f) plot_obsvpred(f, obs = dat$delay.mth, title = "family")))
# inlabru::multiplot(plotlist = plyr::llply(fits.fam, function(f) plot_obsvpred(f$fit, stk = stk, idx = idx, title = "family")))
dev.off()

# ---------------------------------------------------------------------------- #
# Compare village IID family

f.viid <- delay ~ 1 + f(v, model = 'iid',
                   prior = "pc.prec",
                   param = c(10, 0.01)) 

fits.fam2 <- lapply(fam.list, function(fam) inla(formula = f.viid,
                                                data = dat,
                                                family = fam,
                                                control.predictor = list(compute = TRUE),
                                                control.compute = list(dic = TRUE,
                                                                       waic = TRUE,
                                                                       config = TRUE,
                                                                       cpo = TRUE,
                                                                       return.marginals.predictor = TRUE)))

plyr::llply(fits.fam2, function(x) x$waic$waic)
# $Po
# [1] 77429.63
# 
# $NB
# [1] 37529.04
# 
# $NB1
# [1] 37812.7

plyr::llply(fits.fam2, function(x) summary(x))

par(mfrow = c(2,2))
plyr::llply(fits.fam2, function(x) hist(x$cpo$pit, breaks = 30))

disp.check <- plyr::llply(fits.fam2[c(1,2)], function(x) inlatools::dispersion_check(x))
lapply(disp.check, plot)

dist.check <- inlatools::fast_distribution_check(fits.fam2[c(1,2)])
plot(dist.check, scales = "free")

# Plot obs vs fitted

png(here::here(figdir, "compare_villIID_family.png"), height = 900, width = 400)
inlabru::multiplot(plotlist = plyr::llply(fits.fam2,
                                          function(f) plot_obsvpred(f, obs = dat$delay, title = "family")))
# inlabru::multiplot(plotlist = plyr::llply(fits.fam, function(f) plot_obsvpred(f$fit, stk = stk, idx = idx, title = "family")))
dev.off()

# ---------------------------------------------------------------------------- #
# Binomial

f.list.iid <- list(indiv =  gt30 ~ 1 + f(id, model = 'iid',
                                          prior = "pc.prec",
                                          param = c(10, 0.01)),
                   village =  gt30 ~ 1 + f(v, model = 'iid',
                                            prior = "pc.prec",
                                            param = c(10, 0.01)))

fits.bin <- lapply(f.list.iid, function(f) inla(formula = f,
                                                 data = dat,
                                                 family = "binomial",
                                                 control.predictor = list(compute = TRUE),
                                                 control.compute = list(dic = TRUE,
                                                                        waic = TRUE,
                                                                        config = TRUE,
                                                                        cpo = TRUE)))

plyr::llply(fits.bin, function(x) x$waic$waic)
# $indiv
# [1] 5457.177
# 
# $village
# [1] 5401.781

plyr::llply(fits.bin, function(x) summary(x))

par(mfrow = c(2,1))
plyr::llply(fits.bin, function(x) hist(x$cpo$pit, breaks = 30))

png(here::here(figdir, "compare_IID_binomial.png"), height = 600, width = 400)
inlabru::multiplot(plotlist = plyr::llply(fits.bin,
                                          function(f) plot_obsvpred(f, obs = dat$delay_gt30, title = "formula")))
dev.off()

# ---------------------------------------------------------------------------- #
# Less simple model: village level IID random effects

f1 <- y ~ -1 + Intercept + f(v, model = 'iid',
                             prior = "pc.prec",
                             param = c(1, 0.01)) 
f2 <- y ~ -1 + Intercept + f(v, model = 'iid',
                             prior = "pc.prec",
                             param = c(10, 0.01)) 
f3 <- y ~ -1 + Intercept + f(v, model = 'iid',
                             prior = "pc.prec",
                             param = c(50, 0.01)) 

f.list <- list(low = f1, mid = f2, high = f3)

# Fit models
fits.viid <- lapply(f.list, init_inla, data.stack = stk, family = "poisson", cpo = TRUE)


plyr::llply(fits.viid, function(x) x$fit$waic$waic)
# $low
# [1] 77425.9
# 
# $mid
# [1] 77425.64
# 
# $high
# [1] 77425.8

plyr::llply(fits.viid, function(x) summary(x$fit))

par(mfrow = c(3,1))
plyr::llply(fits.viid, function(x) hist(x$fit$cpo$pit, breaks = 30, ylim = c(0,3), prob = T))

# Plot obs vs fitted
inlabru::multiplot(plotlist = plyr::llply(fits.viid, function(f) plot_obsvpred(f$fit, stk, idx = idx)))

################################################################################

