################################################################################
# Description: Univariate regression on covariates and assessment of empirical 
# variogram for spatial dependence
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit/exploratory"
outdir <- "output/exploratory"

dat <- readRDS(here::here("data/analysis","dat_nona.rds"))  %>%
  filter(delay >= 0) %>%
  st_transform(7759)

# Define week/month aggregation
dat <- dat %>%
  mutate(delay_mth = round(delay/30, 0),
         delay_wk = round(delay/7, 0)) 


mesh <- readRDS(here::here("data/analysis","mesh.rds")) 
spde <- readRDS(here::here("data/analysis","spde.rds"))

# ---------------------------------------------------------------------------- #
# Setup data stack

# Try to fit the data to a NB model
params <- MASS::fitdistr(dat$delay, "Negative Binomial")
params

par(mfrow = c(1,2))
plot(qnbinom(ppoints(dat$delay), size=1.11, mu=31.1), sort(dat$delay))
abline(0,1)
# Ok at small values but upper tail much more heavy than an NB
# => over dispersion

params <- MASS::fitdistr(dat$delay, "geometric")
params
plot(qgeom(ppoints(dat$delay), prob = 0.0312), sort(dat$delay))
abline(0,1)
# Geometric actually marginally better it seems 


summary(dat$delay_wk)
summary(dat$delay_mth)

ggplot(dat) +
  geom_bar(aes(x = delay_wk)) +
  labs(x = "Delay (weeks)", y = "Frequency", title = "Delay beyond minimum diagnois criteria of 14 days fever")

ggplot(dat) +
  geom_bar(aes(x = delay_mth)) +
  labs(x = "Delay (months)", y = "Frequency", title = "Delay beyond minimum diagnois criteria of 14 days fever")


par(mfrow = c(1,2))
params <- MASS::fitdistr(dat$delay_mth, "Negative Binomial")
plot(qnbinom(ppoints(dat$delay_mth), size=params$estimate[1], 
             mu=params$estimate[2]), 
     sort(dat$delay_mth))
abline(0,1)

params <- MASS::fitdistr(dat$delay_mth, "geometric")
plot(qgeom(ppoints(dat$delay_mth), prob = params$estimate[1]), 
     sort(dat$delay_mth))
abline(0,1)

# Again geometric quantiles fit better because the distribution allows more weight
# in the upper tail. However, the NB seems to fit the lower end of the distribuion 
# more closely

hist(dat$delay_mth, prob = TRUE, breaks = 20)
lines(0:15 + 0.5, dnbinom(0:15, size = 0.4, mu = 0.5))
hist(dat$delay_mth, prob = TRUE, breaks = 20)
lines(0:15 + 0.5, dgeom(0:15, prob = 0.47))

# ---------------------------------------------------------------------------- #
# Make data stack 

# Covariates of interest
covs <- c("age_s","sex","hiv","prv_tx",
          "marg_caste","occupation",
          "block_endm_2017", "IRS_2017","inc_2017_gt0", 
          "traveltime_s", "traveltime_t_s", "rain",
          "detection")

# Generate the index set for SPDE
indexs <- inla.spde.make.index("s", spde$n.spde)
lengths(indexs)

# Generate a projection matrix A 
coo <- st_coordinates(dat)
A <- inla.spde.make.A(mesh = mesh, loc = coo)

# Define model matrix 
X <- model.matrix(as.formula(paste("~ ",paste(covs, collapse = " + "))), 
                  data = dat)[,-1] 

stk <- inla.stack(
  data = list(y = dat$delay_mth),
  A = list(A, 1, 1, 1),
  effects = list(s = indexs,  # the spatial index,
                 v = dat$v, # village index
                 id = dat$id, # individual index
                 data.frame(  # covariates
                   Intercept = 1, 
                   X)
  )
)

# ---------------------------------------------------------------------------- #
# Simplest model: individual level IID random effects
# Does this capture variation in the observed data?

# Formulae with different prior standard deviations
f1 <- y ~ -1 + Intercept + f(id, model = 'iid',
                             prior = "pc.prec",
                             param = c(1, 0.01)) 
f2 <- y ~ -1 + Intercept + f(id, model = 'iid',
                             prior = "pc.prec",
                             param = c(10, 0.01)) 
f3 <- y ~ -1 + Intercept + f(id, model = 'iid',
                             prior = "pc.prec",
                             param = c(50, 0.01)) 

f.list <- list(low = f1, mid = f2, high = f3)

# Fit models
fits.iid <- lapply(f.list, init_inla, data.stack = stk, family = "nbinomial", cpo = TRUE)

plyr::llply(fits.iid, function(x) x$fit$waic$waic)
# $low
# [1] 12166.45
# 
# $mid
# [1] 12164.79
# 
# $high
# [1] 12164.68

plyr::llply(fits.iid, function(x) summary(x$fit))

par(mfrow = c(3,1))
plyr::llply(fits.iid, function(x) hist(x$fit$cpo$pit, breaks = 30, prob = T))
# Of course indiv IID only model will be bad on PIT as this is a LOO measure and each 
# IID effect doesn't inform on any other observation (so only predict avg via intercept)

# Plot obs vs fitted
inlabru::multiplot(plotlist = plyr::llply(fits.iid, function(f) plot_obsvpred(f$fit, stk, idx = 1:nrow(dat))))

mu_est <- fits.iid$mid$fit$summary.fitted.values$mean[1:nrow(dat)]
prec_est <- fits.iid$mid$fit$summary.fitted.values$mean[1:nrow(dat)]

resids <- (dat$delay_mth - mu_est)/
size_est <- fits.iid$mid$fit$summary.hyperpar$mean[1]

par(mfrow = c(1,1))
plot(dat$delay_mth, resids)
abline(h = 0)

qqplot(qnbinom(ppoints(dat$delay_mth), 
               size = size_est,
               mu = mu_est),
       resids,
       xlab = "Theoretical quantile", ylab = "residuals")
qqline(resids, lty = 2)

plot(qnbinom(ppoints(dat$delay_mth), size = size_est,
           mu = mu_est), 
     sort(dat$delay_mth))
abline(0,1)

# fit.null <- init_inla(y ~ 1, data.stack = stk, family = "nbinomial", cpo = TRUE)
# summary(fit.null$fit)
# 
# mu_est <- fit.null$fit$summary.fitted.values$mean[1:nrow(dat)]
# resids <- dat$delay_mth - mu_est
# size_est <- fit.null$fit$summary.hyperpar$mean[1]
# 
# par(mfrow = c(1,1))
# qqplot(qnbinom(ppoints(dat$delay_mth), 
#                size = size_est,
#                mu = mu_est),
#        resids,
#        xlab = "Theoretical quantile", ylab = "residuals")
# qqline(resids, lty = 2)

# ---------------------------------------------------------------------------- #
# Compare with geometric
fit.iid2 <- init_inla(f2, data.stack = stk, family = "nbinomial2", cpo = TRUE)

fit.iid2$fit$waic$waic
# 14860.08

summary(fit.iid2$fit)
hist(fit.iid2$fit$cpo$pit, breaks = 30)

plot_obsvpred(fit.iid2$fit, stk, idx = 1:nrow(dat))

x <- qgeom(ppoints(dat$delay_mth), 
           prob=fit.iid2$fit$summary.fitted.values[1:nrow(dat),"mean"])
y <- sort(dat$delay_mth)
plot(x, y)
abline(0, 1)

df <- data.frame(x = qgeom(ppoints(dat$delay_mth), 
                           prob=fit.iid2$fit$summary.fitted.values[1:nrow(dat),"mean"]),
                 y = sort(dat$delay_mth))

ggplot(df, aes(x, y)) +
  geom_jitter() +
  geom_abline(0, 1)


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
fits.viid <- lapply(f.list, init_inla, data.stack = stk, family = "nbinomial", cpo = TRUE)

plyr::llply(fits.viid, function(x) x$fit$waic$waic)
# $low
# [1] 12447.59
# 
# $mid
# [1] 12446.08
# 
# $high
# [1] 12445.35

plyr::llply(fits.viid, function(x) summary(x$fit))

par(mfrow = c(3,1))
plyr::llply(fits.viid, function(x) hist(x$fit$cpo$pit, breaks = 30, ylim = c(0,3), prob = T))

# Plot obs vs fitted
inlabru::multiplot(plotlist = plyr::llply(fits.viid, function(f) plot_obsvpred(f$fit, stk, idx = 1:nrow(dat))))

# ---------------------------------------------------------------------------- #
# indiv IID with SPDE

f <- y ~ -1 + Intercept + f(id, model = 'iid',
                             prior = "pc.prec",
                             param = c(10, 0.01)) + f(s, model = spde)

fit.iid.spde <- init_inla(f, data.stack = stk, family = "nbinomial", cpo = TRUE)

fit.iid.spde$fit$waic$waic
# 12316.9

summary(fit.iid.spde$fit)
hist(fit.iid.spde$fit$cpo$pit, breaks = 30, prob = T)
plot_obsvpred(fit.iid.spde$fit, stk, idx = 1:nrow(dat))

# ---------------------------------------------------------------------------- #
# Village iid with SPDE

f <- y ~ -1 + Intercept + f(v, model = 'iid',
                            prior = "pc.prec",
                            param = c(10, 0.01)) + f(s, model = spde)

fit.viid.spde <- init_inla(f, data.stack = stk, family = "nbinomial", cpo = TRUE)

fit.viid.spde$fit$waic$waic
# 12316.9

summary(fit.viid.spde$fit)
hist(fit.viid.spde$fit$cpo$pit, breaks = 30, prob = T)
plot_obsvpred(fit.viid.spde$fit, stk, idx = 1:nrow(dat))

################################################################################
################################################################################
