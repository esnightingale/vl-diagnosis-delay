library(survival)
library(INLA)
library(sf)

dat <- filter(dat, delay >= 0) %>%
  mutate(time = delay/30,
         time.rnd = round(time, 0),
         cens = 1)

hist(dat$time.rnd, prob = TRUE, breaks = 20)
lines(0:15 + 0.5, dnbinom(0:15, size = 0.4, mu = 0.5))

km <- survfit(Surv(delay, cens) ~ 1, dat)
plot(km, conf.int = TRUE)

km <- survfit(Surv(time.rnd, cens) ~ 1, dat)
plot(km, conf.int = TRUE)


km <- survfit(Surv(delay, cens) ~ inc_2017_gt0, dat)
plot(km, conf.int = TRUE, col = 2:1) 
legend('topright', c('Village incidence 2017 = 0', 'Village incidence 2017 > 0'), lty = 1, col = 2:1,
       bty = "n") 


# Define stack

mesh <- readRDS(here::here("data/analysis","mesh.rds")) 
spde <- readRDS(here::here("data/analysis","spde.rds"))
loc <- st_coordinates(dat)
A <- inla.spde.make.A(mesh, loc)

covs <- dat %>%
  st_drop_geometry() %>%
  dplyr::select(age_s, sex, hiv, marg_caste, occupation, detection, prv_tx, traveltime_s,
                inc_2017_gt0, block_endm_2017, rain)

stk.surv <- inla.stack(
  data = list(time = dat$delay, cens = dat$cens), 
  A = list(A, 1, 1, 1), 
  effects = list(s = 1:spde$n.spde,
                 v = dat$v,
                 id = dat$id,
                 data.frame(Intercept = 1, covs)))

# Define formulae
covs2 <- list("age_s","sexMale","hivYes","prv_txYes","marg_casteYes","occupationUnskilled",
              "occupationSkilled", "`occupationSalaried.selfemployed`", "detectionACD",
              "block_endm_2017Endemic", "IRS_2017Yes","inc_2017_gt0Yes",
              "traveltime_s","rain") #"traveltime_t_s", 

# f1 <- Surv(time, cens) ~ -1 + Intercept + 
#   f(id, model = 'iid',
#     prior = "pc.prec",
#     param = c(10, 0.01)) 

f1 <- as.formula(paste0(
  "Surv(time, cens) ~ -1 + Intercept +",
  paste0(covs2, collapse = " + "),
  "+ f(id, model = 'iid', prior = 'pc.prec', param = c(1, 0.01)) + 
   f(s, model = spde)"))

r <- inla(f1, family = "exponential.surv", 
          data = inla.stack.data(stk.surv), 
          control.predictor = list(A = inla.stack.A(stk.surv), compute = TRUE)) 


f2 <- Surv(delay, rep(1, nrow(dat))) ~ -1 + Intercept + f(s, model = spde) 
f3 <- Surv(delay, rep(1, nrow(dat))) ~ -1 + Intercept + f(id, model = 'iid',
                             prior = "pc.prec",
                             param = c(1, 0.01) + f(s, model = spde) 
) + f(s, model = spde) 

f.list <- list(IID = f1, SPDE = f2, IID_SPDE = f3)

# Fit models
fits.base <- lapply(f.list, init_inla, data.stack = stk, family = "nbinomial2", cpo = TRUE)











sinla.delay <- inla.surv(dat$delay,  rep(1, nrow(dat)))

exp.delay <- inla(sinla.delay ~ f(v, model = "iid"),
                  data = dat, 
                  family = "exponential.surv",
                  control.predictor = list(
                    compute = TRUE, link = 1))
summary(exp.delay)

library(parallel)
# Set number of cores to use
options(mc.cores = 4)
# Grid of times
times <- seq(0.01, 300, by = 1)
# Marginal of alpha
alpha.marg <- exp.delay$marginals.fixed[["(Intercept)"]]
# Compute post. marg. of survival function
S.inla <- mclapply(times, function(t){
  S.marg <- inla.tmarginal(function(x) {exp(- exp(x) * t)}, alpha.marg)
  S.stats <- inla.zmarginal(S.marg, silent = TRUE)
  return(unlist(S.stats[c("quant0.025", "mean", "quant0.975")]))
})

S.inla <- as.data.frame(do.call(rbind, S.inla)) 
S.inla$t <- times

km.df <- data.frame(t = km.delay$time,
                    surv = km.delay$surv,
                    lower = km.delay$lower,
                    upper = km.delay$upper)

ggplot() +
  geom_ribbon(data = S.inla, aes(x = t, ymin = quant0.025, ymax = quant0.975),
              alpha = 0.5, fill = "red") +
  geom_line(data = S.inla, aes(x = t,y = mean), col = "red") + 
  geom_ribbon(data = km.df, aes(x = t, ymin = lower, ymax = upper), alpha = 0.5) +
  geom_line(data = km.df, aes(t, surv))


data.frame(obs = km.delay$surv,
           Est = exp.delay$summary.fitted.values[,"mean"],
           low = exp.delay$summary.fitted.values[,"0.025quant"],
           high = exp.delay$summary.fitted.values[,"0.975quant"]) %>%
  ggplot(aes(obs, Est, ymin = low, ymax = high)) +
  geom_point(alpha = 0.3) +
  geom_errorbar(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, col = "grey") +
  geom_smooth() +
  scale_x_continuous(trans = "sqrt") +
  scale_y_continuous(trans = "sqrt")
