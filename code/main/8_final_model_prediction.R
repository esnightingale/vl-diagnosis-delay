################################################################################
# Description: 
# 
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit/excl_lt14"
outdir <- "output/excl_lt14"

dat <- readRDS(here::here("data/analysis/Archive","dat_fit.rds")) %>%
  st_transform(crs = st_crs(7759)) %>%
  filter(delay >= 0)
mesh <- readRDS(here::here("data/analysis/Archive","mesh.rds"))
spde <- readRDS(here::here("data/analysis/Archive","spde.rds"))
stack <- readRDS(here::here("data/analysis/Archive","stack_excl_lt14.rds"))[[1]]

covs <- c("age_s","sex","hiv",
          "marg_caste","occupation",
          "block_endm_2017", "IRS_2017_1","vill_inc_2017_gt0", "prv_tx",
          "traveltime_s", "rain",
          "detection")

blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
  st_transform(7759)

boundary <- sf::st_union(blockmap)
boundary.spdf <- as_Spatial(boundary)

fits.all <- readRDS(here::here(outdir, "fits_covs_all.rds"))
fit.final <- fits.all[["Patient + village awareness"]]

# Fitted values from final fit
# tridx <- inla.stack.index(stack, "train")$data
# df_pred <- data.frame(m = fit.final$summary.fitted.values[tridx, "mean"],
#                       sd = fit.final$summary.fitted.values[tridx, "sd"],
#                       low = fit.final$summary.fitted.values[tridx, "0.025quant"],
#                       med = fit.final$summary.fitted.values[tridx, "0.5quant"],
#                       hi = fit.final$summary.fitted.values[tridx, "0.975quant"]) %>%
#   mutate(Scenario = "Baseline")

#------------------------------------------------------------------------------#
# Predict in scenario where all cases diagnosed with ACD

indexs <- inla.spde.make.index("s", spde$n.spde)

# For training points
coo <- st_coordinates(dat)
A <- inla.spde.make.A(mesh = mesh, loc = coo)

X1 <- model.matrix(as.formula(paste("~ ",paste(covs, collapse = " + "))), 
                   data = dat)[,-1] 

# Training stack
stk.train <- inla.stack(
  tag = "train",
  data = list(y = dat$delay),
  A = list(A, 1, 1),
  effects = list(s = indexs,  # the spatial index,
                 v = dat$v,
                 data.frame(  # covariates
                   Intercept = 1, 
                   X1)))

# Prediction at same points but with detection = ACD for all
dat.pred <- dat %>%
  dplyr::mutate(detection == "ACD")

X2 <- model.matrix(as.formula(paste("~ ",paste(covs, collapse = " + "))), 
                   data = dat.pred)[,-1] 

# Training stack
stk.pred <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(A, 1, 1),
  effects = list(s = indexs,  # the spatial index,
                 v = dat.pred$v,
                 data.frame(  # covariates
                   Intercept = 1, 
                   X2)
  )
)

stk <- inla.stack(stk.train, stk.pred)

# Refit with prediction stack
fit.pred <- inla(fit.final$f,
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

# ---------------------------------------------------------------------------- #
# Extract predicted delays for baseline and 100% ACD scenarios

tidx <- stk$data$index$train
pidx <- stk$data$index$pred

df_fit <- data.frame(m = fit.pred$summary.fitted.values[tidx, "mean"],
                      sd = fit.pred$summary.fitted.values[tidx, "sd"],
                      low = fit.pred$summary.fitted.values[tidx, "0.025quant"],
                      med = fit.pred$summary.fitted.values[tidx, "0.5quant"],
                      hi = fit.pred$summary.fitted.values[tidx, "0.975quant"]) %>%
  mutate(Scenario = "Baseline")

df_pred <- data.frame(m = fit.pred$summary.fitted.values[pidx, "mean"],
                      sd = fit.pred$summary.fitted.values[pidx, "sd"],
                      low = fit.pred$summary.fitted.values[pidx, "0.025quant"],
                      med = fit.pred$summary.fitted.values[pidx, "0.5quant"],
                      hi = fit.pred$summary.fitted.values[pidx, "0.975quant"]) %>%
  mutate(Scenario = "Complete ACD coverage")

df <- bind_cols(df_fit, df_pred)
df <- bind_rows(df_fit, df_pred)

# ---------------------------------------------------------------------------- #
# Compare distributions

ggplot(df) +
  geom_histogram(aes(m), fill = "white", col = "grey") +
  geom_histogram(aes(low), alpha = 0.2, fill = "steelblue") +
  geom_histogram(aes(hi), alpha = 0.2, fill = "indianred") +
  facet_wrap(~Scenario)

ggplot(df, aes(fill = Scenario), position = "identity") +
  geom_histogram(aes(m), alpha = 0.2) 

df %>%
  group_by(Scenario) %>%
  summarise(mean_delay = mean(m),
            median_delay = median(m))

summary(df_fit$m)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8.615  20.657  26.469  29.461  34.824 138.205 
summary(df_pred$m)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 17.87   34.02   39.93   43.53   48.68  176.54 

sum(dat$days_fever > 90)/nrow(dat)
# 0.434327 
sum(df_pred$m > 90)/nrow(df_pred)
# 0.8907778

################################################################################
################################################################################
