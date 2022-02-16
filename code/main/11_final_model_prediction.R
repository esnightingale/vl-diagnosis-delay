################################################################################
# Description: 
# 
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/prediction"
outdir <- "output"

mesh <- readRDS(here::here("data/analysis","mesh.rds"))
spde <- readRDS(here::here("data/analysis","spde.rds"))
dat <- read_data() #readRDS(here::here("data/analysis","dat_nona.rds")) 
coop <- readRDS(here::here("data/analysis","coop.rds"))
stk <- readRDS(here::here("data/analysis","stack.rds")) 

covs <- c("age_s","comorb", 
          "poss_acd",
          "block_endm_2017", "inc_2017_gt0",
          "traveltime_t_s")

boundary <- readRDS(here::here("data","geography","boundary.rds")) %>%
  st_sf() 
boundary.spdf <- as_Spatial(boundary)

fit.final <- readRDS(here::here(outdir, "fits_final.rds"))[["All"]]
# fit.final <- fits[["All"]]

#------------------------------------------------------------------------------#
# Posterior predictive distribution
 zzzzpl-98[p]
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

# Stack for prediction with no ACD
X2 <- X1
X2[,"poss_acdTRUE"] <- 0
stk.pred_noACD <- inla.stack(
  tag = "pred_noACD",
  data = list(y = NA),
  A = list(A, 1, 1),
  effects = list(s = indexs,
                 id = dat$id,
                 data.frame(
                   Intercept = 1,
                   X2)
  )
)

# Stack for prediction with complete ACD
X3 <- X1
X3[,"poss_acdTRUE"] <- 1
stk.pred_fullACD <- inla.stack(
  tag = "pred_fullACD",
  data = list(y = NA),
  A = list(A, 1, 1),
  effects = list(s = indexs,
                 id = dat$id,
                 data.frame(
                   Intercept = 1,
                   X3)
  )
)

stk <- inla.stack(stk.train, stk.pred, stk.pred_noACD, stk.pred_fullACD)
saveRDS(stk, here::here("data/analysis","stack_pred.rds"))

# Refit with prediction stack
fit.pred <- inla(fit.final$f,
            family = "poisson",
            data = inla.stack.data(stk),
            control.predictor = list(
              compute = TRUE, link = 1,
              A = inla.stack.A(stk)),
            control.compute = list(dic = FALSE, 
                                   waic = TRUE, 
                                   config = TRUE,
                                   cpo = FALSE,
                                   return.marginals.predictor = TRUE),
            control.fixed = list(mean = 0, 
                                 prec = 0.1, 
                                 mean.intercept = 0, 
                                 prec.intercept = 0.1),
            verbose = TRUE)

saveRDS(fit.pred, here::here("output","fit_pred.rds"))

#------------------------------------------------------------------------------#
# Plot smooth predictions

# Identify indices which correspond to validation points
index.p <- inla.stack.index(stack = stk, tag = "pred")$data

# Extract summary stats of fitted values at these indices
mean <- fit.pred$summary.fitted.values[index.p, "mean"]
ll <- fit.pred$summary.fitted.values[index.p, "0.025quant"]
ul <- fit.pred$summary.fitted.values[index.p, "0.975quant"]

pred <- data.frame(x = coop[,1], y = coop[,2], Mean = mean, Lower = ll, Upper = ul, iqr = ul - ll) 

ggplot() +
  geom_sf(data = boundary, fill = NA) +
  geom_tile(data = pred, aes(x = x, y = y, fill = Mean)) +
  scale_fill_viridis_c(option = "viridis", direction = -1, end = 0.9) +
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9, 0.85)) +
  geom_sf(data = dat, pch = "+", col = "white", alpha = 0.5) +
  labs(title = "Predicted mean of days fever prior to diagnosis", 
       x = "", y = "", fill = "Delay (days)") 

ggsave(here::here(figdir,"pred_final_wdat.png"), height = 6, width = 9, units = "in")

pred.long <- pred %>%
  tidyr::pivot_longer(-c("x","y","iqr"))

ggplot() +
  geom_sf(data = boundary, fill = NA) +
  geom_tile(data = pred.long, aes(x = x, y = y, fill = value)) +
  scale_fill_viridis_c(option = "viridis", direction = -1, end = 0.9) +
  labs(title = "Fitted mean and 2.5-97.5 quantiles for days fever prior to diagnosis", 
       x = "", y = "", fill = "Delay (days)") +
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) +
  facet_wrap(~name) 

ggsave(here::here(figdir,"pred_final_quants.png"), height = 6, width = 9, units = "in")

#------------------------------------------------------------------------------#
# Exceedance probabilities

# Extract fitted marginals
# Calculate probability of exceeding 30 days from this marginal distribution
pred$excprob <- sapply(fit.pred$marginals.fitted.values[index.p],
                       FUN = function(marg){1 - inla.pmarginal(q = 30, marginal = marg)})
pred <- mutate(pred, 
               excprob_hi = case_when(excprob > 0.5 ~ "greater than 0.5",
                                     excprob <= 0.5 ~ "less/equal to 0.5"),
               excprob_strength = abs(excprob - 0.5))
ggplot() +
  geom_tile(data = pred, aes(x = x, y = y, fill = excprob)) +
  geom_sf(data = boundary, fill = NA) +
  scale_fill_viridis_c(direction = -1, option = "cividis") +
  labs(#title = "Predicted probability of delay exceeding 30 days",
       x = "", y = "", fill = "P(delay > 30)") +
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom" #c(0.9, 0.85) #
        ) -> map_exc
map_exc

# ggsave(here::here(figdir,"excprob30_final.png"), map_exc, height = 7, width = 9, units = "in")

ggplot() +
  geom_tile(data = pred, aes(x = x, y = y, fill = excprob_hi, alpha = excprob_strength)) +
  geom_sf(data = boundary, fill = NA) +
  scale_fill_viridis_d(option = "cividis") +
  # scale_fill_viridis_d(option = "viridis", begin = 0, end = 0.9) + 
  labs(#title = "Predicted probability of delay exceeding 30 days",
       x = "", y = "", fill = "") + #P(delay > 30)
  guides(alpha = "none") +
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position ="bottom"# c(0.9, 0.85) #
        ) -> map_exc2
map_exc2

# ggsave(here::here(figdir,"excprob30_cat_final.png"), map_exc2, height = 7, width = 9, units = "in")

map_exc + map_exc2
ggsave(here::here(figdir,"excprob30_final_combined_alt.png"), height = 7, width = 14, units = "in")


#------------------------------------------------------------------------------#
# Predictions for full/no ACD coverage

## Identify indices which correspond to validation points
index.t <- inla.stack.index(stack = stk, tag = "train")$data
index.p1 <- inla.stack.index(stack = stk, tag = "pred_noACD")$data

# Extract summary stats of fitted values at these indices
m <- fit.pred$summary.fitted.values[index.p1, "mean"]
ll <- fit.pred$summary.fitted.values[index.p1, "0.025quant"]
ul <- fit.pred$summary.fitted.values[index.p1, "0.975quant"]

pred2 <- data.frame(delay = dat$delay, 
                    detection = factor(dat$poss_acd, levels = c(FALSE, TRUE), labels = c("PCD","ACD")),
                    x = coo[,1], y = coo[,2], 
                    pred = m, lower = ll, upper = ul, iqr = ul - ll,
                    fitted = fit.pred$summary.fitted.values[index.t, "mean"]) %>%
  mutate(change = (pred - fitted)*100/fitted)

pred2 %>%
  ggplot() +
  geom_point(aes(fitted, pred, col = detection)) +
  geom_abline(lty = "dashed") +
  scale_fill_viridis_d(option = "plasma") +
  labs(x = "Observed delay", y = "Predicted - no ACD")
ggsave(here::here(figdir, "obsvpred_noacd.png"), height = 6, width = 8, units = "in")

pred2 %>%
  ggplot() +
  geom_point(aes(fitted, pred - delay, col = detection)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_viridis_d(option = "plasma") +
  labs(x = "Observed delay", y = "Predicted increase with no ACD")
ggsave(here::here(figdir, "predchng_noacd.png"), height = 6, width = 8, units = "in")

pred2 %>% 
  dplyr::rename(Fitted = fitted, 
                `Predicted - No ACD` = pred) %>%
  tidyr::pivot_longer(c("Fitted","Predicted - No ACD")) %>% 
  ggplot(aes(x = value, fill = name)) +
  geom_histogram(position = "identity", alpha = 0.3, bins = 50) +
  labs(x = "Delay (days)", y = "Count", fill = "") +
  theme(legend.position = c(0.8,0.4))
ggsave(here::here(figdir, "fit_vs_noacd_hist.png"), height = 4, width = 5, units = "in", dpi = 350)

summary(pred2$fitted)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.262  11.854  16.796  31.523  43.970 493.321 
summary(pred2$pred)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.708  12.284  20.835  34.322  45.361 629.468 

pred2 %>% filter(detection == "ACD") %>% pull(fitted) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.262   9.607  16.427  25.291  30.542 493.321 
pred2 %>% filter(detection == "ACD") %>% pull(pred) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.886  12.258  20.961  32.270  38.971 629.468

summary(pred2$pred - pred2$fitted)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.00028   0.00000   0.00001   2.79918   4.20049 136.14710 

pred2 %>%
  pivot_longer(c("fitted","pred")) %>% 
  mutate(gt30 = (value > 30)) %>% 
  group_by(name) %>% 
  summarise(p_gt30 = mean(gt30))
#   name   p_gt30
# 1 fitted  0.338
# 2 pred    0.351

pred2 %>%
  filter(detection == "ACD") %>%
ggplot() +
  geom_histogram(aes(x = (pred - delay)*100/delay)) + 
  labs(x = "% change from observed")
ggsave(here::here(figdir, "predvsobs_noacd_hist.png"), height = 6, width = 8, units = "in")

pred2 %>% 
  filter(detection == "ACD") %>%
ggplot() +
  geom_sf(data = boundary, fill = NA) +
  geom_point(aes(x = x, y = y, col = (pred - delay)*100/delay), alpha = 0.5) +
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9, 0.85)) +
  scale_color_viridis_c(option = "viridis", direction = -1, end = 0.9, trans = "identity") +
  labs(title = "Predicted change in delay with no ACD", 
       x = "", y = "", col = "% Change") 
ggsave(here::here(figdir, "predvsobs_noacd_map.png"), height = 6, width = 8, units = "in")

# ---------------------------------------------------------------------------- #

## Identify indices which correspond to validation points
index.t <- inla.stack.index(stack = stk, tag = "train")$data
index.p2 <- inla.stack.index(stack = stk, tag = "pred_fullACD")$data

# Extract summary stats of fitted values at these indices
m <- fit.pred$summary.fitted.values[index.p2, "mean"]
ll <- fit.pred$summary.fitted.values[index.p2, "0.025quant"]
ul <- fit.pred$summary.fitted.values[index.p2, "0.975quant"]

pred3 <- data.frame(delay = dat$delay, #
                    detection = factor(dat$poss_acd, levels = c(FALSE,TRUE), labels = c("PCD","ACD")),
                    x = coo[,1], y = coo[,2], 
                    pred = m, lower = ll, upper = ul, iqr = ul - ll,
                    fitted = fit.pred$summary.fitted.values[index.t, "mean"]) %>%
  mutate(change = (pred - fitted)*100/fitted)

pred3 %>%
  ggplot() +
  geom_point(aes(delay, pred, col = detection)) +
  geom_abline(lty = "dashed") +
  scale_fill_viridis_d(option = "plasma") +
  labs(x = "Observed delay", y = "Predicted - complete ACD")
ggsave(here::here(figdir, "obsvpred_fullacd.png"), height = 6, width = 8, units = "in")

pred3 %>%
  ggplot() +
  geom_point(aes(delay, pred - delay, col = detection)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_viridis_d(option = "plasma") +
  labs(x = "Observed delay", y = "Predicted increase with complete ACD")
ggsave(here::here(figdir, "predchng_fullacd.png"), height = 6, width = 8, units = "in")


pred3 %>% 
  dplyr::rename(Fitted = fitted, 
                `Predicted - Complete ACD` = pred) %>%
  tidyr::pivot_longer(c("Fitted","Predicted - Complete ACD")) %>% 
  ggplot(aes(x = value, fill = name)) +
  geom_histogram(position = "identity", alpha = 0.3, bins = 50) +
  labs(x = "Delay (days)", y = "Count", fill = "") +
  theme(legend.position = c(0.8,0.4))
ggsave(here::here(figdir, "fit_vs_fullacd_hist.png"), height = 4, width = 5, units = "in", dpi = 350)

summary(pred3$fitted)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.262  11.854  16.796  31.523  43.970 493.321 
summary(pred3$pred)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.125   9.637  16.329  26.916  35.585 493.321

pred3 %>% filter(detection == "PCD") %>% pull(fitted) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.708  12.290  17.305  35.697  45.648 383.712
pred3 %>% filter(detection == "PCD") %>% pull(pred) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.125   9.641  13.576  28.005  35.812 301.038 

summary(pred3$pred - pred3$fitted)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# -82.67483  -5.66713  -1.76659  -4.60658   0.00000   0.00004 

pred3 %>%
  pivot_longer(c("fitted","pred")) %>% 
  mutate(gt30 = (value > 30)) %>% 
  group_by(name) %>% 
  summarise(p_gt30 = mean(gt30))
# name   p_gt30
# 1 fitted  0.338
# 2 pred    0.293

pred3 %>%
  filter(detection == "PCD") %>%
  ggplot() +
  geom_histogram(aes(x = (pred - delay)*100/delay)) + 
  labs(x = "% change from observed")
ggsave(here::here(figdir, "predvsobs_fullacd_hist.png"), height = 6, width = 8, units = "in")

pred3 %>% 
  filter(detection == "PCD") %>%
  ggplot() +
  geom_sf(data = boundary, fill = NA) +
  geom_point(aes(x = x, y = y, col = (pred - fitted)*100/fitted), alpha = 0.5) +
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9, 0.85)) +
  scale_color_viridis_c(option = "viridis", direction = -1, end = 0.9, trans = "identity") +
  labs(title = "Predicted change in delay with complete ACD", 
       x = "", y = "", col = "% Change") 
ggsave(here::here(figdir, "predchng_fullacd_map.png"), height = 6, width = 8, units = "in")

################################################################################
################################################################################
