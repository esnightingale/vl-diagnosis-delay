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
dat <- read_data()  
coop <- readRDS(here::here("data/analysis","coop.rds"))

covs <- c("age_s","comorb", 
          "poss_acd",
          "block_endm_2017", "inc_2017_gt0",
          "traveltime_t_s")
X <- model.matrix(as.formula(paste("~ ",paste(covs, collapse = " + "))), 
                  data = dat)[,-1] 

boundary <- readRDS(here::here("data","geography","boundary.rds")) %>%
  st_sf() 
boundary.spdf <- as_Spatial(boundary)

blockmap <- readRDS(here::here("data","geography","blockmap.rds")) %>% 
  sf::st_set_crs(7759)

by_block <- readRDS(here::here("data","geography","by_block.rds"))

fit.final <- readRDS(here::here(outdir, "fits_final.rds"))[["All"]]
# fit.final <- fits[["All"]]

#------------------------------------------------------------------------------#
# Final model summary

summary(fit.final$fit)

# Fixed effects:
#   mean    sd 0.025quant 0.5quant 0.975quant   mode kld
# Intercept            3.271 0.075      3.138    3.266      3.436  3.255   0
# age_s                0.118 0.015      0.088    0.118      0.148  0.118   0
# comorb1              0.244 0.080      0.087    0.244      0.402  0.244   0
# poss_acdTRUE        -0.250 0.033     -0.314   -0.250     -0.186 -0.250   0
# block_endm_2017TRUE -0.138 0.054     -0.245   -0.138     -0.032 -0.138   0
# inc_2017_gt0TRUE    -0.078 0.033     -0.142   -0.078     -0.014 -0.078   0
# traveltime_t_s       0.023 0.019     -0.015    0.023      0.061  0.023   0
# 
# Random effects:
#   Name	  Model
# id IID model
# s SPDE2 model
# 
# Model hyperparameters:
#                    mean     sd 0.025quant 0.5quant 0.975quant   mode
# Precision for id  1.165  0.031       1.10    1.166      1.224  1.171
# Range for s      47.259 14.854      26.13   44.520     83.577 39.560
# Stdev for s       0.324  0.048       0.23    0.325      0.416  0.331
# 
# Watanabe-Akaike information criterion (WAIC) ...: 28572.58
# Effective number of parameters .................: 2585.26
# 
# Marginal log-Likelihood:  -18653.69 

#------------------------------------------------------------------------------#
# Posterior predictive distribution

indexs <- inla.spde.make.index("s", spde$n.spde)

coo <- st_coordinates(dat)
A <- inla.spde.make.A(mesh = mesh, loc = coo)

stk.train <- inla.stack(
  tag = "train",
  data = list(y = dat$delay),
  A = list(A, 1, 1),
  effects = list(s = indexs,
                 id = dat$id,
                 data.frame(
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

# Stack for prediction with no ACD
X2 <- X
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
X3 <- X
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

# Stack for prediction with no recent endemicity
# X4 <- X
# X4[,"block_endm_2017"] <- 1
# X4[,"inc_2017_gt0"] <- 1
# stk.pred_nonend <- inla.stack(
#   tag = "pred_nonend",
#   data = list(y = NA),
#   A = list(A, 1, 1),
#   effects = list(s = indexs,
#                  id = dat$id,
#                  data.frame(
#                    Intercept = 1,
#                    X4)
#   )
# )

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

samples <- inla.posterior.sample(1e4, fit.pred)

saveRDS(samples, here::here("output","fit_pred_samples.rds"))

#------------------------------------------------------------------------------#
# Plot smooth predictions

# Identify indices which correspond to validation points
index.p <- inla.stack.index(stack = stk, tag = "pred")$data

# Extract summary stats of fitted values at these indices
m <- fit.pred$summary.fitted.values[index.p, "mean"]
ll <- fit.pred$summary.fitted.values[index.p, "0.025quant"]
ul <- fit.pred$summary.fitted.values[index.p, "0.975quant"]

pred <- data.frame(x = coop[,1], y = coop[,2], Mean = m, Lower = ll, Upper = ul, iqr = ul - ll) 

# dat_in_poly <- sf::st_join(dat, blockmap, join = st_within)
# incl <- which(blockmap$OBJECTID %in% dat_in_poly$OBJECTID)
# blkfill <- rep(NA, nrow(blockmap))
# blkfill[-incl] <- "white"

by_block %>% 
  dplyr::group_by(district) %>% 
  dplyr::summarise(N = sum(N, na.rm = T),
            pop = sum(pop)) %>% 
  ungroup() %>% 
  mutate(inc = N*1e4/pop) -> by_dist

no_inc <- by_dist$district[which(by_dist$inc == 0)]
fillblank <- rep(NA, nrow(blockmap))
fillblank[which(blockmap$district %in% no_inc)] <- "white"

ggplot() +
  geom_tile(data = pred, aes(x = x, y = y, fill = Mean)) + #, alpha = Mean/iqr
  geom_sf(data = blockmap, fill = fillblank) +
  scale_fill_viridis_c(option = "viridis", direction = -1, end = 0.9) +
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9, 0.85)) +
  scale_alpha_continuous(trans = "log2") +
  guides(alpha = "none") +
  # geom_sf(data = dat, pch = "+", col = "white", alpha = 0.5) +
  labs(#title = "Predicted mean of days fever prior to diagnosis", 
       x = "", y = "", fill = "Delay (days)") -> pred_final

pred_final
ggsave(here::here(figdir,"pred_final_dist.png"), pred_final, height = 6, width = 8, units = "in")

ggplot() +
  geom_tile(data = pred, aes(x = x, y = y, fill = iqr/Mean)) +
  geom_sf(data = blockmap, fill = fillblank) +
  scale_fill_viridis_c(option = "viridis", direction = -1, end = 0.9) +
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9, 0.85)) +
  # geom_sf(data = dat, pch = "+", col = "white", alpha = 0.5) +
  labs(#title = "Predicted mean of days fever prior to diagnosis", 
    x = "", y = "", fill = "(Q97.5 - Q2.5)/Mean") 

ggsave(here::here(figdir,"pred_QR_mean.png"), height = 6, width = 8, units = "in")

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
  geom_sf(data = blockmap, fill = fillblank) +
  scale_fill_viridis_c(direction = -1, option = "cividis") +
  labs(#title = "Predicted probability of delay exceeding 30 days",
       x = "", y = "", fill = "P(delay > 30)") +
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom" #c(0.9, 0.85) #
        ) -> map_exc
map_exc

ggsave(here::here(figdir,"excprob30.png"), map_exc, height = 6, width = 8, units = "in")

ggplot() +
  geom_tile(data = pred, aes(x = x, y = y, fill = excprob_hi, alpha = excprob_strength)) +
  geom_sf(data = blockmap, fill = fillblank) +
  scale_fill_viridis_d(option = "cividis") +
  # scale_fill_viridis_d(option = "viridis", begin = 0, end = 0.9) + 
  labs(#title = "Predicted probability of delay exceeding 30 days",
       x = "", y = "", fill = "") + #"Probability of\ndelay > 30 days"
  guides(alpha = "none") +
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9, 0.85) #"bottom"
        ) -> map_exc2
map_exc2

ggsave(here::here(figdir,"excprob30_final.png"), map_exc2, height = 6, width = 8, units = "in")

# map_exc + map_exc2
# ggsave(here::here(figdir,"excprob30_final_combined.png"), height = 7, width = 14, units = "in")

excursions30 <- excursions.inla(
  fit.pred,
  stack = stk,
  tag = "pred",
  method = "EB",
  u = 30,
  u.link = TRUE,
  type = ">",
  n.iter = 10000,
  F.limit = 0,
  verbose = 1,
  seed = 1234
)
saveRDS(excursions30, here::here("output","excursions.rds"))

summary(excursions30)

#Define a fine mesh over region of interest
coo.bnd <- st_coordinates(boundary)[,1:2]
nxy <- c(100, 100)
projgrid <- inla.mesh.projector(mesh, xlim = range(coo.bnd[, 1]), ylim = range(coo.bnd[,2]), dims = nxy)
xy.in <- inout(projgrid$lattice$loc, cbind(coo.bnd[, 1], coo.bnd[, 2]))
submesh <- submesh.grid(matrix(xy.in, nxy[1], nxy[2]), list(loc = projgrid$lattice$loc,
                                                            dims = nxy))

#Interpolate excursion set to the mesh
sets <- continuous(excursions30, submesh, alpha = 0.1)

#Plot excursion set
plot(sets$M["1"], col = "red")

#plot excursion function
proj <- inla.mesh.projector(sets$F.geometry, dims = c(300,200))
image(proj$x, proj$y, inla.mesh.project(proj, field = sets$F))

# exc.set <- data.frame(x = coop[,1], y = coop[,2], exc30$rho

plot(exc30$rho,type="l",
     main="marginal probabilities (black) and excursion function (red)")
lines(exc30$F,col=2)

#------------------------------------------------------------------------------#
# Predictions for full/no ACD coverage

## Identify indices which correspond to validation points
index.t <- inla.stack.index(stack = stk, tag = "train")$data
index.p1 <- inla.stack.index(stack = stk, tag = "pred_fullACD")$data
index.p2 <- inla.stack.index(stack = stk, tag = "pred_noACD")$data

pred_fullACD <- summarise_pred(index.p1)
pred_noACD <- summarise_pred(index.p2)

write.csv(pred_fullACD$tab, here::here("output","tables","tab_fullACD.csv"))
write.csv(pred_noACD$tab, here::here("output","tables","tab_noACD.csv"))

saveRDS(list(fullACD = pred_fullACD, noACD = pred_noACD), 
        here::here(outdir, "summary_fullACD_noACD.rds"))

pred_fullACD$endm %>%
  separate(chg.ci, ", ", into = c("low","high"), convert = TRUE) %>%
  ggplot(aes(x = block_endm_2017, y = chg, ymin = low, ymax = high)) +
    geom_errorbar(width = 0.2) +
    geom_point() +
    labs(x = "", y = "Change in delay (days), total") -> change_tot

pred_fullACD$endm %>%
  separate(chg.pc.ci, ", ", into = c("low","high"), convert = TRUE) %>%
  ggplot(aes(x = block_endm_2017, y = chg.pc, ymin = low, ymax = high)) +  
  geom_errorbar(width = 0.2) +
  geom_point() +
  ylim(c(0,12))+
  labs(x = "", y = "Change in delay (days), per case") -> change_pc

pred_fullACD$endm %>%
  separate(chg.pc.PCD.ci, ", ", into = c("low","high"), convert = TRUE) %>% View()
  ggplot(aes(x = block_endm_2017, y = chg.pc.PCD, ymin = low, ymax = high)) +
  geom_errorbar(width = 0.2) +
  geom_point()+
  ylim(c(0,12)) +
  labs(x = "", y = "Change in delay (days), per PCD case") -> change_pPCD

ggsave(here::here(figdir, "pred_days_change_fullACD.png"), grid.arrange(change_tot, change_pc, change_pPCD, nrow = 1), height = 4, width = 12)

pred_fullACD$endm %>%
  separate(chg.pc.ci, ", ", into = c("low","high"), convert = TRUE) %>%
  mutate(scenario = "100% ACD coverage",
         estimate = chg.pc,
         measure = "per case") %>% 
  dplyr::select(scenario, block_endm_2017, measure, estimate, low, high) -> plotdata1
pred_fullACD$endm %>%
  separate(chg.pc.PCD.ci, ", ", into = c("low","high"), convert = TRUE) %>%
  mutate(scenario = "100% ACD coverage",
         estimate = chg.pc.PCD,
         measure = "per reassigned case") %>% 
  bind_rows(plotdata1)  %>% 
  dplyr::select(scenario, block_endm_2017, measure, estimate, low, high) -> plotdata1

pred_noACD$endm %>%
  separate(chg.pc.ci, ", ", into = c("low","high"), convert = TRUE) %>%
  mutate(scenario = "0% ACD coverage",
         estimate = chg.pc,
         measure = "per case")  %>% 
  dplyr::select(scenario, block_endm_2017, measure, estimate, low, high) -> plotdata2
pred_noACD$endm %>%
  separate(chg.pc.ACD.ci, ", ", into = c("low","high"), convert = TRUE) %>%
  mutate(scenario = "0% ACD coverage",
         estimate = chg.pc.ACD,
         measure = "per reassigned case") %>% 
  dplyr::select(scenario, block_endm_2017, measure, estimate, low, high) %>% 
  bind_rows(plotdata2) %>% 
  bind_rows(plotdata1) -> plotdata

dodge <- position_dodge(width=0.5)  
plotdata %>% 
  ggplot(aes(x = block_endm_2017, y = estimate, ymin = low, ymax = high, lty = measure, colour = scenario)) +  
  geom_errorbar(width = 0.2, position = dodge) +
  geom_point(position = dodge) +
  geom_hline(yintercept = 0, col = "grey") +
  # theme(legend.position = c(0.8,0.95)) +
  labs(x = "", y = "Average change in expected delay (days)", colour = "", lty = "") -> days_change
days_change

ggsave(here::here(figdir, "days_change_ACDcovg.png"), days_change, height = 5, width = 7, dpi = 300)

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

var.pooled <- function(df1,df2,SD1,SD2){
  (df1*SD1^2 + df2*SD2^2)/(df1+df2)
}

SE.diff <- function(var.pool, n1,n2){
  sqrt(var.pool*(1/n1 + 1/n2))
}

dat %>%
  st_drop_geometry() %>%
  mutate(delay_acd = na_if(delay*(poss_acd == TRUE), 0),
         delay_pcd = na_if(delay*(poss_acd != TRUE), 0)) %>%
  group_by(district, block) %>%
  dplyr::summarise(N = n(),
                   med_delay = median(delay),
                   pACD = mean(poss_acd*100),
                   # med_delay_acd = median(delay_acd, na.rm = T),
                   # med_delay_pcd = median(delay_pcd, na.rm = T),
                   mean_delay_acd = mean(delay_acd, na.rm = T),
                   sd_delay_acd = sd(delay_acd, na.rm = T),
                   n_acd = sum(poss_acd == TRUE),
                   mean_delay_pcd = mean(delay_pcd, na.rm = T),
                   sd_delay_pcd = sd(delay_pcd, na.rm = T),
                   n_pcd = sum(poss_acd == FALSE)) %>%
  ungroup() %>% 
  mutate(var.pool = var.pooled(n_acd-1,
                               n_pcd-1,
                               sd_delay_acd,
                               sd_delay_pcd),
         diff_delay_acd = (mean_delay_pcd - mean_delay_acd),
         # diff_med_delay_acd = (med_delay_pcd - med_delay_acd),
         se.diff = SE.diff(var.pool,
                           n1 = n_acd,
                           n2 = n_pcd),
         diff_delay_acd_std = diff_delay_acd/se.diff) -> by_block

by_block <- blockmap %>%
  left_join(by_block, by = c("kamis_master_dist" = "district", "kamis_master_block" = "block")) %>% 
  rowwise() %>%
  dplyr::mutate(pop = median(c_across(`2018`:`2019`)),
                inc = N*1e4/pop) %>%
  ungroup() 

ggplot() +
  geom_sf(data = by_block, aes(geometry = geometry, fill = diff_delay_acd_std)) +
  geom_sf(data = blockmap, fill = NA) +
  scale_fill_viridis_c(option = "plasma", na.value = "white", direction = 1) +
  labs(fill = "ACD impact", caption = "Difference in mean delay between ACD/PCD, divided by Std Err.") + 
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.9,0.85)) -> blk_diff_acd
blk_diff_acd

dat %>%
  st_drop_geometry() %>%
  mutate(delay_acd = na_if(delay*(poss_acd == TRUE), 0),
         delay_pcd = na_if(delay*(poss_acd != TRUE), 0)) %>%
  group_by(block_endm_2017) %>%
  dplyr::summarise(N = n(),
                   mean_delay = mean(delay),
                   pACD = mean(poss_acd*100),
                   mean_delay_acd = mean(delay_acd, na.rm = T),
                   sd_delay_acd = sd(delay_acd, na.rm = T),
                   n_acd = sum(poss_acd == TRUE),
                   mean_delay_pcd = mean(delay_pcd, na.rm = T),
                   sd_delay_pcd = sd(delay_pcd, na.rm = T),
                   n_pcd = sum(poss_acd == FALSE)) %>%
  ungroup() %>% 
  mutate(var.pool = var.pooled(n_acd-1,
                               n_pcd-1,
                               sd_delay_acd,
                               sd_delay_pcd),
         diff_delay_acd = (mean_delay_pcd - mean_delay_acd),
         # diff_med_delay_acd = (med_delay_pcd - med_delay_acd),
         se.diff = SE.diff(var.pool,
                           n1 = n_acd,
                           n2 = n_pcd),
         diff_delay_acd_std = diff_delay_acd/se.diff) -> by_endm

write.csv(t(by_endm), here::here("output","acd_impact_endm.csv"))

ggplot(by_endm, aes(block_endm_2017, diff_delay_acd, 
                    ymin = diff_delay_acd - 1.96*se.diff,
                    ymax = diff_delay_acd + 1.96*se.diff)) +
  geom_errorbar(width = 0.4) +
  geom_point() +
  labs(x = "Block endemic (2017)", y = "Difference in average delay: PCD - ACD")

################################################################################
################################################################################


# pred2 %>%
#   ggplot() +
#   geom_point(aes(fitted, pred, col = detection)) +
#   geom_abline(lty = "dashed") +
#   scale_fill_viridis_d(option = "plasma") +
#   labs(x = "Observed delay", y = "Predicted - no ACD")
# ggsave(here::here(figdir, "obsvpred_noacd.png"), height = 6, width = 8, units = "in")
# 
# pred2 %>%
#   ggplot() +
#   geom_point(aes(fitted, pred - delay, col = detection)) +
#   geom_hline(yintercept = 0, lty = "dashed") +
#   scale_fill_viridis_d(option = "plasma") +
#   labs(x = "Observed delay", y = "Predicted increase with no ACD")
# ggsave(here::here(figdir, "predchng_noacd.png"), height = 6, width = 8, units = "in")
# 
# pred2 %>% 
#   dplyr::rename(Fitted = fitted, 
#                 `Predicted - No ACD` = pred) %>%
#   tidyr::pivot_longer(c("Fitted","Predicted - No ACD")) %>% 
#   ggplot(aes(x = value, fill = name)) +
#   geom_histogram(position = "identity", alpha = 0.3, bins = 50) +
#   labs(x = "Delay (days)", y = "Count", fill = "") +
#   theme(legend.position = c(0.8,0.4))
# ggsave(here::here(figdir, "fit_vs_noacd_hist.png"), height = 4, width = 5, units = "in", dpi = 350)
# 
# summary(pred2$fitted)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 2.309  11.845  16.798  31.523  43.970 493.345 
# summary(pred2$pred)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 2.654  12.276  20.863  34.343  45.384 630.491 
# 
# pred2 %>% filter(detection == "ACD") %>% pull(fitted) %>% summary()
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 2.309   9.635  16.425  25.291  30.540 493.345
# pred2 %>% filter(detection == "ACD") %>% pull(pred) %>% summary()
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 2.95   12.31   20.99   32.32   39.03  630.49 
# 
# summary(pred2$pred - pred2$fitted)
# # Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# # -0.00025   0.00000   0.00001   2.81962   4.23487 137.14586 
# 
# pred2 %>%
#   pivot_longer(c("fitted","pred")) %>% 
#   mutate(gt30 = (value > 30)) %>% 
#   group_by(name) %>% 
#   summarise(p_gt30 = mean(gt30))
# # name   p_gt30
# # 1 fitted  0.338
# # 2 pred    0.351
# 
# pred2 %>%
#   filter(detection == "ACD") %>%
# ggplot() +
#   geom_histogram(aes(x = (pred - delay)*100/delay)) + 
#   labs(x = "% change from observed")
# ggsave(here::here(figdir, "predvsobs_noacd_hist.png"), height = 6, width = 8, units = "in")
# 
# pred2 %>% 
#   filter(detection == "ACD") %>%
# ggplot() +
#   geom_sf(data = boundary, fill = NA) +
#   geom_point(aes(x = x, y = y, col = (pred - delay)*100/delay), alpha = 0.5) +
#   theme(axis.text = element_blank(), panel.grid = element_blank(),
#         legend.position = c(0.9, 0.85)) +
#   scale_color_viridis_c(option = "viridis", direction = -1, end = 0.9, trans = "identity") +
#   labs(title = "Predicted change in delay with no ACD", 
#        x = "", y = "", col = "% Change") 
# ggsave(here::here(figdir, "predvsobs_noacd_map.png"), height = 6, width = 8, units = "in")
# pred3 %>%
#   ggplot() +
#   geom_point(aes(delay, pred, col = detection)) +
#   geom_abline(lty = "dashed") +
#   scale_fill_viridis_d(option = "plasma") +
#   labs(x = "Observed delay", y = "Predicted - complete ACD")
# ggsave(here::here(figdir, "obsvpred_fullacd.png"), height = 6, width = 8, units = "in")
# 
# pred3 %>%
#   ggplot() +
#   geom_point(aes(delay, pred - delay, col = detection)) +
#   geom_hline(yintercept = 0, lty = "dashed") +
#   scale_fill_viridis_d(option = "plasma") +
#   labs(x = "Observed delay", y = "Predicted increase with complete ACD")
# ggsave(here::here(figdir, "predchng_fullacd.png"), height = 6, width = 8, units = "in")
# 
# pred3 %>% 
#   dplyr::rename(Fitted = fitted, 
#                 `Predicted - Complete ACD` = pred) %>%
#   tidyr::pivot_longer(c("Fitted","Predicted - Complete ACD")) %>% 
#   ggplot(aes(x = value, fill = name)) +
#   geom_histogram(position = "identity", alpha = 0.3, bins = 50) +
#   labs(x = "Delay (days)", y = "Count", fill = "") +
#   theme(legend.position = c(0.8,0.4))
# ggsave(here::here(figdir, "fit_vs_fullacd_hist.png"), height = 4, width = 5, units = "in", dpi = 350)
# 
# summary(pred3$fitted)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 2.309  11.845  16.798  31.523  43.970 493.345 
# summary(pred3$pred)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 2.079   9.613  16.325  26.890  35.545 493.345 
# 
# pred3 %>% filter(detection == "PCD") %>% pull(fitted) %>% summary()
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 2.654  12.271  17.293  35.697  45.648 383.763 
# pred3 %>% filter(detection == "PCD") %>% pull(pred) %>% summary()
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 2.078   9.612  13.546  27.961  35.756 300.604
# 
# summary(pred3$pred - pred3$fitted)
# # Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# # -83.15954  -5.69832  -1.76723  -4.63302   0.00000   0.00003
# 
# pred3 %>%
#   pivot_longer(c("fitted","pred")) %>% 
#   mutate(gt30 = (value > 30)) %>% 
#   group_by(name) %>% 
#   summarise(p_gt30 = mean(gt30))
# # name   p_gt30
# # 1 fitted  0.338
# # 2 pred    0.293
# 
# pred3 %>%
#   filter(detection == "PCD") %>%
#   ggplot() +
#   geom_histogram(aes(x = (pred - delay)*100/delay)) + 
#   labs(x = "% change from observed")
# ggsave(here::here(figdir, "predvsobs_fullacd_hist.png"), height = 6, width = 8, units = "in")
# 
# pred3 %>% 
#   filter(detection == "PCD") %>%
#   ggplot() +
#   geom_sf(data = boundary, fill = NA) +
#   geom_point(aes(x = x, y = y, col = (pred - fitted)*100/fitted), alpha = 0.5) +
#   theme(axis.text = element_blank(), panel.grid = element_blank(),
#         legend.position = c(0.9, 0.85)) +
#   scale_color_viridis_c(option = "viridis", direction = -1, end = 0.9, trans = "identity") +
#   labs(title = "Predicted change in delay with complete ACD", 
#        x = "", y = "", col = "% Change") 
# ggsave(here::here(figdir, "predchng_fullacd_map.png"), height = 6, width = 8, units = "in")
