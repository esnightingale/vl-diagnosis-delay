################################################################################
# Description: Univariate regression on covariates and assessment of empirical 
# variogram for spatial dependence
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/descriptive"

dat <- readRDS(here::here("data","analysisdata_individual.rds")) %>%
  dplyr::mutate(vill_inc_2017_gt0 = factor((replace_na(vill_inc_2017,0) > 0), labels = c("No","Yes")),
                vill_inc_2017_t = vill_inc_2017*1e3 + 1e-4,
                IRS_2017_1 = factor((IRS_2017 != 0), labels = c("No","Yes")),
                num_conslt_cat = as.factor(cut(num_conslt, breaks = c(0,1,3,5,8), include.lowest = TRUE)))
# dat.df <- st_drop_geometry(dat) 

dat.fit <- dat %>%
  st_drop_geometry() %>%
  dplyr::select(days_fever, age, sex, hiv, marg_caste, detection, prv_tx_ka, 
                num_conslt, num_conslt_cat, latitude, longitude, traveltime, vill_inc_2017_gt0, vill_inc_2017_t, 
                IRS_2017_1, block_endm_2017, vil_code, patient_id) %>%
  drop_na() 

village <- readRDS(here::here("data","analysisdata_village.rds")) 

# Setup map context
blockmap <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
  st_transform(4326)
boundary <- sf::st_union(blockmap)  
boundary.spdf <- as_Spatial(boundary)

extent <- setNames(st_bbox(blockmap), c("left","bottom","right","top"))
bh_lines <- get_stamenmap(bbox = extent, maptype = "terrain-lines", zoom = 8)

# ---------------------------------------------------------------------------- #
# Visualise the spatial range of the data and pattern of delays

ggplot() +
  geom_sf(data = blockmap, fill = NA) +
  geom_sf(data = st_jitter(dat), aes(col = days_fever), alpha = 0.8, cex = 0.5) +
  scale_colour_viridis_c(trans = "log2") +
  labs(col = "Delay (days)") + 
  theme(axis.text = element_blank(), panel.grid = element_blank()) +
  annotation_scale(location = "br") +
  annotation_north_arrow(location = "tr", pad_x = unit(2, "cm"), pad_y = unit(1, "cm"))

ggsave(here::here(figdir, "fig1.png"), height = 7, width = 9, units = "in")

# ---------------------------------------------------------------------------- #
# Correlation between all covariates

vars <- dat.fit %>%
  dplyr::select(-longitude, -latitude, -vil_code, -patient_id) %>%
  dplyr::mutate(across(everything(), as.numeric))

M <- cor(vars)
corrplot::corrplot(M, method = "color", type = "lower", order = "hclust")

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
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

png(here::here(figdir,"covariate effects", "corrplot.png"), 
    height = 8, width = 8, units = "in", res = 300)
corrplot::corrplot(M, method = "color", type = "lower", order = "hclust", 
         p.mat = p.mat, sig.level = 0.01)
dev.off()

# ---------------------------------------------------------------------------- #
# Fit regression with individual level covariates and investigate residuals

glm.fit.pat <- glm(days_fever ~ age + sex + hiv + marg_caste + detection + 
                     prv_tx_ka + num_conslt,
               family = "poisson",
               data = dat.fit)

summary(glm.fit.pat)
# Note: no evidence of interaction between age and poss ACD

# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -10.784   -3.056   -1.334    0.950   36.822  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)    3.2264918  0.0074025 435.865  < 2e-16 ***
#   age            0.0038424  0.0001226  31.348  < 2e-16 ***
#   sexFemale      0.0164782  0.0046661   3.531 0.000413 ***
#   hivYes         0.3049803  0.0101991  29.903  < 2e-16 ***
#   marg_casteYes  0.0613480  0.0049299  12.444  < 2e-16 ***
#   detectionACD  -0.1973874  0.0048309 -40.859  < 2e-16 ***
#   prv_tx_kaYes   0.0184716  0.0079900   2.312 0.020787 *  
#   num_conslt     0.1584820  0.0014690 107.883  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 97164  on 4321  degrees of freedom
# Residual deviance: 80679  on 4314  degrees of freedom
# AIC: 104190
# 
# 
# With categorical num_consult:
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -11.562   -3.098   -1.462    0.922   38.138  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)          3.4788859  0.0088799 391.769  < 2e-16 ***
#   age                  0.0040203  0.0001226  32.791  < 2e-16 ***
#   sexFemale            0.0159066  0.0046685   3.407 0.000656 ***
#   hivYes               0.2958438  0.0102310  28.916  < 2e-16 ***
#   marg_casteYes        0.0478543  0.0049551   9.658  < 2e-16 ***
#   detectionACD        -0.1989794  0.0048358 -41.147  < 2e-16 ***
#   prv_tx_kaYes        -0.0018002  0.0080283  -0.224 0.822572    
#   num_conslt_cat(1,3]  0.1300913  0.0077991  16.680  < 2e-16 ***
#   num_conslt_cat(3,5]  0.4249002  0.0083179  51.083  < 2e-16 ***
#   num_conslt_cat(5,8]  0.8168850  0.0100145  81.570  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 97164  on 4321  degrees of freedom
# Residual deviance: 81970  on 4312  degrees of freedom
# AIC: 105485

glm.fit.vil <- glm(days_fever ~ log(vill_inc_2017_t,10) + IRS_2017_1 + block_endm_2017 + traveltime,
                   family = "poisson",
                   data = dat.fit)

summary(glm.fit.vil)

# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -8.561  -3.348  -1.995   1.289  38.069  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)               3.8195461  0.0072561 526.392  < 2e-16 ***
#   log(vill_inc_2017_t, 10) -0.0283806  0.0013120 -21.632  < 2e-16 ***
#   IRS_2017_1Yes             0.0222660  0.0058550   3.803 0.000143 ***
#   block_endm_2017Endemic   -0.1501748  0.0048335 -31.070  < 2e-16 ***
#   traveltime               -0.0013195  0.0001883  -7.009 2.41e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 97164  on 4321  degrees of freedom
# Residual deviance: 95146  on 4317  degrees of freedom
# AIC: 118651

glm.fit.all <- glm(days_fever ~ age + sex + hiv + marg_caste + detection + 
                     prv_tx_ka + num_conslt_cat + log(vill_inc_2017_t,10) + 
                     IRS_2017_1 + block_endm_2017 + traveltime,
                   family = "poisson",
                   data = dat.fit)

summary(glm.fit.all)

# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -11.821   -3.039   -1.407    1.021   37.271  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)               3.499e+00  1.128e-02 310.253  < 2e-16 ***
#   age                       4.115e-03  1.231e-04  33.439  < 2e-16 ***
#   sexFemale                 1.620e-02  4.672e-03   3.468 0.000524 ***
#   hivYes                    2.558e-01  1.030e-02  24.836  < 2e-16 ***
#   marg_casteYes             2.762e-02  5.001e-03   5.523 3.34e-08 ***
#   detectionACD             -2.041e-01  4.867e-03 -41.939  < 2e-16 ***
#   prv_tx_kaYes              4.171e-03  8.048e-03   0.518 0.604297    
#   num_conslt_cat(1,3]       1.442e-01  7.840e-03  18.392  < 2e-16 ***
#   num_conslt_cat(3,5]       4.273e-01  8.351e-03  51.172  < 2e-16 ***
#   num_conslt_cat(5,8]       8.090e-01  1.004e-02  80.598  < 2e-16 ***
#   log(vill_inc_2017_t, 10) -1.231e-02  1.328e-03  -9.268  < 2e-16 ***
#   IRS_2017_1Yes             3.145e-02  5.861e-03   5.366 8.07e-08 ***
#   block_endm_2017Endemic   -1.571e-01  4.934e-03 -31.847  < 2e-16 ***
#   traveltime               -2.431e-05  1.905e-04  -0.128 0.898435    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 97164  on 4321  degrees of freedom
# Residual deviance: 80636  on 4308  degrees of freedom
# AIC: 104158

# ---------------------------------------------------------------------------- #
# Explore the relationship between residuals and village covariates

dat.fit$res.pat <- glm.fit.pat$residuals

plot_access <- ggplot(dat.fit, aes(x = traveltime, 
                                   y = res)) + 
  geom_point() +
  labs(x = "Travel time", y = "") +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE) +
  stat_smooth(method = "lm", formula = y ~ x,
              col = "green", lty = "dashed", se = FALSE)

lm.fit.res <- lm(res.pat ~ traveltime, data = dat.fit)
summary(lm.fit.res)

plot_inc <- ggplot(dat.fit, aes(x = vill_inc_2017, 
                                   y = res.pat)) + 
  geom_point() +
  labs(x = "Village incidence, 2017", y = "") +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE) +
  stat_smooth(method = "lm", formula = y ~ x,
              col = "green", lty = "dashed", se = FALSE) 


create_raincloud(dat.fit, "IRS_2017_1", "res.pat", 
                 xlab = "IRS_2017_1", ylab = "Residual", 
                 y_trans = "identity", col_by = "IRS_2017_1") -> plot_irs


create_raincloud(dat.fit, "block_endm_2017", "res.pat", 
                 xlab = "block_endm_2017", ylab = "Residual", 
                 y_trans = "identity", col_by = "block_endm_2017") -> plot_endm

gridExtra::grid.arrange(plot_access,
                        plot_inc,
                        plot_irs,
                        plot_endm,
                        nrow = 2)

# ---------------------------------------------------------------------------- #
# Testing for residual spatial correlation after accounting for patient and 
# village covariates

dat.fit$res.all <- glm.fit.all$residuals

vg <- variogram(res.all~1, data = dat.fit)
plot(vg)

vgmod <- vgm(psill =  0.4, model = "Mat", nugget = 0.3, range = 100)
plot(vg, model = vgmod)

vgfit <- fit.variogram(vg, model = vgmod)    
plot(vg, model = vgfit)

vgfit
# model     psill    range kappa
# 1   Nug 0.3532252   0.0000   0.0
# 2   Mat 0.4729930 209.5673   0.5

png(here::here(figdir, "glmresid_semivariogram.png"), height = 500, width = 600)
plot(vg, model = vgfit)
dev.off()

# ---------------------------------------------------------------------------- #
# Mixed model - check spatial dependence of fitted random effects

dat.std <- standardize::standardize(days_fever ~ age + sex + hiv + marg_caste + detection + 
                                      prv_tx_ka + num_conslt + log(vill_inc_2017_t,10) + 
                                      IRS_2017_1 + block_endm_2017 + traveltime + (1|vil_code),
                                    family = poisson,
                                    data = dat.fit)

dat.std

glm.std <- lme4::lmer(dat.std$formula, dat.std$data)

summary(glm.std)
# Linear mixed model fit by REML ['lmerMod']
# Formula: days_fever ~ age + sex + hiv + marg_caste + detection + prv_tx_ka +      num_conslt + vill_inc_2017 + IRS_2017_1 + block_endm_2017 +  
#   traveltime + (1 | vil_code)
# Data: dat.std$data
# 
# REML criterion at convergence: 43351.4
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.4740 -0.4790 -0.1850  0.2154 10.5286 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# vil_code (Intercept)  334.7   18.29   
# Residual             1030.4   32.10   
# Number of obs: 4336, groups:  vil_code, 2337
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)             52.0900     1.7359  30.008
# age                      3.2314     0.5583   5.788
# sexFemale                0.1319     0.5513   0.239
# hivYes                   8.0337     1.5141   5.306
# marg_casteYes            0.6262     0.6276   0.998
# detectionACD            -4.1498     0.5830  -7.119
# prv_tx_kaYes             0.4744     0.9744   0.487
# num_conslt              11.0837     0.5716  19.390
# vill_inc_2017           -0.5861     0.7246  -0.809
# IRS_2017_1Yes            0.1192     0.7589   0.157
# block_endm_2017Endemic  -3.9773     0.7007  -5.676
# traveltime               0.5087     0.6474   0.786
# 
# Correlation of Fixed Effects:
#   (Intr) age    sexFml hivYes mrg_cY dtcACD prv__Y nm_cns v__201 IRS_20 b__201
# age         -0.107                                                                      
# sexFemale    0.111  0.108                                                               
# hivYes       0.749 -0.094  0.071                                                        
# marg_castYs  0.148  0.092 -0.003  0.042                                                 
# detectinACD  0.149  0.012 -0.020  0.073 -0.075                                          
# prv_tx_kaYs  0.360 -0.067  0.025 -0.138  0.016  0.016                                   
# num_conslt   0.014 -0.062  0.000 -0.014  0.108  0.020  0.051                            
# vll_nc_2017  0.069  0.019 -0.003  0.019 -0.020 -0.018 -0.009  0.027                     
# IRS_2017_1Y -0.202  0.038  0.005  0.023  0.012 -0.026 -0.015  0.017 -0.140              
# blck__2017E  0.193 -0.035  0.003  0.073  0.085  0.046  0.012  0.002 -0.100 -0.170       
# traveltime   0.008  0.013 -0.014  0.009 -0.071 -0.058  0.018  0.043 -0.016 -0.018  0.043



# ---------------------------------------------------------------------------- #
# Empirical variogram and permutation test for spatial dependence

PrevMap::spat.corr.diagnostic(days_fever ~ age + sex + hiv + marg_caste + 
                                detection + prv_tx_ka + num_conslt,
                     coords = ~ longitude + latitude,
                     likelihood = "Poisson",
                     data = dat.std$data)

# With spatial covariates

dat.fit <- dat.fit %>%
  mutate(age_s = scale(age, center = T),
         num_conslt_s = scale(num_conslt, center = T),
         traveltime_s = scale(traveltime, center = T),
         vill_inc_2017_s = scale(vill_inc_2017, center = T))

PrevMap::spat.corr.diagnostic(days_fever ~ age_s + sex + hiv + marg_caste + 
                                detection + prv_tx_ka + num_conslt + 
                                traveltime_s + log(vill_inc_2017_t,10) + IRS_2017_1 + 
                                block_endm_2017,
                              coords = ~ longitude + latitude,
                              likelihood = "Poisson",
                              data = dat.fit)

################################################################################
################################################################################