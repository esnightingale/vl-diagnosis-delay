################################################################################
# Description: Univariate regression on covariates and assessment of empirical 
# variogram for spatial dependence
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/descriptive"

dat <- readRDS(here::here("data","analysisdata_individual.rds")) %>%
  dplyr::mutate(vill_inc_2017_gt0 = factor((replace_na(vill_inc_2017,0) > 1/1000), labels = c("No","Yes")),
                vill_inc_2017_t = vill_inc_2017*1e3 + 1e-4,
                IRS_2017_1 = factor((IRS_2017 != 0), labels = c("No","Yes")),
                num_conslt_cat = as.factor(cut(num_conslt, breaks = c(0,1,3,5,8), include.lowest = TRUE)),
                age_s = as.numeric(scale(age, center = T)),
                traveltime_s = as.numeric(scale(traveltime, center = T)),
                i = row_number())
# dat.df <- st_drop_geometry(dat) 

dat.fit <- dat %>%
  st_drop_geometry() %>%
  dplyr::select(days_fever, age_s, sex, hiv, marg_caste, detection, prv_tx_ka, 
                num_conslt, num_conslt_cat, latitude, longitude, traveltime_s, vill_inc_2017_t, vill_inc_2017_gt0, 
                IRS_2017_1, block_endm_2017, vil_code, i) %>%
  drop_na() 

# ---------------------------------------------------------------------------- #
# Fit regression with individual level covariates and IID random effect on village

glmm.pat <- lme4::glmer(days_fever ~ age_s + sex + hiv + marg_caste + 
                             detection + num_conslt_cat + prv_tx_ka + (1|i),
                        family = "poisson",
                        data = dat.fit)

summary(glmm.pat)

# AIC      BIC   logLik deviance df.resid 
# 57099.2  57156.6 -28540.6  57081.2     4313 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -10.9629  -0.9796  -0.1160   0.3682  28.0756 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# vil_code (Intercept) 0.2789   0.5281  
# Number of obs: 4322, groups:  vil_code, 2329
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    3.246288   0.015060 215.557  < 2e-16 ***
#   age_s          0.058846   0.003459  17.014  < 2e-16 ***
#   sexFemale     -0.014964   0.006964  -2.149   0.0317 *  
#   hivYes         0.144999   0.022321   6.496 8.24e-11 ***
#   marg_casteYes  0.021088   0.008947   2.357   0.0184 *  
#   detectionACD  -0.157695   0.007877 -20.018  < 2e-16 ***
#   num_conslt     0.155835   0.002616  59.561  < 2e-16 ***
#   prv_tx_kaYes  -0.029190   0.012625  -2.312   0.0208 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) age_s  sexFml hivYes mrg_cY dtcACD nm_cns
# age_s       -0.023                                          
# sexFemale   -0.193  0.112                                   
# hivYes      -0.094 -0.065  0.053                            
# marg_castYs -0.212  0.060 -0.026  0.020                     
# detectinACD -0.197  0.011 -0.017  0.073 -0.016              
# num_conslt  -0.552 -0.034 -0.002 -0.008  0.029  0.023       
# prv_tx_kaYs -0.103 -0.069  0.034 -0.084 -0.006  0.048  0.041

# ---------------------------------------------------------------------------- #
# Explore the relationship between residuals and village covariates

dat.fit$res.pat <- glmm.pat$residuals

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
# Testing for residual spatial correlation after accounting for patient covariates

dat.fit$res.pat <- glm.fit.pat$residuals

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
# Empirical variogram and permutation test for spatial dependence

PrevMap::spat.corr.diagnostic(days_fever ~ age + sex + hiv + marg_caste + 
                                detection + prv_tx_ka + num_conslt,
                     coords = ~ longitude + latitude,
                     likelihood = "Poisson",
                     data = dat.std$data)

# ---------------------------------------------------------------------------- #
# Village covariates

glmm.vil <- lme4::glmer(days_fever ~ log(vill_inc_2017_t,10) + IRS_2017_1 + 
                          block_endm_2017 + traveltime_s + (1|patient_id),
                        family = "poisson",
                        data = dat.fit)

summary(glmm.vil)

# AIC      BIC   logLik deviance df.resid 
# 61514.9  61553.1 -30751.4  61502.9     4316 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -10.3479  -1.0916  -0.1650   0.2981  29.0773 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# vil_code (Intercept) 0.3203   0.5659  
# Number of obs: 4322, groups:  vil_code, 2329
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
#  (Intercept)               3.674621   0.033677 109.113  < 2e-16 ***
#   log(vill_inc_2017_t, 10) -0.013767   0.007027  -1.959  0.05010 .  
#   IRS_2017_1Yes            -0.002970   0.028213  -0.105  0.91616    
#   block_endm_2017Endemic   -0.080009   0.026101  -3.065  0.00217 ** 
#   traveltime_s             -0.001783   0.011508  -0.155  0.87688    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) l(__21 IRS_20 b__201
#              l(__2017_,1  0.751                     
#                  IRS_2017_1Y -0.720 -0.352              
#                  blck__2017E -0.283 -0.146 -0.110       
#                  traveltim_s -0.034 -0.045 -0.003  0.043

# With spatial covariates

dat.fit <- dat.fit %>%
  mutate(age_s = scale(age, center = T),
         num_conslt_s = scale(num_conslt, center = T),
         traveltime_s = scale(traveltime, center = T),
         vill_inc_2017_s = scale(vill_inc_2017, center = T))

PrevMap::spat.corr.diagnostic(days_fever ~ age_s + sex + hiv + marg_caste + 
                                detection + prv_tx_ka + num_conslt_s + 
                                traveltime_s + vill_inc_2017_s + IRS_2017_1 + 
                                block_endm_2017,
                              coords = ~ longitude + latitude,
                              likelihood = "Poisson",
                              data = st_drop_geometry(dat.fit))
# ---------------------------------------------------------------------------- #
# With individual and village covariates

glm.fit.all <- lme4::glmer(days_fever ~ age_s + sex + hiv + marg_caste + detection + 
                             prv_tx_ka + num_conslt + log(vill_inc_2017*1000,10) + 
                             IRS_2017_1 + block_endm_2017 + traveltime_s + (1|patient_id),
                           family = "poisson",
                           data = dat.fit)

summary(glm.fit.all)

################################################################################
################################################################################
# # ---------------------------------------------------------------------------- #
# # Mixed model - check spatial dependence of fitted random effects
# 
# dat.std <- standardize::standardize(days_fever ~ age + sex + hiv + marg_caste + detection + 
#                                       prv_tx_ka + num_conslt + vill_inc_2017 + 
#                                       IRS_2017_1 + block_endm_2017 + traveltime + (1|patient_id),
#                                     family = poisson,
#                                     data = dat.fit)
# 
# dat.std
# 
# glm.std <- lme4::lmer(dat.std$formula, dat.std$data)
# 
# summary(glm.std)
# # Linear mixed model fit by REML ['lmerMod']
# # Formula: days_fever ~ age + sex + hiv + marg_caste + detection + prv_tx_ka +      num_conslt + vill_inc_2017 + IRS_2017_1 + block_endm_2017 +  
# #   traveltime + (1 | vil_code)
# # Data: dat.std$data
# # 
# # REML criterion at convergence: 43351.4
# # 
# # Scaled residuals: 
# #   Min      1Q  Median      3Q     Max 
# # -2.4740 -0.4790 -0.1850  0.2154 10.5286 
# # 
# # Random effects:
# #   Groups   Name        Variance Std.Dev.
# # vil_code (Intercept)  334.7   18.29   
# # Residual             1030.4   32.10   
# # Number of obs: 4336, groups:  vil_code, 2337
# # 
# # Fixed effects:
# #   Estimate Std. Error t value
# # (Intercept)             52.0900     1.7359  30.008
# # age                      3.2314     0.5583   5.788
# # sexFemale                0.1319     0.5513   0.239
# # hivYes                   8.0337     1.5141   5.306
# # marg_casteYes            0.6262     0.6276   0.998
# # detectionACD            -4.1498     0.5830  -7.119
# # prv_tx_kaYes             0.4744     0.9744   0.487
# # num_conslt              11.0837     0.5716  19.390
# # vill_inc_2017           -0.5861     0.7246  -0.809
# # IRS_2017_1Yes            0.1192     0.7589   0.157
# # block_endm_2017Endemic  -3.9773     0.7007  -5.676
# # traveltime               0.5087     0.6474   0.786
# # 
# # Correlation of Fixed Effects:
# #   (Intr) age    sexFml hivYes mrg_cY dtcACD prv__Y nm_cns v__201 IRS_20 b__201
# # age         -0.107                                                                      
# # sexFemale    0.111  0.108                                                               
# # hivYes       0.749 -0.094  0.071                                                        
# # marg_castYs  0.148  0.092 -0.003  0.042                                                 
# # detectinACD  0.149  0.012 -0.020  0.073 -0.075                                          
# # prv_tx_kaYs  0.360 -0.067  0.025 -0.138  0.016  0.016                                   
# # num_conslt   0.014 -0.062  0.000 -0.014  0.108  0.020  0.051                            
# # vll_nc_2017  0.069  0.019 -0.003  0.019 -0.020 -0.018 -0.009  0.027                     
# # IRS_2017_1Y -0.202  0.038  0.005  0.023  0.012 -0.026 -0.015  0.017 -0.140              
# # blck__2017E  0.193 -0.035  0.003  0.073  0.085  0.046  0.012  0.002 -0.100 -0.170       
# # traveltime   0.008  0.013 -0.014  0.009 -0.071 -0.058  0.018  0.043 -0.016 -0.018  0.043



################################################################################
################################################################################