################################################################################
# Description: Univariate regression on covariates and assessment of empirical 
# variogram for spatial dependence
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit/exploratory/"
outdir <- "output/exploratory"

dat <- read_data()
mesh <- readRDS(here::here("data/analysis","mesh.rds"))

# Setup map context
boundary <- readRDS(here::here("data","geography","boundary.rds"))

boundary.spdf <- as_Spatial(boundary)

#------------------------------------------------------------------------------#
# Define formulae
make_f_spde <- function(r, s){

    f <- as.formula(sprintf(
      "y ~ -1 + Intercept + f(id, model = 'iid',
                              prior = 'pc.prec', 
                              param = c(10, 0.01)) + 
    f(s, model = inla.spde2.pcmatern(mesh = mesh, 
                                     prior.range = c(%s, 0.01), # P(range < U) = a
                                     prior.sigma = c(%s, 0.01), # P(sigma > U) = a
                                     constr = TRUE))", r, s)
    )
  
  return(f)
}

f.list <- list()
nm <- c()
for (r in c(10, 20, 50, 100, 500)){
  for(s in c(0.1, 1, 5, 10, 50)){
    f <- make_f_spde(r, s)
    f.list <- rlist::list.append(f.list, f)
    nm <- c(nm, paste0("Range ",r, ": SD ", s))
  }
}
names(f.list) <- nm

#------------------------------------------------------------------------------#
# Generate the index set for one of these SPDEs

spde <- inla.spde2.pcmatern(mesh = mesh, 
                    prior.range = c(10, 0.01), # P(range < U) = a
                    prior.sigma = c(1, 0.01), # P(sigma > U) = a
                    constr = TRUE)

indexs <- inla.spde.make.index("s", spde$n.spde)

#------------------------------------------------------------------------------#
# Generate a projection matrix A - projects the spatially continuous GRF from
# the observations to the mesh nodes

# For training points
coo <- st_coordinates(dat)
A <- inla.spde.make.A(mesh = mesh, loc = coo)

dim(A)
# 4271 3674
nrow(dat)

# ---------------------------------------------------------------------------- #
# Define model matrix based on all covariates of interest, removing automatic 
# intercept

# Covariates of interest
covs <- c("age_s","sex","comorb","prv_tx",
          "caste4_r","occ4_cat","poss_acd",
          "block_endm_2017", "IRS_2017","inc_2017_gt0", 
          "traveltime_s", "traveltime_t_s", "rain")

X <- model.matrix(as.formula(paste("~ ",paste(covs, collapse = " + "))), 
                  data = dat)[,-1] 

# ---------------------------------------------------------------------------- #
# Make stack

stk <- inla.stack(
  data = list(y = dat$delay),
  A = list(A, 1, 1),
  effects = list(s = indexs,  # the spatial index,
                 id = dat$id, # observation level index
                 data.frame(
                   Intercept = 1,
                   X) # covariate model matrix
  )
)

saveRDS(stk, here::here("data/analysis","stack.rds"))

# Also make stack for binomial output for later sensitivity analysis
stk_bin <- inla.stack(
  data = list(y = as.numeric(dat$gt30)),
  A = list(A, 1, 1),
  effects = list(s = indexs,  # the spatial index,
                 id = dat$id, # observation level index
                 data.frame(
                   Intercept = 1,
                   X) # covariate model matrix
  )
)

saveRDS(stk_bin, here::here("data/analysis","stack_bin.rds"))

# ---------------------------------------------------------------------------- #
# Fit baseline models with SPDE and IID to explain spatial structure
# Compare priors for SPDE

fits.spde <- lapply(f.list, init_inla, data.stack = stk, family = "poisson")

plyr::llply(fits.spde, function(x) x$fit$waic$waic)
which.min(plyr::llply(fits.spde, function(x) x$fit$waic$waic))
# Range 50: SD 10 

# $`Range 10: SD 0.1`
# [1] 28585.67
# 
# $`Range 10: SD 1`
# [1] 28583.77
# 
# $`Range 10: SD 5`
# [1] 28581.16
# 
# $`Range 10: SD 10`
# [1] 28605.66
# 
# $`Range 10: SD 50`
# [1] 28590.71
# 
# $`Range 20: SD 0.1`
# [1] 28584.25
# 
# $`Range 20: SD 1`
# [1] 28588.29
# 
# $`Range 20: SD 5`
# [1] 28580.57
# 
# $`Range 20: SD 10`
# [1] 28582.07
# 
# $`Range 20: SD 50`
# [1] 28584.16
# 
# $`Range 50: SD 0.1`
# [1] 28587.81
# 
# $`Range 50: SD 1`
# [1] 28582.66
# 
# $`Range 50: SD 5`
# [1] 28585.19
# 
# $`Range 50: SD 10`
# [1] 28576.97
# 
# $`Range 50: SD 50`
# [1] 28591.28
# 
# $`Range 100: SD 0.1`
# [1] 28583.8
# 
# $`Range 100: SD 1`
# [1] 28590.64
# 
# $`Range 100: SD 5`
# [1] 28593.07
# 
# $`Range 100: SD 10`
# [1] 28582.76
# 
# $`Range 100: SD 50`
# [1] 28585.91
# 
# $`Range 500: SD 0.1`
# [1] 28584.66
# 
# $`Range 500: SD 1`
# [1] 28586.44
# 
# $`Range 500: SD 5`
# [1] 28585.5
# 
# $`Range 500: SD 10`
# [1] 28585.21
# 
# $`Range 500: SD 50`
# [1] 28586.36

get_hyperpar <- function(x, nm){
  df <- x$fit$summary.hyperpar[,c(1,3,5)] %>% #
    mutate(prior = nm) %>%
    rownames_to_column("parameter")
}

df <- bind_rows(purrr::map2(fits.spde, names(fits.spde), get_hyperpar)) %>%
  tidyr::separate(prior, into = c("r", "s"), sep = ":", remove = FALSE) %>%
  mutate(r = as.numeric(gsub("Range ","",r)),
         s = as.numeric(gsub("SD ","",s)))

df %>%
  filter(parameter == "Precision for id") %>%
  ggplot(aes(x = r, y = s, col = mean)) +
  geom_jitter() +
  scale_color_viridis_c(trans = "identity") + 
  scale_x_continuous(trans = "log10") + 
  scale_y_continuous(trans = "log10") +
  labs(title = "Precision for id", x = "Prior P[Range < r] = 0.01", y = "Prior P[SD > s] = 0.01") -> p_prec_id

df %>%
  filter(parameter == "Range for s") %>%
  ggplot(aes(x = r, y = s, col = mean)) +
  geom_jitter() +
  scale_color_viridis_c(trans = "log10")+ 
  scale_x_continuous(trans = "log10") + 
  scale_y_continuous(trans = "log10") +
  labs(title = "Range for s", x = "Prior P[Range < r] = 0.01", y = "Prior P[SD > s] = 0.01") -> p_range_s

df %>%
  filter(parameter == "Stdev for s") %>% 
ggplot(aes(x = r, y = s, col = mean)) +
  geom_jitter() +
  scale_color_viridis_c(trans = "identity")+ 
  scale_x_continuous(trans = "log10") + 
  scale_y_continuous(trans = "log10") +
  labs(title = "SD for s", x = "Prior P[Range < r] = 0.01", y = "Prior P[SD > s] = 0.01") -> p_sd_s

png(here::here(figdir, "compare_spde_priors_hyperpars.png"), height = 600, width = 2000, res = 150) 
p_prec_id + p_range_s + p_sd_s
dev.off()

plyr::llply(fits.spde, function(x) x$fit$summary.hyperpar[,c(1,3,5)])
# $`Range 10: SD 0.1`
# mean 0.025quant 0.975quant
# Precision for id  1.1126244  1.0516783  1.1675064
# Range for s      38.8318464 23.8086032 62.8153502
# Stdev for s       0.2874403  0.2290421  0.3538357
# 
# $`Range 10: SD 1`
# mean 0.025quant 0.975quant
# Precision for id  1.1156040  1.0557623  1.1729665
# Range for s      39.6992954 22.5353848 63.6412170
# Stdev for s       0.3404149  0.2666802  0.4489077
# 
# $`Range 10: SD 5`
# mean 0.025quant 0.975quant
# Precision for id  1.1154967  1.0610823  1.1759036
# Range for s      43.1290311 23.9566259 74.8044705
# Stdev for s       0.3460031  0.2736759  0.4312472
# 
# $`Range 10: SD 10`
# mean 0.025quant 0.975quant
# Precision for id  1.0846938  0.9538754  1.1723652
# Range for s      52.6719774 30.8451210 84.4895615
# Stdev for s       0.3616714  0.2619195  0.4911669
# 
# $`Range 10: SD 50`
# mean   0.025quant   0.975quant
# Precision for id   1.02135762 9.660428e-01    1.0703637
# Range for s      892.12623634 2.690338e+02 2224.5666252
# Stdev for s        0.05668075 9.686432e-03    0.1974295
# 
# $`Range 20: SD 0.1`
# mean 0.025quant  0.975quant
# Precision for id  1.0896221  1.0382649   1.1448040
# Range for s      80.0284298 51.6795990 109.8380738
# Stdev for s       0.2651789  0.2164713   0.3281386
# 
# $`Range 20: SD 1`
# mean 0.025quant 0.975quant
# Precision for id  1.1213069  1.0417829  1.1834261
# Range for s      36.9902330 24.3752517 56.2749585
# Stdev for s       0.3386634  0.2502563  0.4187819
# 
# $`Range 20: SD 5`
# mean 0.025quant 0.975quant
# Precision for id  1.1098714  1.0549028  1.1752996
# Range for s      49.1951904 24.7603479 83.5400571
# Stdev for s       0.3562196  0.2731065  0.4627392
# 
# $`Range 20: SD 10`
# mean 0.025quant  0.975quant
# Precision for id   1.0996845   1.046367   1.1550070
# Range for s      144.1660190  77.267190 236.2661270
# Stdev for s        0.5810323   0.382587   0.8250711
# 
# $`Range 20: SD 50`
# mean 0.025quant 0.975quant
# Precision for id   1.095189   1.042938   1.150214
# Range for s      476.529463 353.147495 665.773612
# Stdev for s        1.523797   1.155395   2.029421
# 
# $`Range 50: SD 0.1`
# mean 0.025quant 0.975quant
# Precision for id  1.1027814  1.0355558   1.159082
# Range for s      63.1338034 40.8447112  90.788704
# Stdev for s       0.2953545  0.2111046   0.375904
# 
# $`Range 50: SD 1`
# mean 0.025quant  0.975quant
# Precision for id  1.1003987  1.0480171   1.1615916
# Range for s      78.2786295 31.2608734 142.8338252
# Stdev for s       0.3791411  0.2529505   0.5207186
# 
# $`Range 50: SD 5`
# mean 0.025quant  0.975quant
# Precision for id  1.1104604  1.0498655   1.1650648
# Range for s      80.9400445 50.3876698 128.9066366
# Stdev for s       0.4248773  0.2912753   0.5612059
# 
# $`Range 50: SD 10`
# mean  0.025quant 0.975quant
# Precision for id   1.096920   1.0420826   1.170205
# Range for s      318.099793 168.8113604 549.662656
# Stdev for s        1.108565   0.6087734   1.854747
# 
# $`Range 50: SD 50`
# mean   0.025quant   0.975quant
# Precision for id    1.0235434    0.9689547    1.0699531
# Range for s      4718.2399087 3403.7407883 6539.6219938
# Stdev for s         0.2959752    0.2076895    0.3975926
# 
# $`Range 100: SD 0.1`
# mean 0.025quant  0.975quant
# Precision for id  1.0997272  1.0471058   1.1544839
# Range for s      87.2394807 51.4891252 128.5592177
# Stdev for s       0.3783555  0.2966019   0.4590371
# 
# $`Range 100: SD 1`
# mean  0.025quant  0.975quant
# Precision for id   1.0974764   1.0236118   1.1541867
# Range for s      164.7281969 123.5891945 219.4352470
# Stdev for s        0.6086652   0.4678703   0.7854296
# 
# $`Range 100: SD 5`
# mean 0.025quant  0.975quant
# Precision for id   1.0938053  1.0041862   1.1587744
# Range for s      114.2839682 66.8339774 190.0344913
# Stdev for s        0.4788014  0.3364372   0.6711525
# 
# $`Range 100: SD 10`
# mean 0.025quant  0.975quant
# Precision for id  1.1016555  1.0493779   1.1599057
# Range for s      91.4796235 46.4115830 165.9535711
# Stdev for s       0.4185431  0.2997614   0.5820758
# 
# $`Range 100: SD 50`
# mean  0.025quant  0.975quant
# Precision for id    1.096474    1.042000    1.149680
# Range for s      1455.243811 1035.439341 2246.859352
# Stdev for s         4.529823    3.358466    6.381403
# 
# $`Range 500: SD 0.1`
# mean   0.025quant   0.975quant
# Precision for id 1.019068e+00 9.644123e-01 1.068751e+00
# Range for s      1.300044e+04 9.695584e+03 1.924811e+04
# Stdev for s      1.317468e-02 2.259121e-03 5.969245e-02
# 
# $`Range 500: SD 1`
# mean  0.025quant 0.975quant
# Precision for id   1.094006   1.0692091   1.117649
# Range for s      406.431135 289.1681868 553.050256
# Stdev for s        1.131130   0.6507759   1.833972
# 
# $`Range 500: SD 5`
# mean  0.025quant  0.975quant
# Precision for id   1.092310   1.0397103    1.147099
# Range for s      749.873010 175.6836293 1657.757184
# Stdev for s        2.065357   0.3286333    5.202154
# 
# $`Range 500: SD 10`
# mean 0.025quant  0.975quant
# Precision for id    1.092109   1.039442    1.146756
# Range for s      1551.639457 785.058561 3044.704282
# Stdev for s         4.370104   2.359405    8.211238
# 
# $`Range 500: SD 50`
# mean  0.025quant  0.975quant
# Precision for id    1.094845    1.040025    1.148007
# Range for s      3030.449747 1047.185031 6915.202515
# Stdev for s         8.837375    2.902067   19.867621

# Project fitted SPDEs
p.list <- purrr::imap(fits.spde, function(x, nm) plot_spde_mean(x$fit, nm, limits = c(-1,1))) # limit1 = c(-1,1), limit2 = c(0,1.5))) 
png(here::here(figdir, "compare_spde_priors_new.png"), height = 6000, width = 10000, res = 350) 
gridExtra::grid.arrange(grobs = p.list)
dev.off()

rlist::list.save(fits.spde, here::here(outdir, "fits_compare_priors.rds"))

# Select SPDE prior for modelling
spde <- inla.spde2.pcmatern(mesh = mesh, 
                            prior.range = c(50, 0.01), # P(range < U) = a
                            prior.sigma = c(10, 0.01), # P(sigma > U) = a
                            constr = TRUE)
saveRDS(spde, here::here("data/analysis","spde.rds"))

################################################################################
################################################################################