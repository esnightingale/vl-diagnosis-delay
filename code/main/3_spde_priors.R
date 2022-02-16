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
for (r in c(2, 10, 20, 50, 100, 500)){
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
# 4271 3574
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
# Range 2: SD 5 

# $`Range 2: SD 0.1`
# [1] 28612.76
# 
# $`Range 2: SD 1`
# [1] 28590.37
# 
# $`Range 2: SD 5`
# [1] 28570.42
# 
# $`Range 2: SD 10`
# [1] 28576.42
# 
# $`Range 2: SD 50`
# [1] 28577.87
# 
# $`Range 10: SD 0.1`
# [1] 28585.7
# 
# $`Range 10: SD 1`
# [1] 28582.9
# 
# $`Range 10: SD 5`
# [1] 28581.22
# 
# $`Range 10: SD 10`
# [1] 28586.81
# 
# $`Range 10: SD 50`
# [1] 28584.72
# 
# $`Range 20: SD 0.1`
# [1] 28583.89
# 
# $`Range 20: SD 1`
# [1] 28587.98
# 
# $`Range 20: SD 5`
# [1] 28580.59
# 
# $`Range 20: SD 10`
# [1] 28582.06
# 
# $`Range 20: SD 50`
# [1] 28584.18
# 
# $`Range 50: SD 0.1`
# [1] 28591.46
# 
# $`Range 50: SD 1`
# [1] 28582.67
# 
# $`Range 50: SD 5`
# [1] 28585.85
# 
# $`Range 50: SD 10`
# [1] 28576.9
# 
# $`Range 50: SD 50`
# [1] 28591.27
# 
# $`Range 100: SD 0.1`
# [1] 28581.13
# 
# $`Range 100: SD 1`
# [1] 28584.01
# 
# $`Range 100: SD 5`
# [1] 28593.13
# 
# $`Range 100: SD 10`
# [1] 28582.83
# 
# $`Range 100: SD 50`
# [1] 28585.56
# 
# $`Range 500: SD 0.1`
# [1] 28589.19
# 
# $`Range 500: SD 1`
# [1] 28586.42
# 
# $`Range 500: SD 5`
# [1] 28585.49
# 
# $`Range 500: SD 10`
# [1] 28585.21
# 
# $`Range 500: SD 50`
# [1] 28585.09

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
# $`Range 2: SD 0.1`
# mean 0.025quant 0.975quant
# Precision for id  1.198045   1.148106  1.2394506
# Range for s      39.194501  21.496682 67.8223201
# Stdev for s       0.351706   0.294429  0.4061599
# 
# $`Range 2: SD 1`
# mean 0.025quant 0.975quant
# Precision for id  1.1616047  1.1086969  1.2077516
# Range for s      29.9705983 19.7341518 45.0772375
# Stdev for s       0.3811959  0.3199359  0.4379687
# 
# $`Range 2: SD 5`
# mean 0.025quant 0.975quant
# Precision for id  1.1052089  1.0612900  1.1575003
# Range for s      14.3592746 10.5191355 20.3477236
# Stdev for s       0.3491691  0.2970065  0.4084372
# 
# $`Range 2: SD 10`
# mean 0.025quant 0.975quant
# Precision for id 1.1382216  1.0771021   1.196699
# Range for s      2.8355399  1.6904210   4.761299
# Stdev for s      0.9584598  0.5691429   1.433829
# 
# $`Range 2: SD 50`
# mean 0.025quant 0.975quant
# Precision for id 1.136901   1.079709   1.197136
# Range for s      1.376645   0.925455   2.058752
# Stdev for s      1.854363   1.200303   2.623499
# 
# $`Range 10: SD 0.1`
# mean 0.025quant 0.975quant
# Precision for id  1.1126194  1.0515852  1.1675301
# Range for s      38.8383270 23.8099893 62.8097824
# Stdev for s       0.2874407  0.2290429  0.3538554
# 
# $`Range 10: SD 1`
# mean 0.025quant 0.975quant
# Precision for id  1.1152280  1.0577696   1.173452
# Range for s      40.1711509 22.8806322  65.111731
# Stdev for s       0.3395531  0.2676592   0.439550
# 
# $`Range 10: SD 5`
# mean 0.025quant 0.975quant
# Precision for id  1.1155141  1.0610493   1.175782
# Range for s      43.1170954 23.9597789  74.738587
# Stdev for s       0.3459616  0.2737011   0.431282
# 
# $`Range 10: SD 10`
# mean 0.025quant 0.975quant
# Precision for id  1.1331005  1.0802886  1.1811815
# Range for s      52.7080854 30.8577787 84.5636474
# Stdev for s       0.3617348  0.2618755  0.4912171
# 
# $`Range 10: SD 50`
# mean 0.025quant 0.975quant
# Precision for id  1.1349214  1.0836869  1.1818135
# Range for s      55.6942854 26.4621105 96.6001983
# Stdev for s       0.3708073  0.2581311  0.5546141
# 
# $`Range 20: SD 0.1`
# mean 0.025quant  0.975quant
# Precision for id  1.0894832  1.0383494   1.1454294
# Range for s      81.0917972 52.5103458 111.1611874
# Stdev for s       0.2685568  0.2185472   0.3331646
# 
# $`Range 20: SD 1`
# mean 0.025quant 0.975quant
# Precision for id  1.1216382  1.0431366  1.1832949
# Range for s      36.9594334 24.3626898 56.1853266
# Stdev for s       0.3386045  0.2502041  0.4187178
# 
# $`Range 20: SD 5`
# mean 0.025quant 0.975quant
# Precision for id  1.1098806  1.0549276  1.1752018
# Range for s      49.1872954 24.7526903 83.5341272
# Stdev for s       0.3562171  0.2731213  0.4626895
# 
# $`Range 20: SD 10`
# mean 0.025quant  0.975quant
# Precision for id   1.0997112  1.0464235   1.1550611
# Range for s      144.1633608 77.4750146 235.9913874
# Stdev for s        0.5818093  0.3834856   0.8255426
# 
# $`Range 20: SD 50`
# mean 0.025quant 0.975quant
# Precision for id   1.095134   1.042867   1.150119
# Range for s      475.900679 351.934123 665.435903
# Stdev for s        1.521044   1.151447   2.030308
# 
# $`Range 50: SD 0.1`
# mean 0.025quant 0.975quant
# Precision for id  1.1374070  1.0879995  1.1800335
# Range for s      55.0260679 34.2805598 81.6548054
# Stdev for s       0.2854289  0.2164926  0.3674007
# 
# $`Range 50: SD 1`
# mean 0.025quant  0.975quant
# Precision for id  1.1004081   1.048021   1.1615705
# Range for s      77.5610910  30.351473 142.1524067
# Stdev for s       0.3784106   0.251525   0.5201863
# 
# $`Range 50: SD 5`
# mean 0.025quant 0.975quant
# Precision for id  1.1106171  1.0481159   1.165705
# Range for s      80.7953006 50.3170172 128.666378
# Stdev for s       0.4250843  0.2897382   0.562378
# 
# $`Range 50: SD 10`
# mean  0.025quant 0.975quant
# Precision for id   1.096997   1.0420533   1.170560
# Range for s      318.036048 168.1799147 550.764298
# Stdev for s        1.108299   0.6066277   1.858094
# 
# $`Range 50: SD 50`
# mean   0.025quant   0.975quant
# Precision for id    1.0232649    0.9684924    1.0700990
# Range for s      4758.9183450 3348.7579682 6771.7796756
# Stdev for s         0.2974017    0.1994630    0.4156199
# 
# $`Range 100: SD 0.1`
# mean 0.025quant  0.975quant
# Precision for id  1.0988030  1.0466283   1.1641650
# Range for s      86.1052196 55.7031623 120.8686531
# Stdev for s       0.3771583  0.2939856   0.4601575
# 
# $`Range 100: SD 1`
# mean 0.025quant  0.975quant
# Precision for id   1.0979530  1.0440789   1.1539288
# Range for s      137.3921928 57.6858572 236.0120991
# Stdev for s        0.5123867  0.2053289   0.9006916
# 
# $`Range 100: SD 5`
# mean 0.025quant  0.975quant
# Precision for id   1.093690  1.0038737   1.1587301
# Range for s      114.241921 66.8393485 189.8236313
# Stdev for s        0.478473  0.3363318   0.6706747
# 
# $`Range 100: SD 10`
# mean 0.025quant  0.975quant
# Precision for id  1.1040845  1.0511984   1.1617718
# Range for s      93.7932892 46.0700307 181.9577883
# Stdev for s       0.4202386  0.2981698   0.5695371
# 
# $`Range 100: SD 50`
# mean  0.025quant  0.975quant
# Precision for id    1.092652  1.04045349    1.147508
# Range for s      1069.594780 47.02749867 4049.130422
# Stdev for s         3.353959  0.07254165   14.941932
# 
# $`Range 500: SD 0.1`
# mean   0.025quant   0.975quant
# Precision for id 1.019099e+00 9.685581e-01 1.067129e+00
# Range for s      1.993132e+04 6.692386e+02 9.990379e+04
# Stdev for s      5.759631e-03 1.281536e-04 2.647154e-02
# 
# $`Range 500: SD 1`
# mean  0.025quant 0.975quant
# Precision for id   1.093966   1.0746805   1.111822
# Range for s      403.835389 298.5367829 519.813800
# Stdev for s        1.118522   0.6792547   1.653826
# 
# $`Range 500: SD 5`
# mean  0.025quant  0.975quant
# Precision for id   1.092314   1.0397057    1.147101
# Range for s      752.438647 178.8217222 1656.335547
# Stdev for s        2.073037   0.3349939    5.203954
# 
# $`Range 500: SD 10`
# mean 0.025quant  0.975quant
# Precision for id    1.092112   1.039448    1.146763
# Range for s      1551.978107 784.693961 3046.910126
# Stdev for s         4.371278   2.358620    8.216782
# 
# $`Range 500: SD 50`
# mean  0.025quant  0.975quant
# Precision for id    1.094315    1.041410    1.148883
# Range for s      3384.538577 1099.777098 8547.507494
# Stdev for s         8.893782    4.253304   16.361314

# Project fitted SPDEs
p.list <- purrr::imap(fits.spde, function(x, nm) plot_spde_mean(x$fit, nm, limits = c(-1,1))) # limit1 = c(-1,1), limit2 = c(0,1.5))) 
png(here::here(figdir, "compare_spde_priors_new.png"), height = 6000, width = 10000, res = 350) 
gridExtra::grid.arrange(grobs = p.list) #plyr::llply(fits, function(x) plot_spde(x$fit))
dev.off()

rlist::list.save(fits.spde, here::here(outdir, "fits_compare_priors.rds"))

# Select SPDE prior for modelling
spde <- inla.spde2.pcmatern(mesh = mesh, 
                            prior.range = c(10, 0.01), # P(range < U) = a
                            prior.sigma = c(5, 0.01), # P(sigma > U) = a
                            constr = TRUE)
saveRDS(spde, here::here("data/analysis","spde.rds"))

################################################################################
################################################################################