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
spde <- readRDS(here::here("data/analysis","spde.rds"))
indexs <- inla.spde.make.index("s", spde$n.spde)

# Setup map context
boundary <- readRDS(here::here("data","geography","boundary.rds"))

boundary.spdf <- as_Spatial(boundary)

#------------------------------------------------------------------------------#
# Define formulae
make_f_spde <- function(r, s){

    f <- as.formula(sprintf(
      "y ~ -1 + Intercept + 
      age_s + comorb1 + poss_acdTRUE + 
      block_endm_2017TRUE + inc_2017_gt0TRUE + traveltime_t_s +
      f(id, model = 'iid',
            prior = 'pc.prec', 
            param = c(1, 0.01)) + 
      f(s, model = inla.spde2.pcmatern(mesh = mesh, 
                                       prior.range = c(%s, 0.01), # P(range < U) = a
                                       prior.sigma = c(%s, 0.05), # P(sigma > U) = a
                                       constr = TRUE))", r, s)
    )
  
  return(f)
}

f.list <- list()
nm <- c()
for (r in c(10, 50, 100)){
  for(s in c(0.5, 1, 2, 5)){
    f <- make_f_spde(r, s)
    f.list <- rlist::list.append(f.list, f)
    nm <- c(nm, paste0("Range ",r, ": SD ", s))
  }
}
names(f.list) <- nm


f1 <- y ~ -1 + Intercept + 
  age_s + comorb1 + poss_acdTRUE + 
  block_endm_2017TRUE + inc_2017_gt0TRUE + traveltime_t_s +
  f(id, model = 'iid',
    prior = 'pc.prec', 
    param = c(10, 0.01)) + 
  f(s, model = inla.spde2.pcmatern(mesh = mesh, 
                                   prior.range = c(100, 0.5), # P(range < U) = a
                                   prior.sigma = c(2, 0.05), # P(sigma > U) = a
                                   constr = TRUE))

f2 <- y ~ -1 + Intercept + 
      age_s + comorb1 + poss_acdTRUE + 
      block_endm_2017TRUE + inc_2017_gt0TRUE + traveltime_t_s +
      f(id, model = 'iid',
            prior = 'pc.prec', 
            param = c(10, 0.01)) + 
      f(s, model = inla.spde2.pcmatern(mesh = mesh, 
                                       prior.range = c(100, 0.99), # P(range < U) = a
                                       prior.sigma = c(2, 0.05), # P(sigma > U) = a
                                       constr = TRUE))

f3 <- y ~ -1 + Intercept + 
  age_s + comorb1 + poss_acdTRUE + 
  block_endm_2017TRUE + inc_2017_gt0TRUE + traveltime_t_s +
  f(id, model = 'iid',
    prior = 'pc.prec', 
    param = c(10, 0.01)) + 
  f(s, model = inla.spde2.pcmatern(mesh = mesh, 
                                   prior.range = c(100, 0.99), # P(range < U) = a
                                   prior.sigma = c(5, 0.05), # P(sigma > U) = a
                                   constr = TRUE))

f.list <- rlist::list.append(f.list, `Range 100, p0.5; SD 2` = f1,
                             `Range 100, p0.99; SD 2` = f2, 
                             `Range 100, p0.99; SD 5` = f3)

# ---------------------------------------------------------------------------- #
# Fit baseline models with SPDE and IID to explain spatial structure
# Compare priors for SPDE

fits.spde <- lapply(f.list, init_inla, data.stack = stk, family = "poisson")

plyr::llply(fits.spde, function(x) x$fit$waic$waic)
which.min(plyr::llply(fits.spde, function(x) x$fit$waic$waic))

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
# $`Range 10: SD 0.5`
# mean 0.025quant 0.975quant
# Precision for id  1.223577  1.1080222  1.4376249
# Range for s      37.790133 24.9729611 56.1353291
# Stdev for s       0.311031  0.2420228  0.3851263
# 
# $`Range 10: SD 1`
# mean 0.025quant 0.975quant
# Precision for id  1.1698131  1.1138703  1.2215178
# Range for s      30.3161250 16.2221785 62.6988228
# Stdev for s       0.3194373  0.2420085  0.3827742
# 
# $`Range 10: SD 2`
# mean 0.025quant 0.975quant
# Precision for id  1.1652742  1.1016950  1.2237562
# Range for s      47.2190150 26.1079456 83.4729001
# Stdev for s       0.3238994  0.2293042  0.4166023
# 
# $`Range 10: SD 5`
# mean 0.025quant  0.975quant
# Precision for id  1.1605315  1.1050153   1.2250942
# Range for s      52.8460675 28.9115535 103.9970123
# Stdev for s       0.3226267  0.2473996   0.4193243
# 
# $`Range 50: SD 0.5`
# mean 0.025quant  0.975quant
# Precision for id  1.1521756  1.0849736   1.2106126
# Range for s      79.3138499 48.3717017 123.9009081
# Stdev for s       0.3601556  0.2424272   0.5357295
# 
# $`Range 50: SD 1`
# mean 0.025quant  0.975quant
# Precision for id  1.1551103  1.0899402   1.2132572
# Range for s      76.1976399 46.6969618 123.2796715
# Stdev for s       0.3728776  0.2552521   0.5448813
# 
# $`Range 50: SD 2`
# mean 0.025quant  0.975quant
# Precision for id   1.1476399  1.0930767   1.2099411
# Range for s      135.5379573 72.1276480 225.3464827
# Stdev for s        0.5220168  0.3524045   0.7263679
# 
# $`Range 50: SD 5`
# mean  0.025quant 0.975quant
# Precision for id   1.1457103   1.0894110   1.202354
# Range for s      221.3504729 121.1483168 366.527589
# Stdev for s        0.7223034   0.4510913   1.088742
# 
# $`Range 100: SD 0.5`
# mean 0.025quant  0.975quant
# Precision for id   1.1477930  1.0929678   1.2146598
# Range for s      114.6191699 69.4893399 181.4365429
# Stdev for s        0.4444455  0.3056124   0.6160562
# 
# $`Range 100: SD 1`
# mean 0.025quant  0.975quant
# Precision for id   1.1489675   1.091102   1.2059108
# Range for s      150.8666015  82.779050 249.5180680
# Stdev for s        0.5592905   0.342972   0.8419507
# 
# $`Range 100: SD 2`
# mean  0.025quant 0.975quant
# Precision for id   1.1437319   1.0888321   1.201046
# Range for s      224.7586990 152.1377381 322.098607
# Stdev for s        0.7197768   0.4862455   1.068436
# 
# $`Range 100: SD 5`
# mean 0.025quant  0.975quant
# Precision for id   1.1508212  1.0958703   1.2104867
# Range for s      102.4707514 44.0539313 201.1489868
# Stdev for s        0.4211575  0.2726373   0.5966979
# 
# $`Range 100, p0.5; SD 2`
# mean 0.025quant  0.975quant
# Precision for id   1.1432859  1.0890554   1.2084186
# Range for s      132.3781638 88.1083372 188.9193219
# Stdev for s        0.4853603  0.3409108   0.6975968
# 
# $`Range 100, p0.99; SD 2`
# mean 0.025quant  0.975quant
# Precision for id  1.1550620  1.0934100   1.2124393
# Range for s      70.7525055 39.0805558 113.4361445
# Stdev for s       0.3677498  0.2635034   0.4900774
# 
# $`Range 100, p0.99; SD 5`
# mean 0.025quant  0.975quant
# Precision for id  1.1527112   1.097615   1.2179962
# Range for s      76.4257597  48.239219 114.1823179
# Stdev for s       0.3974952   0.281831   0.5279172


# Project fitted SPDEs
p.list <- purrr::imap(fits.spde, function(x, nm) plot_spde_mean(x$fit, nm, limits = c(-1,1))) # limit1 = c(-1,1), limit2 = c(0,1.5))) 
png(here::here(figdir, "compare_spde_priors.png"), height = 6000, width = 10000, res = 350) 
gridExtra::grid.arrange(grobs = p.list)
dev.off()

rlist::list.save(fits.spde, here::here(outdir, "fits_compare_priors.rds"))

# Select SPDE prior for modelling
spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(10, 0.01), # P(range < U) = a
                            prior.sigma = c(2, 0.05), # P(sigma > U) = a
                            constr = TRUE)
saveRDS(spde, here::here("data/analysis","spde.rds"))

################################################################################
################################################################################