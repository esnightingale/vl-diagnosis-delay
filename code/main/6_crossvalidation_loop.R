################################################################################
# Description: 
# 
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit"
outdir <- "output/cross-validation"

fits <- readRDS(here::here("output", "fits_covs_all.rds"))

data <- readRDS(here::here("data/analysis","dat_nona.rds")) %>%
    st_transform(crs = st_crs(7759))
coop <- readRDS(here::here("data/analysis","coop.rds"))
mesh <- readRDS(here::here("data/analysis","mesh.rds"))
spde <- readRDS(here::here("data/analysis","spde.rds"))
indexs <- inla.spde.make.index("s", spde$n.spde)

# partitions <- readRDS(here::here("data/analysis","stack_partitions.rds"))
# stk.full <- readRDS(here::here("data/analysis","stack.rds"))

# Set up map context
boundary <- readRDS(here::here("data","geography","bihar_block.rds")) %>%
  sf::st_transform(7759) %>%
  sf::st_union()

# Covariates of interest
covs <- c("age_s","sex","hiv",
          "marg_caste","occupation",
          "block_endm_2017", "IRS_2017_1","vill_inc_2017_gt0", "prv_tx",
          "traveltime_s", "rain",
          "detection")

# Number of iterations
M = 10

# Radius of exclusion - 10km
r = 1e4

# Exceedance cutoff
C = 30

# Initialise output - list of data frames which will be different sizes for each
# randomly generated test set
preds_all <- list(data.frame(),data.frame(),data.frame(),data.frame(),
                  data.frame(),data.frame(),data.frame(),data.frame(),
                  data.frame(),data.frame())
mae <- c()
rmse <- c()
  
#------------------------------------------------------------------------------#

# M iterations 
  
  for (m in 1:M){

    # Randomly sample one starting location
    samp <- sample_n(data, 1)
    
    # Define buffer of radius r around sampled observation
    buff <- st_buffer(samp, dist = r)
    
    # Intersect full dataset with buffer
    test.index <- st_intersects(data, buff, sparse = FALSE)
    print(summary(test.index))
  
    dat.train <- data[!test.index,]
    dat.test <- data[test.index,]
    obs.test <- data$delay[test.index]
    
    #-----------------------------------------------------------------------------#
    # DEFINE DATA STACK
  
    # For training points
    coo <- st_coordinates(dat.train)
    A <- inla.spde.make.A(mesh = mesh, loc = coo)
    
    dim(A)
    nrow(dat.train)
    
    # Define model matrix based on all covariates of interest, removing automatic 
    # intercept
    X1 <- model.matrix(as.formula(paste("~ ",paste(covs, collapse = " + "))), 
                       data = dat.train)[,-1] 
    
    # Training stack
    stk.train <- inla.stack(
      tag = "train",
      data = list(y = dat.train$delay),
      A = list(A, 1, 1),
      effects = list(s = indexs,  # the spatial index,
                     v = dat.train$v,
                     data.frame(  # covariates
                       Intercept = 1, 
                       X1)))
  
      # For withheld test points
      coot <- st_coordinates(dat.test)
      At <- inla.spde.make.A(mesh = mesh, loc = coot)
      
      dim(At)
      nrow(dat.test)
      
      X2 <- model.matrix(as.formula(paste("~ ",paste(covs, collapse = " + "))), 
                         data = dat.test)[,-1] 
      
      # Testing stack
      stk.test <- inla.stack(
        tag = "test",
        data = list(y = NA),
        A = list(At, 1, 1),
        effects = list(s = indexs,  # the spatial index,
                       v = dat.test$v,
                       data.frame(  # covariates
                         Intercept = 1, 
                         X2)))
      
      # Stack for smooth prediction from intercept and fitted spatial field (no covariates)
      Ap <- inla.spde.make.A(mesh = mesh, loc = coop)
      
      # Prediction stack
      stk.pred <- inla.stack(
        tag = "pred",
        data = list(y = NA),
        A = list(Ap, 1),
        effects = list(s = indexs,
                       data.frame(
                         Intercept = rep(1, nrow(coop)))
        )
      )
      
      # Combine stacks
      data.stack <- inla.stack(stk.train, stk.test, stk.pred)

#------------------------------------------------------------------------------#
# Refit with sub-sampled training data

    for (i in seq_along(fits)){
      
    refit <- inla(fits[[i]]$f,
                  family = "nbinomial",
                  data = inla.stack.data(data.stack),
                  control.predictor = list(
                    compute = TRUE, link = 1,
                    A = inla.stack.A(data.stack)),
                  control.compute = list(dic = FALSE, 
                                         waic = FALSE, 
                                         config = TRUE,
                                         cpo = FALSE),
                  control.mode = list(result = fits[[i]]$fit, 
                                      restart = TRUE),
                  control.fixed = list(mean = 0, 
                                       prec = 0.1, 
                                       mean.intercept = 0, 
                                       prec.intercept = 0.1),
                  verbose = TRUE)
    
    # Identify test indices in data stack
    idx <- data.stack$data$index$test
    print(length(idx))
    
    # Extract summary stats of fitted values at these indices
    temp <- data.frame(# Match to observed values from the full stack at the specified indices
                       obs = obs.test,
                       ll = refit$summary.fitted.values[idx, "0.025quant"],
                       med = refit$summary.fitted.values[idx, "0.5quant"],
                       mean = refit$summary.fitted.values[idx, "mean"],
                       ul = refit$summary.fitted.values[idx, "0.975quant"],
                       exc.prob = sapply(refit$marginals.fitted.values[idx],
                                         FUN = function(marg){1-inla.pmarginal(q = C, marginal = marg)})) %>%
            dplyr::mutate(abs.err = abs(mean - obs),
                          sq.err = (mean - obs)^2,
                          exc.obs = (obs > C))
    
    # Calculate condtional predictive ordinate - density of the posterior marginal at the observed value
    temp$cpo <- mapply(function(x, m){inla.dmarginal(x = x, marginal = m)},
                       temp$obs,
                       refit$marginals.fitted.values[idx])
    
    preds_all[[i]] <- bind_rows(preds_all[[i]],temp)
  
  }
  
}

names(preds_all) <- names(fits)
saveRDS(preds_all, here::here(outdir, "cv_preds.rds"))

cv_summ <- data.frame(Model = names(fits),
                      MAE.cv = sapply(preds_all, function(x) mean(x$abs.err)),
                      RMSE.cv = sapply(preds_all, function(x) sqrt(mean(x$sq.err))),
                      CPO.cv = sapply(preds_all, function(x) -log(mean(x$cpo))))

cv.out <- list(preds = preds_all, summary = cv_summ)
saveRDS(cv.out, here::here(outdir, "cv_out.rds"))


mod_compare <- read.csv(here::here("output", "mod_compare.csv")) 
tab_random_mav <- read.csv(here::here("output", "tab_random_mav.csv"))

mod_compare <- full_join(mod_compare, tab_random_mav, by = "Model") %>%
  full_join(cv_summ)

write.csv(mod_compare, here::here(outdir,"mod_compare_full.csv"))

# ---------------------------------------------------------------------------- #
# Plot CV predictions against observed

pdf(here::here(figdir,"obs_vs_cvpred.pdf"), height = 7, width = 8)
purrr::imap(preds_all, function(x, nm) plot_preds(x, name = nm))
dev.off()

pdf(here::here(figdir,"obs_vs_cvpred_exc30.pdf"), height = 7, width = 8)
purrr::imap(preds_all, function(x, nm) plot_exc(x, name = nm))
dev.off()

################################################################################
################################################################################