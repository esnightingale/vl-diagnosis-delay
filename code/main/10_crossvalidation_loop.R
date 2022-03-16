################################################################################
# Description: 
# 
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit"
outdir <- "output/cross-validation"

# Final model fits with/without covariates
fits <- readRDS(here::here("output", "fits_final.rds"))[c("All (IID only)", "None", "All")]

# Also read null model with no covariates or spatial field (pure IID)
fit.null <- readRDS(here::here("output/univariate", "fit_null.rds"))

# Combine fits and rename
fits <- rlist::list.prepend(fits, fit.null)
names(fits)[1] <- "Null"
# names(fits) <- c("Null","All (IID only)","None", "All")

data <- read_data() #readRDS(here::here("data/analysis","dat_nona.rds")) 

# coop <- readRDS(here::here("data/analysis","coop.rds"))
mesh <- readRDS(here::here("data/analysis","mesh.rds"))
spde <- readRDS(here::here("data/analysis","spde.rds"))
indexs <- inla.spde.make.index("s", spde$n.spde)

# Set up map context
boundary <- readRDS(here::here("data","geography","boundary.rds")) 

# Covariates of interest
covs <- c("age_s","comorb", "poss_acd",
          "block_endm_2017", "inc_2017_gt0", 
          "traveltime_t_s")

# Number of iterations
M = 10

# Percentage withheld for testing
type = "nonspatial"
test_prop = 0.1

# Radius of exclusion - 20km
# type = "spatial"
# r = 10

# Exceedance cutoff
C = 30

# Initialise output - list of data frames which will be different sizes for each
# randomly generated test set
preds_all <- lapply(1:length(fits), function(x) data.frame())

#------------------------------------------------------------------------------#

# M iterations 
  
for (m in 1:M){
 
    if (type == "spatial"){
      
      # Randomly sample one starting location
      samp <- sample_n(data, 1)
      
      # Define buffer of radius r around sampled observation
      buff <- st_buffer(samp, dist = r)
      
      # Intersect full dataset with buffer
      test.index <- st_intersects(data, buff, sparse = FALSE)
      
    }else if (type == "nonspatial"){
      
      # Randomly sample X% of observations
      test.index <- sample(c(TRUE,FALSE), size = nrow(data), prob = c(test_prop, 1-test_prop), replace = TRUE)

    }
      
    print(summary(test.index))
  
    dat.train <- data[!test.index,]
    dat.test <- data[test.index,]
    obs.test <- data$delay[test.index]
    
    print(summary(dat.train))
    print(summary(dat.test))
    
    #--------------------------------------------------------------------------#
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
                     id = dat.train$id,
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
                       id = dat.test$id,
                       data.frame(  # covariates
                         Intercept = 1, 
                         X2)))
      
      # Combine stacks
      data.stack <- inla.stack(stk.train, stk.test) 

#------------------------------------------------------------------------------#
# Refit with sub-sampled training data

    for (i in seq_along(fits)){
      
    refit <- inla(fits[[i]]$f,
                  family = "poisson",
                  data = inla.stack.data(data.stack),
                  control.predictor = list(
                    compute = TRUE, link = 1,
                    A = inla.stack.A(data.stack)),
                  control.compute = list(dic = FALSE, 
                                         waic = FALSE, 
                                         config = TRUE,
                                         cpo = TRUE,
                                         return.marginals.predictor = TRUE),
                  # control.mode = list(result = fits[[i]]$fit, 
                  #                     restart = TRUE),
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
                       pred = refit$summary.fitted.values[idx, "mean"],
                       ul = refit$summary.fitted.values[idx, "0.975quant"],
                       exc.prob = sapply(refit$marginals.fitted.values[idx],
                                         FUN = function(marg){1-inla.pmarginal(q = C, marginal = marg)})) %>%
            dplyr::mutate(abs.err = abs(pred - obs),
                          sq.err = (pred - obs)^2,
                          exc.obs = obs > C)
    
    # Calculate condtional predictive ordinate - density of the posterior marginal at the observed value
    temp$cpo <- mapply(function(x, m){inla.dmarginal(x = x, marginal = m)},
                       temp$obs,
                       refit$marginals.fitted.values[idx])
    
    preds_all[[i]] <- bind_rows(preds_all[[i]],temp)
  
  }
  
}

names(preds_all) <- names(fits)
# saveRDS(preds_all, here::here(outdir, paste0(type, "_cv_preds.rds")))

cv_summ <- data.frame(Model = names(preds_all),
                      # MAE.cv = sapply(preds_all, function(x) mean(x$abs.err)),
                      MSE.cv = sapply(preds_all, function(x) mean(x$sq.err)),
                      # CPO.cv = sapply(preds_all, function(x) mean(x$cpo)),
                      logs.cv = sapply(preds_all, function(x) -log(mean(x$cpo))),
                      brier.cv = sapply(preds_all, function(x) mean((x$exc.prob - as.numeric(x$exc.obs))^2)))

cv.out <- list(preds = preds_all, summary = cv_summ)
saveRDS(cv.out, here::here(outdir, paste0(type, "_cv_out.rds")))

mod_compare <- read.csv(here::here("output","mod_compare.csv")) 
tab_random_mav <- read.csv(here::here("output","tab_random_mav_refnull.csv"))

mod_compare <- full_join(mod_compare, tab_random_mav, by = "Model") %>%
  full_join(cv_summ)

write.csv(mod_compare, here::here(outdir,paste0(type, "_mod_compare_full.csv")), row.names = FALSE)

# ---------------------------------------------------------------------------- #
# Plot CV predictions against observed

pdf(here::here(figdir,paste0(type, "_obs_vs_cvpred.pdf")), height = 7, width = 8)
purrr::imap(preds_all, function(x, nm) plot_preds(x, name = nm))
dev.off()

pdf(here::here(figdir,paste0(type, "_obs_vs_cvpred_exc30.pdf")), height = 7, width = 8)
purrr::imap(preds_all, function(x, nm) plot_exc(x, name = nm))
dev.off()

################################################################################
################################################################################