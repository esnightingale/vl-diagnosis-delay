################################################################################
# Description: 
# 
# 
################################################################################
################################################################################

set.seed(101)

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
M = 50

# Sample M test points
test.indices <- sample(1:nrow(data), size = M)

# Exclusion radius for spatial
r = 50

# Exceedance cutoff
C = 30

for (type in c("nonspatial", "spatial")){

# Initialise output - list of data frames which will be different sizes for each
# randomly generated test set
preds_all <- lapply(1:length(fits), function(x) data.frame())
n_witheld <- rep(NA, M)

# pdf(here::here(outdir, paste0(type,"_witheld_points.pdf")), height = 4, width = 6)

#------------------------------------------------------------------------------#

# M iterations 
  
for (m in 1:M){
  
  test.index <- test.indices[m]
 
    if (type == "spatial"){
 
      # Define buffer of radius r around sampled observation
      buff <- st_buffer(data[test.index,], dist = r)

      # Intersect full dataset with buffer
      withold <- st_intersects(data, buff, sparse = FALSE)
      print(paste0(sum(withold == TRUE), " observations excluded"))
      
      n_witheld[m] <- sum(withold == TRUE)

      # ggplot() +
      #   geom_sf(data = boundary, aes(geometry = geometry)) +
      #   geom_sf(data = data[!withold,], aes(geometry = geometry), col = "black", cex = 0.7) +
      #   geom_sf(data = data[withold,], aes(geometry = geometry), col = "red", cex = 0.7) +
      #   geom_sf(data = data[test.index,], aes(geometry = geometry), col = "green")
    
    }else if (type == "nonspatial"){
      
      withold <- rep(FALSE, nrow(data))
      withold[test.index] <- TRUE
      
    }
  
    dat.train <- data[!withold,]
    dat.test <- data[test.index,]
    obs.test <- data$delay[test.index]

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
      
      # Single vector so transpose to match dims
      X2 <- t(model.matrix(as.formula(paste("~ ",paste(covs, collapse = " + "))), 
                         data = dat.test)[,-1]) 
      
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

# dev.off()

names(preds_all) <- names(fits)
# saveRDS(preds_all, here::here(outdir, paste0(type, "_cv_preds.rds")))

cv_summ <- data.frame(Model = names(preds_all),
                      # MAE.cv = sapply(preds_all, function(x) mean(x$abs.err)),
                      MSE.cv = sapply(preds_all, function(x) mean(x$sq.err)),
                      # CPO.cv = sapply(preds_all, function(x) mean(x$cpo)),
                      logs.cv = sapply(preds_all, function(x) -log(mean(x$cpo))),
                      brier.cv = sapply(preds_all, function(x) mean((x$exc.prob - as.numeric(x$exc.obs))^2)))

cv.out <- list(preds = preds_all, summary = cv_summ, n_witheld = n_witheld)
saveRDS(cv.out, here::here(outdir, paste0(type, "_cv_out.rds")))

mod_compare <- read.csv(here::here("output","tables","mod_compare.csv")) 
tab_random_mav <- read.csv(here::here("output","tables","tab_random_mav_refnull.csv"))

mod_compare <- full_join(mod_compare, tab_random_mav, by = "Model") %>%
  full_join(cv_summ) %>%
  filter(Model %in% c("Null","All (IID only)","None","All")) %>% 
  mutate(Model.no = c(1, 3, 4, 2),
         dWAIC = WAIC - min(WAIC)) %>% 
  dplyr::select(Model.no, Model, WAIC, dWAIC, MSE, brier, MSE.cv, brier.cv, logs.cv) %>% 
  arrange(Model.no)

write.csv(mod_compare, here::here(outdir,paste0(type, "_mod_compare_full.csv")), row.names = FALSE)

# ---------------------------------------------------------------------------- #
# Plot CV predictions against observed

pdf(here::here(outdir,paste0(type, "_obs_vs_cvpred.pdf")), height = 7, width = 8)
purrr::imap(preds_all, function(x, nm) plot_preds(x, name = nm))
dev.off()

pdf(here::here(outdir,paste0(type, "_obs_vs_cvpred_exc30.pdf")), height = 7, width = 8)
purrr::imap(preds_all, function(x, nm) plot_exc(x, name = nm))
dev.off()

}

################################################################################
################################################################################