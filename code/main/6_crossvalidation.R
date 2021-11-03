################################################################################
# Description: 
# 
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit"
outdir <- "output/cross-validation"

fits.init <- readRDS(here::here("output", "fits_covs_all.rds"))

dat <- readRDS(here::here("data/analysis","dat_nona.rds")) %>%
    st_transform(crs = st_crs(7759))
mesh <- readRDS(here::here("data/analysis","mesh.rds"))
spde <- readRDS(here::here("data/analysis","spde.rds"))
indexs <- inla.spde.make.index("s", spde$n.spde)

partitions <- readRDS(here::here("data/analysis","stack_partitions.rds"))
stk.full <- readRDS(here::here("data/analysis","stack.rds"))

# ---------------------------------------------------------------------------- #

#--------------------------------------#
# Refit with sub-sampled training data #
#--------------------------------------#

refit_validate <- function(res, stack.validate, stack.full){
  
  # Partitioned data stack and index of which obs have been withheld for validation
  stack <- stack.validate$stack
  p.index <- stack.validate$p.index
  
  fit <- inla(res$f,
              family = "nbinomial",
              data = inla.stack.data(stack),
              control.predictor = list(
                compute = TRUE, link = 1,
                A = inla.stack.A(stack)),
              control.compute = list(dic = FALSE, 
                                     waic = FALSE, 
                                     config = TRUE,
                                     cpo = FALSE),
              control.mode = list(result = res$fit, restart = TRUE),
              control.fixed = list(mean = 0, 
                                   prec = 0.1, 
                                   mean.intercept = 0, 
                                   prec.intercept = 0.1),
              verbose = TRUE)
  
  # Identify test indices in data stack
  index <- stack$data$index$test
  print(length(index))
  
  # Extract summary stats of fitted values at these indices
  pred <- data.frame(ll = fit$summary.fitted.values[index, "0.025quant"],
                     med = fit$summary.fitted.values[index, "0.5quant"],
                     ul = fit$summary.fitted.values[index, "0.975quant"],
                     exc.prob = sapply(fit$marginals.fitted.values[index],
                                       FUN = function(marg){1-inla.pmarginal(q = 30, marginal = marg)}),
                     # Match to observed values from the full stack at the specified indices
                     obs = stack.full$stack$data$data$y[p.index]) %>%
          dplyr::mutate(abs.err = abs(med - obs),
                        sq.err = (med - obs)^2)
  
  return(pred)
}

# rf_val <- refit_validate(fits.all$`All (IID only)`, partitions[[1]], stk.full)
# sqrt(mean(rf_val$sq.err))

pred_crossval <- function(res, partitions){
  
  cv.pred <- lapply(partitions, function(x) refit_validate(res, x, stack.full = stk.full))
  
  cv.summ <- data.frame(RMSE = )
  return(cv.pred)
  
}

cv1 <- pred_crossval(fits.init$`All (IID only)`, partitions)

crossval <- plyr::llply(fits.init, refit_inla)

# Starting from initial fits using all data, randomly sub-sample a test set and
# refit with these outcome values set to NA
run_crossval <- function(res.init, M = 10){
  res <- lapply(1:M,
                function(x) {
                  print(x)
                  refit_inla(x)
                }
  )
  return(res)
}


# ---------------------------------------------------------------------------- #

# Fitted spatial field from the null model:
out.field <- INLA::inla.spde2.result(fits.init[[1]]$fit,'s', spde, do.transf=TRUE)
# Range of fitted field
range.out <- INLA::inla.emarginal(function(x) x, out.field$marginals.range.nominal[[1]])
# Define a sensible exclusion radius around removed point
rad <- range.out*0.1 # ~10km

dat %>%
  bind_cols(as.data.frame(st_coordinates(dat))) %>%
  st_drop_geometry() %>%
  rename(y = days_fever) -> data.frame

coords <- data.frame(long = c(rnorm(70), rnorm(30, 3)), lat = rnorm(100))
PA <- rep(c(0, 1), each = 50)
x <- data.frame(x1 = rnorm(100),x2 = c(rnorm(70), rnorm(30, 2)))# x1 no relat., x2 pos. relat.
dataf1 <- sp::SpatialPointsDataFrame(coords = coords, data = data.frame(y = PA, x))
mesh <- INLA::inla.mesh.2d(loc = sp::coordinates(dataf1),max.edge = c(3, 3),cutoff = 1.3)
spde <- INLA::inla.spde2.matern(mesh, alpha=2)#SPDE model is defined
A <- INLA::inla.spde.make.A(mesh, loc = sp::coordinates(dataf1))#projector matrix
dataframe <- data.frame(dataf1) #generate dataframe with response and covariate
modform<-stats::as.formula(paste('y ~ -1+ x1 + x2 + y.intercept + f(spatial.field, model=spde)'))
stk <- INLA::inla.stack(data=list(y=dataframe$y),A=list(A, 1),
                        effects=list(list(spatial.field=1:spde$n.spde),
                                     list(y.intercept=rep(1,length(dataframe$y)),covariate=dataframe[c(-1)])),tag='est')
out <- INLA::inla(modform, family='normal',Ntrials = 1, data=INLA::inla.stack.data(stk, spde=spde),
                  control.predictor = list(A=INLA::inla.stack.A(stk),link=1),
                  control.compute = list( config=TRUE),control.inla = list(int.strategy='eb'))
out.field <- INLA::inla.spde2.result(out,'spatial.field', spde, do.transf=TRUE)
range.out <- INLA::inla.emarginal(function(x) x, out.field$marginals.range.nominal[[1]])

# parameters for the SLOO process
ss <- 1#sample size to process (number of SLOO runs)
rad <- range.out*0.15#define the radius of the spatial buffer surrounding the removed point
mesh <- mesh#use the mesh of the model
dataframe <- dataframe#dataframe with response 'y' and covariates 'x1', 'x2'
dataframe$y <- round(runif(length(dataframe$y), 1, 12))#create positive discrete response
modform <- stats::as.formula(paste('y ~ -1+ y.intercept + x1 + f(spatial.field, model=spde)'))
family <- list('gamma')#one model
ntrials <- rep(round(max(dataframe$y,na.rm=TRUE)*2),length(dataframe$y))# create ntrials for Binomial family
alpha <- 0.05#rmse and mae confidence intervals (1-alpha)



run_sloo <- function(res.init, r) {
  
  sloo <- inlasloo_adapt(data.frame,
                         long = "X",
                         lat = "Y",
                         y = "y",
                         family = "nbinomial",
                         ss = 1,
                         rad = r,
                         modform = res.init$f,
                         mesh = mesh,
                         print = TRUE,
                         plot = TRUE,
                         mae = TRUE,
                         ds = TRUE,
                         sqroot = FALSE,
                         control.mode = list(result = res.init$fit, restart = TRUE),
                         verbose = TRUE)
  
  return(sloo)
  
}

sloo1 <- run_sloo(fits.init$None, r = rad)

fits.sloo <- plyr::llply(fits.all, run_sloo, r = rad)
saveRDS(here::here(outdir, "cross-validation/fits_all_sloo.rds"))

for (m in seq_along(fits.all)) {
  
  sloo <- run_sloo(fits.all[[m]], r = rad)
  saveRDS(here::here(outdir, paste0("cross-validation/sloo_",names(fits.all)[m],".rds")))
  
}



#----------------#
# Spatial LOO CV #
#----------------#
# 
# out.field <- INLA::inla.spde2.result(fits.init$None$fit$fit,'v', spde, do.transf = TRUE)
# range.out <- INLA::inla.emarginal(function(x) x, out.field$marginals.range.nominal[[1]])
# # 0.6886794
# 
# run_sloo <- function(res.init) {
#   
#   sloo <- INLAutils::inlasloo(st_drop_geometry(dat),
#                               long = "longitude",
#                               lat = "latitude",
#                               y = "days_fever",
#                               family = "nbinomial",
#                               ss = 10,
#                               rad = 0.1,
#                               modform = res.init$f,
#                               mesh = mesh,
#                               print = TRUE,
#                               plot = TRUE,
#                               control.mode = list(result = res.init$fit, restart = TRUE))
#   
#   return(sloo)
#   
# }
# 
# fits.sloo <- plyr::llply(covs.list, run_sloo)
# 
# saveRDS(here::here(outdir, "cross-validation/fits_sloo.rds"))
# 
# for (m in seq_along(fits.init)) {
#   
#   sloo <- run_sloo(fits.init[[m]])
#   saveRDS(here::here(outdir, paste0("cross-validation/sloo_",names(fits.init)[m],".rds")))
#   
# }
# 
# In list, tends to crash (parallellise?)
# fits.xval <- plyr::llply(fits.init, run_crossval, M = M)
# saveRDS(fits.xval, here::here(outdir, "fits_xval.rds"))

# In loop, save as go along
# for (model in c("IID",
#                 # "SPDE",
#                 "IID_SPDE")) {
# 
#   saveRDS(run_crossval(fits.init[[model]], M = M),
#           here::here(outdir,
#                      paste0("fits_xval_",model,".rds")))
# }
# 
# ## Pre-define M resampled training/test sets and refit all models with these
# 
# test.percent = 0.2
# lapply(1:M, function(x) sample(1:nrow(dat.fit), floor(nrow(dat.fit)*test.percent)))
# test.idx <- sample(1:nrow(dat.fit), floor(nrow(dat.fit)*test.percent))
# 
# dat.train <- dat.fit
# dat.train$days_fever[test.idx] <- NA
# 
# # spatial.field <- INLA::inla.spde2.result(fit.base.spde$fit,'v', spde, do.transf = TRUE)
# # range.spatial <- INLA::inla.emarginal(function(x) x, spatial.field$marginals.range.nominal[[1]])
# # 
# # f1 <- days_fever ~ 1 + f(v, model = 'iid') 
# # f2 <- days_fever ~ 1 + f(v, model = spde)
# # 
# # INLAutils::inlasloo(st_drop_geometry(dat.fit),
# #                     long = "longitude",
# #                     lat = "latitude",
# #                     y = "days_fever",
# #                     family = "poisson",
# #                     ss = 10,
# #                     rad = range.spatial*0.15,
# #                     modform = list(f1,f2),
# #                     mesh = mesh,
# #                     mae = TRUE,
# #                     ds = TRUE,
# #                     print = TRUE,
# #                     plot = TRUE)
# 
# 
