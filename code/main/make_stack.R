
# Organise the data, projection matrices and fixed/random effects into a "stack".
# Need separate ones for estimation at observed locations and prediction across
# defined grid

# Fixed effects are intercept
# Random effect is spatial GRF

# List A contains 1 to indicate fixed effects are mapped one to one, and a projection
# matrix for the random effects

# Prediction stack has the response set as NA
dat.fit.df <- st_drop_geometry(dat.fit)
dat.val.df <- st_drop_geometry(dat.val)

# stack for estimation stk.e
stk.e <- inla.stack(
  tag = "est", 
  data = list(y = dat.fit$days_fever),
  A = list(A, 1, 1),
  effects = list(s = indexs, i = dat.fit.df$i, 
                 list(data.frame(b0 = 1, 
                                 # patient chars
                                 age_s = dat.fit.df$age_s,
                                 sex = dat.fit.df$sex,
                                 hiv = dat.fit.df$hiv,
                                 marg_caste = dat.fit.df$marg_caste,
                                 detection = dat.fit.df$detection,
                                 prv_tx_ka = dat.fit.df$prv_tx_ka,
                                 # village chars
                                 num_conslt_cat = dat.fit.df$num_conslt_cat,
                                 vill_inc_2017_gt0 = dat.fit.df$vill_inc_2017_gt0,
                                 IRS_2017_1 = dat.fit.df$IRS_2017_1,
                                 block_endm_2017 = dat.fit.df$block_endm_2017,
                                 traveltime_s = dat.fit.df$traveltime_s)))
)


# stack for validation at witheld test points stk.v
stk.v <- inla.stack(
  tag = "val",
  data = list(y = dat.val.df$days_fever),
  A = list(Av, 1, 1),
  effects = list(s = indexs, i = dat.val.df$i, 
                 list(data.frame(b0 = 1, 
                                 # patient chars
                                 age_s = dat.val.df$age_s,
                                 sex = dat.val.df$sex,
                                 hiv = dat.val.df$hiv,
                                 marg_caste = dat.val.df$marg_caste,
                                 detection = dat.val.df$detection,
                                 prv_tx_ka = dat.val.df$prv_tx_ka,
                                 # village chars
                                 num_conslt_cat = dat.val.df$num_conslt_cat,
                                 vill_inc_2017_gt0 = dat.val.df$vill_inc_2017_gt0,
                                 IRS_2017_1 = dat.val.df$IRS_2017_1,
                                 block_endm_2017 = dat.val.df$block_endm_2017,
                                 traveltime_s = dat.val.df$traveltime_s)))
)

# stack for prediction stk.p
stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, Ap),
  effects = list(data.frame(b0 = rep(1, nrow(coop))), s = indexs)
)

# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.v, stk.p)
