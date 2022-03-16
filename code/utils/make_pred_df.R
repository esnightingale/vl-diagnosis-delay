
make_pred_df <- function(sample, index.fit, index.pred, dat){
  
  pred.df <- data.frame(delay = dat$delay, 
                        detection = factor(dat$poss_acd, levels = c(FALSE,TRUE), labels = c("PCD","ACD")),
                        block_endm_2017 = factor(dat$block_endm_2017, levels = c(FALSE,TRUE), labels = c("Non-endemic","Endemic")),
                        fitted = exp(sample$latent[index.fit]),
                        pred = exp(sample$latent[index.pred])) 
  return(pred.df)
  
}