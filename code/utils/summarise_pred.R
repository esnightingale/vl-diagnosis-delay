summarise_pred <- function(index.pred){
  
  samples_pred_df <- lapply(samples, make_pred_df, index.fit = index.t, index.pred = index.pred, dat = dat) 
  
  samples_pred_agg <- lapply(samples_pred_df, agg_pred_df, by = "total") 
  samples_pred_agg_endm <- lapply(samples_pred_df, agg_pred_df, by = "endemic")
  
  tab <- make_tab(samples_pred_agg, samples_pred_agg_endm)
  
  out <- list(tab = tab, total = summary_pred_total, endm = summary_pred_endm)
  return(out)
  
}