summarise_pred <- function(index.pred){
  
  samples_pred_df <- lapply(samples, make_pred_df, index.fit = index.t, index.pred = index.pred, dat = dat) 
  
  samples_pred_agg <- lapply(samples_pred_df, agg_pred_df, by = "total") 
  samples_pred_agg_endm <- lapply(samples_pred_df, agg_pred_df, by = "endemic")
  
  samples_pred_agg %>% 
    bind_rows() %>% 
    summarise(block_endm_2017 = "Total",
              n_cases = median(n_cases),
              n_acd = median(n_acd),
              p_acd = round(median(p_acd),1),
              fitted = round(median(fitted_delay_tot)),
              fitted.ci = paste0(round(quantile(fitted_delay_tot, c(0.01, 0.99))), collapse = ", "),
              fitted.pc = round(median(fitted_delay_tot/n_cases),1),
              fitted.pc.ci = paste0(round(quantile(fitted_delay_tot/n_cases, c(0.01, 0.99)),1), collapse = ", "),
              pred = round(median(pred_delay_tot)),
              pred.ci = paste0(round(quantile(pred_delay_tot, c(0.01, 0.99))), collapse = ", "),
              pred.pc = round(median(pred_delay_tot/n_cases),1),
              pred.pc.ci = paste0(round(quantile(pred_delay_tot/n_cases, c(0.01, 0.99)),1), collapse = ", "),
              chg = round(median(chg_tot)),
              chg.ci = paste0(round(quantile(chg_tot, c(0.01, 0.99))), collapse = ", "),
              chg.pc = round(median(chg_percase),1),
              chg.pc.ci = paste0(round(quantile(chg_percase, c(0.01, 0.99)),1), collapse = ", "),
              chg.pc.PCD = round(median(chg_perPCDcase),1),
              chg.pc.PCD.ci = paste0(round(quantile(chg_perPCDcase, c(0.01, 0.99)),1), collapse = ", "),
              chg.pc.ACD = round(median(chg_perACDcase),1),
              chg.pc.ACD.ci = paste0(round(quantile(chg_perACDcase, c(0.01, 0.99)),1), collapse = ", ")) -> summary_pred_total
  
  samples_pred_agg_endm %>% 
    bind_rows() %>% 
    group_by(block_endm_2017) %>% 
    dplyr::summarise(n_cases = median(n_cases),
                     n_acd = median(n_acd),
                     p_acd = round(median(p_acd),1),
                     fitted = round(median(fitted_delay_tot)),
                     fitted.ci = paste0(round(quantile(fitted_delay_tot, c(0.01, 0.99))), collapse = ", "),
                     fitted.pc = round(median(fitted_delay_tot/n_cases),1),
                     fitted.pc.ci = paste0(round(quantile(fitted_delay_tot/n_cases, c(0.01, 0.99)),1), collapse = ", "),
                     pred = round(median(pred_delay_tot)),
                     pred.ci = paste0(round(quantile(pred_delay_tot, c(0.01, 0.99))), collapse = ", "),
                     pred.pc = round(median(pred_delay_tot/n_cases),1),
                     pred.pc.ci = paste0(round(quantile(pred_delay_tot/n_cases, c(0.01, 0.99)),1), collapse = ", "),
                     chg = round(median(chg_tot)),
                     chg.ci = paste0(round(quantile(chg_tot, c(0.01, 0.99))), collapse = ", "),
                     chg.pc = round(median(chg_percase),1),
                     chg.pc.ci = paste0(round(quantile(chg_percase, c(0.01, 0.99)),1), collapse = ", "),
                     chg.pc.PCD = round(median(chg_perPCDcase),1),
                     chg.pc.PCD.ci = paste0(round(quantile(chg_perPCDcase, c(0.01, 0.99)),1), collapse = ", "),
                     chg.pc.ACD = round(median(chg_perACDcase),1),
                     chg.pc.ACD.ci = paste0(round(quantile(chg_perACDcase, c(0.01, 0.99)),1), collapse = ", ")) -> summary_pred_endm
  
  tab <- make_tab(summary_pred_total, summary_pred_endm)
  
  out <- list(tab = tab, total = summary_pred_total, endm = summary_pred_endm)
  return(out)
  
}