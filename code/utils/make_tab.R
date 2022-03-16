make_tab <- function(samples_pred_agg, samples_pred_agg_endm){
  
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
              saved = round(median(days_saved_tot)),
              saved.ci = paste0(round(quantile(days_saved_tot, c(0.01, 0.99))), collapse = ", "),
              saved.pc = round(median(days_saved_percase),1),
              saved.pc.ci = paste0(round(quantile(days_saved_percase, c(0.01, 0.99)),1), collapse = ", "),
              saved.pc.PCD = round(median(days_saved_perPCDcase),1),
              saved.pc.PCD.ci = paste0(round(quantile(days_saved_perPCDcase, c(0.01, 0.99)),1), collapse = ", "),
              saved.pc.ACD = round(median(days_saved_perACDcase),1),
              saved.pc.ACD.ci = paste0(round(quantile(days_saved_perACDcase, c(0.01, 0.99)),1), collapse = ", ")) -> summary_pred_total
  
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
                     saved = round(median(days_saved_tot)),
                     saved.ci = paste0(round(quantile(days_saved_tot, c(0.01, 0.99))), collapse = ", "),
                     saved.pc = round(median(days_saved_percase),1),
                     saved.pc.ci = paste0(round(quantile(days_saved_percase, c(0.01, 0.99)),1), collapse = ", "),
                     saved.pc.PCD = round(median(days_saved_perPCDcase),1),
                     saved.pc.PCD.ci = paste0(round(quantile(days_saved_perPCDcase, c(0.01, 0.99)),1), collapse = ", "),
                     saved.pc.ACD = round(median(days_saved_perACDcase),1),
                     saved.pc.ACD.ci = paste0(round(quantile(days_saved_perACDcase, c(0.01, 0.99)),1), collapse = ", ")) -> summary_pred_endm
  
  
  full_join(summary_pred_total, summary_pred_endm) %>%
    column_to_rownames("block_endm_2017") %>% 
    mutate(across(ends_with(".ci"), function(x) paste0("[",x,"]"))) %>% 
    unite("fitted", fitted:fitted.ci, sep = " ") %>% 
    unite("fitted.pc", fitted.pc:fitted.pc.ci, sep = " ") %>% 
    unite("pred", pred:pred.ci, sep = " ") %>% 
    unite("pred.pc", pred.pc:pred.pc.ci, sep = " ") %>% 
    unite("saved", saved:saved.ci, sep = " ") %>% 
    unite("saved.pc", saved.pc:saved.pc.ci, sep = " ") %>% 
    unite("saved.pc.PCD", saved.pc.PCD:saved.pc.PCD.ci, sep = " ") %>% 
    unite("saved.pc.ACD", saved.pc.ACD:saved.pc.ACD.ci, sep = " ") %>% 
    sjmisc::rotate_df() -> tab
  
  return(tab)
  
}