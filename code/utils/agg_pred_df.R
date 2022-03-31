
agg_pred_df <- function(pred.df, by){
  
  if (by == "total"){
    pred.df %>%   
      dplyr::summarise(n_cases = n(),
                       n_acd = sum(detection == "ACD"),
                       p_acd = mean(detection == "ACD")*100,
                       fitted_delay_tot = sum(fitted),
                       pred_delay_tot = sum(pred),
                       chg_tot = pred_delay_tot - fitted_delay_tot,
                       chg_perc = chg_tot*100/fitted_delay_tot,
                       chg_percase = chg_tot/n_cases,
                       chg_perPCDcase = na_if(na_if(chg_tot/(n_cases - n_acd), Inf), -Inf),
                       chg_perACDcase = na_if(na_if(chg_tot/n_acd, Inf), -Inf)) -> agg
    
  }else if (by == "endemic"){
    pred.df %>%  
      group_by(block_endm_2017) %>% 
      dplyr::summarise(n_cases = n(),
                       n_acd = sum(detection == "ACD"),
                       p_acd = mean(detection == "ACD")*100,
                       fitted_delay_tot = sum(fitted),
                       pred_delay_tot = sum(pred),
                       chg_tot = pred_delay_tot - fitted_delay_tot,
                       chg_perc = chg_tot*100/fitted_delay_tot,
                       chg_percase = chg_tot/n_cases,
                       chg_perPCDcase = na_if(na_if(chg_tot/(n_cases - n_acd), Inf), -Inf),
                       chg_perACDcase = na_if(na_if(chg_tot/n_acd, Inf), -Inf)) -> agg
    
  } 
  
  return(agg)
}
