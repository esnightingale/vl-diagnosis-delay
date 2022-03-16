
agg_pred_df <- function(pred.df, by){
  
  if (by == "total"){
    pred.df %>%   
      dplyr::summarise(n_cases = n(),
                       n_acd = sum(detection == "ACD"),
                       p_acd = mean(detection == "ACD")*100,
                       fitted_delay_tot = sum(fitted),
                       pred_delay_tot = sum(pred),
                       days_saved_tot = fitted_delay_tot - pred_delay_tot,
                       days_saved_perc = (fitted_delay_tot - pred_delay_tot)*100/fitted_delay_tot,
                       days_saved_percase = (fitted_delay_tot - pred_delay_tot)/n_cases,
                       days_saved_perPCDcase = na_if(na_if((fitted_delay_tot - pred_delay_tot)/(n_cases - n_acd), Inf), -Inf),
                       days_saved_perACDcase = na_if(na_if((fitted_delay_tot - pred_delay_tot)/n_acd, Inf), -Inf)) -> agg
    
  }else if (by == "endemic"){
    pred.df %>%  
      group_by(block_endm_2017) %>% 
      dplyr::summarise(n_cases = n(),
                       n_acd = sum(detection == "ACD"),
                       p_acd = mean(detection == "ACD")*100,
                       fitted_delay_tot = sum(fitted),
                       pred_delay_tot = sum(pred),
                       days_saved_tot = fitted_delay_tot - pred_delay_tot,
                       days_saved_perc = (fitted_delay_tot - pred_delay_tot)*100/fitted_delay_tot,
                       days_saved_percase = (fitted_delay_tot - pred_delay_tot)/n_cases,
                       days_saved_perPCDcase = na_if(na_if((fitted_delay_tot - pred_delay_tot)/(n_cases - n_acd), Inf), -Inf),
                       days_saved_perACDcase = na_if(na_if((fitted_delay_tot - pred_delay_tot)/n_acd, Inf), -Inf)) -> agg
    
  } 
  
  return(agg)
}
