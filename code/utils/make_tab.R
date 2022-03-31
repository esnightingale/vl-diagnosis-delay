make_tab <- function(summary_pred_total, summary_pred_endm){
  
  full_join(summary_pred_total, summary_pred_endm) %>%
    column_to_rownames("block_endm_2017") %>% 
    mutate(across(ends_with(".ci"), function(x) paste0("[",x,"]"))) %>% 
    unite("fitted", fitted:fitted.ci, sep = " ") %>% 
    unite("fitted.pc", fitted.pc:fitted.pc.ci, sep = " ") %>% 
    unite("pred", pred:pred.ci, sep = " ") %>% 
    unite("pred.pc", pred.pc:pred.pc.ci, sep = " ") %>% 
    unite("chg", chg:chg.ci, sep = " ") %>% 
    unite("chg.pc", chg.pc:chg.pc.ci, sep = " ") %>% 
    unite("chg.pc.PCD", chg.pc.PCD:chg.pc.PCD.ci, sep = " ") %>% 
    unite("chg.pc.ACD", chg.pc.ACD:chg.pc.ACD.ci, sep = " ") %>% 
    sjmisc::rotate_df() -> tab
  
  return(tab)
  
}