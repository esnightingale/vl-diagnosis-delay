
read_data <- function(analysis = TRUE){
  
  # Local data folder
  datadir <- "C:/Users/phpuenig/Documents/VL/Data/KAMIS/Clean/ll"
  
  if (analysis == TRUE){
    pat.path <- "analysisdata_pat.rds"
    vill.path <- "analysisdata_vill.rds"
  }else{
    pat.path <- "pat.rds"
    vill.path <- "vill.rds"
  }
  
  dat <- readRDS(file.path(datadir,vill.path)) %>%
    # Matching code from external drive
    full_join(readRDS("E:/vl-diagnosis-delay/data/code.rds")) %>%
    right_join(readRDS(file.path(datadir, pat.path))) 
  
  return(dat)
  
}