
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
    right_join(readRDS("E:/vl-diagnosis-delay/data/decode.rds")) %>%
    right_join(readRDS(file.path(datadir, pat.path))) %>% 
    sf::st_set_crs(7759)
  
  return(dat)
  
}