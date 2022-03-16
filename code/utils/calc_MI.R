calc_MI <- function(values, coo, range){

  S.dist  <-  spdep::dnearneigh(coo, 0, range) 
  lw <- spdep::nb2listw(S.dist, style="W",zero.policy=T) 
  
  MI  <-  moran.mc(values, lw, nsim=599, zero.policy=T) 
  MI
  
  return(plot(MI, main="", las=1))
  
}
