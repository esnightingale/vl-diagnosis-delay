# Summarise fitted SPDE and IID effects over all evaluated points (mean absolute values)

random_mav <- function(fit, name){
  
  # Extract SPDE value at each mesh node, and IID value for each village
  out.spde <- INLA::inla.spde2.result(fit,'s', spde, do.transf=TRUE)
  out.iid <- fit$summary.random$v
  
  # Calculate mean absolute value of SPDE/IID effects across all mesh nodes/villages
  tab <- data.frame(Model = name,
                    mav_iid = mean(abs(out.iid$mean)))
  
  if(!is.null(out.spde$summary.values)){
    tab$mav_spde <- mean(abs(out.spde$summary.values$mean))
  }
  
  return(tab)
  
}
