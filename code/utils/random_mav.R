# Summarise fitted SPDE and IID effects over all evaluated points (mean absolute values)

random_mav <- function(fit, name){
  
  # Extract SPDE value at each mesh node, and IID value for each individual
  out.spde <- INLA::inla.spde2.result(fit,'s', spde, do.transf=TRUE)
  
  # rang <- apply(mesh$loc[, c(1, 2)], 2, range)
  # proj <- inla.mesh.projector(mesh, 
  #                             xlim = rang[, 1], 
  #                             ylim = rang[, 2], 
  #                             dims = c(300, 300))
  # out.spde.proj <- inla.mesh.project(proj, fit$summary.random$s$mean)
  
  out.iid <- fit$summary.random$id
  
  # Calculate mean absolute value of SPDE/IID effects across all mesh nodes/villages
  tab <- data.frame(Model = name)
  
  if(!is.null(out.iid$mean)){
    tab$mav_iid <- mean(abs(out.iid$mean))
    tab$msv_iid <- mean((out.iid$mean)^2)
  }
  
  if(!is.null(out.spde$summary.values)){
    tab$mav_spde <- mean(abs(out.spde$summary.values$mean))
    tab$msv_spde <- mean((out.spde$summary.values$mean)^2)
   }
  
  # if(!is.null(out.spde.proj)){
  #   tab$mav_spde_proj <- mean(abs(out.spde.proj$summary.values$mean))
  #   tab$msv_spde_proj <- mean((out.spde.proj$summary.values$mean)^2)
  # }

  return(tab)
  
}
