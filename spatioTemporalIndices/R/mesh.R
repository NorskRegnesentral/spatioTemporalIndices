#' Create mesh
#'
#' This function creates the mesh to represent the Matern covariance strucutre
#'
#' @param  conf configurations in the model.
#' @details The resolution of the mesh depends on the cutoff configuration. The mesh is constructed based on integration points with high resolution. How much outside of the area is included in the mesh depends on the cbound configuration.
#' @return Returns the mesh used with the SPDE-procedure
#'
#' @export
createMesh <- function(conf){

  cutoff = conf$cutoff
  if(conf$maxEdge[1]<conf$cutoff){
    cutoff = conf$maxEdge[1] #Do not allow maxEdge shorter than cutoff. If maxEdge is provided, no neighboring nodes are distanced further apart. The option to only included cutoff is included for historical reasons.
  }

  utmCRS = paste0("+proj=utm +zone=", conf$zone, " +datum=WGS84 +units=km +no_defs")
  strata_utm <- sf::st_transform(conf$strata, utmCRS)

  buffer = sf::st_buffer(strata_utm$geometry,conf$cbound[1])
  buffer2 = sf::st_buffer(strata_utm$geometry,conf$cbound[2])
  combined <- sf::st_union(buffer)
  combined2 <- sf::st_union(buffer2)

  mesh <- fmesher::fm_mesh_2d(max.edge =conf$maxEdge,
                              boundary = list(combined,combined2),
                              cutoff = conf$cutoff)

  return(list(mesh=mesh, barrier.triangles =NULL))
}



