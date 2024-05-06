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
  maxEdge = c(conf$cutoff,conf$cutoff*4) # Longer distances between nodes outside if inner bounderary

  confPredTmp = list(cellsize = 1000)
  intPoints = constructIntPoints(conf, confPredTmp)$locUTM
  while(dim(intPoints)[1]<5000){#Set up mesh based on a fine grid of integration points
    confPredTmp$cellsize = confPredTmp$cellsize/2
    intPoints = constructIntPoints(conf,confPredTmp)$locUTM
  }
  boundary.loc <- SpatialPoints(as.matrix(intPoints))
  boundary <- list(
    inla.nonconvex.hull(coordinates(boundary.loc), convex  = conf$cbound[1],resolution = 120),
    inla.nonconvex.hull(coordinates(boundary.loc), convex  = conf$cbound[2]))
  mesh <- inla.mesh.2d(boundary=boundary,
                       max.edge=maxEdge,
                       cutoff=conf$cutoff)

  print(paste("Mesh points:",mesh$n))
  plot(mesh)
  return(list(mesh=mesh, barrier.triangles =NULL))
}



