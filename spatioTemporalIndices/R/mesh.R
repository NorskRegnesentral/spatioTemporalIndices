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

  splancs::splancs()#Splancs needed in fmesher::fm_nonconvex_hull_inla
  boundary <- list(
    fmesher::fm_nonconvex_hull_inla(as.matrix(intPoints), convex  = conf$cbound[1],resolution = 120),
    fmesher::fm_nonconvex_hull_inla(as.matrix(intPoints), convex  = conf$cbound[2]))
  mesh <- fmesher::fm_mesh_2d(boundary=boundary,
                              max.edge=maxEdge,
                              cutoff=conf$cutoff)



  return(list(mesh=mesh, barrier.triangles =NULL))
}



