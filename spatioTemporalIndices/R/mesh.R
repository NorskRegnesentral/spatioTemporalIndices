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

  confPredTmp = list(cellsize = 1000)
  intPoints = constructIntPoints(conf, confPredTmp)$locUTM
  while(dim(intPoints)[1]<5000){#Set up mesh based on a fine grid of integration points.
    confPredTmp$cellsize = confPredTmp$cellsize/2
    intPoints = constructIntPoints(conf,confPredTmp)$locUTM
  }

  intPoints = round(intPoints,3) #The mesh can be slightly different across operating system when numbers are not rounded.

  splancs::splancs()#Splancs needed in fmesher::fm_nonconvex_hull_inla
  boundary <- list(
    fmesher::fm_nonconvex_hull_inla(as.matrix(intPoints), convex  = conf$cbound[1],resolution = 120),
    fmesher::fm_nonconvex_hull_inla(as.matrix(intPoints), convex  = conf$cbound[2]))
  mesh <- fmesher::fm_mesh_2d(boundary=boundary,
                              max.edge=conf$maxEdge,
                              cutoff=cutoff)



  return(list(mesh=mesh, barrier.triangles =NULL))
}



