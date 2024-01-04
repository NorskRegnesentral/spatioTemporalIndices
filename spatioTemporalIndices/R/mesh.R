#' Create mesh
#'
#' This function creates the mesh to represent the Matern covariance strucutre
#'
#' @param  conf
#' @return Mesh used
#'
#' @export
createMesh <- function(conf){
  maxEdge = c(conf$cutoff,conf$cutoff*4) # Longer distances between nodes outside if inner bounderary

  confPredTmp = list(cellsize=21)
  intPoints = constructIntPoints(conf, confPredTmp)$locUTM #Mesh is based on integration points
  boundary.loc <- SpatialPoints(as.matrix(intPoints))
#  interior = fmesher::fm_nonconvex_hull(coordinates(boundary.loc), 25)
#  mesh = fmesher::fm_rcdt_2d_inla(coordinates(boundary.loc), refine = list(),interior = interior,  extend = list(offset = conf$offset),
#                                  max.edge=maxEdge,
#                                  cutoff=conf$cutoff)

  boundary <- list(
    inla.nonconvex.hull(coordinates(boundary.loc), 20,resolution = 100),
    inla.nonconvex.hull(coordinates(boundary.loc), conf$cbound))
  mesh <- inla.mesh.2d(boundary=boundary,
                       max.edge=maxEdge,
                       min.angle=c(30, 15),
                       cutoff=conf$cutoff)

  print(paste("Mesh points:",mesh$n))
  return(list(mesh=mesh, barrier.triangles =NULL))
}
