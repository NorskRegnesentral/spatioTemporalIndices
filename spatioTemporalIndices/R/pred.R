
#' constructIntPoints
#' @return
#' @export
#' @examples
#' @return
constructIntPoints<-function(conf,confPred){
  utmCRS = paste0("+proj=utm +zone=", conf$zone," +datum=WGS84 +units=km +no_defs")
  strata_utm <- st_transform(conf$strata,utmCRS)
  points = st_sample(strata_utm,size = round(sum(as.numeric(st_area(conf$strata)))*(1/1.852)^2/1e6/confPred$cellsize,0),type="regular")

  # cross check whether points are inside strata polygons
  pointsSP = suppressWarnings(st_transform(points,st_crs(conf$strata)))
  pointsSP = pointsSP[!st_contains(pointsSP,conf$strata,sparse=FALSE)[,1]]
  points = points[!st_contains(pointsSP,conf$strata,sparse=FALSE)[,1]]
  
  #Define data frame with integration points to be returned
  pointsSP = st_join(st_as_sf(pointsSP),conf$strata)
  
  locUTM = data.frame(st_coordinates(points)) #To be returned
  colnames(locUTM) = c("UTMX", "UTMY")

  idxStrata = pointsSP$id
  
  return(list(locUTM = locUTM, idxStrata = idxStrata))
}
