
#' constructIntPoints
#' @return
#' @export
#' @examples
#' @return
constructIntPoints<-function(conf,confPred){
  utmCRS = CRS(paste0("+proj=utm +zone=", conf$zone," +datum=WGS84 +units=km +no_defs"))
  strata_utm <- spTransform(conf$strata,utmCRS)
  points = sp::makegrid(strata_utm,cellsize = confPred$cellsize)
  pointsSP = sp::SpatialPoints(points,utmCRS)

  pointsSPXY = suppressWarnings(spTransform(pointsSP,conf$strata@proj4string))
  inside = rep(1,dim(pointsSPXY@coords)[1])
  inside[which(is.na(over(pointsSPXY,conf$strata)))] = 0 #over() with lat-lon, some points between strata were not assigned to a strata when using over() with utm
  points = points[which(inside==1),]
  pointsSP = sp::SpatialPoints(points,utmCRS)
  pointsSPXY =  suppressWarnings(spTransform(pointsSP,conf$strata@proj4string))

  #Define data frame with integration points to be returned
  locUTM = data.frame(points[,1],points[,2]) #To be returned
  colnames(locUTM) = c("UTMX", "UTMY")

  idxStrata = rep(-1,dim(locUTM)[1])
  for(i in 1:nrow(conf$strata)){
    insideThis =  which(!is.na(over(pointsSPXY,conf$strata[i,])[,1]))
    idxStrata[insideThis] = i
  }

  return(list(locUTM = locUTM, idxStrata = idxStrata))
}
