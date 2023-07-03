
#' constructIntPoints
#' @return
#' @export
#' @examples
#' @return
constructIntPoints<-function(conf,confPred){
  utmCRS = CRS(paste0("+proj=utm +zone=", conf$zone," +datum=WGS84 +units=km +no_defs"))
  strata_utm <- spTransform(conf$strata,utmCRS)

  #Construct a grid of integration points
  points = sp::makegrid(strata_utm,confPred$nIntPoints)
  pointsSP = sp::SpatialPoints(points,utmCRS)
  inside = rep(1,dim(pointsSP@coords)[1])
  inside[which(is.na(over(pointsSP,strata_utm)))] = 0
  points = points[which(inside==1),]
  pointsSP = sp::SpatialPoints(points,utmCRS)


  #Define data frame with integration points to be returned
  locUTM = data.frame(points[,1],points[,2]) #To be returned
  colnames(locUTM) = c("UTMX", "UTMY")

  idxStrata = rep(-1,dim(locUTM)[1])
  for(i in 1:nrow(conf$strata)){
    insideThis =  which(!is.na(over(pointsSP,strata_utm[i,])[,1]))
    idxStrata[insideThis] = i
  }

  return(list(locUTM = locUTM, idxStrata = idxStrata))
}


