
#' constructIntPoints
#' @return
#' @export
#' @examples
#' @return
constructIntPoints<-function(conf,confPred){
  utmCRS = paste0("+proj=utm +zone=", conf$zone," +datum=WGS84 +units=km +no_defs")
  strata_utm <- st_transform(conf$strata,utmCRS)
  cellength = sqrt(confPred$cellsize)
  points = st_make_grid(strata_utm,cellsize=c(cellength,cellength),what="centers")
  #points = st_sample(strata_utm,size = round(sum(as.numeric(st_area(strata_utm)))*(1/1.852)^2/confPred$cellsize,0),type="regular")
  points = st_as_sf(points)

  #Define data frame with integration points to be returned
  points = st_join(points,strata_utm,left=FALSE)
  
  locUTM = data.frame(st_coordinates(points)) #To be returned
  colnames(locUTM) = c("UTMX", "UTMY")

  # idxStrata = rep(-1,dim(locUTM)[1])
  # for(i in 1:nrow(conf$strata)){
  #     insideThis =  which(!is.na(over(pointsSPXY,test[i,])[,1]))
  #     idxStrata[insideThis] = i
  # }
  # # 
  idxStrata = as.numeric(points$id)
  
  return(list(locUTM = locUTM, idxStrata = idxStrata))
}
