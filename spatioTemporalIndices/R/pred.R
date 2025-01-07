
#' constructIntPoints is used to set up the spatial locations for the integration points.
#' @param conf Configurations for catch-at-length
#' @param confPred Prediction configurations
#' @return A list with integration points. Strata number for each integration point is included.
#' @export
constructIntPoints<-function(conf,confPred){
  utmCRS = paste0("+proj=utm +zone=", conf$zone," +datum=WGS84 +units=km +no_defs")
  strata_utm <- sf::st_transform(conf$strata,utmCRS)
  points = sf::st_make_grid(strata_utm,cellsize=c(confPred$cellsize,confPred$cellsize),what="centers")
  points = sf::st_as_sf(points)

  #Define data frame with integration points to be returned
  points = sf::st_join(points,sf::st_buffer(strata_utm,1),left=FALSE)
  locUTM = data.frame(sf::st_coordinates(points)) #To be returned
  colnames(locUTM) = c("UTMX", "UTMY")
  idxStrata = as.numeric(points$id)

  #Remove possible duplicates
  noDuplicats =!duplicated(locUTM)
  locUTM = locUTM[noDuplicats,]
  idxStrata = idxStrata[noDuplicats]

  return(list(locUTM = locUTM, idxStrata = idxStrata))
}
