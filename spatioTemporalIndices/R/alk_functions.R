##' Set up data
##' @param df raw data
##' @param conf  configurations
##' @details Prepare data
##' @return Data to be provided to TMB
##' @export
setUpData_alk = function(dat_alk, conf_alk,conf_l = NULL){

  #Only use relevant years
  dat_alk$year = as.integer(format(dat_alk$startdatetime, format = "%Y"))
  dat_alk = dat_alk[which(dat_alk$year %in%conf_alk$years),]

  #Use age reading quality if false
  if(conf_alk$readability==0){#Do not utilize age reading quality
    dat_alk$readability[dat_alk$readability==5 | d$readability==6] = 1
  }

  #Truncate age
  dat_alk$ageNotTruncated = dat_alk$age
  dat_alk$age[dat_alk$age>conf_alk$maxAge] = conf_alk$maxAge
  if(length(dat_alk$age<conf_alk$minAge)>0){
    dat_alk$age[dat_alk$age<conf_alk$minAge] = conf_alk$minAge-1
  }
  dat_alk$age = dat_alk$age - min(dat_alk$age)+1 #First age is set to 1
  ageRange = c(min(dat_alk$age), max(dat_alk$age))
  nAge = ageRange[2]-ageRange[1] + 1

  #Bookkeeping, observations within years
  uniqueYears = unique(dat_alk$year)
  idx1 =idx2= rep(0,length(uniqueYears))
  for(y in 1:length(uniqueYears)){
    idx1[y] = min(which(dat_alk$year==uniqueYears[y]))-1
    idx2[y] = max(which(dat_alk$year==uniqueYears[y]))-1
  }

  #Convert to UTM coordinates
  loc = st_as_sf(dat_alk,coords = c("longitude","latitude"),crs="+proj=longlat")
  locUTM = st_coordinates(st_transform(loc,crs=paste0("+proj=utm +zone=", conf_l$zone," +datum=WGS84 +units=km +no_defs")))
  dat_alk$UTMX = locUTM[,1]
  dat_alk$UTMY = locUTM[,2]

  #Set up mesh and spde
  if(conf_alk$meshSimilar){
    mesh = spatioTemporalIndices::createMesh(conf_l)$mesh
  }else{
    stop("Not implemented functionality for different mesh for ALK and catch at length; Must use the same mesh")
  }
  spde <- fmesher::fm_fem(mesh)
  spdeMatricesST = list("M0" = spde$c0, "M1" = spde$g1, "M2" = spde$g2)
  A_list =list()
  for(i in 1:length(uniqueYears)){
    A_list[[i]] = fmesher::fm_basis(mesh, loc = as.matrix(locUTM[which(dat_alk$year == uniqueYears[i]),]))
  }

  #Set up data list
  data = list(age = dat_alk$age,
              ageNotTruncated = dat_alk$ageNotTruncated,
              length = dat_alk$length,
              spdeMatricesST_alk = spdeMatricesST,
              A_alk_list = A_list,
              readability = dat_alk$readability,
              idx1 = idx1,
              idx2 = idx2,
              ageRange = ageRange,
              rwBeta0_alk = conf_alk$rwBeta0,
              maxAge = conf_alk$maxAge,
              minAge = conf_alk$minAge,
              usePCpriorsALK = conf_alk$usePCpriorsALK,
              pcPriorsALKRange = conf_alk$pcPriorsALKRange,
              pcPriorsALKSD = conf_alk$pcPriorsALKSD,
              spatioTemporalALK = conf_alk$spatioTemporal,
              spatialALK = conf_alk$spatial,
              betaLength = conf_alk$betaLength)

  attributes(data)$uniqueYears = uniqueYears
  attributes(data)$loc = loc
  attributes(data)$locUTM = locUTM
  attributes(data)$years = dat_alk$year
  attributes(data)$mesh = mesh

  return(data)
}


##' Set up conf
##' @param years years to include
##' @param maxAge max age
##' @param spatioTemporal include spatio-tempora: 0-no 1-yes 2- yes but without correlation structure in time
##' @param cutoff mesh spesific, min distance between points
##' @param cbound mesh spesific, size of boundary meshSimilar
##' @param meshSimilar If true we apply the same mesh as in the package spatitemporalIndices
##' @param readability If 1: Utilize age reading quality, if 0: Do not utilize age reading quality
##' @details defines the configurations
##' @return Configurations to set up the model
##' @export
defConf_alk = function(years= NULL,minAge,maxAge, spatioTemporal = 0,spatial = 0,rwBeta0 = 1, betaLength = 1,
                       meshSimilar = TRUE,
                       readability = 1,
                       usePCpriorsALK = 0, pcPriorsALKRange = c(300,0.1), pcPriorsALKSD = c(1,0.1)){
  conf = list()
  conf$minAge = minAge
  conf$maxAge = maxAge
  conf$rwBeta0 = rwBeta0
  conf$meshSimilar = meshSimilar
  conf$years = years
  conf$readability = readability
  conf$usePCpriorsALK = usePCpriorsALK
  conf$pcPriorsALKRange = pcPriorsALKRange
  conf$pcPriorsALKSD = pcPriorsALKSD
  conf$spatioTemporal = spatioTemporal
  conf$spatial = spatial
  conf$betaLength= betaLength

  return(conf)
}



##' Set up map variable
##' @param conf Configurations
##' @param par Parameters
##' @details Sets up the map-variable
##' @return Returns the map-variable
##' @export
setMap_alk = function(conf,par){
  map = list()
  map$logKappa_alk = c(1,2)
  map$logSigma_alk = c(1,2)

  if(conf$spatial==0){
    map$xS_alk = as.factor(rep(NA,length(par$xS_alk)))
    map$logKappa_alk[1] = NA
    map$logSigma_alk[1] = NA
  }

  if(conf$spatioTemporal==2){#No use of ar1 structure in time
    map$transRho_alk = as.factor(NA)
  }else if(conf$spatioTemporal==0){#Turn off spatio-temporal structures
    map$transRho_alk = as.factor(NA)
    map$xST_alk = as.factor(rep(NA,length(par$xST_alk)))
    map$logKappa_alk[2] = NA
    map$logSigma_alk[2] = NA
  }

  map$logKappa_alk = as.factor(map$logKappa_alk)
  map$logSigma_alk = as.factor(map$logSigma_alk)

  if(conf$rwBeta0==0){
    map$log_sigma_beta_alk = as.factor(c(NA,NA))
  }else{
    if(conf$betaLength==1){
      map$log_sigma_beta_alk = as.factor(c(1,NA))
    }
  }
  return(map)
}

##' Set up data
##' @param data data
##' @param conf  configurations
##' @details define parameters
##' @return Prameters to be provided to TMB
##' @export
defpar_alk = function(data, conf){

  nAge = data$ageRange[2]-data$ageRange[1] + 1

  xS_alk = array(0.0, dim = c(dim(data$A_alk_list[[1]])[2],nAge-1))
  xST_alk = array(0.0, dim = c(dim(data$A_alk_list[[1]])[2],length(conf$years),nAge-1))

  par = list(beta0_alk = array(0,dim=c(length(conf$years),nAge-1)),
             log_sigma_beta_alk = c(0,0),
             betaLength_alk = rep(0,nAge-1),
             logSigma_alk = c(2,2),
             logKappa_alk = c(-4,-4),
             transRho_alk = 0,
             xS_alk = xS_alk,
             xST_alk = xST_alk)

  if(conf$betaLength == 2){
    par$betaLength_alk = rep(0,(nAge-1)*length(conf$years))
  }

  return(par)
}

