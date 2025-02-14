##' setUpData_alk
##'
##' Set up age data in the format needed by TMB.
##'
##'
##' @param dat_alk age data, a data table with one line for each fish
##' @param conf_alk  configurations for age-at-length model
##' @param conf_l  configurations for length-at-age model
##' @details Prepare data list for the age part. This data-list will be merged with the length data-list later.
##' @return Data to be provided to TMB
##' @export
setUpData_alk = function(dat_alk, conf_alk,conf_l = NULL){


  #Verify that station ID's are unique
  if(min(tapply(dat_alk$startdatetime, dat_alk$station, function(x) length(unique(x))==1))==0){
    stop("Station ID's are not unique. dat_alk$station must be unique for each haul.")
  }

  #Only use relevant years
  dat_alk$year = as.integer(format(dat_alk$startdatetime, format = "%Y"))
  dat_alk = dat_alk[which(dat_alk$year %in%conf_alk$years),]

  #Use age reading quality if false
  if(conf_alk$readability==0){#Do not utilize age reading quality
    dat_alk$readability[dat_alk$readability==5 | dat_alk$readability==6] = 1
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
  loc = sf::st_as_sf(dat_alk,coords = c("longitude","latitude"),crs="+proj=longlat")
  locUTM = sf::st_coordinates(sf::st_transform(loc,crs=paste0("+proj=utm +zone=", conf_l$zone," +datum=WGS84 +units=km +no_defs")))
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


##' defConf_alk
##'
##' Set up the configurations for the age-at-length model
##'
##' @param years Years to include
##' @param minAge minimum age
##' @param maxAge max age
##' @param spatioTemporal include spatio-tempora: 0-no 1-yes 2- yes but without correlation structure in time
##' @param spatial include spatial effect: 0-no 1-yes
##' @param rwBeta0 random walk for intercelt within cohorts:  0-no 1-yes
##' @param betaLength Time varying regression coefficients for length: 1-no, 2-yes
##' @param meshSimilar If TRUE; apply the same mesh as in the package spatitemporalIndices, this is typically the case
##' @param readability If 1: Utilize age reading quality, if 0: Do not utilize age reading quality
##' @param usePCpriorsALK Use pc-priors: 0-no 1-yes. Probably never used; in case it is used, we need set pcPriorsALKRange and pcPriorsALKSD to reasonable values
##' @param pcPriorsALKRange pc-priors for range
##' @param pcPriorsALKSD pc-priors for standard deviations
##' @details Defines the configurations
##' @return Configurations to set up the age-at-length part of the model
##' @export
defConf_alk = function(years= NULL,minAge,maxAge, spatioTemporal = 0,spatial = 0,rwBeta0 = 0, betaLength = 1,
                       meshSimilar = TRUE,
                       readability = 1,
                       usePCpriorsALK = 0, pcPriorsALKRange = c(300,0.1), pcPriorsALKSD = c(1,0.1)){
  conf_alk = list()
  conf_alk$minAge = minAge
  conf_alk$maxAge = maxAge
  conf_alk$rwBeta0 = rwBeta0
  conf_alk$meshSimilar = meshSimilar
  conf_alk$years = years
  conf_alk$readability = readability
  conf_alk$usePCpriorsALK = usePCpriorsALK
  conf_alk$pcPriorsALKRange = pcPriorsALKRange
  conf_alk$pcPriorsALKSD = pcPriorsALKSD
  conf_alk$spatioTemporal = spatioTemporal
  conf_alk$spatial = spatial
  conf_alk$betaLength= betaLength

  return(conf_alk)
}



##' setMap_alk
##'
##' Set up map-variable based on the configurations for the age-at-length model.
##'
##' @param conf_alk Configurations
##' @param par_alk Parameters
##' @details Sets up the map-variable
##' @return Returns the map-variable
##' @export
setMap_alk = function(conf_alk,par_alk){
  map = list()
  map$logKappa_alk = c(1,2)
  map$logSigma_alk = c(1,2)

  if(conf_alk$spatial==0){
    map$xS_alk = as.factor(rep(NA,length(par_alk$xS_alk)))
    map$logKappa_alk[1] = NA
    map$logSigma_alk[1] = NA
  }

  if(conf_alk$spatioTemporal==2){#No use of ar1 structure in time
    map$transRho_alk = as.factor(NA)
  }else if(conf_alk$spatioTemporal==0){#Turn off spatio-temporal structures
    map$transRho_alk = as.factor(NA)
    map$xST_alk = as.factor(rep(NA,length(par_alk$xST_alk)))
    map$logKappa_alk[2] = NA
    map$logSigma_alk[2] = NA
  }

  map$logKappa_alk = as.factor(map$logKappa_alk)
  map$logSigma_alk = as.factor(map$logSigma_alk)

  if(conf_alk$rwBeta0==0){
    map$log_sigma_beta_alk = as.factor(c(NA,NA))
  }else{
    if(conf_alk$betaLength==1){
      map$log_sigma_beta_alk = as.factor(c(1,NA))
    }
  }
  return(map)
}

##' defpar_alk
##'
##' Set up the age-at-length model parameter list used by TMB.
##'
##' @param data_alk data-list for the age-at-length model
##' @param conf_alk  configurations for the age-at-length model
##' @details Set up the age-at-length model parameter list used by TMB.
##' @return The age-at-length model parameter list with initial values.
##' @export
defpar_alk = function(data_alk, conf_alk){

  nAge = data_alk$ageRange[2]-data_alk$ageRange[1] + 1

  xS_alk = array(0.0, dim = c(dim(data_alk$A_alk_list[[1]])[2],nAge-1))
  xST_alk = array(0.0, dim = c(dim(data_alk$A_alk_list[[1]])[2],length(conf_alk$years),nAge-1))

  par = list(beta0_alk = array(0,dim=c(length(conf_alk$years),nAge-1)),
             log_sigma_beta_alk = c(0,0),
             betaLength_alk = rep(0,nAge-1),
             logSigma_alk = c(2,2),
             logKappa_alk = c(-4,-4),
             transRho_alk = 0,
             xS_alk = xS_alk,
             xST_alk = xST_alk)

  if(conf_alk$betaLength == 2){
    par$betaLength_alk = rep(0,(nAge-1)*length(conf_alk$years))
  }

  return(par)
}

