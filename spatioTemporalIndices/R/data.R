

#' setupData
#'
#' Set up data used by model
#'
#' @description
#' @param data A list of the data object returned by getBaseline, where each element of the list is a year
#' @param conf Configurations, see \code{\link{defConf}} for details.
#' @param confPred Prediction configurations
#' @return
#' @export
#' @examples
setupData = function(dat_l,conf_l,confPred){

  #Remove too short lengths
  dat_l = dat_l[dat_l$lengthGroup>=conf_l$minLength,]

  #Remove too long lengths
  for(id in unique(dat_l$station)){
    index =which(dat_l$station==id & dat_l$lengthGroup>=conf_l$maxLength)
    if(conf_l$plusGroup==1){
      dat_l$catch[index[1]] = sum(dat_l$catch[index])
    }
  }
  if(length(which(dat_l$lengthGroup> conf_l$maxLength)>0)){
    dat_l = dat_l[-which(dat_l$lengthGroup> conf_l$maxLength),]
  }

  # Add strata polygon if none was provided
  if(is.null(conf_l$strata)) {
     hullpol = st_convex_hull(st_union(st_as_sf(dat_l,coords=c("longitude","latitude"),crs="+proj=longlat")))
     conf_l$strata = st_as_sf(hullpol)
     conf_l$strata$id = 1
     conf_l$zone = floor((mean(st_bbox(conf_l$strata)[c(1,3)]) + 180) / 6) + 1
     conf_l$stratasystem = "data"
     conf_l$strata_number = nrow(conf_l$strata)
     print("Strata polygon created from survey data.")
  }

  #Convert to UTM coordinates
  #loc = data.frame(dat_l$longitude,dat_l$latitude)
  loc = st_as_sf(dat_l,coords = c("longitude","latitude"),crs="+proj=longlat")
  locUTM = st_coordinates(st_transform(loc,crs=paste0("+proj=utm +zone=", conf_l$zone," +datum=WGS84 +units=km +no_defs")))

  dat_l$UTMX = locUTM[,1]
  dat_l$UTMY = locUTM[,2]

  dat_l$year = as.integer(format(dat_l$startdatetime, format = "%Y"))

  #Remove observations in years not used
  dat_l = dat_l[dat_l$year %in% conf_l$years, ]

  #Set up structure used for the SPDE-procedure
  meshS=createMesh(conf_l)$mesh
  spdeS = inla.spde2.matern(meshS, alpha=2)
#  spdeS <- fmesher::fm_fem(meshS)
  spdeMatricesS = spdeS$param.inla[c("M0","M1","M2")]
#  spdeMatricesS = list("M0" = spdeS$c0, "M1" = spdeS$g1, "M2" = spdeS$g2)

  A_ListS=list(rep(1,length(conf_l$years)))
  plot(meshS)
  meshST=meshS; spdeST= spdeS; spdeMatricesST = spdeMatricesS;  A_ListST = A_ListS; #Apply same setup for spatial and spatio-temporal part

  singleHauls = which(dat_l$lengthGroup==conf_l$maxLength)

  dist = dat_l$distance[singleHauls]
  yearObs = dat_l$year[singleHauls]
  station = dat_l$station[singleHauls]
  locObs = cbind(dat_l$UTMX[singleHauls],dat_l$UTMY[singleHauls])
  locObsLatLon= cbind(dat_l$longitude[singleHauls],dat_l$latitude[singleHauls])
  depth = dat_l$depth[singleHauls]
  date =  format(dat_l$startdatetime,format="%d;%m;%Y")[singleHauls]

  fishObs=dat_l$catch

  DateTime = dat_l$startdatetime[singleHauls]
  lat = dat_l$latitude[singleHauls]; lon = dat_l$longitude[singleHauls]
  sunAlt=altOfSun(min=as.numeric(format(DateTime,format="%M")),
                  hour=as.numeric(format(DateTime,format="%H")),
                  day=as.numeric(format(DateTime,format="%d")),
                  month=as.numeric(format(DateTime,format="%m")),
                  lat=lat, lon=lon)$alt.of.sun


  sunAlt2=altOfSun(min=ifelse(as.numeric(format(DateTime,format="%M"))>59.8,59.9,
                              as.numeric(format(DateTime,format="%M"))+.1),
                   hour=as.numeric(format(DateTime,format="%H")),
                   day=as.numeric(format(DateTime,format="%d")),
                   month=as.numeric(format(DateTime,format="%m")),
                   lat=lat, lon=lon)$alt.of.sun
  sunSetting=ifelse(sunAlt>=sunAlt2,1,0)

  counter=1
  nStationsEachYear=rep(0,length(conf_l$years))
  for(y in conf_l$years){
    loc=cbind(dat_l$UTMX[dat_l$year==y & dat_l$lengthGroup==conf_l$maxLength],dat_l$UTMY[dat_l$year==y& dat_l$lengthGroup==conf_l$maxLength])

    nStationsEachYear[counter]=dim(loc)[1]
    A_ListS[[counter]]=inla.spde.make.A(meshS,loc)
    A_ListST[[counter]]=inla.spde.make.A(meshST,loc)
#    A_ListS[[counter]]= fmesher::fm_basis(meshS, loc)
#    A_ListST[[counter]]= fmesher::fm_basis(meshST, loc)

    counter = counter + 1
  }

  #If depth is missing:
  missingDepth=which(is.na(depth))
  if(sum(missingDepth>0)>0){
    for(i in 1:length(missingDepth)){
      distances=spDistsN1(locObs,locObs[missingDepth[i],])
      minimumDistanceIndex=which(order(distances)==2)
      depth[missingDepth[i]]=depth[minimumDistanceIndex]
    }
  }
  depth[depth<conf_l$minDepth]=conf_l$minDepth
  depth[depth>conf_l$maxDepth]=conf_l$maxDepth


  fishObsMatrix = t(matrix(fishObs,nrow = length(conf_l$lengthGroups)))#length(unique(dat_l$Station)))
  colnames(fishObsMatrix) = as.character(seq(min(conf_l$lengthGroups),max(conf_l$lengthGroups),by = conf_l$dLength))
  rownames(fishObsMatrix) = paste0(yearObs,"-",station)

  predMatrix = fishObsMatrix *0 #If an element is 1, it is not included in the likelihood
  predMatrix[substring(rownames(predMatrix),1,4) %in% conf_l$skipYears] = 1

  #Define the data and parameters given to TMB----------
  if(conf_l$sunAlt[1] <2){ #Include only one Fourier basis
    sunAltFormula = as.formula(paste( "fishObsMatrix ~ sin+cos"))
    sunAltFormulaIntRep = as.formula(paste( "sunAlt ~ sin+cos"))
  }else{
    stop("To many basis in Fourier approximation of sun heigth effect.")
  }


  sunAltTrans = sunAlt*0
  altMax = rep(-999,length(sunAlt));altMin = rep(999,length(sunAlt))
  for(t in c(seq(0,3,by = 0.1),seq(9,13,by = 0.1), seq(21,24,by = 0.1))){
    hour = floor(t)
    min = (t-hour)*60
    alt = altOfSun(min=min, hour=hour,
                 day=as.numeric(format(DateTime,format="%d")),
                 month=as.numeric(format(DateTime,format="%m")),
                 lat=lat, lon=lon)$alt.of.sun
    altMax[which(altMax<alt)] = alt[which(altMax<alt)]
    altMin[which(altMin>alt)] = alt[which(altMin>alt)]
  }

  sunAltTrans[which(sunSetting==0)] = (sunAlt[which(sunSetting==0)] - altMin[which(sunSetting==0)])/(-altMin[which(sunSetting==0)] +altMax[which(sunSetting==0)]) /2
  sunAltTrans[which(sunSetting!=0)] = 1-(sunAlt[which(sunSetting!=0)] - altMin[which(sunSetting!=0)])/(-altMin[which(sunSetting!=0)] +altMax[which(sunSetting!=0)]) /2
  sunAltTrans[which(sunAltTrans>1)]=1
  sunAltTrans[which(sunAltTrans<0)]=0

  gamSetup_sunAlt=mgcv::gam(sunAltFormula,
                            data=data.frame(haulID=rownames(fishObsMatrix),sin= sin(sunAltTrans *2*pi),
                                            cos=cos(sunAltTrans*2*pi),
                                            fishObsMatrix=fishObsMatrix[,2]),fit=FALSE)

  gamSetup_depth=mgcv::gam(fishObsMatrix~s(depth,bs="cs",k = conf_l$splineDepth[1]),
                           data=data.frame(haulID=rownames(fishObsMatrix),depth=depth,
                                           fishObsMatrix=fishObsMatrix[,2]),fit=FALSE)

  X_sunAlt = gamSetup_sunAlt$X[,-1]
  X_depth = gamSetup_depth$X[,-1]
  S_depth=as(gamSetup_depth$smooth[[1]]$S[[1]], "TsparseMatrix")
  Sdim=nrow(S_depth)

  gamSetup_sunAltReport = gam(sunAltFormulaIntRep,
                              data=data.frame(sunAlt=seq(0,1,by = 0.025),sin=sin(seq(0,1,by = 0.025) *2*pi),
                                              cos = cos(seq(0,1,by = 0.025)*2*pi)),fit=FALSE)
  X_sunAltReport = gamSetup_sunAltReport$X[,-1]

  #Maximum sun height is used in index calculation
  maxSunAlt = 0.5
  gamSetup_sunAltIntegrate = gam(sunAltFormulaIntRep,
                                 data=data.frame(sunAlt=rep(maxSunAlt,2),sin=sin(rep(maxSunAlt,2) *2*pi),
                                                 cos = cos(rep(maxSunAlt,2)*2*pi)),fit=FALSE)
  X_sunAltIntegrate = gamSetup_sunAltIntegrate$X[,-1]

  X_depthReport = PredictMat(gamSetup_depth$smooth[[1]],data = data.frame(depth=as.numeric(seq(min(depth),max(depth),length.out = 40))))

  #Weighting used when smoothing length dimension in latent effect
  weigthLength = conf_l$lengthGroupsReduced*0
  nCollaps = sum(conf_l$lengthGroupsReduced==1)
  counter = 1
  current = 1
  for(i in 1:length(weigthLength)){
    if(current==conf_l$lengthGroupsReduced[i]){
      current = current +1
      weigthLength[i] = 1
      counter = 1
    }else{
      weigthLength[i] = 1-counter/nCollaps
      counter = counter +1
    }
  }

  #Extract areas of each strata
  strata_utm <- st_transform(conf_l$strata,crs=paste0("+proj=utm +zone=",conf_l$zone," +datum=WGS84"))
  areas <- as.numeric(st_area(strata_utm))*(1/1.852)^2/1000000

  #Include observations as a vector (needed for keep functionality)
  obsVector = apply(fishObsMatrix, 1, rbind)
  idxStart = seq(0,length(fishObsMatrix), by = dim(fishObsMatrix)[2])

  data <- list(A_ListS = A_ListS,
               A_ListST = A_ListST,
               dist = dist,
               fishObsMatrix = fishObsMatrix,
               obsVector = fishObs,#obsVector,
               idxStart = idxStart,
               numberOfLengthGroups=length(conf_l$lengthGroups),
               spdeMatricesS = spdeMatricesS,
               spdeMatricesST = spdeMatricesST,
               nStationsEachYear=nStationsEachYear,
               predMatrix = predMatrix,
               X_depth=X_depth,
               S_depth=S_depth,
               Sdim=Sdim,
               X_sunAlt = X_sunAlt,
               X_sunAltReport = X_sunAltReport,
               X_depthReport=X_depthReport,
               lengthGroupsReduced = conf_l$lengthGroupsReduced-1,
               weigthLength = weigthLength,
               areas = areas,
               nBasisSunAlt = max(1,conf_l$sunAlt[1]), #Currently only implemented with one Fourier basis
               obsModel = conf_l$obsModel,
               simulateProcedure = conf_l$simulateProcedure,
               rwBeta0 = conf_l$rwBeta0,
               applyALK = conf_l$applyALK,
               lengthGroups = conf_l$lengthGroups,
               dL = conf_l$dLength,
               trawlWidth = conf_l$trawlWidth,
               lowMemory = conf_l$lowMemory)

  #Configurations used on C-side
  data$spatial=conf_l$spatial
  data$spatioTemporal=conf_l$spatioTemporal
  data$splineDepth=conf_l$splineDepth[2]
  data$useNugget=conf_l$nugget
  data$pcPriorsRange = conf_l$pcPriorRange
  data$pcPriorsSD = conf_l$pcPriorsd
  data$usePCpriors = conf_l$usePcPriors
  data$zeroInflated = conf_l$zeroInflated
  data$sunEffect = conf_l$sunAlt[2]
  data$strataReport = conf_l$strataReport

  attributes(data)$year = yearObs
  attributes(data)$meshS = meshS
  attributes(data)$meshST = meshST
  attributes(data)$locObs = locObs
  attributes(data)$locObsLatLon = locObsLatLon
  attributes(data)$depth = depth #Used when constructing integration points
  attributes(data)$X_sunAltIntegrate = X_sunAltIntegrate #Used when constructing integration points
  attributes(data)$sunAltTrans = sunAltTrans

  #Include integration points
  data = includeIntPoints(data,conf_l,confPred, gamSetup_depth)

  data$selectedStratas = 1:nrow(conf_l$strata)
  if(conf_l$applyALK ==0){
    data = includeDummyAge(data)
  }

  return(data)
}

#' includeIntPoints
#' @param data data in model
#' @param conf_l configurations in model
#' @param confPred prediction configurations
#' @param gamSetup_depth gam setup for depth, needed to set up the spline for integration points
#' @return
#' @export
#' @examples
#' @return
includeIntPoints<-function(data,conf_l,confPred, gamSetup_depth){
  points = constructIntPoints(conf_l,confPred)
  ApredS = inla.spde.make.A(attributes(data)$meshS,loc = as.matrix(points$locUTM))
  ApredST = inla.spde.make.A(attributes(data)$meshST,loc = as.matrix(points$locUTM))
#  ApredS = fmesher::fm_basis(attributes(data)$meshS,loc = as.matrix(points$locUTM))
#  ApredST = fmesher::fm_basis(attributes(data)$meshST,loc = as.matrix(points$locUTM))

  data$ApredS = ApredS
  data$ApredST = ApredST
  data$idxStrata = points$idxStrata

  data$xInt = points$locUTM[,1]
  data$yInt = points$locUTM[,2]

  #Find depth covariate
  if(!is.null(confPred$Depth)) {
    if(grepl(".nc",confPred$Depth)) {
      tryCatch({
        b <- marmap::readGEBCO.bathy(confPred$Depth,res=25)
        bf <- marmap::fortify.bathy(b)
        bf$z <- -1*bf$z
        bf <- subset(bf,z < conf_l$maxDepth & z > conf_l$minDepth)

        bf = st_as_sf(bf,coords=c("x","y"),crs="+proj=longlat")
        bfUTM = st_transform(bf,crs=paste0("+proj=utm +zone=", conf_l$zone," +datum=WGS84 +units=km +no_defs"))

        intPoints = st_as_sf(points$locUTM,coords=c("UTMX","UTMY"),crs=paste0("+proj=utm +zone=", conf_l$zone," +datum=WGS84 +units=km +no_defs"))
        intPoints= st_join(intPoints,bfUTM,join=st_nearest_feature)

        depthGEBCO = intPoints$z

        depthGEBCO[depthGEBCO<conf_l$minDepth]=conf_l$minDepth
        depthGEBCO[depthGEBCO>conf_l$maxDepth]=conf_l$maxDepth

        X_depth = PredictMat(gamSetup_depth$smooth[[1]],data = data.frame(depth=depthGEBCO))
        data$X_depth_int = X_depth },
      error = function(e) {
        confPred$Depth=NULL
        print("No depth data loaded, using depth information from survey stations.") })
    }
    # with depth data from NOAA database
    else if(confPred$Depth == "NOAA") {
      tryCatch({
        b = getNOAA.bathy(lon1 = floor(min(attributes(data)$locObsLatLon[,1])), lon2 = ceiling(max(attributes(data)$locObsLatLon[,1])),
                          lat1 = floor(min(attributes(data)$locObsLatLon[,2])), lat2 = ceiling(max(attributes(data)$locObsLatLon[,2])),
                          resolution = 3)
        bf = fortify.bathy(b)

        bf = st_as_sf(bf,coords=c("x","y"),crs="+proj=longlat")
        bfUTM = st_transform(bf,crs=paste0("+proj=utm +zone=", conf_l$zone," +datum=WGS84 +units=km +no_defs"))

        intPoints = st_as_sf(points$locUTM,coords=c("UTMX","UTMY"),crs=paste0("+proj=utm +zone=", conf_l$zone," +datum=WGS84 +units=km +no_defs"))
        intPoints= st_join(intPoints,bfUTM,join=st_nearest_feature)

        depthNOAA = -intPoints$z

        depthNOAA[depthNOAA<conf_l$minDepth]=conf_l$minDepth
        depthNOAA[depthNOAA>conf_l$maxDepth]=conf_l$maxDepth

        X_depth = PredictMat(gamSetup_depth$smooth[[1]],data = data.frame(depth=depthNOAA[minDist]))

        data$X_depth_int = X_depth },
        error = function(e) {
          confPred$Depth=NULL
          print("Not able to access NOAA database, using depth information from survey stations.") })
    } else {
      confPred$Depth=NULL
    }
  }
  if(is.null(confPred$Depth)) {
    obs = st_as_sf(data.frame(UTMX=attributes(data)$locObs[,1],UTMY=attributes(data)$locObs[,2],depth=attributes(data)$depth),coords=c("UTMX","UTMY"),crs=paste0("+proj=utm +zone=", conf_l$zone," +datum=WGS84 +units=km +no_defs"))
    intPoints = st_as_sf(points$locUTM,coords=c("UTMX","UTMY"),crs=paste0("+proj=utm +zone=", conf_l$zone," +datum=WGS84 +units=km +no_defs"))

    intPoints= st_join(intPoints,obs,join=st_nearest_feature)

    depthDATA = intPoints$depth

    X_depth = PredictMat(gamSetup_depth$smooth[[1]],data = data.frame(depth=depthDATA))
    data$X_depth_int = X_depth
    print("No depth information for prediction provided, using depth from observations.")
  }

  #Set sun altitude given at maximum
  X_sunAltIntegrateTmp = attributes(data)$X_sunAltIntegrate
  X_sunAltIntegrate = matrix(0,dim(points$locUTM)[1],dim(X_sunAltIntegrateTmp)[2])
  for(i in 1:dim(points$locUTM)[1]){
    X_sunAltIntegrate[i,1:dim(X_sunAltIntegrateTmp)[2]] = X_sunAltIntegrateTmp[1,1:dim(X_sunAltIntegrateTmp)[2]]
  }
  data$X_sunAltIntegrate = X_sunAltIntegrate
  return(data)
}



#' combineLengthALK combines data needed for length and age needed in the joint model for catch-at-length and the age-at-length model
#' @param dat_l data in catch-at-length model
#' @param par_l parameters in catch-at-length model
#' @param map_l map-variable in catch-at-length model
#' @param dat_a data in age-at-length model
#' @param par_a parameters in age-at-length model
#' @param map_a map-variable in age-at-length model
#' @param conf_l configurations in catch-at-length model
#' @param confPred prediction configurations
#' @return A list with combined data for catch-at-length and age-at-length needed by the joint model.
#' @export
#' @examples
#' @return
combineLengthALK<-function(dat_l, par_l,map_l,dat_a, par_a,map_a,conf_l,confPred){

  dat_joint = c(dat_l,dat_a)
  par_joint = c(par_l,par_a)
  map_joint = c(map_l,map_a)

  points = constructIntPoints(conf_l,confPred)
  Apred_alk = inla.spde.make.A(attributes(dat_a)$mesh,loc = as.matrix(points$locUTM))
#  Apred_alk = fmesher::fm_basis(attributes(dat_a)$mesh,loc = as.matrix(points$locUTM))

  dat_joint$Apred_alk = Apred_alk
  attributes(dat_joint)$year = attributes(dat_l)$year
  attributes(dat_joint)$meshS = attributes(dat_l)$meshS
  attributes(dat_joint)$meshST = attributes(dat_l)$meshST
  attributes(dat_joint)$locObs = attributes(dat_l)$locObs
  attributes(dat_joint)$locObsLatLon = attributes(dat_l)$locObsLatLon
  attributes(dat_joint)$depth = attributes(dat_l)$depth #Used when constructing integration points
  attributes(dat_joint)$X_sunAltIntegrate = attributes(dat_l)$X_sunAltIntegrate #Used when constructing integration points
  attributes(dat_joint)$sunAltTrans = attributes(dat_l)$sunAltTrans

  return(list(dat_joint = dat_joint, par_joint = par_joint,map_joint = map_joint))
}

