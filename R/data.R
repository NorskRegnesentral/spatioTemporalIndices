

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
setupData = function(dataLength,conf,confPred){

  #Remove too short lengths
  dataLength = dataLength[dataLength$lengthGroup>=conf$minLength,]

  #Remove too long lengths
  for(id in unique(dataLength$station)){
    index =which(dataLength$station==id & dataLength$lengthGroup>=conf$maxLength)
    if(conf$plusGroup==1){
      dataLength$catch[index[1]] = sum(dataLength$catch[index])
    }
  }
  if(length(which(dataLength$lengthGroup> conf$maxLength)>0)){
    dataLength = dataLength[-which(dataLength$lengthGroup> conf$maxLength),]
  }

  #Convert to UTM coordinates
  loc = data.frame(dataLength$longitude,dataLength$latitude)
  names(loc) = c("X","Y")
  attr(loc, "projection") = "LL"
  attr(loc, "zone") = conf$zone
  locUTM <- PBSmapping::convUL(loc)
  colnames(locUTM) = c("UTMX", "UTMY")

  dataLength$UTMX = locUTM$UTMX
  dataLength$UTMY = locUTM$UTMY

  dataLength$year = as.integer(format(dataLength$startdatetime, format = "%Y"))

  #Remove observations in years not used
  dataLength = dataLength[dataLength$year %in% conf$years, ]

  #Set up structure used for the SPDE-procedure
  meshS=createMesh(conf)$mesh
  spdeS = inla.spde2.matern(meshS, alpha=2)
  spdeMatricesS = spdeS$param.inla[c("M0","M1","M2")]
  A_ListS=list(rep(1,length(conf$years)))
  plot(meshS)
  meshST=meshS; spdeST= spdeS; spdeMatricesST = spdeMatricesS;  A_ListST = A_ListS; #Apply same setup for spatial and spatio-temporal part

  singleHauls = which(dataLength$lengthGroup==conf$maxLength)

  dist = dataLength$distance[singleHauls]
  yearObs = dataLength$year[singleHauls]
  station = dataLength$station[singleHauls]
  locObs = cbind(dataLength$UTMX[singleHauls],dataLength$UTMY[singleHauls])
  locObsLatLon= cbind(dataLength$longitude[singleHauls],dataLength$latitude[singleHauls])
  depth = dataLength$depth[singleHauls]
  date =  format(dataLength$startdatetime,format="%d;%m;%Y")[singleHauls]

  fishObs=dataLength$catch

  DateTime = dataLength$startdatetime[singleHauls]
  lat = dataLength$latitude[singleHauls]; lon = dataLength$longitude[singleHauls]
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
  nStationsEachYear=rep(0,length(conf$years))
  for(y in conf$years){
    loc=cbind(dataLength$UTMX[dataLength$year==y & dataLength$lengthGroup==conf$maxLength],dataLength$UTMY[dataLength$year==y& dataLength$lengthGroup==conf$maxLength])

    nStationsEachYear[counter]=dim(loc)[1]
    A_ListS[[counter]]=inla.spde.make.A(meshS,loc)
    A_ListST[[counter]]=inla.spde.make.A(meshST,loc)
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
  depth[depth<conf$minDepth]=conf$minDepth
  depth[depth>conf$maxDepth]=conf$maxDepth


  fishObsMatrix = t(matrix(fishObs,nrow = length(conf$lengthGroups)))#length(unique(dataLength$Station)))
  colnames(fishObsMatrix) = as.character(seq(min(conf$lengthGroups),max(conf$lengthGroups),by = conf$dLength))
  rownames(fishObsMatrix) = paste0(yearObs,"-",station)

  predMatrix = fishObsMatrix *0 #If an element is 1, it is not included in the likelihood
  predMatrix[substring(rownames(predMatrix),1,4) %in% conf$skipYears] = 1

  #Define the data and parameters given to TMB----------
  if(conf$sunAlt[1] <2){ #Include only one Fourier basis
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

  gamSetup_depth=mgcv::gam(fishObsMatrix~s(depth,bs="cs",k = conf$splineDepth[1]),
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

  X_depthReport = PredictMat(gamSetup_depth$smooth[[1]],data = data.frame(depth=as.numeric(seq(min(depth),max(depth),length.out = 20))))

  #Weighting used when smoothing length dimension in latent effect
  weigthLength = conf$lengthGroupsReduced*0
  nCollaps = sum(conf$lengthGroupsReduced==1)
  counter = 1
  current = 1
  for(i in 1:length(weigthLength)){
    if(current==conf$lengthGroupsReduced[i]){
      current = current +1
      weigthLength[i] = 1
      counter = 1
    }else{
      weigthLength[i] = 1-counter/nCollaps
      counter = counter +1
    }
  }

  #Extract areas of each strata
  strata_utm <- spTransform(conf$strata,CRS(paste0("+proj=utm +zone=",conf$zone," +datum=WGS84")))
  areas <- raster::area(strata_utm)*(1/1.852)^2/1000000

  #Include observations as a vector (needed for keep functionality)
  obsVector = apply(fishObsMatrix, 1, rbind)
  idxStart = seq(0,length(fishObsMatrix), by = dim(fishObsMatrix)[2])

  data <- list(A_ListS = A_ListS,
               A_ListST = A_ListST,
               dist = dist,
               fishObsMatrix = fishObsMatrix,
               obsVector = fishObs,#obsVector,
               idxStart = idxStart,
               numberOfLengthGroups=length(conf$lengthGroups),
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
               lengthGroupsReduced = conf$lengthGroupsReduced-1,
               weigthLength = weigthLength,
               areas = areas,
               nBasisSunAlt = max(1,conf$sunAlt[1]), #Currently only implemented with one Fourier basis
               obsModel = conf$obsModel,
               simulateProcedure = conf$simulateProcedure,
               rwBeta0 = conf$rwBeta0,
               applyALK = conf$applyALK,
               lengthGroups = conf$lengthGroups,
               dL = conf$dLength,
               trawlWidth = conf$trawlWidth,
               lowMemory = conf$lowMemory)

  #Configurations used on C-side
  data$spatial=conf$spatial
  data$spatioTemporal=conf$spatioTemporal
  data$splineDepth=conf$splineDepth[2]
  data$useNugget=conf$nugget
  data$pcPriorsRange = conf$pcPriorRange
  data$pcPriorsSD = conf$pcPriorsd
  data$usePCpriors = conf$usePcPriors
  data$zeroInflated = conf$zeroInflated
  data$sunEffect = conf$sunAlt[2]
  data$strataReport = conf$strataReport

  attributes(data)$year = yearObs
  attributes(data)$meshS = meshS
  attributes(data)$meshST = meshST
  attributes(data)$locObs = locObs
  attributes(data)$locObsLatLon = locObsLatLon
  attributes(data)$depth = depth #Used when constructing integration points
  attributes(data)$X_sunAltIntegrate = X_sunAltIntegrate #Used when constructing integration points
  attributes(data)$sunAltTrans = sunAltTrans

  #Include integration points
  data = includeIntPoints(data,conf,confPred, gamSetup_depth)

  data$selectedStratas = 1:nrow(conf$strata)
  if(conf$applyALK ==0){
    data = includeDummyAge(data)
  }

  return(data)
}

#' includeIntPoints
#' @return
#' @export
#' @examples
#' @return
includeIntPoints<-function(data,conf,confPred, gamSetup_depth){
  points = constructIntPoints(conf,confPred)
  ApredS = inla.spde.make.A(attributes(data)$meshS,loc = as.matrix(points$locUTM))
  ApredST = inla.spde.make.A(attributes(data)$meshST,loc = as.matrix(points$locUTM))

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
        bf <- marmap::fortify.bathy(bathy)
        bf$z <- -1*bf$z
        bf <- subset(bf,z < conf$maxDepth & z > conf$minDepth)

        loc = data.frame(bf$x,bf$y)
        names(loc) = c("X","Y")
        attr(loc, "projection") = "LL"
        attr(loc, "zone") = conf$zone
        locUTM <- PBSmapping::convUL(loc)
        colnames(locUTM) = c("UTMX", "UTMY")

        obs = SpatialPoints(locUTM,CRS(paste0("+proj=utm +zone=", conf$zone," +datum=WGS84 +units=km +no_defs")))

        intPoints = SpatialPoints(points$locUTM,CRS(paste0("+proj=utm +zone=", conf$zone," +datum=WGS84 +units=km +no_defs")))

        dist = gDistance(obs,intPoints, byid=T)
        minDist <- apply(dist, 1, function(x) order(x, decreasing=F)[1])

        depthGEBCO = -bf$z

        depthGEBCO[depthGEBCO<conf$minDepth]=conf$minDepth
        depthGEBCO[depthGEBCO>conf$maxDepth]=conf$maxDepth

        X_depth = PredictMat(gamSetup_depth$smooth[[1]],data = data.frame(depth=depthNOAA[minDist]))
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

        loc = data.frame(bf$x,bf$y)
        names(loc) = c("X","Y")
        attr(loc, "projection") = "LL"
        attr(loc, "zone") = conf$zone
        locUTM <- PBSmapping::convUL(loc)
        colnames(locUTM) = c("UTMX", "UTMY")

        obs = SpatialPoints(locUTM,CRS(paste0("+proj=utm +zone=", conf$zone," +datum=WGS84 +units=km +no_defs")))

        intPoints = SpatialPoints(points$locUTM,CRS(paste0("+proj=utm +zone=", conf$zone," +datum=WGS84 +units=km +no_defs")))

        dist = gDistance(obs,intPoints, byid=T)
        minDist <- apply(dist, 1, function(x) order(x, decreasing=F)[1])

        depthNOAA = -bf$z

        depthNOAA[depthNOAA<conf$minDepth]=conf$minDepth
        depthNOAA[depthNOAA>conf$maxDepth]=conf$maxDepth

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
    obs = SpatialPoints(attributes(data)$locObs,CRS(paste0("+proj=utm +zone=", conf$zone," +datum=WGS84 +units=km +no_defs")))
    intPoints = SpatialPoints(points$locUTM,CRS(paste0("+proj=utm +zone=", conf$zone," +datum=WGS84 +units=km +no_defs")))
    dist = gDistance(obs,intPoints, byid=T)
    minDist <- apply(dist, 1, function(x) order(x, decreasing=F)[1])

    xDepth = matrix(0,dim(points$locUTM)[1],dim(data$X_depth)[1])
    xDepth = data$X_depth[minDist,]
    data$X_depth_int = xDepth
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



#' combineLengthALK
#' @return
#' @export
#' @examples
#' @return
combineLengthALK<-function(dat_l, par_l,map_l,dat_a, par_a,map_a,conf,confPred){

  dat_joint = c(dat_l,dat_a)
  par_joint = c(par_l,par_a)
  map_joint = c(map_l,map_a)

  points = constructIntPoints(conf,confPred)
  Apred_alk = inla.spde.make.A(attributes(dat_a)$mesh,loc = as.matrix(points$locUTM))
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

