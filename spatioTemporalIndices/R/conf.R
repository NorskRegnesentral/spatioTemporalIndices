
##' defConf
##'
##' Configurations used for catch-at-length. The object returned is
##' used in \code{\link{setupData}} and in \code{\link{setPar}}.
##'
##' @param years A numeric vector specifying the range of years to include. All years between the earliest and latest values in the vector are used. For example, if `years = c(2018, 2020)`, the data for 2018, 2019, and 2020 is included. Use the `skipYears` parameter to exclude specific years within this range.
##' @param skipYears Years to skip in the estimation within the range of `years`.
##' @param spatial If 0-no spatial effect, 1-Include spatial effect
##' @param spatioTemporal If 0-no spatial-temporal effect, 1-Include spatial-temporal effect, 2-Include spatial-temporal effect with independence between years. Note: Setting this to 1 makes inference time consuming and include correlation between years.
##' @param nugget If 0-no nugget effect, 1-Include nugget effect. Note that this nugget effect includes correlation between length groups within the haul.
##' @param splineDepth Vector with depth configurations, first element is the k-variable in spline (6 means 6 basis function). Second element: 0: no depth effect, 1: depth effect(length independent), 2:length dependent depth effect.
##' @param sunAlt Vector with sun altitude configurations, first element is the number of basis functions. Second element: 0: No sun effect, 1:sun effect (length independent), 2: length dependent sun effect.
##' @param maxLength Maximum length, used when defining the pluss group
##' @param minLength Minimum length, used when defining the minus group (probably never used)
##' @param dLength Length of each length group bin
##' @param reduceLength The resolution in length dimension used in latent effect. If 1: all length groups are included as in latent effect. If 2: latent dimension for length group is defined as every second length group... Including all length groups in latent effect often result in high memory usage and computation time.
##' @param cutoff Cutoff used in mesh, a smaller value results in a more dense mesh.
##' @param cbound boundary configurations for the mesh. The first element is how far to extend the inner boundary, and the second is for the outer boundary.
##' @param pcPriorRange pc-priors for spatial and spatio-temporal effect (probably never used).
##' @param pcPriorsd pc-priors for spatial and spatio-temporal effect (probably never used).
##' @param usePcPriors If 0-No pc-priors used, 1- use pc-priors (probably never used).
##' @param zeroInflated If 0-No zero inflation, 1- use zero inflation (NB!! not fully implemented).
##' @param stratasystem List containing dsn and layer of shapefile
##' @param obsModel Observation model: 1- negative binomial (not fully implemented), 2-Poisson (most relevant when including a nugget effect), 3-Tweedie (not fully implemented)
##' @param rwBeta0 Include intercept of catch-at-length as a random walk? 0-No, 1-Yes
##' @param applyALK Calculate indices-at-age by using the age-at-length model: 0-No , 1-yes
##' @param mapRhoL How to couple correlation parameters for length in space, space-time and nugget; For experimental use
##' @param simulateProcedure How to do simulations; For experimental use
##' @param minDepth Minimum depth considered
##' @param maxDepth Maximum depth considered
##' @param trawlWidth Swept width of trawl; Used when defining the unit when comparing with other estimation methods.
##' @param plusGroup Include plus group? 1:yes, 0: No
##' @param strataReport ADREPORT index in each strata? 1:yes, 0: No. NB!: Currently not working in combination with ALK
##' @param lowMemory 0-No, 1-Yes. If 1: Uses less memory by not reporting uncertainties of splines and more.
##' @details This function sets up the configurations for the catch-at-length model
##' @export
defConf <- function(years, skipYears=NULL,spatial = 1,spatioTemporal = 0,nugget = 1,splineDepth=c(6,0),sunAlt=c(1,0),
                    maxLength = NULL,minLength = NULL,dLength = 1, reduceLength = 3,
                    cutoff = 100,cbound = c(18,130),
                    pcPriorRange = c(100,0.1),pcPriorsd = c(1,0.1), usePcPriors = 0, zeroInflated = 0,
                    stratasystem = list(),
                    obsModel = 2,rwBeta0 = 1,applyALK = 0,
                    mapRhoL = c(0,1,2), simulateProcedure = 1,
                    minDepth=50,maxDepth=600,trawlWidth=1,
                    plusGroup = 1, strataReport = 0, lowMemory = 0){
  conf= list()
  conf$years = min(years):max(years)
  conf$skipYears = skipYears
  # Numeric = use directly; Data = use input data to determine
  if(is.numeric(minLength)){
    conf$minLength = minLength
  }else{
    stop("Need minLength")
  }

  if(is.numeric(maxLength)){
    conf$maxLength=maxLength
  }else{
    stop("Need maxLength")
  }

  conf$dLength = dLength
  conf$lengthGroups = seq(conf$minLength,conf$maxLength,by = conf$dLength)

  if(length(stratasystem)>0) {
    strata = sf::st_read(dsn = stratasystem[[1]], layer = stratasystem[[2]],quiet = TRUE) # dsn path must be adjusted to folder location
    conf$strata = sf::st_transform(strata,crs="+proj=longlat")
    conf$strata$id = row.names(conf$strata)
    conf$zone = floor((mean(sf::st_bbox(conf$strata)[c(1,3)]) + 180) / 6) + 1
    conf$stratasystem = stratasystem
    conf$strata_number = nrow(conf$strata)
  } else {
    conf$strata=NULL;conf$zone=NULL;conf$stratasystem=NULL;conf$strata_number=NULL
  }
  conf$minDepth = minDepth
  conf$maxDepth = maxDepth
  conf$sunAlt = sunAlt
  conf$splineDepth=splineDepth
  conf$spatial = spatial
  conf$spatioTemporal = spatioTemporal
  conf$cutoff = cutoff
  conf$cbound = cbound
  lTmp = rep(1:100,each = reduceLength)
  conf$lengthGroupsReduced = lTmp[1:length(conf$lengthGroups)]
  conf$reduceLength = reduceLength
  conf$pcPriorRange = pcPriorRange
  conf$pcPriorsd= pcPriorsd
  conf$usePcPriors = usePcPriors
  conf$zeroInflated = zeroInflated
  if(zeroInflated==1){
    warning("Zero-inflation not fully implemented")
  }
  conf$obsModel = obsModel
  conf$nugget = nugget
  conf$mapRhoL = mapRhoL
  conf$simulateProcedure = simulateProcedure
  conf$rwBeta0 = rwBeta0
  conf$applyALK = applyALK
  conf$trawlWidth = trawlWidth
  conf$plusGroup =plusGroup
  conf$strataReport = strataReport
  conf$lowMemory = lowMemory
  return(conf)
}

##' defConfPred
##'
##'Configurations used for prediction.
##'
##' @param conf Configurations used when fitting the model
##' @param cellsize provide distance between integration points
##' @param Depth if "NOAA": use NOAA data base for estimating depth in integration points; if GEBCO (.nc) file: use file for estimating depth in integration
##' @details This function sets up the configurations for the predictions, i.e. the index-at-length and index-at-age.
##' @export
defConfPred <- function(conf,cellsize=20,Depth="Data"){
  confPred = list()
  confPred$Strata=1:conf$strata_number
  confPred$cellsize = cellsize
  confPred$Depth = Depth
  return(confPred)
}

##' setMap
##'
##'\code{setMap} defines the map-argument used in MakeADFun.
##'
##' @param par Parameters included
##' @param conf Configurations
##' @details This function sets up the map-list used to couple and simplify the model with the map-argument used by MakeADFun in TMB
##' @export
setMap <- function(par, conf){
  map= list()
  map$log_kappa = as.factor(c(0,1))
  map$tan_rho_l = conf$mapRhoL #Default use same parameter for length correlation

  if(conf$rwBeta0==0){
    map$log_sigma_beta0 = as.factor(rep(NA,length(par$log_sigma_beta0)))
  }

  if(conf$nugget[1] ==0){
    map$log_sigma = as.factor(c(0,1,NA))
    map$nugget = as.factor(rep(NA,length(par$nugget)))
    map$tan_rho_l[3] = NA
  }else{
    map$log_sigma = as.factor(c(0,1,2))
  }


  if(conf$sunAlt[2]==0){
    map$betaSun = as.factor(rep(NA,length(par$betaSun))) #Not use time in day
  }else if(conf$sunAlt[2]==1){
    tmp = 0:(conf$sunAlt[1]*2-1)
    map$betaSun = as.factor(c(tmp,tmp)) #Not use length dependent time in day
  }
  if(length(conf$years)==1){
    map$tan_rho_t = as.factor(NA)
    if(conf$spatial==1 & conf$spatioTemporal==1){
      warning("Using both spatial and spatio-temporal with only one year of data will overparametrize the model, turn of spatio-temporal contribution")
      conf$spatioTemporal=0;
    }
  }
  if(conf$splineDepth[2]==0){
    map$betaDepth=as.factor(rep(NA,length(par$betaDepth)))
    map$log_lambda=as.factor(c(NA,NA))
  }else if(conf$splineDepth[2]==1){
    tmp = 0:(length(par$betaDepth)/2-1)
    map$log_lambda=as.factor(c(0,0))
    map$betaDepth = as.factor(c(tmp,tmp)) #Not use length dependent depth effect
  }else if(conf$splineDepth[2]==2){
    #map$log_lambda=as.factor(c(0,0))#Use same lambda in both depth splines
  }


  if(conf$spatial ==0){
    map$xS = as.factor(rep(NA,length(par$xS)))
    if(conf$nugget[1] ==0){
      map$log_sigma= as.factor(c(NA,0,NA))
    }else{
      map$log_sigma= as.factor(c(NA,0,1))
    }
    map$log_kappa= as.factor(c(NA,0))
    map$tan_rho_l[1]= NA
    map$tan_rho_l[2:3]  = map$tan_rho_l[2:3]-1
  }
  if(conf$spatioTemporal ==2){
    map$tan_rho_t = as.factor(NA)
  }
  if(conf$spatioTemporal ==0){
    map$xST = as.factor(rep(NA,length(par$xST)))
    map$tan_rho_t = as.factor(NA)
    if(conf$spatial ==0){
      map$log_kappa= as.factor(c(NA,NA))
      if(conf$nugget[1] ==0){
        map$log_sigma= as.factor(c(NA,NA,NA))
        map$tan_rho_l[3]= NA
      }else{
        map$log_sigma= as.factor(c(NA,NA,1))
        map$tan_rho_l[3]= 0
      }
      map$tan_rho_l[1:2]= c(NA,NA)

    }else{
      map$log_kappa= as.factor(c(0,NA))
      if(conf$nugget[1] ==0){
        map$log_sigma= as.factor(c(0,NA,NA))
        map$tan_rho_l[3] =   NA
      }else{
        map$log_sigma= as.factor(c(0,NA,1))
        map$tan_rho_l[3] =   map$tan_rho_l[3]-1
      }
      map$tan_rho_l[2] =  NA
    }
  }

  if(conf$obsModel==2){
    map$logSize = as.factor(NA)
  }

  if(conf$zeroInflated==0){
    map$delta_z = as.factor(rep(NA,length(par$delta_z)))
  }else if(conf$zeroInflated==1){
    map$delta_z = as.factor(c(0,1,NA))
  }else if(conf$zeroInflated==2){
    map$delta_z = as.factor(c(0,NA,NA))
  }else{
    stop("No such zero inflation included")
  }

  map$tan_rho_l = as.factor(map$tan_rho_l)

  if(conf$applyALK==0){
    map$beta0_alk = as.factor(rep(NA, length(par$beta0_alk)))
    map$log_sigma_beta_alk = as.factor(rep(NA, length(par$log_sigma_beta_alk)))
    map$betaLength_alk = as.factor(rep(NA, length(par$betaLength_alk)))
    map$logSigma_alk = as.factor(rep(NA, length(par$logSigma_alk)))
    map$logKappa_alk = as.factor(rep(NA, length(par$logKappa_alk)))
    map$transRho_alk = as.factor(rep(NA, length(par$transRho_alk)))
    map$xST_alk = as.factor(rep(NA, length(par$xST_alk)))
    map$xS_alk = as.factor(rep(NA, length(par$xS_alk)))
  }

  if(length(conf$skipYears)>0) {
    if(conf$rwBeta0 ==0){ #Free beta0, and needs to not estimate these in years without data
      skipped = which(conf$years==conf$skipYears)
      map$beta0 = matrix(data=c(1:length(par$beta0)),ncol=ncol(par$beta0),nrow=nrow(par$beta0))
      map$beta0[skipped,] <- NA
      map$beta0 = as.factor(map$beta0)
    }
  }

  return(map)
}

