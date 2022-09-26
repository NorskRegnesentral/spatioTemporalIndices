
##' defConf
##'
##' Configurations used. The object returned is
##' used in \code{\link{setupData}} and in  \code{\link{setPar}}.
##'
##' @param years Years of data included
##' @param skipYears Years to skip in the estimation
##' @param spatial If 0-no spatial effect, 1-Include spatial effect
##' @param spatioTemporal If 0-no spatial-temporal effect, 1-Include spatial-temporal effect, 2-Include spatial-temporal effect with independence between years
##' @param nugget If 0-no nugget effect, 1-Include nugget effect
##' @param splineDepth Vector with depth configurations, first element is the k-variable in spline (6 means 6 basis function). Second element: 0- no depth effect, 1-depth effectk, 2-length dependent depth effect.
##' @param sunAlt Vector with sun altitude configurations, first element is the numberof basis functions. Second element: 0-sun effect, 1-length dependent sun effect.
##' @param maxLength Maximum length, used when defining the pluss group
##' @param minLength Minimum length, used when defining the minus group (probably never used)
##' @param reduceLength The resolution in length dimension used in latent effect
##' @param cutoff Cutoff used in mesh, first element is for spatial effect, second for spatio-temporal effect, third is 0 if same mesh is used for both spatial and spatio-temporal effect
##' @param cbound cbound used in mesh, first element is for spatial effect, second for spatio-temporal effect
##' @param pcPriorRange pc-priors for spatial and spatio-temporal effect
##' @param pcPriorsd pc-priors for spatial and spatio-temporal effect
##' @param usePcPriors If 0-No pc-priors used, 1- use pc-priors (probably never used)
##' @param zeroInflated If 0-No zero inflation, 1- use zero inflation (NB!! not fully implemented)
##' @param stratasystem List containing dsn and layer of shapefile
##' @param minDepth Minimum depth considered
##' @param maxDepth Maximum depth considered
##' @param trawlWidth Swept width of trawl
##' @param plusGroup Include plus goup? 1:yes, 0: No
##' @details
##' @export
defConf <- function(years, skipYears=NULL,spatial = 1,spatioTemporal = 0,nugget = 1,splineDepth=c(6,1),sunAlt=c(1,0),
                    maxLength = NULL,dLength = 1, minLength = NULL,reduceLength = 3,
                    cutoff = 100, cbound = 200,
                    pcPriorRange = c(100,0.1),pcPriorsd = c(1,0.1), usePcPriors = 0, zeroInflated = 0,
                    obsModel = 2,rwBeta0 = 1,applyALK = 0,
                    mapRhoL = c(0,1,2), simulateProcedure = 1,
                    stratasystem = list(),
                    minDepth=50,maxDepth=600,trawlWidth=1,
                    plusGroup = 1){
  conf= list()
  conf$years = years
  conf$skipYears = skipYears
  # Numeric = use directly; Data = use input data to determine
  if(is.numeric(minLength)) conf$minLength = minLength
  if(!is.numeric(minLength)) conf$minLength = ifelse(floor(min(ld$IndividualTotalLength,na.rm=T)) %% dLength == 0,
                     floor(min(ld$IndividualTotalLength,na.rm=T)),floor(min(ld$IndividualTotalLength,na.rm=T))-(dLength-1))

  if(is.numeric(maxLength)) conf$maxLength=maxLength
  if(!is.numeric(maxLength)) {
    print("Using maximum observed length as maxLength")
    conf$maxLength = ceiling(max(ld$IndividualTotalLength)) }
  conf$dLength = dLength
  conf$lengthGroups = seq(conf$minLength,conf$maxLength,by = conf$dLength)

  if(length(stratasystem)>0) {
    strata = readOGR(dsn = stratasystem[[1]], layer = stratasystem[[2]]) # dsn path must be adjusted to folder location
    conf$strata = spTransform(strata,CRS("+proj=longlat"))
    conf$zone = floor((mean(bbox(conf$strata)[1,]) + 180) / 6) + 1
    conf$stratasystem = stratasystem
    conf$strata_number = nrow(conf$strata)
  } else {
    conf$strata=NULL;conf$zone=NULL;conf$stratasystem=NULL;conf$strata_number=NULL
  }
  conf$minDepth = minDepth
  conf$maxDepth = maxDepth

  #Covariates
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

  conf$obsModel = obsModel
  conf$nugget = nugget
  conf$mapRhoL = mapRhoL
  conf$simulateProcedure = simulateProcedure
  conf$rwBeta0 = rwBeta0
  conf$applyALK = applyALK

  conf$trawlWidth = trawlWidth

  conf$plusGroup =plusGroup
  conf$smartStart = NULL #TODO
  return(conf)
}

##' defConfPred
##'
##'Configurations used for prediction.
##'
##' @param conf Configurations used when fitting the model
##' @param nIntPoints provide number of integration points
##' @param Depth if "NOAA": use NOAA data base for estimating depth in integration points
##' @details
##' @export
defConfPred <- function(conf,nIntPoints= 4000,Depth="Data"){
  confPred = list()
  confPred$Strata=1:conf$strata_number
  confPred$nIntPoints = nIntPoints
  confPred$Depth = Depth
  return(confPred)
}

##' setMap
##'
##'\code{setMap} defines the map-argument used in MakeADFun.
##'
##' @param par Parameters included
##' @param conf Configurations
##' @details
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
    map$nuggetIndex = as.factor(rep(NA,length(par$nuggetIndex)))
    map$tan_rho_l[3] = NA
  }else{
    map$log_sigma = as.factor(c(0,1,2))
  }


  if(conf$sunAlt[1]==0){
    map$betaSun = as.factor(rep(NA,length(par$betaSun))) #Not use time in day
  }else if(conf$sunAlt[2]==0){
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
  }else if(conf$splineDepth[2]==1){
    map$log_lambda=as.factor(c(0,0))#Use same lambda in both depth splines
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
  }else{
    warning("Map-functionality not yet implemented for selected zero-inflation procedure.")
  }

  map$tan_rho_l = as.factor(map$tan_rho_l)

  if(conf$applyALK==0){
    map$beta0_alk = as.factor(rep(NA, length(par$beta0_alk)))
    map$log_sigma_beta0_alk = as.factor(NA)
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

