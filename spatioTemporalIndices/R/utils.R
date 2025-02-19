##' loglik
##'
##' @param object The fitted model retuned by fitModel
##' @param ... Extra argument
##' @export
logLik.stim<-function(object, ...){
  ret<- -object$opt$objective

  sunCov = 0
  if(object$conf_l$sunAlt[2]==1)sunCov = object$conf_l$sunAlt[1]*2
  if(object$conf_l$sunAlt[2]==2)sunCov = object$conf_l$sunAlt[1]*4

  beta0 = 0
  if(object$conf_l$rwBeta0==0){
    beta0 = beta0+ length(object$conf_l$lengthGroups)*length(object$conf_l$years)
  }else{
    beta0 = beta0+ length(object$conf_l$lengthGroups)
  }

  betaLengthALK = 0
  if(object$conf_l$applyALK!=0){
    if(object$conf_alk$rwBeta0==0){
      beta0 = beta0+ length(object$conf_alk$minAge:object$conf_alk$maxAge )*length(object$conf_l$years)
    }else{
      beta0 = beta0+ length(object$conf_alk$minAge:object$conf_alk$maxAge) + length(object$conf_l$years) -1
    }
    betaLengthALK = length(object$pl$betaLength_alk)

  }
  attr(ret,"df")<-length(object$opt$par) + beta0 + sunCov + betaLengthALK
  class(ret)<-"logLik"
  ret
}


##' Print stim object
##' @method print stim
##' @param  x fitted model
##' @param  ... Extra argument
##' @details Print log-likelihood and the main convergence criteria
##' @export
print.stim<-function(x,...){
  cat("STIM model: log likelihood is", logLik.stim(x),"Convergence", ifelse(0==x$opt$convergence, "OK\n", "failed\n"))
}

# ##' Simulate from a STIM object
# ##' @param  run
# ##' @param  nsim
# ##' @param  seed
# ##' @details Simulate observations from STIM as explained in Breiviek et al. (2020).
# ##' @export
#simulate.stim<-function(run, nsim=1, seed=NULL, ...){
#  if(!is.null(seed)) set.seed(seed)
#  pl <- as.list(run$rep,"Est")
#  est <- unlist(pl)
#  ret <- replicate(nsim,
#                   c(run$data[names(run$data)!="fishObsMatrix"],#all the old data
#                     run$obj$simulate()["fishObsMatrix"])#simulated observations
#                   , simplify=FALSE)
#
#  ret
#}


##' Find sun height and sun rise. Thanks to Espen Johnsen who provided the code.
##' @param min Minutes in hour
##' @param hour Hour in day
##' @param day Day of month
##' @param month Month of year
##' @param lat Latitude of station
##' @param lon Longitude of station
##' @export
altOfSun <- function(min, hour, day, month,lat, lon){
  # altitude of sun
  UTC <- hour+min/60
  CET <- (UTC + 1) %% 24
  dayadd <- cumsum(c(0,31,28,31,30,31,30,31,31,30,31,30,31))
  cumday <- day + dayadd[month]
  K1 <- (lon - 15 - 0.4083 * sin(0.0172 * (cumday-80))
         - 1.7958 * cos(0.0172 * (cumday-80))
         + 2.4875 * sin(0.0344 * (cumday-80)))
  SST <- ((CET*15) + K1) / (180/pi)
  dkl <- asin(0.3979 * sin((0.0172 * (cumday - 80))
                           + 0.03346 * (sin(0.0172 * cumday) - 0.98112)))
  Brq <- lat/(180/pi)
  sinush <- (sin(dkl)*sin(Brq)) - (cos(dkl)*cos(Brq)*cos(SST))
  alt.of.sun <- asin(sinush) * (180/pi)

  # time when altitude of sun = asun.0
  asun.0 <- 0
  K2 <- (sin(dkl)*sin(Brq) - sin(asun.0/(180/pi))) / (cos(dkl)*cos(Brq))
  K2[K2 < (-1)] <- -1        # polar night
  K2[K2 > ( 1)] <-  1         # midnight sun
  SST0 <- acos(K2)
  CET0 <- (SST0 * (180/pi) - K1) / 15
  UTC0 <- (CET0 - 1) + 24*(CET0 < 1)
  sun.rise <- UTC0%%24
  list(alt.of.sun=alt.of.sun, sun.rise=sun.rise)
}

##' Includes dummy variables for age
##' @param data Data for catch-at-length model internally used by TMB
##' @details This is just a dummy function to include empty variables in age-at-length model.
##' @export
includeDummyAge = function(data){

  data$age = c(1,2,3)
  data$ageNotTruncated = c(1,2,3)
  data$length = c(1,2,3)
  data$readability = c(1,2,3)
  data$ageRange = c(1,2,3)
  data$idx1 = c(1,2,3)
  data$idx2 = c(1,2,3)
  data$spdeMatricesST_alk = data$spdeMatricesS
  data$A_alk_list = data$A_list
  data$Apred_alk = data$ApredS
 # data$lengthGroups = c(1,2,3)
  data$dL = c(1,2,3)
  data$rwBeta0_alk = 0
  data$maxAge = 99
  data$minAge = 0
  data$pcPriorsALKRange = c(200,0.1)
  data$pcPriorsALKSD = c(1,0.1)
  data$usePCpriorsALK = 0
  data$spatialALK = 0
  data$spatioTemporalALK = 0
  data$betaLength = 0

  return(data)
}

