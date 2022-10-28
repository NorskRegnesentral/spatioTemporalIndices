## Return AIC, used by stats::AIC()
##'
##' @param object
##' @details
##' @export
logLik.stim<-function(object, ...){
  ret<- -object$opt$objective

  sunCov = 0
  if(object$conf_l$sunAlt[2]==0)sunCov = object$conf_l$sunAlt[1]*2
  if(object$conf_l$sunAlt[2]==1)sunCov = object$conf_l$sunAlt[1]*4

  beta0 = 0
  if(object$conf_l$rwBeta0==0){
    beta0 = beta0+ length(object$conf_l$lengthGroups)*length(object$conf_l$years)
  }

  if(object$conf_l$applyALK!=0){
    if(object$conf_alk$rwBeta0==0){
      beta0 = beta0+ length(object$conf_alk$minAge:object$conf_alk$maxAge )*length(object$conf_l$years)
    }
  }
  attr(ret,"df")<-length(object$opt$par) + beta0 + sunCov #NB, beta_0 included in the inner optimization in TMB
  class(ret)<-"logLik"
  ret
}


##' Print stim object
##' @method print stim
##' @param  x
##' @details Print log-likelihood and the main convergence criteria
##' @export
print.stim<-function(x, ...){
  cat("STIM model: log likelihood is", logLik.stim(x,...),"Convergence", ifelse(0==x$opt$convergence, "OK\n", "failed\n"))
}

##' Simulate from a STIM object
##' @param  run
##' @param  nsim
##' @param  seed
##' @details Simulate observations from STIM as explained in Breiviek et al. (2020).
##' @export
simulate.stim<-function(run, nsim=1, seed=NULL, ...){
  if(!is.null(seed)) set.seed(seed)
  pl <- as.list(run$rep,"Est")
  est <- unlist(pl)
  ret <- replicate(nsim,
                   c(run$data[names(run$data)!="fishObsMatrix"],#all the old data
                     run$obj$simulate()["fishObsMatrix"])#simulated observations
                   , simplify=FALSE)

  ret
}



##' Print parameters
##' @param  run stim object
##' @details
##' @export
partable<-function(run, ...){

  pl = as.list(run$rep,"Est")
  plsd = as.list(run$rep,"Std. Error")

  rho_l =  2/(1 + exp(-2 * pl$tan_rho_l)) - 1
  rho_l_U =  2/(1 + exp(-2 * (pl$tan_rho_l + 1.96*plsd$tan_rho_l))) - 1
  rho_l_L =  2/(1 + exp(-2 * (pl$tan_rho_l - 1.96*plsd$tan_rho_l))) - 1
  rho_t =  2/(1 + exp(-2 * pl$tan_rho_t)) - 1
  rho_t_U =  2/(1 + exp(-2 * (pl$tan_rho_t+ 1.96*plsd$tan_rho_t))) - 1
  rho_t_L =  2/(1 + exp(-2 * (pl$tan_rho_t- 1.96*plsd$tan_rho_t))) - 1

  rho_l[1:2] = rho_l[1:2]^(1/run$conf_l$reduceLength) #Accommodates for reduced length dimension
  rho_l_U[1:2] = rho_l_U[1:2]^(1/run$conf_l$reduceLength)
  rho_l_L[1:2] = rho_l_L[1:2]^(1/run$conf_l$reduceLength)

  rho_l = round(rho_l,2)
  rho_l_U = round(rho_l_U,2)
  rho_l_L = round(rho_l_L,2)
  rho_t = round(rho_t,2)
  rho_t_U = round(rho_t_U,2)
  rho_t_L = round(rho_t_L,2)


  sigma = round(exp(pl$log_sigma),2)
  sigmaU = round(exp(pl$log_sigma + 1.96*plsd$log_sigma),2)
  sigmaL = round(exp(pl$log_sigma - 1.96*plsd$log_sigma),2)

  kappa = round(exp(pl$log_kappa),4)
  kappaU = round(exp(pl$log_kappa +  1.96*plsd$log_kappa),4)
  kappaL = round(exp(pl$log_kappa -  1.96*plsd$log_kappa),4)

  sigmaRW = c(exp(pl$log_sigma_beta0), exp(pl$log_sigma_beta0 + 1.96*plsd$log_sigma_beta0*c(-1,1)))
  sigmaRW_alk = c(exp(pl$log_sigma_beta0_alk), exp(pl$log_sigma_beta0_alk + 1.96*plsd$log_sigma_beta0_alk*c(-1,1)))

  sigma_alkS = c(exp(pl$logSigma_alk)[1],exp(pl$logSigma_alk[1] + 1.96*plsd$logSigma_alk[1]*c(-1,1)))
  sigma_alkST = c(exp(pl$logSigma_alk)[2],exp(pl$logSigma_alk[2] + 1.96*plsd$logSigma_alk[2]*c(-1,1)))

  kappa_alkS = c(exp(pl$logKappa_alk)[1],exp(pl$logKappa_alk[1] + 1.96*plsd$logKappa_alk[1]*c(-1,1)))
  kappa_alkST = c(exp(pl$logKappa_alk)[2],exp(pl$logKappa_alk[2] + 1.96*plsd$logKappa_alk[2]*c(-1,1)))


  sigmaRW = round(sigmaRW,2);
  sigmaRW_alk = round(sigmaRW_alk,2);
  sigma_alkS = round(sigma_alkS,2);
  sigma_alkST = round(sigma_alkST,2);
  kappa_alkS = round(kappa_alkS,4);
  kappa_alkST = round(kappa_alkST,4);


  if(run$conf_l$spatial==0){
    sigma[1] = "-";sigmaL[1] = "-";sigmaU[1] = "-"
    kappa[1] = "-";kappaL[1] = "-";kappaU[1] = "-"
    rho_l[1] = "-";rho_l_L[1] = "-";rho_l_U[1] = "-"
  }
  if(run$conf_l$spatioTemporal==0){
    sigma[2] = "-";sigmaL[2] = "-";sigmaU[2] = "-"
    kappa[2] = "-";kappaL[2] = "-";kappaU[2] = "-"
    rho_l[2] = "-";rho_l_L[2] = "-";rho_l_U[2] = "-"
    rho_t = "-";rho_t_L = "-";rho_t_U = "-"
  }
  if(run$conf_alk$spatial==0){
    sigma_alkS = rep("-",3);
    kappa_alkS = rep("-",3);
  }
  if(run$conf_alk$spatioTemporal==0){
    sigma_alkST = rep("-",3);
    kappa_alkST = rep("-",3);
  }


  par = matrix(0,15,3)
  par[1,1] = rho_l[1]; par[1,2:3] = c(rho_l_L[1],rho_l_U[1])
  par[2,1] = rho_l[2]; par[2,2:3] = c(rho_l_L[2],rho_l_U[2])
  par[3,1] = rho_l[3]; par[3,2:3] = c(rho_l_L[3],rho_l_U[3])
  par[4,1] = rho_t; par[4,2:3] = c(rho_t_L,rho_t_U)

  par[5,1] = sigma[1]; par[5,2:3] = c(sigmaL[1],sigmaU[1])
  par[6,1] = sigma[2]; par[6,2:3] = c(sigmaL[2],sigmaU[2])
  par[7,1] = sigma[3]; par[7,2:3] = c(sigmaL[3],sigmaU[3])
  par[8,1] = kappa[1]; par[8,2:3] = c(kappaL[1],kappaU[1])
  par[9,1] = kappa[2]; par[9,2:3] = c(kappaL[2],kappaU[2])

  par[10,] = sigmaRW
  par[11,] = sigmaRW_alk
  par[12,] = sigma_alkS
  par[13,] = sigma_alkST
  par[14,] = kappa_alkS
  par[15,] = kappa_alkST

  rownames(par) = c("rho_l_S","rho_l_ST","rho_l_nugg","rho_t","sigmaS", "sigmaST","sigmaNugget", "kappaS", "kappaST", "sigma_beta0_RW",
                    "sigmaRW_alk","sigmaS_alk","sigmaST_alk", "kappaS_alk", "kappaST_alk")
  colnames(par) = c("MLE", "0.025Percentile", "0.975Percentile")

  return(par)
}


##' Find sun heigth and sun rise. Thanks to Espen Johnsen who provided the code.
##' @param
##' @details
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
##' @param
##' @details
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

  return(data)
}

