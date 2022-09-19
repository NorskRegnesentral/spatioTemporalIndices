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

  attr(ret,"df")<-length(object$opt$par) + length(object$conf_l$lengthGroups)*length(object$conf_l$years) + sunCov #NB, beta_0 included in the inner optimization in TMB
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

  betaSun=pl$betaSun[1:2]
  betaSun_L=betaSun - 1.96*plsd$betaSun[1:2]
  betaSun_U=betaSun + 1.96*plsd$betaSun[1:2]

  rho_l =  2/(1 + exp(-2 * pl$tan_rho_l)) - 1
  rho_l_U =  2/(1 + exp(-2 * (pl$tan_rho_l + 1.96*plsd$tan_rho_l))) - 1
  rho_l_L =  2/(1 + exp(-2 * (pl$tan_rho_l - 1.96*plsd$tan_rho_l))) - 1
  rho_t =  2/(1 + exp(-2 * pl$tan_rho_t)) - 1
  rho_t_U =  2/(1 + exp(-2 * (pl$tan_rho_t+ 1.96*plsd$tan_rho_t))) - 1
  rho_t_L =  2/(1 + exp(-2 * (pl$tan_rho_t- 1.96*plsd$tan_rho_t))) - 1

  rho_l[1:2] = rho_l[1:2]^(1/run$conf_l$reduceLength) #Accommodates for reduced length dimension
  rho_l_U[1:2] = rho_l_U[1:2]^(1/run$conf_l$reduceLength)
  rho_l_L[1:2] = rho_l_L[1:2]^(1/run$conf_l$reduceLength)

  sigma = exp(pl$log_sigma)
  sigmaU = exp(pl$log_sigma + 1.96*plsd$log_sigma)
  sigmaL = exp(pl$log_sigma - 1.96*plsd$log_sigma)

  kappa = exp(pl$log_kappa)
  kappaU = exp(pl$log_kappa +  1.96*plsd$log_kappa)
  kappaL = exp(pl$log_kappa -  1.96*plsd$log_kappa)

  nu = exp(pl$logSize)
  nuU = exp(pl$logSize+  1.96*plsd$logSize)
  nuL = exp(pl$logSize -  1.96*plsd$logSize)

  lambda = 1/exp(pl$log_lambda)[1]
  lambdaU = 1/exp(pl$log_lambda +  1.96*plsd$log_lambda)[1]
  lambdaL = 1/exp(pl$log_lambda -  1.96*plsd$log_lambda)[1]

  sigmaRW = c(exp(pl$log_sigma_beta0), exp(pl$log_sigma_beta0 + 1.96*plsd$log_sigma_beta0*c(-1,1)))

  sigmaRW_alk = c(exp(pl$log_sigma_beta0_alk), exp(pl$log_sigma_beta0_alk + 1.96*plsd$log_sigma_beta0_alk*c(-1,1)))
  sigma_alk = c(exp(pl$logSigma_alk), exp(pl$logSigma_alk + 1.96*plsd$logSigma_alk*c(-1,1)))
  kappa_alk = c(exp(pl$logKappa_alk), exp(pl$logKappa_alk + 1.96*plsd$logKappa_alk*c(-1,1)))


  par = matrix(0,16,3)
  par[1,1] = rho_l[1]; par[1,2:3] = c(rho_l_L[1],rho_l_U[1])
  par[2,1] = rho_l[2]; par[2,2:3] = c(rho_l_L[2],rho_l_U[2])
  par[3,1] = rho_l[3]; par[3,2:3] = c(rho_l_L[3],rho_l_U[3])
  par[4,1] = rho_t; par[4,2:3] = c(rho_t_L,rho_t_U)

  par[5,1] = sigma[1]; par[5,2:3] = c(sigmaL[1],sigmaU[1])
  par[6,1] = sigma[2]; par[6,2:3] = c(sigmaL[2],sigmaU[2])
  par[7,1] = sigma[3]; par[7,2:3] = c(sigmaL[3],sigmaU[3])

  par[8,1] = kappa[1]; par[8,2:3] = c(kappaL[1],kappaU[1])
  par[9,1] = kappa[2]; par[9,2:3] = c(kappaL[2],kappaU[2])
  par[10,1] = lambda; par[10,2:3] = c(lambdaL,lambdaU)

  par[11,1] = betaSun[1]; par[11,2:3] = c(betaSun_L[1],betaSun_U[1])
  par[12,1] = betaSun[2]; par[12,2:3] = c(betaSun_L[2],betaSun_U[2])

  par[13,] = sigmaRW;
  par[14,] = sigmaRW_alk;
  par[15,] = sigma_alk;
  par[16,] = kappa_alk;

  rownames(par) = c("rho_l_S","rho_l_ST","rho_l_nugg","rho_t","sigmaS", "sigmaST","sigmaNugget", "kappaS", "kappaST","lambda", "betaSun1", "betaSun2", "sigma_beta0_RW",
                    "sigmaRW_alk","sigmaST_alk", "kappaST_alk")
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
  data$lengthGroups = c(1,2,3)
  data$dL = c(1,2,3)
  data$rwBeta0_alk = 0
  data$maxAge = 99
  data$minAge = 0


  return(data)
}

