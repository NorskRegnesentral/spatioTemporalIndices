##' RTMB objective function of STIM
##'
##' Joint negative log-likelihood of the catch-at-length model and (if \code{dat$applyALK==1}) the age-at-length model,
##' including prediction and reporting of the indices. This is the R/RTMB translation of the TMB C++ template
##' (src/spatioTemporalIndices.cpp and inst/include/*.hpp) in the package spatioTemporalIndices.
##'
##' @param par Parameter list, see \code{\link{setPar}} and \code{\link{defpar_alk}}
##' @param dat Data list returned by \code{\link{setupData}} (combined with \code{\link{combineLengthALK}} if age data is used)
##' @return The joint negative log-likelihood
##' @details Used as \code{RTMB::MakeADFun(function(p) nllSTIM(p, dat), par, ...)}, see \code{\link{makeObj}}.
##' @import RTMB
##' @export
nllSTIM <- function(par, dat){

  nll <- nllIndex(par, dat)

  ALKres <- NULL
  if(dat$applyALK==1){
    ALKres <- nllALK(par, dat)
    nll <- nll + ALKres$nll
  }

  indexPred(par, dat)

  if(dat$applyALK==1){
    ALK_hauls <- ALKhauls(par, dat) #ALKs at haul locations
    RTMB::REPORT(ALK_hauls)
  }

  return(nll)
}

##' Construct the RTMB model object
##'
##' @param data Data list
##' @param par Parameter list
##' @param random Names of random effects
##' @param profile Names of parameters to profile out
##' @param map Map-list
##' @param silent Passed to RTMB::MakeADFun
##' @return An RTMB (TMB-compatible) model object
##' @export
makeObj <- function(data, par, random, profile, map, silent){
  RTMB::MakeADFun(function(p) nllSTIM(p, data), par, random = random, profile = profile, map = map, silent = silent)
}


##' Free memory used by the AD tapes of a fitted model
##'
##' Replaces \code{TMB::FreeADFun}, which must not be used on RTMB objects (the tapes are then freed again by their
##' finalizers, which crashes R). The references to the tapes are dropped, so that they are freed at the next garbage collection.
##' After this, \code{obj$fn} and \code{obj$gr} can no longer be used, while \code{obj$report()} still works.
##' @param obj Model object returned by \code{\link{makeObj}}
##' @return NULL
##' @export
freeObj <- function(obj){
  obj$env$ADFun <- NULL
  obj$env$ADGrad <- NULL
  if(!is.null(obj$env$spHess)){
    hessEnv <- environment(obj$env$spHess)
    if(exists("ADHess", envir = hessEnv, inherits = FALSE)) assign("ADHess", NULL, envir = hessEnv)
  }
  invisible(NULL)
}


#Helpers----------------------------------------------------------------------

#Sparse matrix product; returns a base matrix when evaluated with numeric parameters (e.g. obj$report()) instead of a Matrix-class object
mm <- function(A, x){
  r <- A %*% x
  if(methods::is(r, "Matrix")) r <- as.matrix(r)
  r
}

#Precision matrix of the SPDE approximation (alpha = 2), see Lindgren et al. (2011); equivalent to R_inla::Q_spde in TMB
Qspde <- function(spde, kappa){
  k2 <- kappa*kappa
  k2*k2*methods::as(spde$M0, "CsparseMatrix") + 2*k2*methods::as(spde$M1, "CsparseMatrix") + methods::as(spde$M2, "CsparseMatrix")
}

#Log densities used with RTMB::dseparable; equivalent to GMRF(), AR1() and an iid GMRF in TMB
dGMRF <- function(Q) function(x) RTMB::dgmrf(x, 0, Q, log = TRUE)
dAR1 <- function(rho) function(x) RTMB::dautoreg(x, phi = rho, log = TRUE)
dIID <- function(x) sum(RTMB::dnorm(x, 0, 1, log = TRUE))

#Contribution of the pc-prior for range and standard deviation
pcPrior <- function(kappa, sigma, pcRange, pcSD){
  d <- 2
  R <- -log(pcRange[2])*pcRange[1]^(d/2)
  S <- -log(pcSD[2])/pcSD[1]
  rhoP <- sqrt(8)/kappa
  log(d/2*R*S*rhoP^(-1-d/2)*exp(-R*rhoP^(-d/2) - S*sigma))
}

#Matrix (nLengthGroups x nLatentLengths) mapping latent length dimension to the length groups (used when reducing length dimension)
lengthMap <- function(dat, nLatent){
  nL <- dat$numberOfLengthGroups
  lgr <- dat$lengthGroupsReduced #zero based
  reduced <- isTRUE(lgr[1]==lgr[2])
  Lw <- base::matrix(0, nL, nLatent)
  for(l in 1:nL){
    if(reduced & (dat$weigthLength[l]<0.9999)){
      Lw[l, lgr[l]+1] <- dat$weigthLength[l]
      Lw[l, lgr[l]+2] <- 1-dat$weigthLength[l]
    }else{
      Lw[l, lgr[l]+1] <- 1
    }
  }
  Lw
}

#Coefficient in convex combination of length dependent effects
convexW <- function(nL) (nL - (1:nL))/(nL - 1)

#Reshape array x (d1 x nY x d3) to a d1 x (nY*d3) matrix; column y + nY*(k-1) holds x[,y,k]
flatten3 <- function(x){
  d <- dim(x)
  dim(x) <- c(d[1], d[2]*d[3])
  x
}

#Columns of flatten3(x) belonging to year y
colsYear <- function(y, nY, d3) y + nY*(0:(d3-1))

#Continuation-ratio age-length key: lp is (n x (nAges-1)) linear predictor, returns n x nAges probabilities
contRatio <- function(lp){
  nA <- ncol(lp) + 1
  P <- vector("list", nA)
  probLess <- 0
  for(a in 1:nA){
    if(a < nA){
      P[[a]] <- 1/(1 + exp(-lp[, a]))*(1 - probLess)
      probLess <- probLess + P[[a]]
    }else{
      P[[a]] <- 1 - probLess
    }
  }
  do.call(cbind, P)
}


#Catch-at-length likelihood (index.hpp)--------------------------------------
nllIndex <- function(par, dat){
  nll <- 0
  nY <- length(dat$nStationsEachYear)
  nL <- dat$numberOfLengthGroups
  nH <- sum(dat$nStationsEachYear)

  if(dat$rwBeta0==1 & nY > 1){
    sigma_beta0 <- exp(par$log_sigma_beta0[1])
    nll <- nll - sum(RTMB::dnorm(par$beta0[-1, , drop = FALSE], par$beta0[-nY, , drop = FALSE], sigma_beta0, log = TRUE))
  }

  sigma <- exp(par$log_sigma)
  kappa <- exp(par$log_kappa)
  rho_t <- 2/(1 + exp(-2*par$tan_rho_t[1])) - 1
  rho_l <- 2/(1 + exp(-2*par$tan_rho_l)) - 1

  Q_s <- Qspde(dat$spdeMatricesS, kappa[1])
  Q_st <- Qspde(dat$spdeMatricesST, kappa[2])
  scaleS <- 1/((4*pi)*kappa[1]*kappa[1])
  scaleST <- 1/((4*pi)*kappa[2]*kappa[2])

  if(dat$spatial==1){
    if(dat$usePCpriors==1){
      nll <- nll - pcPrior(kappa[1], sigma[1], dat$pcPriorsRange, dat$pcPriorsSD)
    }
    nll <- nll - RTMB::dseparable(dGMRF(Q_s), dAR1(rho_l[1]))(par$xS)
  }
  if(dat$spatioTemporal != 0){
    if(dat$usePCpriors==1){
      nll <- nll - pcPrior(kappa[2], sigma[2], dat$pcPriorsRange, dat$pcPriorsSD)
    }
    nll <- nll - RTMB::dseparable(dGMRF(Q_st), dAR1(rho_t), dAR1(rho_l[2]))(par$xST)
  }
  if(dat$useNugget==1){
    nll <- nll - RTMB::dseparable(dAR1(rho_l[3]), dIID)(par$nugget)
  }

  #p-spline for depth
  nSpline <- length(par$betaDepth)/2
  parDepth1 <- par$betaDepth[1:nSpline]
  parDepth2 <- par$betaDepth[nSpline + 1:nSpline]
  lambda <- exp(par$log_lambda)
  if(dat$splineDepth >= 1){
    nll <- nll - (0.5*dat$Sdim*par$log_lambda[1] - 0.5*lambda[1]*sum(parDepth1*mm(dat$S_depth, parDepth1)))
  }
  if(dat$splineDepth == 2){ #Two length dependent splines
    nll <- nll - (0.5*dat$Sdim*par$log_lambda[2] - 0.5*lambda[2]*sum(parDepth2*mm(dat$S_depth, parDepth2)))
  }

  #fourier basis for sun altitude
  nB2 <- dat$nBasisSunAlt*2
  betaSunLow <- par$betaSun[1:nB2]
  betaSunHigh <- par$betaSun[nB2 + 1:nB2]

  #Linear predictor (hauls x length groups)
  w <- convexW(nL)
  lowEffect <- dat$X_sunAlt %*% betaSunLow + dat$X_depth %*% parDepth1
  highEffect <- dat$X_sunAlt %*% betaSunHigh + dat$X_depth %*% parDepth2
  covariates <- lowEffect %*% t(w) + highEffect %*% t(1 - w)

  yearHaul <- rep(1:nY, dat$nStationsEachYear)
  Ymat <- outer(yearHaul, 1:nY, "==")*1
  intercept <- Ymat %*% par$beta0

  nLatent <- dim(par$xS)[2]
  Lw <- lengthMap(dat, nLatent)
  A_S <- do.call(rbind, dat$A_ListS)
  A_ST <- Matrix::bdiag(dat$A_ListST)
  deltaS <- mm(A_S, par$xS) %*% t(Lw)
  xST <- par$xST
  dim(xST) <- c(dim(xST)[1]*dim(xST)[2], dim(xST)[3])
  deltaST <- mm(A_ST, xST) %*% t(Lw)

  eta <- intercept + covariates +
    deltaS*sigma[1]/sqrt(scaleS) +
    deltaST*sigma[2]/sqrt(scaleST) +
    log(dat$dist) +
    t(par$nugget)*sigma[3]

  #Poisson likelihood, only obsModel==2 is implemented
  obs <- t(base::matrix(dat$obsVector, nrow = nL))
  use <- dat$predMatrix == 0
  nll <- nll - sum((obs*eta - exp(eta) - lgamma(obs + 1))[use])

  return(nll)
}


#Age-at-length likelihood (alk.hpp)------------------------------------------
nllALK <- function(par, dat){
  nll <- 0
  nAges <- dat$ageRange[2] - dat$ageRange[1] + 1
  minAge <- dat$ageRange[1]
  nY <- length(dat$nStationsEachYear)
  nYalk <- length(dat$idx1)

  sigma_beta_alk <- exp(par$log_sigma_beta_alk)
  if(dat$rwBeta0_alk==1){
    if(nYalk > 1 & nAges > 2){
      nll <- nll - sum(RTMB::dnorm(par$beta0_alk[2:nYalk, 2:(nAges-1), drop = FALSE],
                                   par$beta0_alk[1:(nYalk-1), 1:(nAges-2), drop = FALSE], sigma_beta_alk[1], log = TRUE))
    }
    nll <- nll - RTMB::dnorm(par$beta0_alk[1, nAges-1], par$beta0_alk[2, nAges-1], sigma_beta_alk[1], log = TRUE) #random walk for oldest age in first year
    nll <- nll - RTMB::dnorm(par$beta0_alk[nYalk, 1], par$beta0_alk[nYalk-1, 1], sigma_beta_alk[1], log = TRUE) #random walk for youngest age in last year
  }

  if(dat$betaLength==1 & nYalk > 1){
    bL <- par$betaLength_alk
    dim(bL) <- c(nAges-1, length(bL)/(nAges-1))
    nll <- nll - sum(RTMB::dnorm(bL[, 2:nYalk, drop = FALSE], bL[, 1:(nYalk-1), drop = FALSE], sigma_beta_alk[2], log = TRUE))
  }

  sigma <- exp(par$logSigma_alk)
  kappa <- exp(par$logKappa_alk)
  rho <- 2/(1 + exp(-par$transRho_alk[1])) - 1

  QS <- Qspde(dat$spdeMatricesST_alk, kappa[1])
  QST <- Qspde(dat$spdeMatricesST_alk, kappa[2])

  if(dat$spatialALK != 0){
    nll <- nll - RTMB::dseparable(dGMRF(QS), dIID)(par$xS_alk)
  }
  if(dat$spatioTemporalALK != 0){
    nll <- nll - RTMB::dseparable(dGMRF(QST), dAR1(rho), dIID)(par$xST_alk)
  }

  if(dat$usePCpriorsALK==1){
    if(dat$spatialALK != 0){
      nll <- nll - pcPrior(kappa[1], sigma[1], dat$pcPriorsALKRange, dat$pcPriorsALKSD)
    }
    if(dat$spatioTemporalALK != 0){
      nll <- nll - pcPrior(kappa[2], sigma[2], dat$pcPriorsALKRange, dat$pcPriorsALKSD)
    }
  }

  scaleS <- 1/((4*pi)*kappa[1]*kappa[1])
  scaleST <- 1/((4*pi)*kappa[2]*kappa[2])

  #Linear predictor (fish x (nAges-1)), observations are sorted by year
  xST_alk <- flatten3(par$xST_alk)
  linPred <- vector("list", nY)
  for(y in 1:nY){
    rows <- (dat$idx1[y]:dat$idx2[y]) + 1
    A <- dat$A_alk_list[[y]]
    bL <- betaLengthYear(par$betaLength_alk, dat$betaLength, y, nAges)
    linPred[[y]] <- base::matrix(1, length(rows), 1) %*% par$beta0_alk[y, , drop = FALSE] +
      base::matrix(dat$length[rows], ncol = 1) %*% bL +
      mm(A, par$xS_alk)/sqrt(scaleS)*sigma[1] +
      mm(A, xST_alk[, colsYear(y, nY, nAges-1), drop = FALSE])/sqrt(scaleST)*sigma[2]
  }
  ALK <- contRatio(do.call(rbind, linPred))

  #Likelihood, accounting for age reading quality
  W <- readabilityWeights(dat, nAges, minAge)
  used <- base::rowSums(W) > 0 #W is data; base version keeps it numeric
  nll <- nll - sum(log(rowSums(ALK*W)[used]))

  RTMB::REPORT(ALK)
  return(list(nll = nll))
}

#Regression coefficients for length (1 x (nAges-1)) in year y
betaLengthYear <- function(betaLength_alk, betaLength, y, nAges){
  if(betaLength==0){
    b <- betaLength_alk[1:(nAges-1)]
  }else{
    b <- betaLength_alk[(y-1)*(nAges-1) + 1:(nAges-1)]
  }
  dim(b) <- c(1, nAges-1)
  b
}

#Matrix (fish x ages) with 1 for the ages each observation is compatible with, see readability in setUpData_alk
readabilityWeights <- function(dat, nAges, minAge){
  nObs <- length(dat$age)
  W <- base::matrix(0, nObs, nAges)
  k <- dat$age - minAge + 1
  for(s in 1:nObs){
    r <- dat$readability[s]
    if(r==1){
      W[s, k[s]] <- 1
    }else if(r==5){
      W[s, k[s]] <- 1
      if(dat$ageNotTruncated[s] < dat$maxAge) W[s, k[s]+1] <- 1
    }else if(r==6){
      W[s, k[s]] <- 1
      if((dat$ageNotTruncated[s] <= dat$maxAge) & (dat$ageNotTruncated[s] >= dat$minAge) & (k[s] > 1)) W[s, k[s]-1] <- 1
    }
  }
  W
}

#Linear predictor of ALK (n x (nAges-1)) for given length, spatial projector and year
alkLinPred <- function(par, dat, Aproj, y, len, nAges, nY){
  sigma <- exp(par$logSigma_alk)
  kappa <- exp(par$logKappa_alk)
  scaleS <- 1/((4*pi)*kappa[1]*kappa[1])
  scaleST <- 1/((4*pi)*kappa[2]*kappa[2])
  bL <- betaLengthYear(par$betaLength_alk, dat$betaLength, y, nAges)
  space <- mm(Aproj, par$xS_alk)/sqrt(scaleS)*sigma[1] +
    mm(Aproj, flatten3(par$xST_alk)[, colsYear(y, nY, nAges-1), drop = FALSE])/sqrt(scaleST)*sigma[2]
  space + base::matrix(1, nrow(Aproj), 1) %*% (par$beta0_alk[y, , drop = FALSE] + len*bL)
}

#ALKs at the haul locations (ALKhauls in alk.hpp). Array: length group x age x haul x year
ALKhauls <- function(par, dat){
  nAges <- dat$ageRange[2] - dat$ageRange[1] + 1
  nY <- length(dat$nStationsEachYear)
  nL <- dat$numberOfLengthGroups
  ALK_obs <- RTMB::AD(array(0, c(nL, nAges, max(dat$nStationsEachYear), nY))) #AD from the start so that "[<-" dispatches to the advector method
  for(y in 1:nY){
    nS <- dat$nStationsEachYear[y]
    for(l in 1:nL){
      lp <- alkLinPred(par, dat, dat$A_ListST[[y]], y, dat$lengthGroups[l] + 0.5*dat$dL[1], nAges, nY)
      ALK_obs[l, , 1:nS, y] <- t(contRatio(lp))
    }
  }
  ALK_obs
}


#Prediction of indices (indexPred.hpp)---------------------------------------
indexPred <- function(par, dat){
  nL <- dat$numberOfLengthGroups
  nY <- length(dat$nStationsEachYear)
  nInt <- length(dat$xInt)
  sigma <- exp(par$log_sigma)
  kappa <- exp(par$log_kappa)
  scaleS <- 1/((4*pi)*kappa[1]*kappa[1])
  scaleST <- 1/((4*pi)*kappa[2]*kappa[2])

  nSpline <- length(par$betaDepth)/2
  parDepth1 <- par$betaDepth[1:nSpline]
  parDepth2 <- par$betaDepth[nSpline + 1:nSpline]
  nB2 <- dat$nBasisSunAlt*2
  betaSunLow <- par$betaSun[1:nB2]
  betaSunHigh <- par$betaSun[nB2 + 1:nB2]

  #Covariate effects at integration points (nInt x nL), sun altitude at its maximum
  w <- convexW(nL)
  sunLow <- (dat$X_sunAltIntegrate %*% betaSunLow)[1, 1]
  sunHigh <- (dat$X_sunAltIntegrate %*% betaSunHigh)[1, 1]
  covariates <- (dat$X_depth_int %*% parDepth1 + sunLow) %*% t(w) + (dat$X_depth_int %*% parDepth2 + sunHigh) %*% t(1 - w)

  nLatent <- dim(par$xS)[2]
  Lw <- lengthMap(dat, nLatent)
  deltaPredS <- mm(dat$ApredS, par$xS) %*% t(Lw)*(sigma[1]/sqrt(scaleS))
  xSTpred <- mm(dat$ApredST, flatten3(par$xST))

  #Strata bookkeeping
  nStrata <- max(dat$idxStrata)
  Sidx <- outer(dat$idxStrata, 1:nStrata, "==")*1
  nIntStrata <- base::colSums(Sidx)
  strataScale <- dat$areas[1:nStrata]/nIntStrata
  totalScale <- sum(dat$areas)/nInt

  lengthIndexDetailed <- RTMB::AD(array(0, c(nY, nL, nInt))) #AD from the start so that "[<-" dispatches to the advector method
  lengthIndex <- RTMB::AD(base::matrix(0, nY, nL))
  lengthIndexStrata <- RTMB::AD(array(0, c(nY, nL, nStrata)))

  applyALK <- dat$applyALK==1
  nAges <- 2 #Dummy-number
  if(applyALK){
    nAges <- dat$ageRange[2] - dat$ageRange[1] + 1
    ALK_int <- RTMB::AD(array(0, c(nL, nAges, nInt, nY)))
  }
  ageIndex <- RTMB::AD(base::matrix(0, nY, nAges))
  ageIndexDetailed <- RTMB::AD(array(0, c(nY, nAges, nInt)))
  ageIndexStrata <- RTMB::AD(array(0, c(nY, nAges, nStrata)))

  for(y in 1:nY){
    deltaPredST <- xSTpred[, colsYear(y, nY, nLatent), drop = FALSE] %*% t(Lw)*(sigma[2]/sqrt(scaleST))
    beta0y <- par$beta0[y, , drop = FALSE]
    D <- exp(base::matrix(1, nInt, 1) %*% beta0y + covariates + deltaPredS + deltaPredST)*dat$trawlWidth #nInt x nL

    lengthIndexDetailed[y, , ] <- t(D)
    lengthIndex[y, ] <- colSums(D)*totalScale
    lengthIndexStrata[y, , ] <- (t(D) %*% Sidx)*(base::matrix(1, nL, 1) %*% t(strataScale))

    if(applyALK){
      ageDetailed <- 0
      for(l in 1:nL){
        lp <- alkLinPred(par, dat, dat$Apred_alk, y, dat$lengthGroups[l] + 0.5*dat$dL[1], nAges, nY)
        P <- contRatio(lp) #nInt x nAges
        ALK_int[l, , , y] <- t(P)
        ageDetailed <- ageDetailed + D[, l]*P
      }
      ageIndexDetailed[y, , ] <- t(ageDetailed)
      ageIndex[y, ] <- colSums(ageDetailed)*totalScale
      ageIndexStrata[y, , ] <- (t(ageDetailed) %*% Sidx)*(base::matrix(1, nAges, 1) %*% t(strataScale))
    }
  }

  logLengthIndex <- log(lengthIndex)
  logLengthIndexStrata <- log(lengthIndexStrata)
  logLengthIndexTotal <- log(rowSums(lengthIndex))

  RTMB::REPORT(logLengthIndex)
  RTMB::REPORT(logLengthIndexStrata)
  RTMB::REPORT(logLengthIndexTotal)
  RTMB::REPORT(lengthIndexDetailed) #For plotting

  if((dat$strataReport==1) & (dat$applyALK==0)){
    RTMB::ADREPORT(logLengthIndexStrata)
  }
  if(dat$applyALK==0){#Do not adreport length if calculating index per length to reduce computation time.
    RTMB::ADREPORT(logLengthIndex)
    if(dat$lowMemory==0){
      RTMB::ADREPORT(logLengthIndexTotal)
    }
  }

  if(applyALK){
    logAgeIndex <- log(ageIndex)
    logAgeIndexStrata <- log(ageIndexStrata)
    logAgeIndexTotal <- log(rowSums(ageIndex))
    RTMB::ADREPORT(logAgeIndex)
    RTMB::REPORT(logAgeIndex)
    if(dat$lowMemory==0){
      RTMB::ADREPORT(logAgeIndexTotal)
    }else{
      RTMB::REPORT(logAgeIndexTotal)
    }
    RTMB::REPORT(ageIndexDetailed) #For plotting
    RTMB::REPORT(logAgeIndexStrata)
    RTMB::REPORT(ALK_int)
    if(dat$strataReport==1){#Optionally sdreport index-at-age in strata, this may result in high memory consumption
      RTMB::ADREPORT(logAgeIndexStrata)
    }
  }

  if(dat$lowMemory==0){
    #Report covariates
    fourierReportLow <- (dat$X_sunAltReport %*% betaSunLow)[, 1]
    RTMB::ADREPORT(fourierReportLow)
    if(dat$sunEffect==2){
      fourierReportHigh <- (dat$X_sunAltReport %*% betaSunHigh)[, 1]
      RTMB::ADREPORT(fourierReportHigh)
    }
    depthReport1 <- (dat$X_depthReport %*% parDepth1)[, 1]
    RTMB::ADREPORT(depthReport1)
    if(dat$splineDepth==2){ #Two length dependent splines
      depthReport2 <- (dat$X_depthReport %*% parDepth2)[, 1]
      RTMB::ADREPORT(depthReport2)
    }
  }
  invisible(NULL)
}
