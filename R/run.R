#' Run model
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats nlminb
#' @param dat_l Data frame with length data
#' @param conf Configurations for length part of model
#' @param confPred Configurations for predictions
#' @param dat_alk Data frame with age data
#' @param conf_alk Configurations for ALK model.
#' @param ... Parameters sendt to sdreport.
#' @useDynLib spatioTemporalIndices
#' @return A fitted stim object
#' @details
#' @export
fitModel<-function(dat_l,conf_l,confPred,dat_alk = NULL, conf_alk = NULL,parSet = NULL,runModel = TRUE, mapSet = NULL,...){

  tryCatch({
      setMKLthreads(1) #not profiting much by using more cores
    },
    error=function(cond) {
      message("MKL library not used, couputation time may be reduced by using MLK library")
    }
  )

  print("Set up length data")

  data_l = setupData(dat_l,conf_l,confPred)
  par = setPar(data_l,conf_l)
  map = setMap(par,conf_l)

  if(conf_l$applyALK==1){
    print("Set up age data")

    if(is.null(dat_alk))stop("Needs age data")
    #Set up data
    conf_alk$meshSimilar = TRUE #Apply same mesh as used for length
    conf_alk$zone = conf_l$zone #Apply same UTM center
    conf_alk$years= conf_l$years #Apply same year range
    data_alk = setUpData_alk(dat_alk,conf_alk,conf_l)

    #Define parameters
    par_alk = defpar_alk(data_alk,conf_alk)
    map_alk = setMap_alk(conf_alk,par_alk)

    #Combine length and ALK input data
    joint = combineLengthALK(data_l,par,map,data_alk,par_alk,map_alk,conf_l,confPred)

    data = joint$dat_joint
    par = joint$par_joint
    map = joint$map_joint
  } else {
    data=data_l
  }
  print("Start inference")


  if(!is.null(parSet)){par = parSet}
  if(!is.null(mapSet)){map = mapSet}

  random = c("xS","xST", "nugget","betaDepth")
  profile = NULL
  if(conf_l$sunAlt[2]!=0){
    profile = c("betaSun")
  }
  if(conf_l$rwBeta0==1){
    random=c(random,"beta0")
  }else{
    profile = c(profile,"beta0")
  }
  if(conf_l$applyALK ==1){
    random = c(random, "xS_alk","xST_alk")
    profile = c(profile,"betaLength_alk")
    if(data$rwBeta0_alk==1){
      random=c(random,"beta0_alk")
    }else{
      profile = c(profile,"beta0_alk")
    }
    if(conf_l$rwBeta0==0& data$rwBeta0_alk==0 &conf_l$nugget==0 & conf_l$spatial==0 & conf_l$spatioTemporal==0 & conf_l$splineDepth[2]==0 & conf_l$sunAlt[2]==0){
      profile = NULL#Need one parameter that is not profiled
    }
  }else{
    if(conf_l$rwBeta0==0 &conf_l$nugget==0 & conf_l$spatial==0 & conf_l$spatioTemporal==0 & conf_l$splineDepth[2]==0 & conf_l$sunAlt[2]==0){
      profile = NULL#Need one parameter that is not profiled
    }
  }

  obj <- MakeADFun(data, par, random=random,profile = profile, DLL="spatioTemporalIndices",map = map)

  if(runModel){
    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  control = list(trace = 1,iter.max = 1000, eval.max = 1000))
    rep <- sdreport(obj,...)
    pl = as.list(rep,"Est")
    plSd = as.list(rep,"Std")

    rl = as.list(rep,"Est", report = TRUE)
    rlSd = as.list(rep,"Std", report = TRUE)

    FreeADFun(obj)#Free memory from C-side
    toReturn = list(obj = obj,opt = opt,rep = rep,conf_l = conf_l,confPred = confPred,conf_alk = conf_alk,data = data,map = map,par = par,dat_l = dat_l,dat_alk = dat_alk,
                    pl = pl, plSd = plSd, rl = rl, rlSd = rlSd)
  }else{
    toReturn = list(obj = obj,conf_l = conf_l,confPred = confPred,conf_alk = conf_alk,data = data,map = map,par = par,dat_l = dat_l,dat_alk = dat_alk)
  }

  class(toReturn) = "stim"
  return(toReturn)
}




#' jit
#' @param run The result of running fitModel
#' @param njit Number of jitter runs
#' @param ncores Number of cores to use
#' @param sd Standard deviation of noise to start values
#' @details
#' @return
#' @export
#'
jit<-function(run,njit,ncores = 1,sd = 0.1){

  #Construct jitter
  data_l = setupData(run$dat_l,run$conf_l,run$confPred)
  par = setPar(data_l,run$conf_l)
  if(run$conf_l$applyALK==1){
    dat_alk = setUpData_alk(run$dat_alk,run$conf_alk, run$conf_l)
    par_alk = defpar_alk(dat_alk,run$conf_alk)
    par = c(par,par_alk)
  }
  parOriginal = par

  pOriginal <- unlist(parOriginal)
  par <- lapply(1:njit, function(i)relist(pOriginal+rnorm(length(pOriginal),sd=sd), run$par))

  #Set to spatio-temporal parameters to arrays as needed in the current implementation


  for(i in 1:length(par)){
    scaleS = 1/((4*3.14159265)*exp(par[[i]]$logKappa)*exp(par[[i]]$logKappa))
    scaleS_alk = 1/((4*3.14159265)*exp(par[[i]]$logKappa_alk)*exp(par[[i]]$logKappa_alk))

    par[[i]]$xS = array(par[[i]]$xS, dim = dim(parOriginal$xS))* sqrt(scaleS[1])
    par[[i]]$xST = array(par[[i]]$xST, dim = dim(parOriginal$xST))* sqrt(scaleS[2])
    par[[i]]$xST_alk = array(par[[i]]$xST_alk, dim = dim(parOriginal$xST_alk))* sqrt(scaleS_alk[1])
    par[[i]]$xS_alk = array(par[[i]]$xS_alk, dim = dim(parOriginal$xS_alk))* sqrt(scaleS_alk[2])
  }

  #Set those who are mapped as NA to its initial value
  par = lapply(par, function(p){
    for(i in 1:length(run$map)){
      j = which(names(p) == names(run$map)[i])
      p[[j]][which(is.na(run$map[[i]]))] = parOriginal[[j]][which(is.na(run$map[[i]]))]
    }
    p
  })

  if(run$conf_l$applyALK==1){
    for(i in 1:njit){
      tmp = sum(exp(par[[i]]$betaLength_alk[-1]))
      for(a in 2:length(parOriginal$betaLength_alk)){
        #Very often run into  convergence issues when these are not decreasing. Reasonable they are decreasing because fish grows with age.
        par[[i]]$betaLength_alk[a] = par[[i]]$betaLength_alk[a-1] - exp(par[[i]]$betaLength_alk[a])/tmp
      }
    }
  }

  if(ncores ==1){
    runs=lapply(par,function(f) fitModel(run$dat_l,conf_l = run$conf_l,confPred = run$confPred, dat_alk = run$dat_alk,conf_alk = run$conf_alk,parPrior = f))
  }else{
    stop("not working with several cores, TODO")
    cl <-makeCluster(min(njit,ncores),outfile = "")
    clusterExport(cl,varlist=c("run","par"),envir=environment())
    clusterEvalQ(cl, library("spatioTemporalIndices"))
    runs=parLapply(cl,par,function(f) fitModel(run$dat_l,conf_l = run$conf_l,confPred = run$confPred, dat_alk = run$dat_alk,conf_alk = run$conf_alk,parPrior = f))
    stopCluster(cl)
  }

  attributes(runs)$runOriginal = run

  p <- lapply(runs, function(f){
    unlist(as.list(f$rep,"Est"))
  })
  plList <- lapply(runs, function(f){
    as.list(f$rep,"Est")
  })

  #Max diff of indices and log-likelihood
  rl = as.list(run$rep, "Est", report = TRUE)

  if(run$conf_l$applyALK==1){
    maxDiffIndices = lapply(runs,function(f){
      rlJ = as.list(f$rep, "Est", report = TRUE)
      max(abs(rlJ$logAgeIndex- rl$logAgeIndex))
    })
  }else{
    maxDiffIndices = lapply(runs,function(f){
      rlJ = as.list(f$rep, "Est", report = TRUE)
      max(abs(rlJ$logLengthIndex- rl$logLengthIndex))
    })
  }



  maxDiffLogLik = lapply(runs,function(f){
    max(abs(f$opt$objective- run$opt$objective))
  })

  maxDiffSun = lapply(runs,function(f){
    rlJ = as.list(f$rep, "Est", report = TRUE)
    max(abs(rlJ$fourierReportLow- rl$fourierReportLow))
  })

  maxDiffDepth= lapply(runs,function(f){
    rlJ = as.list(f$rep, "Est", report = TRUE)
    max(abs(rlJ$depthReport1- rl$depthReport1))
  })


  pl = as.list(run$rep,"Est")
  scaleS = 1/((4*3.14159265)*exp(pl$logKappa)*exp(pl$logKappa))
  pl$xS = pl$xS/scaleS[1]* exp(pl$log_sigma[1])
  pl$xST = pl$xST/scaleS[2]* exp(pl$log_sigma[2])
  pl$nugget = pl$nugget* exp(pl$log_sigma[3])

  scaleS_alk = 1/((4*3.14159265)*exp(pl$logKappa_alk)*exp(pl$logKappa_alk))
  pl$xS_alk = pl$xS_alk/scaleS_alk[1]* exp(pl$logSigma_alk[1])
  pl$xST_alk = pl$xST_alk/scaleS_alk[2]* exp(pl$logSigma_alk[2])

  diffPar = lapply(runs, function(f) {
    plJit = as.list(f$rep, "Est")
    scaleS = 1/((4*3.14159265)*exp(plJit$logKappa)*exp(plJit$logKappa))
    plJit$xS = plJit$xS/scaleS[1]* exp(plJit$log_sigma[1])
    plJit$xST = plJit$xST/scaleS[2]* exp(plJit$log_sigma[2])
    plJit$nugget = plJit$nugget* exp(plJit$log_sigma[3])

    scaleS_alk = 1/((4*3.14159265)*exp(plJit$logKappa_alk)*exp(plJit$logKappa_alk))
    plJit$xS_alk = plJit$xS_alk/scaleS_alk[1]* exp(plJit$logSigma_alk[1])
    plJit$xST_alk = plJit$xST_alk/scaleS_alk[2]* exp(plJit$logSigma_alk[2])

    mpd = mapply(function(g,h){
      max(abs(g-h))
    } , plJit, pl)
  })

  maxMat = t(matrix(unlist(diffPar), nrow = length(diffPar[[1]])))
  maxVec = apply(maxMat,2, max)

  maxVecAll = c(maxVec, max(unlist(maxDiffIndices)),
                    max(unlist(maxDiffLogLik)),
                    max(unlist(maxDiffSun)),
                    max(unlist(maxDiffDepth)))
  names(maxVecAll) = c(names(pl), "Index", "logLik", "Sun","Depth")
  colnames(maxMat) = names(pl)

  return(list(runs = runs, maxVecAll = maxVecAll,maxMat = maxMat))
}




#' retroSTIM
#' @param run The result of running fitModel
#' @param nyears The number of years to remove sequentially
#' @param years Default NULL, will overwrite nyears and select which years to remove sequentially
#' @param ncores Number of cores to use.
#' @details
#' @return
#' @export
#'
retroSTIM = function(run,nyears,years = NULL,ncores = 1){
  if(is.null(years)){
    if(ncores>1){
      cl <- makeCluster(ncores) #set up nodes
      on.exit(stopCluster(cl)) #shut it down
      clusterExport(cl, varlist=c("run","confPred"), envir=environment()) #TODO: by some reason confPred is needed here. Look more into that.
      lib.ver <- dirname(path.package("spatioTemporalIndices"))
      clusterExport(cl, varlist="lib.ver", envir=environment())
      clusterEvalQ(cl, {library(spatioTemporalIndices ,lib.loc=lib.ver)
        library(spatioTemporalALK ,lib.loc=lib.ver)})
      ret = parLapply(cl,1:nyears, function(y,dat_l = run$dat_l,dat_alk = run$dat_alk,
                                            conf_l = run$conf_l,conf_alk = run$conf_alk,confPred= run$confPred){
        maxYear = max(conf_l$years)
        yy_l = as.numeric(format(dat_l$startdatetime,"%Y"))
        yy_alk = as.numeric(format(dat_alk$startdatetime,"%Y"))
        dat_l = dat_l[which(yy_l<= (maxYear-y)),]
        dat_alk = dat_alk[which(yy_alk<= (maxYear-y)),]
        conf_l$years = conf_l$years[conf_l$years <= ((maxYear-y))]
        runTmp = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk)
        runTmp})
    }else{
      ret = lapply(1:nyears, function(y,dat_l = run$dat_l,dat_alk = run$dat_alk,
                                      conf_l = run$conf_l,conf_alk = run$conf_alk,confPred= run$confPred){
        maxYear = max(conf_l$years)
        yy_l = as.numeric(format(dat_l$startdatetime,"%Y"))
        yy_alk = as.numeric(format(dat_alk$startdatetime,"%Y"))
        dat_l = dat_l[which(yy_l<= (maxYear-y)),]
        dat_alk = dat_alk[which(yy_alk<= (maxYear-y)),]
        conf_l$years = conf_l$years[conf_l$years <= ((maxYear-y))]
        runTmp = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk)
        runTmp})
    }
  }else{
    if(ncores >1){
      cl <- makeCluster(ncores) #set up nodes
      on.exit(stopCluster(cl)) #shut it down
      clusterExport(cl, varlist=c("run","confPred"), envir=environment())
      lib.ver <- dirname(path.package("spatioTemporalIndices"))
      clusterExport(cl, varlist="lib.ver", envir=environment())
      clusterEvalQ(cl, {library(spatioTemporalIndices ,lib.loc=lib.ver)
        library(spatioTemporalALK ,lib.loc=lib.ver)})
      ret = parLapply(cl, years, function(year,dat_l = run$dat_l,dat_alk = run$dat_alk,
                                          conf_l = run$conf_l,conf_alk = run$conf_alk,confPred= run$confPred){
        yy_l = as.numeric(format(dat_l$startdatetime,"%Y"))
        yy_alk = as.numeric(format(dat_alk$startdatetime,"%Y"))
        dat_l = dat_l[which(yy_l< year),]
        dat_alk = dat_alk[which(yy_alk< year),]
        conf_l$years = conf_l$years[conf_l$years <year]
        runTmp = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk)
        runTmp})
    }else{
      ret = lapply(years, function(year,dat_l = run$dat_l,dat_alk = run$dat_alk,
                                   conf_l = run$conf_l,conf_alk = run$conf_alk,confPred= run$confPred){
        yy_l = as.numeric(format(dat_l$startdatetime,"%Y"))
        yy_alk = as.numeric(format(dat_alk$startdatetime,"%Y"))
        dat_l = dat_l[which(yy_l)< year,]
        dat_alk = dat_alk[which(yy_alk)< year,]
        conf_l$years = conf_l$years[conf_l$years <year]
        runTmp = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk)
        runTmp})
    }
  }
  attributes(ret)$run = run
  class(ret) = "stimRetro"
  ret
}


