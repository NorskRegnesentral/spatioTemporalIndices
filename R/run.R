#' Run model
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats nlminb
#' @param dat_l Data frame with length data
#' @param conf Configurations for length part of model
#' @param confPred Configurations for predictions
#' @param dat_alk Data frame with age data
#' @param conf_alk Configurations for ALK model
#' @param low Lower bounds in the optimization, typically not used.
#' @param up Upper bounds in the optimization, typically not used.
#' @param ... Parameters sendt to sdreport.
#' @useDynLib spatioTemporalIndices
#' @return A fitted stim object
#' @details
#' @export
fitModel<-function(dat_l,conf_l,confPred,dat_alk = NULL, conf_alk = NULL,parPrior = NULL,  low = list(),up = list(),...){

  tryCatch(
    {
      #Will fail if not using MKL (efficient linear algebra library)
      setMKLthreads(1) #not profiting much by using more cores
    },
    error=function(cond) {
      message("MKL library not used, couputation time may be reduced by using MLK library")
      return(NA)
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

  if(!is.null(parPrior)){par = parPrior}


  if(conf_l$applyALK ==1){
    if(conf_l$rwBeta0==1){
      if(data$rwBeta0_alk==1){
        obj <- MakeADFun(data, par, random=c("xST_alk","xS","xST","betaDepth", "nugget","nuggetIndex","beta0","beta0_alk"),profile = c("betaSun","betaLength_alk"), DLL="spatioTemporalIndices",map = map)
      }else{
        obj <- MakeADFun(data, par, random=c("xST_alk","xS","xST","betaDepth", "nugget","nuggetIndex","beta0"),profile = c("betaSun","betaLength_alk","beta0_alk"), DLL="spatioTemporalIndices",map = map)
      }
    }else{
      if(data$rwBeta0_alk==1){
        obj <- MakeADFun(data, par, random=c("xST_alk","xS","xST","betaDepth", "nugget","nuggetIndex","beta0_alk"),profile = c("beta0","betaSun","betaLength_alk"), DLL="spatioTemporalIndices",map = map)
      }else{
        obj <- MakeADFun(data, par, random=c("xST_alk","xS","xST","betaDepth", "nugget","nuggetIndex"),profile = c("beta0","betaSun","beta0_alk","betaLength_alk"), DLL="spatioTemporalIndices",map = map)
      }
    }
  }else{
    if(conf_l$rwBeta0==1){
      if(conf_l$sunAlt[1]==1){
        obj <- MakeADFun(data, par, random=c("xS","xST","betaDepth", "nugget","nuggetIndex","beta0"),profile = c("betaSun"), DLL="spatioTemporalIndices",map = map)
      }else{#Profile needs not to be mapped to only constants, TODO: make this part neater
        obj <- MakeADFun(data, par, random=c("xS","xST","betaDepth", "nugget","nuggetIndex","beta0"), DLL="spatioTemporalIndices",map = map)
      }
    }else{
      obj <- MakeADFun(data, par, random=c("xS","xST","betaDepth", "nugget","nuggetIndex"),profile = c("beta0","betaSun"), DLL="spatioTemporalIndices",map = map)
    }
  }


  lower<-rep(-Inf,length(obj$par))
  upper<-rep(Inf,length(obj$par))
  for(nn in names(low)) lower[names(obj$par)==nn]=low[[nn]]
  for(nn in names(up)) upper[names(obj$par)==nn]=up[[nn]]

  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(trace = 1,iter.max = 1000, eval.max = 1000),
                lower = lower, upper = upper)

  rep <- sdreport(obj,...)
  pl = as.list(rep,"Est")
  plSd = as.list(rep,"Std")

  rl = as.list(rep,"Est", report = TRUE)
  rlSd = as.list(rep,"Std", report = TRUE)


  toReturn = list(obj = obj,opt = opt,rep = rep,conf_l = conf_l,confPred = confPred,conf_alk = conf_alk,data = data,map = map,par = par,dat_l = dat_l,dat_alk = dat_alk,
                  pl = pl, plSd = plSd, rl = rl, rlSd = rlSd)
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
jit<-function(run,njit,ncores = 1,sd = 0.2){

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
    par[[i]]$xS = array(par[[i]]$xS, dim = dim(parOriginal$xS))
    par[[i]]$xST = array(par[[i]]$xST, dim = dim(parOriginal$xST))
    par[[i]]$xST_alk = array(par[[i]]$xST_alk, dim = dim(parOriginal$xST_alk))
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
    warning("Jitter somethimes fails because of betaLength_alk")
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

  #Max diff of indices and log-likelihood
  rl = as.list(run$rep, "Est", report = TRUE)

  maxDiffIndices = lapply(runs,function(f){
    rlJ = as.list(f$rep, "Est", report = TRUE)
    max(abs(rlJ$logLengthIndex- rl$logLengthIndex))
  })

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
  pl$xS = pl$xS* exp(pl$log_sigma[1])
  pl$xST = pl$xST* exp(pl$log_sigma[2])
  pl$nugget = pl$nugget* exp(pl$log_sigma[3])


  diffPar = lapply(runs, function(f) {
    plJit = as.list(f$rep, "Est")
    plJit$xS = plJit$xS* exp(plJit$log_sigma[1])
    plJit$xST = plJit$xST* exp(plJit$log_sigma[2])
    plJit$nugget = plJit$nugget* exp(plJit$log_sigma[3])
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


  return(list(runs = runs, maxVecAll = maxVecAll,maxMat = maxMat))
}


