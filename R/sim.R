#' Simulation study
#' @param run Run retruned by fitModel
#' @param nsim Number of simulations
#' @return A runSim object containtin all runs based on simulated data
#' @details
#' @export
simStudy = function(run,nsim = 5){
  simData = list()
  simRuns = list()
  for(i in 1:nsim){
    simData[[i]] = run$obj$simulate()
  }

  for(i in 1:nsim){
    simRuns[[i]] = fitModelSim(run,simData[[i]])
  }

  class(simRuns)  = "stimSim"
  attributes(simRuns)$run = run
  attributes(simRuns)$simData = simData
  return(simRuns)
}

#' Simulation study, run models
#' @param run Run retruned by fitModel
#' @param simData One set of simulated data
#' @return Returns a STIM object fitted with the simulated data
#' @details
#' @export
fitModelSim<-function(run,simData){
  data = run$data
  map = run$map

  data$obsVector = simData$obsVector
  data$age =simData$age
  data$readability = rep(1,length(data$readability)) #NB
  par = run$par
  if(conf_l$applyALK ==1){
    if(conf_l$rwBeta0==1){
      if(data$rwBeta0_alk==1){
        obj <- MakeADFun(data, par, random=c("xS_alk","xST_alk","xS","xST","betaDepth", "nugget","beta0","beta0_alk"),profile = c("betaSun","betaLength_alk"), DLL="spatioTemporalIndices",map = map)
      }else{
        obj <- MakeADFun(data, par, random=c("xS_alk","xST_alk","xS","xST","betaDepth", "nugget","beta0"),profile = c("betaSun","betaLength_alk","beta0_alk"), DLL="spatioTemporalIndices",map = map)
      }
    }else{
      if(data$rwBeta0_alk==1){
        obj <- MakeADFun(data, par, random=c("xS_alk","xST_alk","xS","xST","betaDepth", "nugget","beta0_alk"),profile = c("beta0","betaSun","betaLength_alk"), DLL="spatioTemporalIndices",map = map)
      }else{
        obj <- MakeADFun(data, par, random=c("xS_alk","xST_alk","xS","xST","betaDepth", "nugget"),profile = c("beta0","betaSun","beta0_alk","betaLength_alk"), DLL="spatioTemporalIndices",map = map)
      }
    }
  }else{
    if(conf_l$rwBeta0==1){
      if(conf_l$sunAlt[2]>0){
        obj <- MakeADFun(data, par, random=c("xS","xST","betaDepth", "nugget","beta0"),profile = c("betaSun"), DLL="spatioTemporalIndices",map = map)
      }else{
        obj <- MakeADFun(data, par, random=c("xS","xST","betaDepth", "nugget","beta0"), DLL="spatioTemporalIndices",map = map)
      }
    }else{
      obj <- MakeADFun(data, par, random=c("xS","xST","betaDepth", "nugget"),profile = c("beta0","betaSun"), DLL="spatioTemporalIndices",map = map)
    }
  }


  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(trace = 1,iter.max = 1000, eval.max = 1000))
  rep <- sdreport(obj)
  pl = as.list(rep,"Est")
  plSd = as.list(rep,"Std")

  rl = as.list(rep,"Est", report = TRUE)
  rlSd = as.list(rep,"Std", report = TRUE)

  FreeADFun(obj)#Free memory from C-side

  toReturn = list(obj = obj,opt = opt,rep = rep,conf_l = conf_l,confPred = confPred,conf_alk = conf_alk,data = data,map = map,par = par,dat_l = dat_l,dat_alk = dat_alk,
                  pl = pl, plSd = plSd, rl = rl, rlSd = rlSd)
  class(toReturn) = "stim"

  return(toReturn)
}


