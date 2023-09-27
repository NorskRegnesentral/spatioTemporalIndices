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
    simRuns[[i]] = tryCatch(
      {
        fitModelSim(run,simData[[i]])
      },
      error = function(e){
        NA
      }
    )
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
  data$readability = rep(1,length(data$readability)) #NB!!! Do not simulate readability
  par = run$par
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


