##' Posterior mean. Only do the bias correction for spatial effect for catch-at-length
##' @param run The fitted model
##' @details Return posterior mean of the index (using the bias-correction flag in TMB::sdreport for the catch-at-length spatial and spatio-temporal random effects only)
##' @export
posteriorMean = function(run){

  if(run$conf_l$applyALK ==0){
    #Remove random effect in ALK to reduce memory consumption
    pl = run$pl
    map = run$map

    if(run$conf_l$rwBeta0==1){#No bias-correction for intercept
      map$beta0 = as.factor(pl$beta0*NA)
    }

    if(run$conf_l$splineDepth[2]!=0){#No bias-correction for splines
      map$betaDepth = as.factor(pl$betaDepth*NA)
    }

    dat_l = run$dat_l
    conf_l = run$conf_l
    confPred = run$confPred
    runMAP = fitModel(dat_l,conf_l, confPred,mapSet = map,parSet = pl,runModel = FALSE)


    #Do bias-correct one year at the time to reduce memory consuption
    split = apply(run$obj$env$ADreportIndex()$logLengthIndex,1,function(f) list(f))
    split = lapply(split,function(f)f[[1]])

    sdrep <- TMB::sdreport(runMAP$obj,bias.correct = TRUE,
                           bias.correct.control = list(split = split,sd = FALSE),
                           skip.delta.method = TRUE,
                           getReportCovariance = FALSE, ignore.parm.uncertainty = TRUE)

  }else{
    #Remove random effect in ALK to reduce memory consumption
    pl = run$pl
    map = run$map

    if(run$conf_l$rwBeta0==1){#No bias-correction for intercept
      map$beta0 = as.factor(pl$beta0*NA)
    }

    if(run$conf_l$splineDepth[2]!=0){#No bias-correction for splines
      map$betaDepth = as.factor(pl$betaDepth*NA)
    }

    if(run$conf_alk$rwBeta0==1){
      map$beta0_alk = as.factor(pl$beta0_alk*NA)
    }
    map$xS_alk = as.factor(pl$xS_alk*NA)
    map$xST_alk = as.factor(pl$xST_alk*NA)

    dat_l = run$dat_l
    dat_alk = run$dat_alk
    conf_l = run$conf_l
    conf_alk = run$conf_alk
    confPred = run$confPred
    runMAPalk = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,mapSet = map,parSet = pl,runModel = FALSE)

    #Do bias-correct one year at the time to reduce memory consuption
    split = apply(run$obj$env$ADreportIndex()$logAgeIndex,1,function(f) list(f))
    split = lapply(split,function(f)f[[1]])

    sdrep <- TMB::sdreport(runMAPalk$obj,bias.correct = TRUE,
                           bias.correct.control = list(split = split,sd = FALSE),
                           skip.delta.method = TRUE,
                           getReportCovariance = FALSE, ignore.parm.uncertainty = TRUE)
  }

  return(sdrep)
}
