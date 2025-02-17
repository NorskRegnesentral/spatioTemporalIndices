suppressMessages(library(spatioTemporalIndices))

dat_l = readRDS("NEAhadLengthAgePosteriorMean/haddock2018-2020_length_ex_rus_reduced.rds")
dat_alk = readRDS("NEAhadLengthAgePosteriorMean/haddock2018-2020_age_ex_rus_reduced.rds")

conf_l = defConf(years = 2019:2020, # years to use, use all years with data by default
                 maxLength = 60, # Numeric = use directly; NULL = use input data to determine
                 minLength = 20, # Numeric = use directly; NULL = use input data to determine
                 spatioTemporal =2 ,
                 spatial =0,
                 rwBeta0 = 1,
                 sunAlt = c(1,0),
                 splineDepth = c(6,1),
                 cbound = c(18,130),
                 dLength = 5,
                 reduceLength = 3,
                 stratasystem = list(dsn="NEAhadLengthAgePosteriorMean/strata", layer = "Vintertoktet_nye_strata"),
                 minDepth=150,maxDepth=400,
                 applyALK = 1,
                 cutoff =100)


#Define configurations age part
conf_alk = defConf_alk(maxAge = 8,
                       minAge = 3,
                       spatioTemporal = 2,
                       spatial =0,
                       rwBeta0 = 0)

confPred = defConfPred(conf=conf_l,Depth="DATA",cellsize = 400)

# run model
runALK = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE)

sdrepALK_bc = posteriorMean(runALK)

rlALK_bc = as.list(sdrepALK_bc,what = "Est. (bias.correct)",report = TRUE)


# run model
conf_l$applyALK = 0
run = fitModel(dat_l,conf_l, confPred,ignore.parm.uncertainty = TRUE,silent = TRUE)

sdrep_bc = posteriorMean(run)

rl_bc = as.list(sdrep_bc,what = "Est. (bias.correct)",report = TRUE)

resultsOut = list(logAgeIndex = rlALK_bc$logAgeIndex,
                  logLengthIndex = rl_bc$logLengthIndex)

load("NEAhadLengthAgePosteriorMean/resultsExp.RData")
expect_equal(resultsOut$logAgeIndex, resultsExp$logAgeIndexExp,tolerance = 1e-3)
expect_equal(resultsOut$logLengthIndex, resultsExp$logLengthIndexExp,tolerance = 1e-3)

if(FALSE){
  resultsExp = list(logAgeIndexExp = rlALK_bc$logAgeIndex,
                    logLengthIndexExp = rl_bc$logLengthIndex)
  save(resultsExp,file = "NEAhadLengthAgePosteriorMean/resultsExp.RData")

}

