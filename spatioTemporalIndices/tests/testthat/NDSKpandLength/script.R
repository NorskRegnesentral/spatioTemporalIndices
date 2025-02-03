suppressMessages(library(spatioTemporalIndices))

dat_l = readRDS("NDSKpandLength/NDSKpand2018-2020_length.rds")

conf_l = defConf(years = 2018:2020,
               dLength = 1, # length intervall in mm
               maxLength = 30, # Numeric = use directly; NULL = use input data to determine
               minLength = 7, # Numeric = use directly; NULL = use input data to determine
               cutoff = 30,
               spatioTemporal = 1,
               spatial = 0,
               sunAlt=c(0,0),
               reduceLength = 3,
               rwBeta0 = 0,
               stratasystem = list(dsn="NDSKpandLength/strata", layer = "shrimp_areas_NSSK"),
               minDepth=50,maxDepth=600,
               trawlWidth=11.7,
               applyALK=0,
               strataReport=0)

confPred = defConfPred(conf=conf_l,Depth="DATA",cellsize = 20)

# run model
run = fitModel(dat_l,conf_l,confPred,ignore.parm.uncertainty = TRUE)

objectiveExp = run$opt$objective
rlIndex = round(run$rl$logLengthIndex,5)
rlIndexSd = round(run$rlSd$logLengthIndex,5)
par = round(run$opt$par,5)

resultsOut = list(objectiveExp = objectiveExp,
                  rlIndex = rlIndex,
                  rlIndexSd = rlIndexSd,
                  par = par)


load("NDSKpandLength/resultsExp.RData")
expect_equal(resultsOut$objectiveExp, resultsExp$objectiveExp,tolerance = 1e-4)
expect_equal(resultsOut$rlIndex, resultsExp$rlIndex,tolerance = 1e-2)
expect_equal(resultsOut$rlIndexSd, resultsExp$rlIndexSd,tolerance = 1e-2)
expect_equal(resultsOut$par, resultsExp$par,tolerance = 1e-2)

runTwoStage = fitModel(dat_l,conf_l,confPred,twoStage = TRUE,ignore.parm.uncertainty = TRUE)
expect_equal(runTwoStage$opt$objective, resultsExp$objectiveExp,tolerance = 1e-4)


if(FALSE){
  resultsExp = list(objectiveExp = objectiveExp,
                    rlIndex = rlIndex,
                    rlIndexSd = rlIndexSd,
                    par = par)
  save(resultsExp,file = "NDSKpandLength/resultsExp.RData")
}

