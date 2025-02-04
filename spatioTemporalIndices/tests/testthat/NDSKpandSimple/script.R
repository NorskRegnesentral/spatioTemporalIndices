suppressMessages(library(spatioTemporalIndices))

dat_l = readRDS("NDSKpandLength/NDSKpand2018-2020_length.rds")

conf_l = defConf(years = 2020:2020,
               dLength = 1, # length intervall in mm
               maxLength = 30, # Numeric = use directly; NULL = use input data to determine
               minLength = 7, # Numeric = use directly; NULL = use input data to determine
               cutoff = 20,
               spatioTemporal = 0,
               nugget = 1,
               spatial = 0,
               sunAlt=c(1,1),
               splineDepth = c(6,1),
               reduceLength = 3,
               rwBeta0 = 1,
               stratasystem = list(dsn="NDSKpandLength/strata", layer = "shrimp_areas_NSSK"),
               minDepth=50,maxDepth=600,
               trawlWidth=11.7,
               applyALK=0,
               strataReport=0)

confPred = defConfPred(conf=conf_l,Depth="DATA",cellsize = 20)

# run model
run = fitModel(dat_l,conf_l,confPred,ignore.parm.uncertainty = TRUE)

conf_l$nugget = 0
run2 = fitModel(dat_l,conf_l,confPred,ignore.parm.uncertainty = TRUE)

resultsOut = list(AIC = AIC(run,run2))

load("NDSKpandSimple/resultsExp.RData")
expect_equal(resultsOut$objectiveExp, resultsExp$objectiveExp,tolerance = 1e-4)



if(FALSE){
  resultsExp = list(AIC = AIC(run,run2))
  save(resultsExp,file = "NDSKpandSimple/resultsExp.RData")
}

