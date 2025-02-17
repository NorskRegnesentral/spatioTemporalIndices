suppressMessages(library(spatioTemporalIndices))

dat_l = readRDS("NDSKpandLength/NDSKpand2018-2020_length.rds")
dat_l$station <- paste(format(dat_l$startdatetime,"%Y"),dat_l$station, sep = "_")


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
run = fitModel(dat_l,conf_l,confPred,ignore.parm.uncertainty = TRUE,silent = TRUE)




objective = run$opt$objective
rlIndex = run$rl$logLengthIndex
rlIndexSd = run$rlSd$logLengthIndex
par = run$opt$par

resultsOut = list(objective = objective,
                  rlIndex = rlIndex,
                  rlIndexSd = rlIndexSd,
                  par = par)

#Verify all indices-at-length and parameters are as expected
load("NDSKpandLength/resultsExp.RData")
expect_equal(resultsOut$objective, resultsExp$objective,tolerance = 1e-4)
expect_equal(resultsOut$rlIndex, resultsExp$rlIndex,tolerance = 1e-3)
expect_equal(resultsOut$rlIndexSd, resultsExp$rlIndexSd,tolerance = 1e-3)
expect_equal(resultsOut$par, resultsExp$par,tolerance = 1e-3)

#Verify that the two-stage approach leads to the same objective
runTwoStage = fitModel(dat_l,conf_l,confPred,twoStage = TRUE,ignore.parm.uncertainty = TRUE,silent = TRUE)
expect_equal(runTwoStage$opt$objective, resultsExp$objective,tolerance = 1e-4)

#Verify skip years
conf_l$skipYears = 2019
run_skip = fitModel(dat_l,conf_l,confPred,ignore.parm.uncertainty = TRUE,silent = TRUE)
expect_equal(run_skip$opt$objective, resultsExp$objectiveSkip,tolerance = 1e-4)


#Verify that save covaraince structures are as expected
write_covariance_matrices(run,"NDSKpandLength/yearlyCov.rds")
cov = readRDS("NDSKpandLength/yearlyCov.rds")
covExp = readRDS("NDSKpandLength/yearlyCovExp.rds")
expect_equal(cov,
             covExp,tolerance = 1e-3)
#TODO: check "write_indices_ICES_format" for indices-at-length.


if(FALSE){
  resultsExp = list(objective = objective,
                    rlIndex = rlIndex,
                    rlIndexSd = rlIndexSd,
                    par = par,
                    objectiveSkip = run_skip$opt$objective)
  save(resultsExp,file = "NDSKpandLength/resultsExp.RData")
  write_covariance_matrices(run,"NDSKpandLength/yearlyCovExp.rds")

}

