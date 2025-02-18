suppressMessages(library(spatioTemporalIndices))

dat_l = readRDS("NEAhadLengthSimple/haddock2018-2020_length_ex_rus_reduced.rds")

conf_l = defConf(years = 2018:2020, # years to use, use all years with data by default
                 maxLength = 65, # Numeric = use directly; NULL = use input data to determine
                 minLength = 20, # Numeric = use directly; NULL = use input data to determine
                 spatioTemporal =0 ,
                 spatial =0,
                 rwBeta0 = 1,
                 sunAlt = c(1,2),
                 splineDepth = c(6,2),
                 cbound = c(18,130),
                 dLength = 5,
                 reduceLength = 3,
                 stratasystem = list(dsn="NEAhadLengthSimple/strata", layer = "Vintertoktet_nye_strata"),
                 minDepth=150,maxDepth=400,
                 applyALK = 0,
                 cutoff =100)


confPred = defConfPred(conf=conf_l,Depth="DATA",cellsize = 50)


###Length dependent covariates:
runLenghtDepCov = fitModel(dat_l,conf_l,confPred,ignore.parm.uncertainty = TRUE,silent = TRUE)


resultsOut = list()
resultsOut$objective = runLenghtDepCov$opt$objective
resultsOut$rlIndexLenghtDepCov = runLenghtDepCov$rl$logLengthIndex


load("NEAhadLengthSimple/resultsExp.RData")
expect_equal(resultsOut$objective, resultsExp$objective,tolerance = 1e-4)
expect_equal(resultsOut$rlIndexLenghtDepCov, resultsExp$rlIndexLenghtDepCov,tolerance = 1e-3)

test_that("Plot runs without error", {
  expect_silent(plotResults(runLenghtDepCov, what = "sunAlt"))
  expect_silent(plotResults(runLenghtDepCov, what = "depth"))
})

if(FALSE){
  resultsExp = list()
  resultsExp$objective = runLenghtDepCov$opt$objective
  resultsExp$rlIndexLenghtDepCov = runLenghtDepCov$rl$logLengthIndex
  save(resultsExp,file = "NEAhadLengthSimple/resultsExp.RData")
}
