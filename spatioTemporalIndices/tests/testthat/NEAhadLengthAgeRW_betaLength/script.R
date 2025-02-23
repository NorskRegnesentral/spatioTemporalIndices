suppressMessages(library(spatioTemporalIndices))

dat_l = readRDS("NEAhadLengthAgeRW_betaLength/haddock2018-2020_length_ex_rus_reduced.rds")
dat_alk = readRDS("NEAhadLengthAgeRW_betaLength/haddock2018-2020_age_ex_rus_reduced.rds")

conf_l = defConf(years = 2018:2020, # years to use, use all years with data by default
                 maxLength = 60, # Numeric = use directly; NULL = use input data to determine
                 minLength = 20, # Numeric = use directly; NULL = use input data to determine
                 spatioTemporal =0 ,
                 spatial =0,
                 rwBeta0 = 1,
                 dLength = 5,
                 reduceLength = 3,
                 stratasystem = list(dsn="NEAhadLengthAgeRW_betaLength/strata", layer = "Vintertoktet_nye_strata"),
                 applyALK = 1,
                 cutoff =100)


#Define configurations age part
conf_alk = defConf_alk(maxAge = 8,
                       minAge = 3,
                       spatioTemporal = 0,
                       spatial =0,
                       rwBeta0 = 1,
                       betaLength = 1)

confPred = defConfPred(conf=conf_l,Depth="DATA",cellsize = 200)

# run model
run = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE)


#no RW for beta0 and no RW for beta_length
conf_alk$rwBeta0 = 0
conf_alk$betaLength = 2
run2 = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE)

resultsOut = list(objective = run$opt$objective,
                  objective2 = run2$opt$objective,
                  logAgeIndex = run$rl$logAgeIndex,
                  logAgeIndex2 = run2$rl$logAgeIndex)
resultsOut$AIC = AIC(run,run2)

load("NEAhadLengthAgeRW_betaLength/resultsExp.RData")

expect_equal(resultsOut$logAgeIndex, resultsExp$logAgeIndexExp,tolerance = 1e-4)
expect_equal(resultsOut$objective, resultsExp$objectiveExp,tolerance = 1e-4)
expect_equal(resultsOut$logAgeIndex2, resultsExp$logAgeIndex2Exp,tolerance = 1e-4)
expect_equal(resultsOut$objective2, resultsExp$objective2Exp,tolerance = 1e-4)
expect_equal(resultsOut$AIC, resultsExp$AICExp,tolerance = 1e-4)

test_that("Plot runs without error", {
  expect_silent(plotResults(run, what = "ALK", year = 2020))
})


if(FALSE){
  resultsExp = list(objectiveExp = run$opt$objective,
                    objective2Exp = run2$opt$objective,
                    logAgeIndexExp = run$rl$logAgeIndex,
                    logAgeIndex2Exp = run2$rl$logAgeIndex,
                    AICExp = AIC(run,run2))
  save(resultsExp,file = "NEAhadLengthAgeRW_betaLength/resultsExp.RData")
}


