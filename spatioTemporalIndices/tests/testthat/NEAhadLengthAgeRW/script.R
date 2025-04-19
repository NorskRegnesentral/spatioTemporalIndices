suppressMessages(library(spatioTemporalIndices))

dat_l = readRDS("NEAhadLengthAgeRW/haddock2018-2020_length_ex_rus_reduced.rds")
dat_alk = readRDS("NEAhadLengthAgeRW/haddock2018-2020_age_ex_rus_reduced.rds")

conf_l = defConf(years = 2018:2020, # years to use, use all years with data by default
                 maxLength = 60, # Numeric = use directly; NULL = use input data to determine
                 minLength = 20, # Numeric = use directly; NULL = use input data to determine
                 spatioTemporal =0 ,
                 spatial =0,
                 rwBeta0 = 0,
                 dLength = 5,
                 reduceLength = 3,
                 stratasystem = list(dsn="NEAhadLengthAgeRW/strata", layer = "Vintertoktet_nye_strata"),
                 applyALK = 1,
                 cutoff =100)


#Define configurations age part
conf_alk = defConf_alk(maxAge = 8,
                       minAge = 3,
                       spatioTemporal = 0,
                       spatial =0,
                       rwBeta0 = 1)

confPred = defConfPred(conf=conf_l,Depth="DATA",cellsize = 200)

# run model with random walk beta0 in ALK
run = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE)

run2 = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE)

#RW for beta0 and beta_length in ALK
conf_alk$rwBeta0 = 1
conf_alk$betaLength = 1
run2 = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE)

#no RW for beta0 and beta_length in ALK
conf_alk$rwBeta0 = 0
conf_alk$betaLength = 2
run3 = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE)

resultsOut = list(objective = run$opt$objective,
                  objective2 = run2$opt$objective,
                  objective3 = run3$opt$objective,
                  logAgeIndex = run$rl$logAgeIndex,
                  logAgeIndex2 = run2$rl$logAgeIndex,
                  logAgeIndex3 = run3$rl$logAgeIndex)
resultsOut$AIC = AIC(run,run2,run3)



load("NEAhadLengthAgeRW/resultsExp.RData")

expect_equal(resultsOut$logAgeIndex, resultsExp$logAgeIndexExp,tolerance = 1e-4)
expect_equal(resultsOut$objective, resultsExp$objectiveExp,tolerance = 1e-4)
expect_equal(resultsOut$logAgeIndex2, resultsExp$logAgeIndex2Exp,tolerance = 1e-4)
expect_equal(resultsOut$objective2, resultsExp$objective2Exp,tolerance = 1e-4)
expect_equal(resultsOut$logAgeIndex3, resultsExp$logAgeIndex3Exp,tolerance = 1e-4)
expect_equal(resultsOut$objective3, resultsExp$objective3Exp,tolerance = 1e-4)
expect_equal(resultsOut$AIC, resultsExp$AICExp,tolerance = 1e-4)

test_that("Plot runs without error", {
  expect_silent(plotResults(run2, what = "ALK", year = 2020))
})




if(FALSE){
  resultsExp = list(objectiveExp = run$opt$objective,
                    objective2Exp = run2$opt$objective,
                    objective3Exp = run3$opt$objective,
                    logAgeIndexExp = run$rl$logAgeIndex,
                    logAgeIndex2Exp = run2$rl$logAgeIndex,
                    logAgeIndex3Exp = run3$rl$logAgeIndex)
  resultsExp$AICExp = AIC(run,run2,run3)
  save(resultsExp,file = "NEAhadLengthAgeRW/resultsExp.RData")
}


