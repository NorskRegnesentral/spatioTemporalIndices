suppressMessages(library(spatioTemporalIndices))

dat_l = readRDS("NEAhadLengthAge/haddock2018-2020_length_ex_rus_reduced.rds")
dat_alk = readRDS("NEAhadLengthAge/haddock2018-2020_age_ex_rus_reduced.rds")

conf_l = defConf(years = 2018:2020, # years to use, use all years with data by default
                 maxLength = 60, # Numeric = use directly; NULL = use input data to determine
                 minLength = 20, # Numeric = use directly; NULL = use input data to determine
                 spatioTemporal =0 ,
                 spatial =1,
                 rwBeta0 = 1,
                 sunAlt = c(1,0),
                 splineDepth = c(6,0),
                 cbound = c(18,130),
                 dLength = 5,
                 reduceLength = 3,
                 stratasystem = list(dsn="NEAhadLengthAge/strata", layer = "Vintertoktet_nye_strata"),
                 minDepth=150,maxDepth=400,
                 applyALK = 1,
                 cutoff =100)


#Define configurations age part
conf_alk = defConf_alk(maxAge = 8,
                       minAge = 3,
                       spatioTemporal = 2,
                       spatial =1,
                       rwBeta0 = 0)

confPred = defConfPred(conf=conf_l,Depth="DATA",cellsize = 50)

# run model
run = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE)

#Two stage age
run_twoStage = fitModel(dat_l,conf_l, twoStage = TRUE, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE)
expect_equal(run$rl$logAgeIndex, run_twoStage$rl$logAgeIndex,tolerance = 1e-3)


objectiveExp = run$opt$objective
rlIndex = round(run$rl$logAgeIndex,5)
rlIndexSd = round(run$rlSd$logAgeIndex,5)
par = round(run$opt$par,5)

resultsOut = list(objectiveExp = objectiveExp,
                  rlIndex = rlIndex,
                  rlIndexSd = rlIndexSd,
                  par = par)


#Verify all indices and parameters are as expected
load("NEAhadLengthAge/resultsExp.RData")
expect_equal(resultsOut$objectiveExp, resultsExp$objectiveExp,tolerance = 1e-4)
expect_equal(resultsOut$rlIndex, resultsExp$rlIndex,tolerance = 1e-2)
expect_equal(resultsOut$rlIndexSd, resultsExp$rlIndexSd,tolerance = 1e-2)
expect_equal(resultsOut$par, resultsExp$par,tolerance = 1e-2)

#Verify that save indices on ICES-format are as expected
write_indices_ICES_format(run,file = "NEAhadLengthAge/indexFile.dat", name = "nameOfSurvey",digits = 0)
write_indices_ICES_format(run,file = "NEAhadLengthAge/indexFileVar.dat",variance = TRUE, name = "nameOfSurvey",digits = 2)
expect_equal(readLines("NEAhadLengthAge/indexFile.dat"),
             readLines("NEAhadLengthAge/indexFileExp.dat"))
expect_equal(readLines("NEAhadLengthAge/indexFileVar.dat"),
             readLines("NEAhadLengthAge/indexFileVarExp.dat"))

#Verify that save covaraince structures are as expected
write_covariance_matrices(run,"NEAhadLengthAge/yearlyCov.rds")
cov = readRDS("NEAhadLengthAge/yearlyCov.rds")
covExp = readRDS("NEAhadLengthAge/yearlyCovExp.rds")
expect_equal(cov,
             covExp,tolerance = 1e-2)


#Reduce complexity for time efficency
conf_alk$spatioTemporal = 0
conf_alk$spatial = 0
runSimpler = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE)


#Test simulation
set.seed(1)
sim = simStudy(runSimpler,nsim = 1)
objectiveSimExp = sim[[1]]$opt$objective
resultsOut$objectiveSimExp = objectiveSimExp
expect_equal(resultsOut$objectiveSimExp, resultsExp$objectiveSimExp,tolerance = 1e-4)

#Test jitter
set.seed(1)
jj = jit(runSimpler,njit = 1)
resultsJitter = jj$maxVecAll
load("NEAhadLengthAge/resultsJitterExp.RData")
expect_equal(resultsJitter, resultsJitterExp,tolerance = 1e-4)

#Test retro
ret = retroSTIM(runSimpler,nyears = 2)
resultsRetro = ret[[2]]$opt$objective
load("NEAhadLengthAge/resultsRetroExp.RData")
expect_equal(resultsRetro, resultsRetroExp,tolerance = 1e-4)

#Test no errors in plots
test_that("Plot runs without error", {
  expect_silent(plotResults(run, what = "ALK", year = 2020))
  expect_silent(plotResults(run, what = "ALK", year = 2020, lon_lat = c(13,78)))
  expect_silent(plotResults(run, what = "space",year = 2020, age = 5))
  expect_silent(plotResults(run, what = "space", year = 2020, length = 40))
  expect_silent(plotResults(run, what = "variance"))
  expect_silent(plotResults(run, what = "correlation"))
})


if(FALSE){
  resultsExp = list(objectiveExp = objectiveExp,
                    rlIndex = rlIndex,
                    rlIndexSd = rlIndexSd,
                    par = par,
                    objectiveSimExp = objectiveSimExp)
  save(resultsExp,file = "NEAhadLengthAge/resultsExp.RData")

  write_indices_ICES_format(run,file = "NEAhadLengthAge/indexFileExp.dat", name = "nameOfSurvey",digits = 0)
  write_indices_ICES_format(run,file = "NEAhadLengthAge/indexFileVarExp.dat",variance = TRUE, name = "nameOfSurvey",digits = 2)
  write_covariance_matrices(run,"NEAhadLengthAge/yearlyCovExp.rds")

  resultsJitterExp = resultsJitter
  save(resultsJitterExp,file = "NEAhadLengthAge/resultsJitterExp.RData")

  resultsRetroExp = resultsRetro
  save(resultsRetroExp,file = "NEAhadLengthAge/resultsRetroExp.RData")

}
