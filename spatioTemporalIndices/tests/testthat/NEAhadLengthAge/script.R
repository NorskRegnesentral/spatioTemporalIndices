suppressMessages(library(spatioTemporalIndices))

dat_l = readRDS("NEAhadLengthAge/haddock2018-2020_length_ex_rus.rds")
dat_alk = readRDS("NEAhadLengthAge/haddock2018-2020_age_ex_rus.rds")

conf_l = defConf(years = 2018:2020, # years to use, use all years with data by default
                 maxLength = 60, # Numeric = use directly; NULL = use input data to determine
                 minLength = 20, # Numeric = use directly; NULL = use input data to determine
                 spatioTemporal =0 ,
                 spatial =1,
                 rwBeta0 = 1,
                 sunAlt = c(1,2),
                 splineDepth = c(6,2),
                 cbound = c(18,130),
                 dLength = 5,
                 reduceLength = 3,
                 stratasystem = list(dsn="NEAhadLengthAge/strata", layer = "Vintertoktet_nye_strata"),
                 minDepth=150,maxDepth=400,
                 applyALK = 1,
                 cutoff =100)


#Define configurations age part
conf_alk = defConf_alk(maxAge = 10,
                       minAge = 3,
                       spatioTemporal = 2,
                       spatial =1,
                       rwBeta0 = 1)

confPred = defConfPred(conf=conf_l,Depth="DATA",cellsize = 20)

# run model
run = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE)


objectiveExp = run$opt$objective
rlIndex = round(run$rl$logAgeIndex,5)
rlIndexSd = round(run$rlSd$logAgeIndex,5)
par = round(run$opt$par,5)

resultsOut = list(objectiveExp = objectiveExp,
                  rlIndex = rlIndex,
                  rlIndexSd = rlIndexSd,
                  par = par)


load("NEAhadLengthAge/resultsExp.RData")
expect_equal(resultsOut$objectiveExp, resultsExp$objectiveExp,tolerance = 1e-4)
expect_equal(resultsOut$rlIndex, resultsExp$rlIndex,tolerance = 1e-2)
expect_equal(resultsOut$rlIndexSd, resultsExp$rlIndexSd,tolerance = 1e-2)
expect_equal(resultsOut$par, resultsExp$par,tolerance = 1e-2)



saveIndex(run,file = "testthat.txt", folder = "NEAhadLengthAge/")

#Verify that save indices and standard deviations are not changed
expect_equal(read.table("NEAhadLengthAge/testthat.txt"),
             read.table("NEAhadLengthAge/testthatExp.txt"),tolerance = 1e-2)
expect_equal(read.table("NEAhadLengthAge/sdtestthat.txt"),
             read.table("NEAhadLengthAge/sdtestthatExp.txt"),tolerance = 1e-2)

#Verify that saved correlation structures are not changed
load("NEAhadLengthAge/cov_testthatExp.Rda")
covYearsExp = covYears
load("NEAhadLengthAge/cov_testthat.Rda")
expect_equal(covYearsExp,
             covYears,tolerance = 1e-2)


#test simulation
set.seed(1)
sim = simStudy(run,nsim = 1)
objectiveSimExp = sim[[1]]$opt$objective
resultsOut$objectiveSimExp = objectiveSimExp
expect_equal(resultsOut$objectiveSimExp, resultsExp$objectiveSimExp,tolerance = 1e-4)

#test jitter
set.seed(1)
jj = jit(run,njit = 1)
resultsJitter = jj$maxVecAll
load("NEAhadLengthAge/resultsJitterExp.RData")
expect_equal(resultsJitter, resultsJitterExp,tolerance = 1e-4)

#Test retro
ret = retroSTIM(run,nyears = 1)
resultsRetro = ret[[1]]$opt$objective
load("NEAhadLengthAge/resultsRetroExp.RData")
expect_equal(resultsRetro, resultsRetroExp,tolerance = 1e-4)

if(FALSE){
  resultsExp = list(objectiveExp = objectiveExp,
                    rlIndex = rlIndex,
                    rlIndexSd = rlIndexSd,
                    par = par,
                    objectiveSimExp = objectiveSimExp)
  save(resultsExp,file = "NEAhadLengthAge/resultsExp.RData")
  saveIndex(run,file = "testthatExp.txt", folder = "NEAhadLengthAge/")

  resultsJitterExp = resultsJitter
  save(resultsJitterExp,file = "NEAhadLengthAge/resultsJitterExp.RData")

  resultsRetroExp = resultsRetro
  save(resultsRetroExp,file = "NEAhadLengthAge/resultsRetroExp.RData")

}
