suppressMessages(library(spatioTemporalIndices))

dat_l = readRDS("NEAhadLengthAgeAR1/haddock2018-2020_length_ex_rus_reduced.rds")
dat_alk = readRDS("NEAhadLengthAgeAR1/haddock2018-2020_age_ex_rus_reduced.rds")

conf_l = defConf(years = 2018:2020, # years to use, use all years with data by default
                 maxLength = 60, # Numeric = use directly; NULL = use input data to determine
                 minLength = 20, # Numeric = use directly; NULL = use input data to determine
                 spatioTemporal =0 ,
                 spatial =0,
                 rwBeta0 = 0,
                 dLength = 5,
                 reduceLength = 3,
                 stratasystem = list(dsn="NEAhadLengthAgeAR1/strata", layer = "Vintertoktet_nye_strata"),
                 applyALK = 1,
                 cutoff =100)


#Define configurations age part
conf_alk = defConf_alk(maxAge = 8,
                       minAge = 3,
                       spatioTemporal = 0,
                       spatial =0,
                       rwBeta0 = 1)

confPred = defConfPred(conf=conf_l,Depth="DATA",cellsize = 200)

# run model
run = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE)


resultsOut = list(objective = run$opt$objective,
                  logAgeIndex = run$rl$logAgeIndex)

load("NEAhadLengthAgeAR1/resultsExp.RData")

expect_equal(resultsOut$logAgeIndex, resultsExp$logAgeIndexExp,tolerance = 1e-3)
expect_equal(resultsOut$objective, resultsExp$objectiveExp,tolerance = 1e-4)

if(FALSE){
  resultsExp = list(objectiveExp = run$opt$objective,
                    logAgeIndexExp = run$rl$logAgeIndex)
  save(resultsExp,file = "NEAhadLengthAgeAR1/resultsExp.RData")

}


