suppressMessages(library(spatioTemporalIndices))

dat_l = readRDS("NEAhadTechnical/haddock2018-2020_length_ex_rus_reduced.rds")

conf_l = defConf(years = 2018:2020, # years to use, use all years with data by default
                 maxLength = 65, # Numeric = use directly; NULL = use input data to determine
                 minLength = 20, # Numeric = use directly; NULL = use input data to determine
                 spatial =1,
                 cbound = c(18,130),
                 dLength = 5,
                 reduceLength = 3,
                 stratasystem = list(dsn="NEAhadTechnical/strata", layer = "Vintertoktet_nye_strata"),
                 applyALK = 0,
                 cutoff =10)

confPred = defConfPred(conf=conf_l,Depth="DATA",cellsize = 100)

###Check number of mesh points is not changed:
runLenghtDepCov = fitModel(dat_l,conf_l,confPred,ignore.parm.uncertainty = TRUE,silent = TRUE,runModel = FALSE)
mesh = attributes(runLenghtDepCov$data)$meshS

resultsOut = list()
resultsOut$mesh_n = mesh$n

load("NEAhadTechnical/resultsExp.RData")
expect_equal(resultsOut$mesh_n, resultsExp$mesh_n,tolerance = 1e-10)

if(FALSE){
  resultsExp = list()
  resultsExp$mesh_n = mesh$n
  save(resultsExp,file = "NEAhadTechnical/resultsExp.RData")
}
