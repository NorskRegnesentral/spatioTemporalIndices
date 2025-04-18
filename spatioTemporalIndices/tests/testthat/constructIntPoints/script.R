suppressMessages(library(spatioTemporalIndices))

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
                 stratasystem = list(dsn="constructIntPoints/strata", layer = "Vintertoktet_nye_strata"),
                 minDepth=150,maxDepth=400,
                 applyALK = 1,
                 cutoff =100)

confPred = defConfPred(conf=conf_l,Depth="DATA",cellsize = 50)


#Unit test for constructIntPoints
intPoints = constructIntPoints(conf_l,confPred)
resultsOut =intPoints

load("constructIntPoints/resultsExp.RData")
expect_equal(resultsOut$locUTM, resultsExp$locUTM,tolerance = 1e-6)
expect_equal(resultsOut$idxStrata, resultsExp$idxStrata,tolerance = 1e-6)

if(FALSE){
  resultsExp =intPoints 
  save(resultsExp,file = "constructIntPoints/resultsExp.RData")
}
