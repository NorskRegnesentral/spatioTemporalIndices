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
                 stratasystem = list(dsn="createMesh/strata", layer = "Vintertoktet_nye_strata"),
                 minDepth=150,maxDepth=400,
                 applyALK = 1,
                 cutoff =100)


mesh = createMesh(conf_l)$mesh

resultsOut =list(n = mesh$n, meanX = mean(mesh$loc[,1]), meanY = mean(mesh$loc[,2]))
load("createMesh/resultsExp.RData")

expect_equal(resultsOut$n, resultsExp$n,tolerance = 1e-6)
expect_equal(resultsOut$meanX, resultsExp$meanX,tolerance = 1e-6)
expect_equal(resultsOut$meanY, resultsExp$meanY,tolerance = 1e-6)


if(FALSE){
  resultsExp =list(n = mesh$n, meanX = mean(mesh$loc[,1]), meanY = mean(mesh$loc[,2]))
  save(resultsExp,file = "createMesh/resultsExp.RData")
}
