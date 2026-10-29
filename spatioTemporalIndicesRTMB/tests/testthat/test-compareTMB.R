#Verify that the RTMB implementation reproduces the TMB implementation (package spatioTemporalIndices)
test_that("RTMB and TMB implementations give identical fits", {
  skip_if_not_installed("spatioTemporalIndices")
  T_ <- asNamespace("spatioTemporalIndices")
  R_ <- asNamespace("spatioTemporalIndicesRTMB")

  dat_l = readRDS("NEAhadLengthSimple/haddock2018-2020_length_ex_rus_reduced.rds")
  fit = function(ns){
    conf_l = ns$defConf(years = 2018:2020, maxLength = 65, minLength = 20, spatioTemporal = 1, spatial = 1, rwBeta0 = 1,
                        sunAlt = c(1,2), splineDepth = c(6,2), dLength = 5, reduceLength = 3,
                        stratasystem = list(dsn="NEAhadLengthSimple/strata", layer = "Vintertoktet_nye_strata"),
                        minDepth=150, maxDepth=400, applyALK = 0, cutoff = 100, strataReport = 1)
    confPred = ns$defConfPred(conf=conf_l, Depth="DATA", cellsize = 50)
    ns$fitModel(dat_l, conf_l, confPred, silent = TRUE)
  }
  runT = fit(T_)
  runR = fit(R_)

  expect_equal(runR$opt$objective, runT$opt$objective, tolerance = 1e-8)
  expect_equal(runR$opt$par, runT$opt$par, tolerance = 1e-6)
  expect_identical(names(runR$rep$value), names(runT$rep$value))
  expect_equal(runR$rl, runT$rl, tolerance = 1e-6)
  expect_equal(runR$rlSd, runT$rlSd, tolerance = 1e-6)
  repT = runT$obj$report()
  repR = runR$obj$report()
  for(n in names(repT)) expect_equal(repR[[n]], repT[[n]], tolerance = 1e-6, info = n)
})
