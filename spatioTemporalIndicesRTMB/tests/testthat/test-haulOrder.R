#The fit must not depend on the order of the hauls in the length data
test_that("fit is independent of the order of hauls", {
  suppressMessages(library(spatioTemporalIndicesRTMB))
  dat_l = readRDS("NEAhadLengthSimple/haddock2018-2020_length_ex_rus_reduced.rds")
  set.seed(1)
  st = unique(dat_l$station)
  shuffled = dat_l[order(match(dat_l$station, sample(st)), dat_l$lengthGroup), ]
  expect_true(is.unsorted(as.integer(format(shuffled$startdatetime, "%Y")))) #Hauls are no longer sorted by year

  conf_l = defConf(years = 2018:2020, maxLength = 65, minLength = 20, spatioTemporal = 0, spatial = 1, rwBeta0 = 1,
                   sunAlt = c(1,1), splineDepth = c(6,1), dLength = 5, reduceLength = 3,
                   stratasystem = list(dsn="NEAhadLengthSimple/strata", layer = "Vintertoktet_nye_strata"),
                   minDepth=150, maxDepth=400, applyALK = 0, cutoff = 100)
  confPred = defConfPred(conf=conf_l, Depth="DATA", cellsize = 50)
  runSorted = fitModel(dat_l, conf_l, confPred, silent = TRUE, ignore.parm.uncertainty = TRUE)
  runShuffled = fitModel(shuffled, conf_l, confPred, silent = TRUE, ignore.parm.uncertainty = TRUE)
  expect_equal(runShuffled$opt$objective, runSorted$opt$objective, tolerance = 1e-8)
  expect_equal(runShuffled$rl$logLengthIndex, runSorted$rl$logLengthIndex, tolerance = 1e-6)
})
