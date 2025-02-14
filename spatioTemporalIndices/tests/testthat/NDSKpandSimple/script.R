suppressMessages(library(spatioTemporalIndices))

dat_l = readRDS("NDSKpandSimple/NDSKpand2018-2020_length.rds")
dat_l$station <- paste(format(dat_l$startdatetime,"%Y"),dat_l$station, sep = "_")

conf_l = defConf(years = 2020:2020,
               dLength = 1, # length intervall in mm
               maxLength = 25, # Numeric = use directly; NULL = use input data to determine
               minLength = 7, # Numeric = use directly; NULL = use input data to determine
               cutoff = 20,
               spatioTemporal = 0,
               nugget = 1,
               spatial = 0,
               sunAlt=c(1,1),
               splineDepth = c(6,1),
               reduceLength = 3,
               rwBeta0 = 1,
               stratasystem = list(dsn="NDSKpandSimple/strata/", layer = "shrimp_areas_NSSK"),
               minDepth=50,maxDepth=600,
               trawlWidth=11.7,
               applyALK=0,
               strataReport=0)

confPred = defConfPred(conf=conf_l,Depth="NDSKpandSimple/gebco_2023_NDSK.nc",cellsize = 50)

####Covariates
runCovariates = fitModel(dat_l,conf_l,confPred,ignore.parm.uncertainty = TRUE,silent = TRUE)

###Length dependent covariates:
conf_l$sunAlt=c(1,2)
conf_l$splineDepth = c(6,2)
confPred = defConfPred(conf=conf_l,Depth="NDSKpandSimple/gebco_2023_NDSK.nc",cellsize = 100)
runLenghtDepCov = fitModel(dat_l,conf_l,confPred,ignore.parm.uncertainty = TRUE,silent = TRUE)

#no Nugget
conf_l$sunAlt=c(1,0)
conf_lsplineDepth = c(6,0)
conf_l$nugget = 0
runNoNugget = fitModel(dat_l,conf_l,confPred,ignore.parm.uncertainty = TRUE,silent = TRUE)

#missing depth
conf_l$sunAlt=c(1,0)
conf_lsplineDepth = c(6,1)
conf_l$nugget = 1
dat_l$depth[which(dat_l$station==dat_l$station[1])] = NA
runMissingDepth = fitModel(dat_l,conf_l,confPred,ignore.parm.uncertainty = TRUE,silent = TRUE)

#Test no reduced latent space
conf_l = defConf(years = 2018:2020,
                 dLength = 1, # length intervall in mm
                 maxLength = 20, # Numeric = use directly; NULL = use input data to determine
                 minLength = 15, # Numeric = use directly; NULL = use input data to determine
                 cutoff = 20,
                 spatioTemporal = 2,
                 nugget = 1,
                 spatial = 1,
                 reduceLength = 1,
                 rwBeta0 = 0,
                 stratasystem = list(dsn="NDSKpandSimple/strata/", layer = "shrimp_areas_NSSK"),
                 trawlWidth=11.7,
                 applyALK=0)
runNoReducedSpace = fitModel(dat_l,conf_l,confPred,ignore.parm.uncertainty = TRUE,silent = FALSE)


resultsOut = list(AIC = AIC(runCovariates,runNoNugget,runLenghtDepCov,runMissingDepth,runNoReducedSpace))
resultsOut$rlIndexCovariates = runCovariates$rl$logLengthIndex
resultsOut$rlIndexLenghtDepCov = runLenghtDepCov$rl$logLengthIndex
resultsOut$rlIndexNoNugget = runNoNugget$rl$logLengthIndex
resultsOut$rlIndexMissingDepth = runMissingDepth$rl$logLengthIndex
resultsOut$rlIndexNoReducedSpace = runNoReducedSpace$rl$logLengthIndex


load("NDSKpandSimple/resultsExp.RData")
expect_equal(resultsOut$AIC, resultsExp$AIC,tolerance = 1e-4)
expect_equal(resultsOut$rlIndexCovariates, resultsExp$rlIndexCovariates,tolerance = 1e-2)
expect_equal(resultsOut$rlIndexLenghtDepCov, resultsExp$rlIndexLenghtDepCov,tolerance = 1e-2)
expect_equal(resultsOut$rlIndexNoNugget, resultsExp$rlIndexNoNugget,tolerance = 1e-2)
expect_equal(resultsOut$rlIndexMissingDepth, resultsExp$rlIndexMissingDepth,tolerance = 1e-2)
expect_equal(resultsOut$rlIndexNoReducedSpace, resultsExp$rlIndexNoReducedSpace,tolerance = 1e-2)

test_that("Plot runs without error", {
  expect_silent(plotResults(runLenghtDepCov, what = "sunAlt"))
  expect_silent(plotResults(runLenghtDepCov, what = "depth"))
})


if(FALSE){
  resultsExp = list(AIC = AIC(runCovariates,runNoNugget,runLenghtDepCov,runMissingDepth,runNoReducedSpace))
  resultsExp$rlIndexCovariates = runCovariates$rl$logLengthIndex
  resultsExp$rlIndexLenghtDepCov = runLenghtDepCov$rl$logLengthIndex
  resultsExp$rlIndexNoNugget = runNoNugget$rl$logLengthIndex
  resultsExp$rlIndexMissingDepth = runMissingDepth$rl$logLengthIndex
  resultsExp$rlIndexNoReducedSpace = runNoReducedSpace$rl$logLengthIndex
  save(resultsExp,file = "NDSKpandSimple/resultsExp.RData")
}

