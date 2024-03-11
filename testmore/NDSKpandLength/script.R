suppressMessages(library(spatioTemporalIndices))

dat_l = readRDS("NDSKpand2018-2020_length.rds")

conf_l = defConf(years = 2018:2020,
               dLength = 1, # length intervall in mm
               maxLength = 30, # Numeric = use directly; NULL = use input data to determine
               minLength = 7, # Numeric = use directly; NULL = use input data to determine
               cutoff = 30,
               spatioTemporal = 1,
               spatial = 0,
               sunAlt=c(0,0),
               reduceLength = 3,
               rwBeta0 = 0,
               stratasystem = list(dsn="strata", layer = "shrimp_areas_NSSK"),
               minDepth=50,maxDepth=600,
               trawlWidth=11.7,
               applyALK=0,
               strataReport=0)

confPred = defConfPred(conf=conf_l,Depth="DATA",cellsize = 20)

# run model
start_time <- Sys.time()
run = fitModel(dat_l,conf_l,confPred,ignore.parm.uncertainty = TRUE)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed


cat(run$opt$objective,"\n", file="res.out")
#cat(run$opt$objective,"\n", file="res.EXP")
apply(round(run$rl$logLengthIndex,3),1, function(f)cat(f,"\n", file="res.out", append = TRUE))
#apply(round(run$rl$logAgeIndex,3),1, function(f)cat(f,"\n", file="res.EXP", append = TRUE))



