suppressMessages(library(spatioTemporalIndices))
suppressMessages(library(spatioTemporalALK))


dat_l = readRDS("haddock2018-2020_length_ex_rus.rds")
dat_alk = readRDS("haddock2018-2020_age_ex_rus.rds")

conf_l = defConf(years = 2018:2020, # years to use, use all years with data by default
                 maxLength = 60, # Numeric = use directly; NULL = use input data to determine
                 minLength = 20, # Numeric = use directly; NULL = use input data to determine
                 spatioTemporal =0 ,
                 spatial =1,
                 rwBeta0 = 1,
                 sunAlt = c(1,2),
                 splineDepth = c(6,2),
                 dLength = 5,
                 reduceLength = 3,
                 stratasystem = list(dsn="strata", layer = "Vintertoktet_nye_strata"),
                 minDepth=150,maxDepth=400,
                 applyALK = 1,
                 cutoff =120, cbound = 130)


#Define configurations age part
conf_alk = defConf_alk(maxAge = 10,
                       minAge = 3,
                       spatioTemporal = 2,
                       spatial =1,
                       rwBeta0 = 1)

confPred = defConfPred(conf=conf_l,Depth="DATA",nInt=5000)

# run model
start_time <- Sys.time()
run = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed


cat(run$opt$objective,"\n", file="res.out")
#cat(run$opt$objective,"\n", file="res.EXP")
apply(round(run$rl$logAgeIndex,3),1, function(f)cat(f,"\n", file="res.out", append = TRUE))
#apply(round(run$rl$logAgeIndex,3),1, function(f)cat(f,"\n", file="res.EXP", append = TRUE))



