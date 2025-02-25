suppressMessages(library(spatioTemporalIndices))

dat_l = readRDS("NEAhadLengthAgeErrors/haddock2018-2020_length_ex_rus_reduced.rds")
dat_alk = readRDS("NEAhadLengthAgeErrors/haddock2018-2020_age_ex_rus_reduced.rds")

conf_l = defConf(years = 2018:2020, # years to use, use all years with data by default
                 maxLength = 60, # Numeric = use directly; NULL = use input data to determine
                 minLength = 20, # Numeric = use directly; NULL = use input data to determine
                 spatioTemporal =0 ,
                 spatial =0,
                 rwBeta0 = 0,
                 dLength = 5,
                 reduceLength = 3,
                 stratasystem = list(dsn="NEAhadLengthAgeErrors/strata", layer = "Vintertoktet_nye_strata"),
                 applyALK = 1,
                 cutoff =100)


#Define configurations age part
conf_alk = defConf_alk(maxAge = 8,
                       minAge = 3,
                       spatioTemporal = 0,
                       spatial =0,
                       rwBeta0 = 1)

confPred = defConfPred(conf=conf_l,Depth="DATA",cellsize = 200)

#Not unique haul ID's
dat_l$station[1] =dat_l$station[100]
expect_error(fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE))
dat_l = readRDS("NEAhadLengthAgeErrors/haddock2018-2020_length_ex_rus_reduced.rds")

#Wrong conf_l$minLength
conf_l$minLength = -5
expect_error(fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE))
conf_l$minLength = 20

#Not increasing length groups
dat_l$lengthGroup[5] = dat_l$lengthGroup[6]
expect_error(fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE))
dat_l = readRDS("NEAhadLengthAgeErrors/haddock2018-2020_length_ex_rus_reduced.rds")

#Missing length group in a haul
dat_l = dat_l[-3,]
expect_error(fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE))
dat_l = readRDS("NEAhadLengthAgeErrors/haddock2018-2020_length_ex_rus_reduced.rds")

#Not unique length group widths
dat_l$lengthGroup[10] = dat_l$lengthGroup[10] + 1
expect_error(fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE))
dat_l = readRDS("NEAhadLengthAgeErrors/haddock2018-2020_length_ex_rus_reduced.rds")

#A catch is negative
dat_l$catch[100] = -1
expect_error(fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE))
dat_l = readRDS("NEAhadLengthAgeErrors/haddock2018-2020_length_ex_rus_reduced.rds")

#A catch is NA
dat_l$catch[100] = NA
expect_error(fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE))
dat_l = readRDS("NEAhadLengthAgeErrors/haddock2018-2020_length_ex_rus_reduced.rds")

#Conf years outside of year range in data
conf_l$years = 2018:2021
expect_error(fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE))
conf_l$years = 2018:2020

#Not unique hauls in dat_alk
dat_alk$station[1] =dat_alk$station[100]
expect_error(fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE))
dat_alk = readRDS("NEAhadLengthAgeErrors/haddock2018-2020_age_ex_rus_reduced.rds")

#min age less than age in dat_alk
dat_alk = dat_alk[-which(dat_alk$age<5),]
expect_error(fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE))
dat_alk = readRDS("NEAhadLengthAgeErrors/haddock2018-2020_age_ex_rus_reduced.rds")

#Max age older than age in dat_alk
conf_alk$maxAge = 999
expect_error(fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE,silent = TRUE))
conf_alk$maxAge = 8
