## spatioTemporalIndices
Git page for the R-package `spatioTemporalIndices`. 

This model uses catch-at-length and age-at-length observations to construct indices-at-length and indices-at-age, along with yearly covariance matrices that include uncertainty in both age-at-length and catch-at-length.

### Installation

The model is installed by typing 

```R
devtools::install_github("NorskRegnesentral/spatioTemporalIndices")
```

# Quick example

Here is a quick example of how to generate indices-at-age with associated covariance structures. For the full R code to run a similar example, we refer to the folder  `testmore/NEAhadLengthAge`. We use Norwegian haddock observations from the Barents Sea winter survey as our data.


### Data

The length data must be in the format of a data frame with the following columns: haul ID, length group, time, distance trawled, latitude, longitude, and the number of fish caught. Note that each row represents one observed length group in a haul.

```R
station   lengthGroup   startdatetime   distance  latitude longitude   catch
idHaul1   5         2018-02-02 11:10:46    0.89     73.34    18.13     0
idHaul1   10        2018-02-02 11:10:46    0.89     73.34    18.13     20
idHaul1   15        2018-02-02 11:10:46    0.89     73.34    18.13     52
idHaul1   20        2018-02-02 11:10:46    0.89     73.34    18.13     22
```

The age-at-length data must be in the format of a data frame with the following columns: haul ID, time, latitude, longitude, length of fish, and readability. Note that each row represents one observed fish. The station ID needs to match the ID given in the length data above.


```R
station   startdatetime       latitude longitude length readability
idHaul1   2018-02-02 11:10:46  73.34    18.13     32           1
idHaul1   2018-02-02 11:10:46  73.34    18.13     28           1
idHaul1   2018-02-02 11:10:46  73.34    18.13     17           1
idHaul1   2018-02-02 11:10:46  73.34    18.13     54           1
```

### Confgurations

Set up configurations for catch-at-length model:

```R
conf_l = defConf(years = 2018:2020, # years to use, 
                 maxLength = 75, 
                 minLength = 20, 
                 spatioTemporal =0 ,
                 spatial =1,
                 stratasystem = list(dsn="strata", layer = "Vintertoktet_nye_strata"),
                 applyALK = 1)
```

Set up configurations for age-at-length model. Note that the package `spatioTemporalALK` needs to be installed (https://github.com/NorskRegnesentral/spatioTemporalALK). 

```R
conf_alk = defConf_alk(maxAge = 10,
                       minAge = 3,
                       spatioTemporal = 2,
                       spatial =1)
```

For documentation of the configurations, see `?defConf` and `?defConf_alk`.


Set up prediction configurations:

```R
confPred = defConfPred(conf=conf_l,cellsize = 20)
```

### Fit model

Fit the model 
```R
run = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk)
```


### Extract indices and covariance structures
The indices and their associated standard deviations can be accessed in the list of reported quantities:
```R
run$rl$logAgeIndex
run$rlSd$logAgeIndex
```

The indices and corresponding covariance structures can be saved by

```R
saveIndex(run,file = "index.txt", folder = "")
```

This will save the files `index.txt` and `cov_index.Rda`, containing the indices and a list with all yearly covariance matrices.


### Use index and covaraince structures in assessment

For the use of the indices and covariance structures in the state space assessment model SAM, we refer to the SAM help file at  http://www.nielsensweb.org/configurations.html.
