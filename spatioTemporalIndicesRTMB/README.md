# spatioTemporalIndicesRTMB

RTMB implementation of the spatio-temporal index model STIM ([Breivik et al. 2021](https://doi.org/10.1093/icesjms/fsab073), [Breivik et al. 2024](https://doi.org/10.1093/icesjms/fsae079)).

The package has the same user interface as `spatioTemporalIndices` (`defConf`, `defConf_alk`, `defConfPred`, `fitModel`, `plotResults`, `write_indices_ICES_format`, ...). The difference is that the objective function (likelihood of the catch-at-length and age-at-length models, and prediction of the indices) is written in R with [RTMB](https://github.com/kaskr/RTMB) instead of as a TMB C++ template. No compilation is needed, and the model can be read, modified and debugged in R.

## Where the model lives

| TMB (`spatioTemporalIndices`) | RTMB (this package), `R/objective.R` |
|---|---|
| `src/spatioTemporalIndices.cpp` | `nllSTIM()` |
| `inst/include/index.hpp` (`nllIndex`) | `nllIndex()` |
| `inst/include/alk.hpp` (`nllALK`, `ALK`, `ALKhauls`) | `nllALK()`, `alkLinPred()`, `ALKhauls()` |
| `inst/include/indexPred.hpp` (`indexPred`) | `indexPred()` |
| `TMB::MakeADFun(data, par, DLL = ...)` | `makeObj(data, par, ...)` → `RTMB::MakeADFun(function(p) nllSTIM(p, data), par, ...)` |

Translation of the main TMB building blocks:

| TMB C++ | RTMB |
|---|---|
| `R_inla::Q_spde(spde, kappa)` | `Qspde()`: `kappa^4*M0 + 2*kappa^2*M1 + M2` |
| `SEPARABLE(f, g)(x)` (g on inner, f on outer dimension) | `RTMB::dseparable(g, f)(x)` (first function on first dimension) |
| `GMRF(Q)`, `AR1(rho)` | `RTMB::dgmrf(x, 0, Q, log = TRUE)`, `RTMB::dautoreg(x, phi = rho, log = TRUE)` |
| loops over years, hauls and length groups | sparse projector matrices and matrix products |

## Agreement with the TMB version

Fitted with both packages on the test data (`fitModel(..., newtonsteps = 1)`), maximum absolute differences:

| Configuration | objective | parameters | log-indices | SD of log-indices |
|---|---|---|---|---|
| length only, length-dependent sun and depth effects | 1.5e-11 | 2.6e-13 | 4.3e-14 | 3.4e-12 |
| spatial + AR1 spatio-temporal, strata report | 7.3e-12 | 1.3e-11 | 8.7e-12 | 1.6e-10 |
| length + ALK (spatial and spatio-temporal ALK) | 7.3e-12 | 2.8e-12 | 8.9e-13 | 1.8e-11 |
| length + ALK with random walks, AR1 spatio-temporal ALK, time-varying length effect, length-dependent covariates, strata report | 2.2e-11 | 3.4e-10 | 3.5e-11 | 3.7e-11 |

Note on the optimization path: with identical fixed-effect estimates (difference ~1e-11), the random-effect modes, and hence the indices, can differ by ~3e-7 (relative) depending on the sequence of function evaluations. This happens within TMB itself (e.g. `newtonsteps = 1` vs `2` on the haddock test data) and is far below the uncertainty of the indices. The test of the ICES index files (written as rounded integers) therefore compares the numbers with a relative tolerance of 1e-6 instead of character by character.

REPORTed quantities used for plotting (`lengthIndexDetailed`, `ageIndexDetailed`, `ALK_int`, `ALK`, `ALK_hauls`) agree to < 1e-9. Computation times were similar to TMB. `tests/testthat/test-compareTMB.R` repeats this comparison when `spatioTemporalIndices` is installed; the other tests check the RTMB fits against the expected values stored from the TMB version.

## Not (yet) implemented

- Simulation (`simStudy()`): the TMB version uses `SIMULATE` blocks in C++; `simStudy()` stops with an error here.
- One-step-ahead residuals (`keep`/`OBS`) are not used by the package and are not implemented.
- As in the TMB version, only the Poisson observation model (`obsModel = 2`) is available.

## Implementation notes for developers

- The package imports RTMB so that RTMB's generics (`colSums`, `rowSums`, `%*%`, densities, ...) dispatch on AD types. Use `base::matrix`, `base::rowSums`, etc. for data-only objects, since RTMB's versions return AD objects during taping.
- Arrays that are filled by assignment inside the objective are created with `RTMB::AD()`, so that `[<-` dispatches to the AD method.
- `obj$report()` evaluates the objective with numeric parameters; sparse products go through `mm()` so that they return base matrices in that case.
- Use `RTMB::sdreport()`, not `TMB::sdreport()`, on RTMB objects; otherwise ADREPORTed quantities are missing.
- Do not call `TMB::FreeADFun()` on RTMB objects: R crashes at the next garbage collection (the tapes are freed twice). `freeObj()` drops the references to the tapes instead.
- `obj$env$ADreportIndex()` does not work after `freeObj()`; positions of ADREPORTed quantities are taken from `names(run$rep$value)` instead (see `posteriorMean()`).
