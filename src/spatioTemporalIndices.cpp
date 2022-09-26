#define TMB_LIB_INIT R_init_spatioTemporalIndices

#include <TMB.hpp>
#include <cmath>
using namespace tmbutils;

#include "../inst/include/define.hpp"
#include "../inst/include/alk.hpp"
#include "../inst/include/index.hpp"
#include "../inst/include/indexPred.hpp"


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()
  using namespace density;

  //Load data and parameters
  dataSet<Type> dat;
  DATA_MATRIX(fishObsMatrix); dat.fishObsMatrix = fishObsMatrix; //Observations in a matrix
  DATA_VECTOR(dist); dat.dist = dist; // Distance of each haul
  DATA_STRUCT(A_ListS, LOSM_t); //TODO
  DATA_STRUCT(A_ListST, LOSM_t);
  DATA_STRUCT(spdeMatricesS,spde_t); //Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  DATA_STRUCT(spdeMatricesST,spde_t); //Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  DATA_IVECTOR(nStationsEachYear) dat.nStationsEachYear = nStationsEachYear;
  DATA_INTEGER(numberOfLengthGroups); dat.numberOfLengthGroups = numberOfLengthGroups;
  DATA_SPARSE_MATRIX(S_depth); dat.S_depth = S_depth;
  DATA_INTEGER(Sdim); dat.Sdim = Sdim;
  DATA_IVECTOR(lengthGroupsReduced); dat.lengthGroupsReduced = lengthGroupsReduced;//used when reducing length dimension
  DATA_MATRIX(predMatrix);  dat.predMatrix = predMatrix;//used to leave out data from likelihood
  DATA_MATRIX(X_sunAlt); dat.X_sunAlt = X_sunAlt;
  DATA_MATRIX(X_sunAltReport); dat.X_sunAltReport = X_sunAltReport;
  DATA_MATRIX(X_sunAltIntegrate); dat.X_sunAltIntegrate = X_sunAltIntegrate;
  DATA_MATRIX(X_depth); dat.X_depth = X_depth;
  DATA_MATRIX(X_depthReport); dat.X_depthReport = X_depthReport;
  DATA_MATRIX(X_depth_int); dat.X_depth_int = X_depth_int;
  DATA_VECTOR(weigthLength);  dat.weigthLength = weigthLength;//wieght used when reducing length dimension
  DATA_INTEGER(zeroInflated); dat.zeroInflated = zeroInflated;//if !=0, zero inflated model is used, work in progress
  DATA_INTEGER(nBasisSunAlt); dat.nBasisSunAlt = nBasisSunAlt;
  DATA_INTEGER(obsModel); dat.obsModel = obsModel;
  DATA_VECTOR(obsVector); dat.obsVector = obsVector;//Observations in a vector. TODO: remove the observation matrix.
  DATA_IVECTOR(idxStart); dat.idxStart = idxStart;//Index helper when refering to the observation vector
  DATA_VECTOR_INDICATOR(keep,obsVector);
  DATA_VECTOR(xInt);  dat.xInt = xInt;//Used to calcualte COF
  DATA_VECTOR(yInt);  dat.yInt = yInt;//Used to calcualte COF
  DATA_SPARSE_MATRIX(ApredS); dat.ApredS = ApredS;
  DATA_SPARSE_MATRIX(ApredST); dat.ApredST = ApredST;
  DATA_INTEGER(spatial);  dat.spatial = spatial;//Configuration spatial, include if 1
  DATA_INTEGER(spatioTemporal);  dat.spatioTemporal = spatioTemporal;//Configuration spatio-temporal, include if 1
  DATA_INTEGER(splineDepth);  dat.splineDepth = splineDepth;//Configuration spatio-temporal, include if 1
  DATA_INTEGER(useNugget);  dat.useNugget = useNugget;//Configuration useNugget, include if 1
  DATA_INTEGER(simulateProcedure);  dat.simulateProcedure = simulateProcedure;
  DATA_VECTOR(pcPriorsRange); dat.pcPriorsRange = pcPriorsRange;
  DATA_VECTOR(pcPriorsSD); dat.pcPriorsSD = pcPriorsSD;
  DATA_INTEGER(usePCpriors); dat.usePCpriors = usePCpriors;
  DATA_INTEGER(rwBeta0); dat.rwBeta0 = rwBeta0;
  DATA_IVECTOR(idxStrata); dat.idxStrata = idxStrata;
  DATA_VECTOR(areas); dat.areas = areas;


  //ALK stuff----------------
  DATA_IVECTOR(age); dat.age = age;//The response
  DATA_IVECTOR(ageNotTruncated); dat.ageNotTruncated = ageNotTruncated;
  DATA_VECTOR(length); dat.length = length;//Covariate
  DATA_IVECTOR(readability); dat.readability = readability;
  DATA_IVECTOR(ageRange);dat.ageRange = ageRange;
  DATA_IVECTOR(idx1); dat.idx1 = idx1;//Bookkeeping
  DATA_IVECTOR(idx2); dat.idx2 = idx2;
  DATA_STRUCT(spdeMatricesST_alk,spde_t); //TODO: Include in dataset
  DATA_STRUCT(A_alk_list, LOSM_t);
  DATA_SPARSE_MATRIX(Apred_alk);dat.Apred_alk = Apred_alk;
  DATA_VECTOR(lengthGroups); dat.lengthGroups = lengthGroups;//used when applying ALK
  DATA_VECTOR(dL);dat.dL = dL;
  DATA_INTEGER(applyALK);dat.applyALK = applyALK;
  DATA_INTEGER(rwBeta0_alk); dat.rwBeta0_alk = rwBeta0_alk;//Indicator if including AR1-structure on beta0
  DATA_INTEGER(maxAge); dat.maxAge = maxAge;
  DATA_INTEGER(minAge); dat.minAge = minAge;
  DATA_VECTOR(pcPriorsALKRange); dat.pcPriorsALKRange = pcPriorsALKRange;
  DATA_VECTOR(pcPriorsALKSD); dat.pcPriorsALKSD = pcPriorsALKSD;
  DATA_INTEGER(usePCpriorsALK); dat.usePCpriorsALK = usePCpriorsALK;

  //-----------------------------

  paraSet<Type> par;
  PARAMETER_MATRIX(beta0); par.beta0 = beta0;
  PARAMETER_VECTOR(betaSun); par.betaSun = betaSun;
  PARAMETER_VECTOR(betaDepth); par.betaDepth = betaDepth;
  PARAMETER_VECTOR(log_lambda); par.log_lambda = log_lambda;
  PARAMETER_VECTOR(log_sigma); par.log_sigma = log_sigma;
  PARAMETER_VECTOR(log_kappa); par.log_kappa = log_kappa;
  PARAMETER(logSize); par.logSize = logSize;
  PARAMETER_VECTOR(tan_rho_t); par.tan_rho_t = tan_rho_t;
  PARAMETER_VECTOR(tan_rho_l); par.tan_rho_l = tan_rho_l;
  PARAMETER_VECTOR(delta_z); par.delta_z = delta_z;
  PARAMETER_ARRAY(xS); par.xS = xS;
  PARAMETER_ARRAY(xST); par.xST = xST;
  PARAMETER_ARRAY(nugget); par.nugget = nugget;
  PARAMETER_VECTOR(log_sigma_beta0); par.log_sigma_beta0 = log_sigma_beta0;

  //ALK stuff-------------------
  PARAMETER_MATRIX(beta0_alk); par.beta0_alk = beta0_alk;//Intercepts
  PARAMETER_VECTOR(log_sigma_beta0_alk);par.log_sigma_beta0_alk = log_sigma_beta0_alk;
  PARAMETER_VECTOR(betaLength_alk); par.betaLength_alk = betaLength_alk;//Regression parameters
  PARAMETER_VECTOR(logSigma_alk);par.logSigma_alk = logSigma_alk;
  PARAMETER_VECTOR(logKappa_alk);par.logKappa_alk = logKappa_alk;
  PARAMETER_VECTOR(transRho_alk);par.transRho_alk = transRho_alk;
  PARAMETER_ARRAY(xS_alk); par.xS_alk = xS_alk;// Ordering: xST.col(age).col(year).col(spatial)
  PARAMETER_ARRAY(xST_alk); par.xST_alk = xST_alk;// Ordering: xST.col(age).col(year).col(spatial)
  //----------------------------


  Type nll = 0.0;
  //Likelihood contribution from length part
  nll += nllIndex(dat,par,spdeMatricesS,spdeMatricesST,A_ListS,A_ListST, keep);

  //Likelihood contribution from ALK part
  if(dat.applyALK==1){
    nll += nllALK(dat,par,spdeMatricesST_alk,A_alk_list);
  }

  //Predict and report indices
  indexPred(dat,par, this);

  array<Type> ALK_hauls;
  if(dat.applyALK==1){
    ALK_hauls = ALKhauls(dat,par,A_ListST); //ALKs at for integration points
    REPORT(ALK_hauls);
  }


  return nll;
}
