#define TMB_LIB_INIT R_init_spatioTemporalIndices

#include <TMB.hpp>
#include <cmath>
using namespace tmbutils;

#include "../inst/include/define.hpp"
#include "../inst/include/alk.hpp"
#include "../inst/include/index.hpp"


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
  PARAMETER_ARRAY(xST_alk); par.xST_alk = xST_alk;// Ordering: xST.col(age).col(year).col(spatial)
  //----------------------------

  Type nll = 0.0;

  nll += nllIndex(dat,par,spdeMatricesS,spdeMatricesST,A_ListS,A_ListST, keep);

  array<Type> ALK_int;
  if(dat.applyALK==1){
    nll += nllALK(dat,par,spdeMatricesST_alk,A_alk_list);
    ALK_int = ALK(dat,par); //Calcualate ALKs at for integration points
  }

  //Calculate indices-----------------------------------
  int nYears = dat.nStationsEachYear.size();
  int nInt = dat.xInt.size();
  vector<Type> sigma = exp(par.log_sigma);
  vector<Type> kappa = exp(par.log_kappa);

  Type scaleS = Type(1)/((4*3.14159265)*kappa[0]*kappa[0]); //No effect on results, but needed for interpreting the sigma^2 parameter as marginal variance. See section 2.1 in Lindgren (2011)
  Type scaleST = Type(1)/((4*3.14159265)*kappa[1]*kappa[1]); //No effect on results, but needed for interpreting the sigma^2 parameter as marginal variance


  array<Type>  lengthIndex(nYears,numberOfLengthGroups);
  array<Type>  lengthIndexDetailed(nYears,numberOfLengthGroups, nInt);
  lengthIndex.setZero();
  lengthIndexDetailed.setZero();

  vector<Type> deltaPredST;
  vector<Type> deltaPredS;

  int nSplineDepth = par.betaDepth.size()/2;
  vector<Type> parDepth1(nSplineDepth);
  vector<Type> parDepth2(nSplineDepth);
  for(int i =0; i<nSplineDepth; ++i){
    parDepth1(i) = par.betaDepth(i);
    parDepth2(i) = par.betaDepth(i + nSplineDepth);
  }

  vector<Type> betaSunLow(dat.nBasisSunAlt*2);
  vector<Type> betaSunHigh(dat.nBasisSunAlt*2);
  for(int i =0; i<(dat.nBasisSunAlt*2); ++i){
    betaSunLow(i) = par.betaSun(i);
    betaSunHigh(i) = par.betaSun(i + dat.nBasisSunAlt*2);
  }

  vector<Type> depthEffectInt1=X_depth_int*parDepth1;
  vector<Type> depthEffectInt2=X_depth_int*parDepth2;
  vector<Type> timeInDayEffectIntLow = X_sunAltIntegrate*betaSunLow;
  vector<Type> timeInDayEffectIntHigh = X_sunAltIntegrate*betaSunHigh;

  for(int y=0; y<nYears; ++y){
    for(int l =0; l<numberOfLengthGroups; ++l){
      Type covariatesConvexW = (numberOfLengthGroups-l-1)/(numberOfLengthGroups-1);
      if(lengthGroupsReduced(0)==lengthGroupsReduced(1)){
        deltaPredS = weigthLength(l)*ApredS * xS.col(lengthGroupsReduced(l)).matrix()+
          (1-weigthLength(l))*ApredS * xS.col(lengthGroupsReduced(l)+1).matrix();
        deltaPredST = weigthLength(l)*ApredST * xST.col(lengthGroupsReduced(l)).col(y).matrix()+
          (1-weigthLength(l))*ApredST * xST.col(lengthGroupsReduced(l)+1).col(y).matrix();
      }else{
        deltaPredS = ApredS * xS.col(lengthGroupsReduced(l)).matrix();
        deltaPredST = ApredST * xST.col(lengthGroupsReduced(l)).col(y).matrix();
      }
      for(int i =0; i<nInt; ++i){
        lengthIndexDetailed(y,l,i) =  exp(beta0.row(y)(l) +
          covariatesConvexW*timeInDayEffectIntLow(0) + (1-covariatesConvexW)*timeInDayEffectIntHigh(0)+
          covariatesConvexW*depthEffectInt1(i) + (1-covariatesConvexW)*depthEffectInt2(i)+
          deltaPredS(i)/sqrt(scaleS)*sigma(0)+
          deltaPredST(i)/sqrt(scaleST)*sigma(1));
      }
    }
  }

  int nAges = 2;//Dummy-number
  if(dat.applyALK==1){
    nAges = dat.ageRange(1) - dat.ageRange(0) + 1;
  }
  array<Type> ageIndexDetailed(nYears,nAges, nInt);
  array<Type> ageIndex(nYears,nAges);
  ageIndexDetailed.setZero();
  ageIndex.setZero();

  int nStrata = 0;
  for(int i=0; i<idxStrata.size(); ++i){
    if(nStrata<idxStrata(i)){
      nStrata = idxStrata(i);
    }
  }
  array<Type> lengthIndexStrata(nYears,numberOfLengthGroups,nStrata);
  lengthIndexStrata.setZero();

  vector<Type> nIntStrata(nStrata);
  nIntStrata.setZero();
  for(int i=0; i<nInt; ++i){
    nIntStrata(idxStrata(i)-1) = nIntStrata(idxStrata(i)-1) +1;
  }


  for(int y=0; y<nYears; ++y){
    for(int l = 0; l<numberOfLengthGroups; ++l){
      for(int i=0; i<nInt; ++i){
        if(dat.applyALK==1){
          for(int a = 0; a<nAges; ++a){
            ageIndex(y,a) = ageIndex(y,a) + lengthIndexDetailed(y,l,i)* ALK_int(l,a,i,y);
            ageIndexDetailed(y,a,i) = ageIndexDetailed(y,a,i) +  lengthIndexDetailed(y,l,i)*ALK_int(l,a,i,y);
          }
        }
        lengthIndex(y,l) = lengthIndex(y,l) + lengthIndexDetailed(y,l,i);
        lengthIndexStrata(y,l, idxStrata(i)-1) = lengthIndexStrata(y,l, idxStrata(i)-1) + lengthIndexDetailed(y,l,i);
      }
    }
  }

  for(int y=0; y<nYears; ++y){
    for(int l = 0; l<numberOfLengthGroups; ++l){
      for(int s=0; s<nStrata; ++s){
        lengthIndexStrata(y,l,s) = lengthIndexStrata(y,l,s)*areas(s)/nIntStrata(s);
      }
      lengthIndex(y,l) = lengthIndex(y,l) *areas.sum()/nInt;
    }
    if(dat.applyALK==1){
      for(int a = 0; a<nAges; ++a){
        ageIndex(y,a) = ageIndex(y,a) *areas.sum()/nInt;
      }
    }
  }
  //--------------------------------------

  array<Type> logLengthIndex(nYears,numberOfLengthGroups);
  array<Type> logLengthIndexStrata(nYears,numberOfLengthGroups,nStrata);
  logLengthIndex= log(lengthIndex);
  logLengthIndexStrata = log(lengthIndexStrata);

  REPORT(lengthIndex);
  REPORT(lengthIndexStrata);

  if(dat.applyALK==0){//Do not adreport length if calculating index per length to reduce computation time.
    ADREPORT(logLengthIndex);
  }


  array<Type> logAgeIndex(nYears,nAges);
  if(dat.applyALK==1){
    logAgeIndex = log(ageIndex);
    ADREPORT(logAgeIndex);
  }


  //Report COG
  if(dat.applyALK==1){
    matrix<Type> COGAgeX(nYears,nAges);
    matrix<Type> COGAgeY(nYears,nAges);
    COGAgeX.setZero();
    COGAgeY.setZero();

    Type COGageTotal=0;

    for(int a = 0; a<nAges; ++a){
      for(int y = 0; y<nYears; ++y){
        COGageTotal=0;
        for(int i =0; i<nInt; ++i){
          COGAgeX(y,a) += ageIndexDetailed(y,a,i)*dat.xInt(i);
          COGAgeY(y,a) += ageIndexDetailed(y,a,i)*dat.yInt(i);
          COGageTotal += ageIndexDetailed(y,a,i);
        }
        COGAgeX(y,a) = COGAgeX(y,a)/COGageTotal;
        COGAgeY(y,a) = COGAgeY(y,a)/COGageTotal;
      }
    }

    REPORT(ageIndexDetailed);//For plotting
//    ADREPORT(COGAgeX);
//    ADREPORT(COGAgeY);
    REPORT(ALK_int);
  }

  //For plotting
  REPORT(lengthIndexDetailed);


  //Report covariates
  vector<Type> fourierReportLow = X_sunAltReport*betaSunLow;
//  vector<Type> fourierReportHigh = X_sunAltReport*betaSunHigh;
  ADREPORT(fourierReportLow);
//  ADREPORT(fourierReportHigh);

  vector<Type> depthReport1 =X_depthReport*parDepth1;
//  vector<Type> depthReport2 =X_depthReport*parDepth2;
  ADREPORT(depthReport1);
//  ADREPORT(depthReport2);

  return nll;
}
