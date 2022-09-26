

#define ADREPORT_F(name,F) F->reportvector.push(name,#name);

#define REPORT_F(name,F)					                               \
if(isDouble<Type>::value && F->current_parallel_region<0) {	\
  Rf_defineVar(Rf_install(#name),					                      \
               PROTECT(asSEXP(name)),F->report);			         \
  UNPROTECT(1);						                                       \
}



/* List of sparse matrices */
template<class Type>
struct LOSM_t : vector<SparseMatrix<Type> > {
  LOSM_t(SEXP x){  /* x = List passed from R */
(*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      if(!isValidSparseMatrix(sm))
        error("Not a sparse matrix");
      (*this)(i) = asSparseMatrix<Type>(sm);
    }
  }
};


template <class Type>
struct dataSet{
  vector<int> age;
  vector<int> ageNotTruncated;
  vector<Type> length;
  vector<int> readability;
  vector<int> ageRange;
  vector<int> idx1;
  vector<int> idx2;
  matrix<Type> fishObsMatrix;
  vector<Type> dist;
  vector<int> nStationsEachYear;
  int numberOfLengthGroups;
  SparseMatrix<Type> S_depth;
  int Sdim;
  vector<int> lengthGroupsReduced;
  vector<Type> lengthGroups;
  vector<Type> dL;
  matrix<Type> predMatrix;
  matrix<Type> X_sunAlt;
  matrix<Type> X_sunAltReport;
  matrix<Type> X_sunAltIntegrate;
  matrix<Type> X_depth;
  matrix<Type> X_depthReport;
  matrix<Type> X_depth_int;
  vector<Type> weigthLength;
  int zeroInflated;
  int nBasisSunAlt;
  int obsModel;
  vector<Type> obsVector;
  vector<int> idxStart;
  vector<Type> xInt;
  vector<Type> yInt;
  SparseMatrix<Type> ApredS;
  SparseMatrix<Type> ApredST;
  SparseMatrix<Type> Apred_alk;
  int spatial;
  int spatioTemporal;
  int splineDepth;
  int useNugget;
  vector<Type> pcPriorsRange;
  vector<Type> pcPriorsSD;
  int usePCpriors;
  int usePCpriorsALK;
  vector<Type> pcPriorsALKRange;
  vector<Type> pcPriorsALKSD;
  int simulateProcedure;
  int rwBeta0;
  int maxAge;
  int minAge;
  int applyALK;
  vector<int> idxStrata;
  vector<Type> areas;
  int spatioTemporalALK;
  int spatialALK;


  int rwBeta0_alk;

};


template <class Type>
struct paraSet{
  matrix<Type> beta0_alk;
  vector<Type> log_sigma_beta0_alk;
  vector<Type> betaLength_alk;
  vector<Type> logSigma_alk;
  vector<Type> logKappa_alk;
  vector<Type> transRho_alk;
  array<Type> xST_alk;
  array<Type> xS_alk;

  matrix<Type> beta0;
  vector<Type> log_sigma_beta0;
  vector<Type> betaSun;
  vector<Type> betaDepth;
  vector<Type> log_lambda;
  vector<Type> log_sigma;
  vector<Type> log_kappa;
  Type logSize;
  vector<Type> tan_rho_t;
  vector<Type> tan_rho_l;
  vector<Type> delta_z;
  array<Type> xS;
  array<Type> xST;
  array<Type> nugget;

};
