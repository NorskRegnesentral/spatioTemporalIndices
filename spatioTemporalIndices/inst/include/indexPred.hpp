using namespace Eigen; //Utilize sparse structures

template <class Type>
  void indexPred(dataSet<Type> dat, paraSet<Type> par, objective_function<Type> *of){

    array<Type> ALK_int;
    if(dat.applyALK==1){
      ALK_int = ALK(dat,par); //ALKs at for integration points
    }
    int numberOfLengthGroups =dat.numberOfLengthGroups;

    //Calculate indices-----------------------------------
    int nYears = dat.nStationsEachYear.size();
    int nInt = dat.xInt.size();
    vector<Type> sigma = exp(par.log_sigma);
    vector<Type> kappa = exp(par.log_kappa);

    Type scaleS = Type(1)/((4*M_PI)*kappa(0)*kappa(0)); //No effect on results, but needed for interpreting the sigma^2 parameter as marginal variance. See section 2.1 in Lindgren (2011)
    Type scaleST = Type(1)/((4*M_PI)*kappa(1)*kappa(1)); //No effect on results, but needed for interpreting the sigma^2 parameter as marginal variance

    array<Type>  lengthIndex(nYears,numberOfLengthGroups);
    array<Type>  lengthIndexDetailed(nYears,numberOfLengthGroups, nInt);
    vector<Type>  lengthIndexTotal(nYears);
    lengthIndex.setZero();
    lengthIndexDetailed.setZero();
    lengthIndexTotal.setZero();

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

    vector<Type> depthEffectInt1=dat.X_depth_int*parDepth1;
    vector<Type> depthEffectInt2=dat.X_depth_int*parDepth2;
    vector<Type> timeInDayEffectIntLow = dat.X_sunAltIntegrate*betaSunLow;
    vector<Type> timeInDayEffectIntHigh = dat.X_sunAltIntegrate*betaSunHigh;

    for(int y=0; y<nYears; ++y){
      for(int l =0; l<numberOfLengthGroups; ++l){
        Type covariatesConvexW = (numberOfLengthGroups-l-1)/(numberOfLengthGroups-1);
        if(dat.lengthGroupsReduced(0)==dat.lengthGroupsReduced(1)){
          deltaPredS = dat.weigthLength(l)*dat.ApredS * par.xS.col(dat.lengthGroupsReduced(l)).matrix()+
            (1-dat.weigthLength(l))*dat.ApredS * par.xS.col(dat.lengthGroupsReduced(l)+1).matrix();
          deltaPredST = dat.weigthLength(l)*dat.ApredST * par.xST.col(dat.lengthGroupsReduced(l)).col(y).matrix()+
            (1-dat.weigthLength(l))*dat.ApredST * par.xST.col(dat.lengthGroupsReduced(l)+1).col(y).matrix();
        }else{
          deltaPredS = dat.ApredS * par.xS.col(dat.lengthGroupsReduced(l)).matrix();
          deltaPredST = dat.ApredST * par.xST.col(dat.lengthGroupsReduced(l)).col(y).matrix();
        }
        for(int i =0; i<nInt; ++i){
          lengthIndexDetailed(y,l,i) =  exp(par.beta0.row(y)(l) +
            covariatesConvexW*timeInDayEffectIntLow(0) + (1-covariatesConvexW)*timeInDayEffectIntHigh(0)+
            covariatesConvexW*depthEffectInt1(i) + (1-covariatesConvexW)*depthEffectInt2(i)+
            deltaPredS(i)/sqrt(scaleS)*sigma(0)+
            deltaPredST(i)/sqrt(scaleST)*sigma(1));
          lengthIndexDetailed(y,l,i) = lengthIndexDetailed(y,l,i)*dat.trawlWidth;
        }
      }
    }

//Zero inflation not fully implemented
//    Type  muZero;
//    Type  pZero;
//    if(dat.zeroInflated !=0){
//      for(int y=0; y<nYears; ++y){
//        for(int l =0; l<numberOfLengthGroups; ++l){
//          for(int i =0; i<nInt; ++i){
//            muZero = exp(par.delta_z(0) +
//              par.delta_z(1)*log(lengthIndexDetailed(y,l,i)));
//            pZero = dpois(Type(0), muZero);
//            lengthIndexDetailed(y,l,i) = lengthIndexDetailed(y,l,i)*(1-pZero);
//          }
//        }
//      }
//    }


    int nAges = 2;//Dummy-number
    if(dat.applyALK==1){
      nAges = dat.ageRange(1) - dat.ageRange(0) + 1;
    }
    array<Type> ageIndex(nYears,nAges);
    array<Type> ageIndexDetailed(nYears,nAges, nInt);
    vector<Type>  ageIndexTotal(nYears);
    ageIndex.setZero();
    ageIndexDetailed.setZero();
    ageIndexTotal.setZero();

    int nStrata = 0;
    for(int i=0; i<dat.idxStrata.size(); ++i){
      if(nStrata<dat.idxStrata(i)){
        nStrata = dat.idxStrata(i);
      }
    }
    array<Type> lengthIndexStrata(nYears,dat.numberOfLengthGroups,nStrata);
    lengthIndexStrata.setZero();

    vector<Type> nIntStrata(nStrata);
    nIntStrata.setZero();
    for(int i=0; i<nInt; ++i){
      nIntStrata(dat.idxStrata(i)-1) = nIntStrata(dat.idxStrata(i)-1) +1;
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
          lengthIndexStrata(y,l, dat.idxStrata(i)-1) = lengthIndexStrata(y,l, dat.idxStrata(i)-1) + lengthIndexDetailed(y,l,i);
        }
      }
    }

    for(int y=0; y<nYears; ++y){
      for(int l = 0; l<numberOfLengthGroups; ++l){
        for(int s=0; s<nStrata; ++s){
          lengthIndexStrata(y,l,s) = lengthIndexStrata(y,l,s)*dat.areas(s)/nIntStrata(s);
        }
        lengthIndex(y,l) = lengthIndex(y,l) *dat.areas.sum()/nInt;
        lengthIndexTotal(y) = lengthIndexTotal(y) + lengthIndex(y,l);
      }
      if(dat.applyALK==1){
        for(int a = 0; a<nAges; ++a){
          ageIndex(y,a) = ageIndex(y,a) *dat.areas.sum()/nInt;
          ageIndexTotal(y) = ageIndexTotal(y) + ageIndex(y,a);
        }
      }
    }
    //--------------------------------------

    array<Type> logLengthIndex(nYears,numberOfLengthGroups);
    array<Type> logLengthIndexStrata(nYears,numberOfLengthGroups,nStrata);
    vector<Type>  logLengthIndexTotal(nYears);
    logLengthIndex= log(lengthIndex);
    logLengthIndexStrata = log(lengthIndexStrata);
    logLengthIndexTotal = log(lengthIndexTotal);

    REPORT_F(logLengthIndex, of);
    REPORT_F(logLengthIndexStrata, of);
    REPORT_F(logLengthIndexTotal, of);
    REPORT_F(lengthIndexDetailed, of);//For plotting

    if((dat.strataReport==1) & (dat.applyALK==0)){
      ADREPORT_F(logLengthIndexStrata, of);
    }

    if(dat.applyALK==0){//Do not adreport length if calculating index per length to reduce computation time.
      ADREPORT_F(logLengthIndex, of);
      if(dat.lowMemory==0){
        ADREPORT_F(logLengthIndexTotal, of);
      }else{
        REPORT_F(logLengthIndexTotal, of);
      }
    }

    array<Type> logAgeIndex(nYears,nAges);
    vector<Type>  logAgeIndexTotal(nYears);
    if(dat.applyALK==1){
      logAgeIndex = log(ageIndex);
      logAgeIndexTotal = log(ageIndexTotal);
      ADREPORT_F(logAgeIndex, of);
      REPORT_F(logAgeIndex, of);
      if(dat.lowMemory==0){
        ADREPORT_F(logAgeIndexTotal, of);
      }else{
        REPORT_F(logAgeIndexTotal, of);
      }
      REPORT_F(ageIndexDetailed, of);//For plotting
      REPORT_F(ALK_int, of);
    }

    if(dat.lowMemory==0){
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
        //ADREPORT_F(COGAgeX, of);
        //ADREPORT_F(COGAgeY, of);
      }

    //Report covariates
      vector<Type> fourierReportLow = dat.X_sunAltReport*betaSunLow;
      ADREPORT_F(fourierReportLow, of);

      if(dat.sunEffect==2){
        vector<Type> fourierReportHigh = dat.X_sunAltReport*betaSunHigh;
        ADREPORT_F(fourierReportHigh, of);
      }

      vector<Type> depthReport1 =dat.X_depthReport*parDepth1;
      ADREPORT_F(depthReport1, of);
      if(dat.splineDepth ==2){ //Two length dependent splines
        vector<Type> depthReport2 =dat.X_depthReport*parDepth2;
        ADREPORT_F(depthReport2, of);
      }
    }
  }
