using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()
using namespace density; //use GMRF
using namespace Eigen; //Utilize sparse structures

template <class Type>
  Type nllIndex(dataSet<Type> dat, paraSet<Type> par, spde_t<Type> spdeMatricesS, spde_t<Type> spdeMatricesST ,
                LOSM_t<Type> A_ListS, LOSM_t<Type> A_ListST, data_indicator<vector<Type>,Type> keep, objective_function<Type> *of){

    Type nll = 0;

    if(dat.rwBeta0 ==1){
      Type sigma_beta0 = exp(par.log_sigma_beta0(0));
      for(int y=1;y<dat.nStationsEachYear.size();y++){
        for(int l=0; l<dat.numberOfLengthGroups; ++l){
          nll -= dnorm(par.beta0(y,l),par.beta0(y-1,l),sigma_beta0,true);
        }
      }
    }else{
      //The intercepts are free model parameters
    }
    //Transform parameters
    vector<Type> sigma = exp(par.log_sigma);
    vector<Type> kappa = exp(par.log_kappa);
//    Type size= exp(par.logSize);

//Only poisson distributed catch implemented currently
//    Type scale; //If apply dgamma
//    Type shape; //If apply dgamma
//    Type pTweedie;

    vector<Type> rho_t(1);//Time correlation parameter
    rho_t(0)=Type(2)/(Type(1) + exp(-Type(2) * par.tan_rho_t(0))) - Type(1);

    vector<Type> rho_l(3);//Length correlation parameter
    rho_l(0) = Type(2)/(Type(1) + exp(-Type(2) * par.tan_rho_l(0))) - Type(1);
    rho_l(1) = Type(2)/(Type(1) + exp(-Type(2) * par.tan_rho_l(1))) - Type(1);
    rho_l(2) = Type(2)/(Type(1) + exp(-Type(2) * par.tan_rho_l(2))) - Type(1);

    //Construct Q in spatial dimension
    SparseMatrix<Type> Q_s = Q_spde(spdeMatricesS,kappa(0));
    SparseMatrix<Type> Q_st = Q_spde(spdeMatricesST,kappa(1));

    //Latent spatio-terporal effects
    Type scaleS = Type(1)/((4*M_PI)*kappa[0]*kappa[0]); //No effect on results, but needed for interpreting the sigma^2 parameter as marginal variance. See section 2.1 in Lindgren (2011)
    Type scaleST = Type(1)/((4*M_PI)*kappa(1)*kappa(1)); //No effect on results, but needed for interpreting the sigma^2 parameter as marginal variance
    Type d = 2; //Part of spatial pc-prior
    Type rhoP;
    Type R = -log(dat.pcPriorsRange(1))*pow(dat.pcPriorsRange(0),d/2);
    Type S = -log(dat.pcPriorsSD(1))/dat.pcPriorsSD(0);
    if(dat.spatial==1){
      if(dat.usePCpriors==1){
        rhoP = sqrt(8)/kappa(0);
        nll -= log( d/2 * R *S * pow(rhoP,(-1-d/2))* exp(-R* pow(rhoP,(-d/2)) -S* sigma(0)) ); //pc-prior contribution
      }
      nll += SEPARABLE(AR1(rho_l(0)),GMRF(Q_s))(par.xS);
      if(dat.simulateProcedure==1){
        SIMULATE_F(of){
          SEPARABLE(AR1(rho_l(0)),GMRF(Q_s)).simulate(par.xS);
        }
      }
    }
    if(dat.spatioTemporal != 0){
      if(dat.usePCpriors==1){
        rhoP = sqrt(8)/kappa(1);
        nll -= log( d/2 * R *S * pow(rhoP,(-1-d/2))* exp(-R* pow(rhoP,(-d/2)) -S* sigma(1)) ); //pc-prior contribution
      }
      nll += SEPARABLE(AR1(rho_l(1)),SEPARABLE(AR1(rho_t(0)),GMRF(Q_st)))(par.xST);
      if(dat.simulateProcedure==1){
        SIMULATE_F(of){
            SEPARABLE(AR1(rho_l(1)),SEPARABLE(AR1(rho_t(0)),GMRF(Q_st))).simulate(par.xST);
        }
      }
    }

    if(dat.useNugget==1){
      int nHaul = dat.nStationsEachYear.sum();
      SparseMatrix<Type> Q_nuggetIID(nHaul,nHaul);
      for(int i = 0; i< nHaul; ++i){
        Q_nuggetIID.coeffRef(i,i)=1;
      }
      nll += SEPARABLE(GMRF(Q_nuggetIID),AR1(rho_l(2)))(par.nugget);
      SIMULATE_F(of){
          SEPARABLE(GMRF(Q_nuggetIID),AR1(rho_l(2))).simulate(par.nugget);
      }
    }

    //p-spline
    vector<Type> lambda = exp(par.log_lambda);
    int nSplineDepth = par.betaDepth.size()/2;
    vector<Type> parDepth1(nSplineDepth);
    vector<Type> parDepth2(nSplineDepth);
    for(int i =0; i<nSplineDepth; ++i){
      parDepth1(i) = par.betaDepth(i);
      parDepth2(i) = par.betaDepth(i + nSplineDepth);
    }
    if(dat.splineDepth ==1){
      nll -= Type(0.5)*dat.Sdim*par.log_lambda(0) - 0.5*lambda(0)*GMRF(dat.S_depth).Quadform(parDepth1);
    }else if(dat.splineDepth ==2){ //Two length dependent splines
      nll -= Type(0.5)*dat.Sdim*par.log_lambda(0) - 0.5*lambda(0)*GMRF(dat.S_depth).Quadform(parDepth1);
      nll -= Type(0.5)*dat.Sdim*par.log_lambda(1) - 0.5*lambda(1)*GMRF(dat.S_depth).Quadform(parDepth2);
    }

    //fourier
    vector<Type> betaSunLow(dat.nBasisSunAlt*2);
    vector<Type> betaSunHigh(dat.nBasisSunAlt*2);
    for(int i =0; i<(dat.nBasisSunAlt*2); ++i){
      betaSunLow(i) = par.betaSun(i);
      betaSunHigh(i) = par.betaSun(i + dat.nBasisSunAlt*2);
    }


//    Type log_var_minus_mu; //Needed if applying negative binomial observation model
    Type covariatesConvexW; //Coefficient in convex combination of length dependent effects
    matrix<Type> mu(dat.fishObsMatrix.rows(),dat.fishObsMatrix.cols());
    vector<Type> deltaS, deltaS2, deltaST, deltaST2; //Smoothing in length dimension
    matrix<Type> deltaMatrixS(dat.numberOfLengthGroups,999), deltaMatrixST(dat.numberOfLengthGroups,999);//Smoothed effect in length dimension
    vector<Type> timeInDayLow = dat.X_sunAlt*betaSunLow;
    vector<Type> timeInDayHigh = dat.X_sunAlt*betaSunHigh;
    vector<Type> depthEffect1 = dat.X_depth*parDepth1;
    vector<Type> depthEffect2 = dat.X_depth*parDepth2;

    int counter=0;
    for(int y=0;y<dat.nStationsEachYear.size();y++){
      SparseMatrix<Type> As = A_ListS(y);
      SparseMatrix<Type> Ast = A_ListST(y);
      //Latent effects contribution
      for(int l=0; l <dat.numberOfLengthGroups;++l){
        deltaS = As * par.xS.col(dat.lengthGroupsReduced(l)).matrix();
        deltaST = Ast * par.xST.col(dat.lengthGroupsReduced(l)).col(y).matrix();
        if((dat.lengthGroupsReduced(0)==dat.lengthGroupsReduced(1)) & (dat.weigthLength(l)<0.9999)){ //Length dimension is reduced and there exists longer length in dat.lengthGroupsReduced
          deltaS2 = As * par.xS.col(dat.lengthGroupsReduced(l)+1).matrix();
          deltaST2 = Ast * par.xST.col(dat.lengthGroupsReduced(l)+1).col(y).matrix();
          for(int s=0; s<dat.nStationsEachYear(y);++s){
            deltaMatrixS(l,s) = dat.weigthLength(l)*deltaS(s) + (1-dat.weigthLength(l))*deltaS2(s);
            deltaMatrixST(l,s) = dat.weigthLength(l)*deltaST(s) + (1-dat.weigthLength(l))*deltaST2(s);
          }
        }else{
          for(int s=0; s<dat.nStationsEachYear(y);++s){
            deltaMatrixS(l,s) = deltaS(s);
            deltaMatrixST(l,s) = deltaST(s);
          }
        }
      }
//      Type muZero = 0;//Used if applying zero inflation
//      Type pZero = 0;//Used if applying zero inflation
//      Type pPos = 0;//Used if applying zero inflation
      for(int s=0; s<dat.nStationsEachYear(y);++s){
        for(int l=0; l <dat.numberOfLengthGroups;++l){
          covariatesConvexW = (dat.numberOfLengthGroups-l-1)/(dat.numberOfLengthGroups-1);
          mu(counter,l)= exp( par.beta0.row(y)(l)+
            covariatesConvexW*timeInDayLow(counter) + (1-covariatesConvexW)*timeInDayHigh(counter) +
            deltaMatrixS.row(l)(s)*sigma(0)/sqrt(scaleS)+
            deltaMatrixST.row(l)(s)*sigma(1)/sqrt(scaleST)+
            log(dat.dist(counter))+
            covariatesConvexW*depthEffect1(counter) + (1-covariatesConvexW)*depthEffect2(counter)+
            par.nugget.col(counter)(l)*sigma(2));
          if(dat.predMatrix(counter,l)==0){
            if(dat.zeroInflated !=0){//Include additional zero-probability //Under development
//              muZero = exp(par.delta_z(0) +
//                par.delta_z(1)* par.beta0.row(y)(l) +
//                par.delta_z(1)* (covariatesConvexW*timeInDayLow(counter) + (1-covariatesConvexW)*timeInDayHigh(counter)) +
//                par.delta_z(1)* deltaMatrixS.row(l)(s)*sigma(0)/sqrt(scaleS)+
//                par.delta_z(1)* deltaMatrixST.row(l)(s)*sigma(1)/sqrt(scaleST)+
//                par.delta_z(1)* (covariatesConvexW*depthEffect1(counter) + (1-covariatesConvexW)*depthEffect2(counter))+
//                par.delta_z(1)* par.nugget.col(counter)(l)*sigma(2)+
//                par.delta_z(2)* log(dat.dist(counter)));
//              pZero = dpois(Type(0), muZero,true);
//              if(dat.obsModel==1){
//                log_var_minus_mu=log(mu(counter,l)*mu(counter,l)*size);
//                pPos = dnbinom_robust(dat.obsVector(dat.idxStart(counter) +l), log(mu(counter,l)),log_var_minus_mu,true)  + logspace_sub(Type(0),pZero);
//              }else if (dat.obsModel==2){
//                pPos = dpois(dat.obsVector(dat.idxStart(counter) +l), mu(counter,l),true)  + logspace_sub(Type(0),pZero);
//              }else if (dat.obsModel==3){
//                //Gamma distributed response
//                if(dat.obsVector(dat.idxStart(counter) +l) >0){
//                  scale = size;
//                  shape = mu(counter,l)/scale;
//                  pPos = dgamma(dat.obsVector(dat.idxStart(counter) +l), shape,scale,true)  + logspace_sub(Type(0),pZero);
//                }else{
//                  //Not needed, gamma cant be exactly zero
//                }
//              }else{
//                //Not implemented
//                pPos = 0;
//              }
//
//              if(dat.fishObsMatrix(counter,l)==0){
//                if(dat.obsModel==3){
//                  nll -=keep(dat.idxStart(counter) +l)*pZero;
//                }else{
//                  nll -=keep(dat.idxStart(counter) +l)*logspace_add(pZero,pPos);
//                  }
//              }else{
//                nll -=keep(dat.idxStart(counter) +l)*pPos;
//              }
            }else{
              if(dat.obsModel==1){
//                nll -= keep(dat.idxStart(counter) +l)*dnbinom_robust(dat.obsVector(dat.idxStart(counter) +l), log(mu(counter,l)),log_var_minus_mu,true);
                //TODO simulation
              }else if(dat.obsModel==2){
                //Currenlty only obsModel==2 implemented
                nll -= keep(dat.idxStart(counter) +l)*dpois(dat.obsVector(dat.idxStart(counter) +l), mu(counter,l),true);
                SIMULATE_F(of){
                    dat.obsVector(dat.idxStart(counter) +l) = rpois(mu(counter,l));
                }
              }else if (dat.obsModel==3){
//                //Tweedie distributed response
//                scale = size;
//                shape = mu(counter,l)/scale;
//                pTweedie =  Type(1)/(Type(1) + exp(-Type(2) * par.tweedieP(0))) + Type(1);
//                //TODO simulation
//                nll -= keep(dat.idxStart(counter) +l)*dtweedie(dat.obsVector(dat.idxStart(counter) +l), mu(counter,l),size,pTweedie,true);
              }
            }
          }
        }
        counter++;
      }
    }
    //---------------------------------------------

    SIMULATE_F(of){
      vector<Type> obsVector=dat.obsVector;
      REPORT_F(obsVector,of);
      vector<Type> nugget=par.nugget;
      REPORT_F(nugget,of);
      array<Type> xS=par.xS;
      array<Type> xST=par.xST;
      REPORT_F(xS,of);
      REPORT_F(xST,of);
    }

    return(nll);
  }
