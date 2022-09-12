#' Initialie values for all model parameters and random effects.
#' @param data Data
#' @param conf Configurations
#' @details
#' @return a list containing initial values for all model parameters and random effects in the model.
#' @export
setPar <- function(data,conf){

  if(conf$lengthGroupsReduced[1] == conf$lengthGroupsReduced[2]){ #Latent effect when reducing length dimension
    xS = array(0.0, dim = c(dim(data$A_ListS[[1]])[2],NumberOfLengthGroups=max(conf$lengthGroupsReduced)+1))
    xST = array(0.0, dim = c(dim(data$A_ListST[[1]])[2],length(conf$years),NumberOfLengthGroups=max(conf$lengthGroupsReduced)+1))
  }else{
    xS = array(0.0, dim = c(dim(data$A_ListS[[1]])[2],NumberOfLengthGroups=max(conf$lengthGroupsReduced)))
    xST = array(0.0, dim = c(dim(data$A_ListST[[1]])[2],length(conf$years),NumberOfLengthGroups=max(conf$lengthGroupsReduced)))
  }

  nugget = array(0.0, dim = c(dim(data$fishObsMatrix)[2],dim(data$fishObsMatrix)[1]))


  nuggetIndex = array(0.0, dim = c(dim(data$fishObsMatrix)[2],length(data$xInt)))

  beta0 = array(0,dim=c(length(conf$years), length(conf$lengthGroups)))

  parameters <- list(beta0 = beta0,
                     betaSun =rep(0,max(1,conf$sunAlt[1])*4),
                     betaDepth = rep(0,data$Sdim*2),
                     log_lambda =c(5,5),
                     log_sigma =c(0,0,0),
                     log_kappa =c(-5,-5),
                     logSize =0.2,
                     tan_rho_t =0,
                     tan_rho_l =rep(0,3),
                     delta_z = rep(0,9),
                     xS = xS,
                     xST  = xST,
                     nugget = nugget,
                     nuggetIndex = nuggetIndex,
                     log_sigma_beta0 = 0)

  if(conf$applyALK==0){
    parameters = includeDummyPar(parameters)
  }
  return(parameters)
}

#' Includes dummy-parameters if ALK is not applied.
#' @param par Parameters
#' @details
#' @return ...
#' @export
includeDummyPar <- function(par){

  par$beta0_alk = matrix(rep(0,4),nrow = 2)
  par$log_sigma_beta0_alk = 0
  par$betaLength_alk = c(1,2,3)
  par$logSigma_alk = c(1,2,3)
  par$logKappa_alk = c(1,2,3)
  par$transRho_alk = c(1,2,3)
  par$xST_alk = array(0.0, dim = c(2,2))

  return(par)
}


