##' Save indices with assosiated variance and covariance
##' @param  run stim object
##' @param  file name of file with index.
##' @param  folder Path to the folder the index is saved.
##' @details Saves the esitmated index with assosiated variance estimates and covariance structures in a format that can be directly utilized by SAM.
##' @export
saveIndex = function(run,file,folder = ""){

  rl = as.list(run$rep,what = "Est", report = TRUE)
  rlSd = as.list(run$rep, "Std", report = TRUE)

  nYears = length(run$conf_l$years)
  if(run$conf_l$applyALK==1){
    ageLogIndex = rl$logAgeIndex
    ageLogIndexSd = rlSd$logAgeIndex

    ageIndex = round(as.data.frame(exp(ageLogIndex)),3)
    ageLogIndexSd = round(as.data.frame(ageLogIndexSd),3)

    ageRange = run$data$ageRange + (run$conf_alk$maxAge- run$data$ageRange[2])
    rownames(ageIndex) = run$conf_l$years
    colnames(ageIndex) = ageRange[1]:ageRange[2]
    rownames(ageLogIndexSd) = run$conf_l$years
    colnames(ageLogIndexSd) = ageRange[1]:ageRange[2]


    #Extract yearly covariance matrices
    id = which(names(run$rep$value)=="logAgeIndex")
    cov = run$rep$cov[id,id]
    covYears = list()
    for(i in 1:length(run$conf_l$years)){
      id = i + (0:(run$data$ageRange[2]-1))*length(run$conf_l$years)
      covYears[[i]] = cov[id,id]
      rownames(covYears[[i]]) = ageRange[1]:ageRange[2]
      colnames(covYears[[i]]) = ageRange[1]:ageRange[2]
      covYears[[i]] = covYears[[i]]
    }

  }else{
    ageLogIndex = rl$logLengthIndex
    ageLogIndexSd = rlSd$logLengthIndex
    ageIndex = round(as.data.frame(exp(ageLogIndex)),3)
    ageLogIndexSd = round(as.data.frame(ageLogIndexSd),3)

    #Extract yearly covariance matrices
    id = which(names(run$rep$value)=="logLengthIndex")
    cov = run$rep$cov[id,id]

    covYears = list()
    for(i in 1:length(run$conf_l$years)){
      id = i + (0:(run$data$numberOfLengthGroups-1))*length(run$conf_l$years)
      covYears[[i]] = cov[id,id]
      covYears[[i]] = covYears[[i]]
    }
  }

  write.table(ageIndex, file = paste0(folder,file)) #Save index
  write.table(ageLogIndexSd, file = paste0(folder,"sd", file)) #Save sd of index

  fileNameCov = paste0(folder,"cov_",strsplit(file,'[.]')[[1]][1],".Rda")
  save(covYears,file = fileNameCov) #Can be directly utilized by SAM
}
