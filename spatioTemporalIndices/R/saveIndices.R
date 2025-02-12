##' Save indices with assosiated variance and covariance
##' @param  run stim object
##' @param  file name of file with index.
##' @param  folder Path to the folder the index is saved.
##' @param  sdrep_bc Output of TMB::sdreport(..., bias.correct = TRUE)
##' @details Saves the esitmated index with assosiated variance estimates and covariance structures in a format that can be directly utilized by SAM.
##' @importFrom utils write.table
##' @export
saveIndex = function(run,file,folder = "",sdrep_bc = NULL){

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




##' Save indices in standard ICES format
##' @param  run stim object
##' @param  name name of index.
##' @param  file Name of file. Must end with .dat
##' @param  variance Save variances of log-indices?
##' @param  sdrep_bc sdreport from TMB with bias.correct
##' @param  digits number of decimal places in the output file
##' @details Saves the estimated index or associated variance of log-index in standard ICES format
##' @importFrom utils write.table
##' @export
write_indices_ICES_format = function(run,file,name,sdrep_bc = NULL,variance = FALSE, digits = 1){

  if(run$conf_l$applyALK==0) stop("Not yet implemented for indices-at-length. See run$rl and run$rlSd")

  if(!is.null(sdrep_bc)) stop("Writing posterior mean estiamtes on standard ICES format not yet implemented yet")

  if(variance){
    name = paste0(name, " (Variance of log-observations)")
  }

  time_in_year <- as.numeric(format(run$dat_l$startdatetime, "%j"))/365
  time_survey = quantile(time_in_year, c(0.25,0.75))
  header_lines <- c(name,
                    paste(min(run$conf_l$years), max(run$conf_l$years)),
                    paste("1 1", round(time_survey[1],2), round(time_survey[2],2)),
                    paste(run$conf_alk$minAge, run$conf_alk$maxAge))

  writeLines(header_lines, file)

  indices = round(exp(run$rl$logAgeIndex)[,-1],digits) #Remove youngest age that we needed when calculating the ALK
  logIndicesVar = round((run$rlSd$logAgeIndex)[,-1]^2,digits)
  logIndicesVar = signif(logIndicesVar,digits = digits)

  if(variance){
    write.table(logIndicesVar, file, append = TRUE, col.names = FALSE,
                quote = FALSE, sep = " ", eol = "\n", dec = ".",
                row.names = rep(1, nrow(indices)))
  }else{
    write.table(indices, file, append = TRUE, col.names = FALSE,
                quote = FALSE, sep = " ", eol = "\n", dec = ".",
                row.names = rep(1, nrow(indices)))
  }

}

##' Save indices in standard ICES format
##' @param  run stim object
##' @param  file Name of the file. Must end with .rds
##' @details Saves a list with estimated yearly covariance matrices of the log-index.
##' @export
write_covriance_matrices = function(run, file){
  if(run$conf_l$applyALK==1){
    ageRange = run$data$ageRange + (run$conf_alk$maxAge- run$data$ageRange[2])
    id = which(names(run$rep$value)=="logAgeIndex")
    cov = run$rep$cov[id,id]
    covYears = list()
    for(i in 1:length(run$conf_l$years)){
      id = i + (0:(run$data$ageRange[2]-1))*length(run$conf_l$years)
      covYears[[i]] = cov[id,id]
      rownames(covYears[[i]]) = ageRange[1]:ageRange[2]
      colnames(covYears[[i]]) = ageRange[1]:ageRange[2]
      covYears[[i]] = covYears[[i]][-1,]
      covYears[[i]] = covYears[[i]][,-1]
    }
  }else{
    id = which(names(run$rep$value)=="logLengthIndex")
    cov = run$rep$cov[id,id]

    covYears = list()
    for(i in 1:length(run$conf_l$years)){
      id = i + (0:(run$data$numberOfLengthGroups-1))*length(run$conf_l$years)
      covYears[[i]] = cov[id,id]
      covYears[[i]] = covYears[[i]]
    }
  }

  names(covYears) = run$conf_l$years
  saveRDS(covYears,file = file)
}
