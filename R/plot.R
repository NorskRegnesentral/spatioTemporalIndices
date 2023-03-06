

#' plotResults Plot results
#' @param run fitted object returned by \code{\link{fitModel}}
#' @param what what to plot, options in first element "sunAlt", "depth", "S", "ST" and "SST". If spatial effect is plotted, year and length must be provided
#'
#' @export
plotResults  <- function(run,what=NULL, legend = FALSE){
  pl = as.list(run$rep,what = "Est")
  if(what[1] == "sunAlt"){
    plotSunAlt(run)
  }else if(what[1] == "depth"){

    depthSpline = run$rep$value[names(run$rep$value)=="depthReport1"]
    sdDepthSpline<-run$rep$sd[names(run$rep$value)=="depthReport1"]
    depth=seq(from=min(attributes(run$data)$depth),to=max(attributes(run$data)$depth),length.out = length(depthSpline))

    plot(depth, depthSpline, lty=1,type = 'l',ylim = c(min(depthSpline - 1.96*sdDepthSpline)-0.2,max(depthSpline + 1.96*sdDepthSpline)+0.2),
         ylab = "Effect in linear predictor",main = "Depth effect",col="red",xlab = "Depth (m)",
         cex.main = 1.7,cex.lab = 1.5, cex = 1.5)
    lines(depth, depthSpline - 1.96*sdDepthSpline, lty=2,col="red")
    lines(depth, depthSpline + 1.96*sdDepthSpline, lty=2,col="red")

    if(run$conf_l$splineDepth[2]==2){
      depthSpline2 = run$rep$value[names(run$rep$value)=="depthReport2"]
      sdDepthSpline2<-run$rep$sd[names(run$rep$value)=="depthReport2"]
      lines(depth, depthSpline2,col = 'blue')
      lines(depth, depthSpline2 - 1.96*sdDepthSpline2, lty=2,col="blue")
      lines(depth, depthSpline2 + 1.96*sdDepthSpline2, lty=2,col="blue")

      minLength = min(run$conf_l$lengthGroups)
      maxLength = max(run$conf_l$lengthGroups)
      legend(legend=c(paste0(minLength, " cm"), paste0(maxLength, " cm")),col=c("red","blue"),"topright",lty=1,cex =1.5)
    }

    abline(h = 0)
    abline(v=100,lty=3)
    abline(v=200,lty=3)
    abline(v=300,lty=3)
    abline(v=400,lty=3)
  }else if(what[1] == "S"){
    mesh = attributes(run$data)$meshS
    randomEffects=run$rep$par.random[names(run$rep$par.random)=="xS"]
    l = run$conf_l$lengthGroupsReduced[which(run$conf_l$lengthGroups == as.numeric(what[3]))]
    indexStart=mesh$n *(l-1)+1
    indexEnd=indexStart+mesh$n-1
    indexStart2 = mesh$n *l+1
    indexEnd2=indexStart2+mesh$n-1

    kappa = exp(run$pl$log_kappa)[1]
    scalingConstant = exp(run$rep$par.fixed[which(names(run$rep$par.fixed)=="log_sigma")])[1]
    scalingConstant2 = sqrt(1/((4*3.14159265)*kappa*kappa))
    scale = scalingConstant/scalingConstant2

    lD = which(run$conf_l$lengthGroups == as.numeric(what[3]))
    spatialE = (run$data$weigthLength[lD]*randomEffects[indexStart:indexEnd] + (1-run$data$weigthLength[lD])*randomEffects[indexStart2:indexEnd2] )*scale

    year = as.numeric(what[2])
    yearPosition = year-min(run$conf_l$years)+1
    beta0 = summary(run$rep)[which(rownames(summary(run$rep))=="beta0")]
    ll = which(run$conf_l$lengthGroups==as.numeric(what[3]))
    beta0This = beta0[(ll-1)*length(run$conf_l$years) + yearPosition]

    proj = inla.mesh.projector(mesh)
    latentFieldMAP = beta0This + spatialE

    image(proj$x,proj$y, inla.mesh.project(proj, latentFieldMAP),col =  colorRampPalette(c("white","yellow", "red"))(12),
               xlab = '', ylab = "",
#               zlim = c(range[1],range[2]),
               main = paste0("Spatial effect for length group: ", what[3]),
               cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1,xaxt='n',yaxt='n')

    contour(proj$x, proj$y,inla.mesh.project(proj, latentFieldMAP) ,add = T,labcex  = 1,cex = 1)
    points(attributes(run$data)$locObs[attributes(run$data)$year == as.numeric(what[2]),],cex = 0.01,col = 'blue')

    #Convert map to UTM coordinates-------------------------------------------------------------------------
    newmap <- maps::map("world", c("Norway","Sweden","Finland","Russia"),fill = TRUE,plot = FALSE, col = "transparent")
    mapTmp = data.frame(newmap$x,newmap$y)
    mapTmp[which(is.na(mapTmp[,1])),] = 3.141592 #Need something different from NA
    names(mapTmp) = c("X","Y")
    attr(mapTmp, "projection") = "LL"
    attr(mapTmp, "zone") = conf_l$zone
    ddpcr::quiet(mapTmp <- PBSmapping::convUL(mapTmp))
    colnames(mapTmp) = c("UTMX", "UTMY")
    mapTmp[which(is.na(newmap$x)),] = NA
    polygon(mapTmp,col = 'lightgrey')
    #----------------------------------------------------------------------------------------------------

    if(legend == TRUE){
      image.plot(proj$x,proj$y, inla.mesh.project(proj, latentFieldMAP),col =  colorRampPalette(c("white","yellow", "red"))(12),
                 add = TRUE,legend.width = 5,legend.only = TRUE, cex = 1.6,
                 zlim = c(range[1],range[2]), axis.args = list(cex.axis = 1.3))
    }


  }else if(what[1] == "ST"){

    mesh = attributes(run$data)$meshST
    randomEffects=run$rep$par.random[names(run$rep$par.random)=="xST"]
    year = as.numeric(what[2])
    yearPosition = year-min(run$conf_l$years)+1
    l = run$conf_l$lengthGroupsReduced[which(run$conf_l$lengthGroups == as.numeric(what[3]))]
    indexStart=length(run$conf_l$years)*mesh$n *(l-1)+(yearPosition-1)*mesh$n+1
    indexEnd=indexStart+mesh$n-1

    indexStart2=length(run$conf_l$years)*mesh$n *(l)+(yearPosition-1)*mesh$n+1
    indexEnd2=indexStart2+mesh$n-1

    kappa = exp(run$pl$log_kappa)[2]
    scalingConstant = exp(run$rep$par.fixed[which(names(run$rep$par.fixed)=="log_sigma")])[2]
    scalingConstant2 = sqrt(1/((4*3.14159265)*kappa*kappa))
    scale = scalingConstant/scalingConstant2

    lD = which(run$conf_l$lengthGroups == as.numeric(what[3]))
    ST = (run$data$weigthLength[lD]*randomEffects[indexStart:indexEnd] + (1-run$data$weigthLength[lD])*randomEffects[indexStart2:indexEnd2] )*scale

    beta0 = summary(run$rep)[which(rownames(summary(run$rep))=="beta0")]
    ll = which(run$conf_l$lengthGroups==as.numeric(what[3]))
    beta0This = beta0[(ll-1)*length(run$conf_l$years) + yearPosition]
    proj = inla.mesh.projector(mesh)
    latentFieldMAP = beta0This + ST

    image(proj$x,proj$y, inla.mesh.project(proj, latentFieldMAP),
          col =  colorRampPalette(c("white","yellow", "red"))(12),
#          xlab = '', ylab = "",
        #  zlim = c(range[1],range[2]),
          main = "",
          xlim = c(floor(min(proj$x)),ceiling(max(proj$x))),
          ylim = c(floor(min(proj$y)),ceiling(max(proj$y))), # offset removed
          cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1,xaxt = 'n',yaxt = 'n' )
#    if(addX)axis(1, at=seq(0,10000,by = 300), labels=seq(0,10000,by = 300), col.axis="black", las=2)
#    if(addY)axis(2, at=seq(0,10000,by = 300), labels=seq(0,10000,by = 300), col.axis="black", las=2)
    contour(proj$x, proj$y,inla.mesh.project(proj, latentFieldMAP) ,add = T,labcex  = 1,cex = 1)
    cex = rowMeans(run$data$fishObsMatrix[which(attributes(run$data)$year == year), (lD-1):(lD+1)])
    points(attributes(run$data)$locObs[attributes(run$data)$year == year,],cex = log(cex+2),col = 'blue')


    #Convert map to UTM coordinates-------------------------------------------------------------------------
    newmap <- maps::map("world", c("Norway","Sweden","Finland","Russia"),fill = TRUE,plot = FALSE, col = "transparent")
    mapTmp = data.frame(newmap$x,newmap$y)
    mapTmp[which(is.na(mapTmp[,1])),] = 3.141592 #Need something different from NA
    names(mapTmp) = c("X","Y")
    attr(mapTmp, "projection") = "LL"
    attr(mapTmp, "zone") = conf_l$zone
    ddpcr::quiet(mapTmp <- PBSmapping::convUL(mapTmp))
    colnames(mapTmp) = c("UTMX", "UTMY")
    mapTmp[which(is.na(newmap$x)),] = NA
    polygon(mapTmp,col = 'lightgrey')
    #----------------------------------------------------------------------------------------------------

    title(main=what[2],line = -1,cex.main=1.2,outer=FALSE)

    if(legend == TRUE){

      image.plot(proj$x,proj$y, inla.mesh.project(proj, latentFieldMAP),col =  colorRampPalette(c("white","yellow", "red"))(12),
            xlab = '', ylab = "",
        #    zlim = c(range[1],range[2]),
            main = "",
            xlim = c(floor(min(proj$x)),ceiling(max(proj$x))),
            ylim = c(floor(min(proj$y)),ceiling(max(proj$y))), # offset removed
            cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1,xaxt = 'n',yaxt = 'n',
            smallplot = c(0.1, 0.15, .1, .9), axes=F)

    }

    title(main=paste0("Spatial effect for length group: ", what[3]),line = 0,cex.main=1.5,outer=TRUE)

    }else if(what[1]=="SST"){
      mesh = attributes(run$data)$meshST
      randomEffects=run$rep$par.random[names(run$rep$par.random)=="xST"]
      year = as.numeric(what[2])
      yearPosition = year-min(run$conf_l$years)+1
      l = run$conf_l$lengthGroupsReduced[which(run$conf_l$lengthGroups == as.numeric(what[3]))]
      indexStart=length(run$conf_l$years)*mesh$n *(l-1)+(yearPosition-1)*mesh$n+1
      indexEnd=indexStart+mesh$n-1
      indexStart2=length(run$conf_l$years)*mesh$n *(l)+(yearPosition-1)*mesh$n+1
      indexEnd2=indexStart2+mesh$n-1


      kappa = exp(run$pl$log_kappa)[2]
      scalingConstant = exp(run$rep$par.fixed[which(names(run$rep$par.fixed)=="log_sigma")])[2]
      scalingConstant2 = sqrt(1/((4*3.14159265)*kappa*kappa))
      scale = scalingConstant/scalingConstant2

      lD = which(run$conf_l$lengthGroups == as.numeric(what[3]))
      ST = (run$data$weigthLength[lD]*randomEffects[indexStart:indexEnd] + (1-run$data$weigthLength[lD])*randomEffects[indexStart2:indexEnd2] )*scale


      meshS = attributes(run$data)$meshS
      randomEffectsS=run$rep$par.random[names(run$rep$par.random)=="xS"]
      indexStartS=meshS$n *(l-1)+1
      indexEndS=indexStartS+meshS$n-1
      indexStartS2=meshS$n *(l)+1
      indexEndS2=indexStartS2+meshS$n-1

      kappa = exp(run$pl$log_kappa)[1]
      scalingConstant = exp(run$rep$par.fixed[which(names(run$rep$par.fixed)=="log_sigma")])[1]
      scalingConstant2 = sqrt(1/((4*3.14159265)*kappa*kappa))
      scale = scalingConstant/scalingConstant2


      S = (run$data$weigthLength[lD]*randomEffects[indexStartS:indexEndS] + (1-run$data$weigthLength[lD])*randomEffectsS[indexStartS2:indexEndS2] )*scale

      beta0 = pl$beta0
      ll = which(run$conf_l$lengthGroups==as.numeric(what[3]))
      beta0This = beta0[yearPosition,ll]
      proj = inla.mesh.projector(mesh)
      latentFieldMAP =beta0This + S +  ST
      image(proj$x,proj$y, inla.mesh.project(proj, latentFieldMAP),col =  colorRampPalette(c("white","yellow", "red"))(12),
            xlab = '', ylab = "",
   #         zlim = c(range[1],range[2]),
            main = paste0("Year ", what[2]),
            cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1,xaxt='n',yaxt='n')

      contour(proj$x, proj$y,inla.mesh.project(proj, latentFieldMAP) ,add = T,labcex  = 1,cex = 1)
      cex = rowMeans(run$data$fishObsMatrix[which(attributes(run$data)$year == year), (lD-1):(lD+1)])
      points(attributes(run$data)$locObs[attributes(run$data)$year == year,],cex = log(cex+2)/2,col = 'blue')

      #Convert map to UTM coordinates-------------------------------------------------------------------------
      newmap <- maps::map("world", c("Norway","Sweden","Finland","Russia"),fill = TRUE,plot = FALSE, col = "transparent")
      mapTmp = data.frame(newmap$x,newmap$y)
      mapTmp[which(is.na(mapTmp[,1])),] = 3.141592 #Need something different from NA
      names(mapTmp) = c("X","Y")
      attr(mapTmp, "projection") = "LL"
      attr(mapTmp, "zone") = conf_l$zone
      ddpcr::quiet(mapTmp <- PBSmapping::convUL(mapTmp))
      colnames(mapTmp) = c("UTMX", "UTMY")
      mapTmp[which(is.na(newmap$x)),] = NA
      polygon(mapTmp,col = 'lightgrey')
      #----------------------------------------------------------------------------------------------------
      if(legend == TRUE){
        image.plot(proj$x,proj$y, inla.mesh.project(proj, latentFieldMAP),col =  colorRampPalette(c("white","yellow", "red"))(12),
                   xlab = '', ylab = "",
                  # zlim = c(range[1],range[2]),
                   main = "",
                   xlim = c(floor(min(proj$x)),ceiling(max(proj$x))),
                   ylim = c(floor(min(proj$y)),ceiling(max(proj$y))), # offset removed
                   cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1,xaxt = 'n',yaxt = 'n',
                   smallplot = c(0.1, 0.15, .1, .9), axes=F)
      }

    }else if(is.null(what)){
      print("Unknown plotting procedure")
    }
}

#' plotTimeofDay
#'
#' Plot time in day effect
#' @param run
#' @export
plotSunAlt<-function(run){

  sunAltEffectLower=  run$rep$value[which(names(run$rep$value)=="fourierReportLow")]
  sunAltEffectLowerL=  sunAltEffectLower - 1.96*run$rep$sd[which(names(run$rep$value)=="fourierReportLow")]
  sunAltEffectLowerU=  sunAltEffectLower + 1.96*run$rep$sd[which(names(run$rep$value)=="fourierReportLow")]

  sunAltEffectUpper=  run$rep$value[which(names(run$rep$value)=="fourierReportHigh")]
  sunAltEffectUpperL=  sunAltEffectUpper - 1.96*run$rep$sd[which(names(run$rep$value)=="fourierReportHigh")]
  sunAltEffectUpperU=  sunAltEffectUpper + 1.96*run$rep$sd[which(names(run$rep$value)=="fourierReportHigh")]

  plot(sunAltEffectLower,
       ylim=c(min(sunAltEffectLowerL,sunAltEffectUpperL),max(sunAltEffectLowerU,sunAltEffectUpperU)),
       main="Sun altitude effect",type="l",col="red",xaxt="n",xlab = "Sun altitude", ylab  = "Effect in linear predictor",
       cex.main = 1.7,cex.lab = 1.5, cex = 1.5)
  abline(a=0,b=0)
  lines(sunAltEffectLowerL,main="Lower",col="red",lty = 2)
  lines(sunAltEffectLowerU,main="Lower",col="red",lty = 2)

  if(run$conf_l$sunAlt[2] ==2){
    lines(sunAltEffectUpper,col="blue")
    lines(sunAltEffectUpperL,main="Lower",col="blue",lty = 2)
    lines(sunAltEffectUpperU,main="Lower",col="blue",lty = 2)
  }


  tmp = length(sunAltEffectLower)
  x = c(1,tmp/4,tmp/2,tmp*3/4,tmp)
  text = c("Lowest (morning)", "", "Heighest", "", "Lowest (evening)")
  axis(1, at=x,labels=text)
  minLength = min(run$conf_l$lengthGroups)
  maxLength = max(run$conf_l$lengthGroups)
  if(run$conf_l$sunAlt[2] ==2){
    legend(legend=c(paste0(minLength, " cm"), paste0(maxLength, " cm")),col=c("red","blue"),"topright",lty=1,cex =1.5)
  }

  abline(v=1,lty=3)
  abline(v=tmp/4,lty=3)
  abline(v=tmp/2,lty=3)
  abline(v=tmp*3/4,lty=3)
  abline(v=tmp,lty=3)

}
