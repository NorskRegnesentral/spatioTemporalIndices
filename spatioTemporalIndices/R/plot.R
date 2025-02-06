

#' plotResults Plot results
#' @param run fitted object returned by \code{\link{fitModel}}
#' @param what What to plot. Options in first element: "sunAlt", "depth", "ALK", or "space". If "ALK" is selected, year also must be provided, e.g., what = c("ALK",2020). If a spatial plot is selected; year and length or age must also provided. Eg.; what = c("space", 2020,5,"age") or what = c("space", 2020,50,"length").
#' @param xlim optional xlim sent to fields::image.plot
#' @param ylim optional xlim sent to fields::image.plot
#' @param zlim optional zlim sent to fields::image.plot
#' @importFrom graphics abline legend lines
#' @export
plotResults  <- function(run,what=NULL, xlim = NULL, ylim = NULL, zlim = NULL){
  pl = as.list(run$rep,what = "Est")
  if(what[1] == "sunAlt"){
    plotSunAlt(run)
  }else if(what[1] == "depth"){

    depthSpline = run$rep$value[names(run$rep$value)=="depthReport1"]
    sdDepthSpline<-run$rep$sd[names(run$rep$value)=="depthReport1"]
    depth=seq(from=min(attributes(run$data)$depth),to=max(attributes(run$data)$depth),length.out = length(depthSpline))

    plot(depth, depthSpline, lty=1,type = 'l',ylim = c(min(depthSpline - 1.96*sdDepthSpline),max(depthSpline + 1.96*sdDepthSpline)),
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
    for(i in 1:20){
      abline(v=i*50,lty=3)
    }
  }else if(what[1]=="space"){

    xInt = run$data$xInt
    yInt = run$data$yInt

    if(is.null(xlim)) xlim = range(xInt) + c(-diff(range(xInt))/10,diff(range(xInt))/10)
    if(is.null(ylim)) ylim = range(yInt) + c(-diff(range(yInt))/10,diff(range(yInt))/10)

    yearOriginal = what[2]
    year = which(run$conf_l$years == yearOriginal)

    if(what[4] == "length"){
      report = run$obj$report()$lengthIndexDetailed
      lengthOriginal = as.numeric(what[3])
      length = max(which(run$conf_l$lengthGroups < lengthOriginal))
      if(is.null(zlim)) zlim = c(0,log(max(report[year,length,])))
    }else{
      report = run$obj$report()$ageIndexDetailed
      ageOriginal = as.numeric(what[3])
      age = ageOriginal - run$conf_alk$minAge + 2
      if(is.null(zlim)) zlim = c(0,log(max(report[year,age,])))
    }

    col = grDevices::colorRampPalette(c("#FFFFCC", "#FFEDA0", "#FEB24C", "#F03B20", "#BD0026"))(20)

    x = sort(unique(run$data$xInt))
    y = sort(unique(run$data$yInt))
    z = matrix(NA,nrow = length(x), ncol = length(y))
    xC = 1
    for(xx in x){
      yC = 1
      for(yy in y){
        if(length(which(run$data$xInt ==xx & run$data$yInt ==yy))>0){
          if(what[4] == "length"){
            z[xC,yC] = log(report[year,length,which(run$data$xInt ==xx & run$data$yInt ==yy)])
          }else{
            z[xC,yC] = log(report[year,age,which(run$data$xInt ==xx & run$data$yInt ==yy)])
          }
        }
        yC = yC+1
      }
      xC = xC +1
    }

    if(what[4] == "length"){
      fields::image.plot(x,y, z,col =  col,
            cex = 1.6,zlim = zlim, xlim = xlim,ylim = ylim,
            xlab = "Eastern direction (km)", ylab = "Northern direction (km)",
            main = paste0("Log CPUE at length ", lengthOriginal, " in year ", yearOriginal ),
            cex.lab = 1.4,cex.main = 1.5)
    }else{
      fields::image.plot(x,y, z,col =  col,
            cex = 1.6,zlim = zlim, xlim = xlim,ylim = ylim,
            xlab = "Eastern direction (km)", ylab = "Northern direction (km)",
            main = paste0("Log CPUE at age ", ageOriginal, " in year ",yearOriginal ),
            cex.lab = 1.4,cex.main = 1.5)
    }

  }else if(what[1] == "ALK"){

    lengthInt = seq(min(run$data$length), max(run$data$length), length.out = 200)
    nAges = run$data$ageRange[2]- run$data$ageRange[1] + 1

    col <- grDevices::rainbow(nAges, start = 0, end = 0.95)

    linPredMatrix = matrix(0,length(lengthInt), nAges-1)
    for(a in 1:(nAges-1)){
      linPredMatrix[,a] = run$pl$beta0_alk[which(run$conf_alk$years == what[2]), a] + run$pl$betaLength_alk[a]*lengthInt
    }
    ALK = matrix(0,length(lengthInt), nAges)
    for(a in 1:nAges){
      probLess = rep(0,dim(ALK)[1])
      if(a>1){
        for(b in 1:(a-1)){
          tmp = ALK[,b];
          probLess = probLess + tmp;
        }
      }
      if(a <(nAges)){
        tmp2 = linPredMatrix[,a];
        ALK[,a] = stats::plogis(tmp2)*(1-probLess);
      }else{
        ALK[,nAges] = (1-probLess);
      }
    }

    plot(lengthInt,ALK[,1], type = "l", ylab = "Probability", xlab = "Length",
         ylim = c(0,1), main = paste0("Spatially averaged ALK in  year ",what[2]),
         col = col[1],lwd = 3, cex.main = 1.5,cex.lab = 1.4, cex.axis =1.2)
    for(l in 2:nAges){
      lines(lengthInt,ALK[,l],col = col[l],lwd = 3)
    }
    #Include labels
    ages = run$conf_alk$minAge:run$conf_alk$maxAge
    for(i in 0:length(ages)){
      graphics::text(x = graphics::par("usr")[1] + 0.02 * diff(graphics::par("usr")[1:2]),  # 2% from the left boundary
           y = i / (length(ages) + 1) + 0.05,
           labels = paste("a =", i + ages[1] - 1),
           col = col[i + 1],
           cex = 1.3,
           adj = 0)  # Left alignment
    }
  }
}

#' plotTimeofDay
#'
#' Plot time in day effect
#' @param run Fitted model returned by, \code{\link{fitModel}}
#' @importFrom graphics abline axis legend lines
#' @export
plotSunAlt<-function(run){

  sunAltEffectLower=  run$rep$value[which(names(run$rep$value)=="fourierReportLow")]
  sunAltEffectLowerL=  sunAltEffectLower - 1.96*run$rep$sd[which(names(run$rep$value)=="fourierReportLow")]
  sunAltEffectLowerU=  sunAltEffectLower + 1.96*run$rep$sd[which(names(run$rep$value)=="fourierReportLow")]

  sunAltEffectUpper=  run$rep$value[which(names(run$rep$value)=="fourierReportHigh")]
  sunAltEffectUpperL=  sunAltEffectUpper - 1.96*run$rep$sd[which(names(run$rep$value)=="fourierReportHigh")]
  sunAltEffectUpperU=  sunAltEffectUpper + 1.96*run$rep$sd[which(names(run$rep$value)=="fourierReportHigh")]

  pi = 3.1415
  tmp = seq(0,2*pi, length.out = length(sunAltEffectLower))
  plot(tmp,sunAltEffectLower,
       ylim=c(min(sunAltEffectLowerL,sunAltEffectUpperL),max(sunAltEffectLowerU,sunAltEffectUpperU)),
       main="Sun altitude effect",type="l",col="red",xaxt="n",xlab = "Sun altitude", ylab  = "Effect in linear predictor",
       cex.main = 1.7,cex.lab = 1.5, cex = 1.5)
  abline(a=0,b=0)
  lines(tmp,sunAltEffectLowerL,main="Lower",col="red",lty = 2)
  lines(tmp,sunAltEffectLowerU,main="Lower",col="red",lty = 2)

  if(run$conf_l$sunAlt[2] ==2){
    lines(tmp,sunAltEffectUpper,col="blue")
    lines(tmp,sunAltEffectUpperL,main="Lower",col="blue",lty = 2)
    lines(tmp,sunAltEffectUpperU,main="Lower",col="blue",lty = 2)
  }

  x = c(0,0.5*pi,pi,pi*3/2,2*pi)
  text = c("Lowest (morning)", "", "Highest", "", "Lowest (evening)")
  axis(1, at=x,labels=text)
  minLength = min(run$conf_l$lengthGroups)
  maxLength = max(run$conf_l$lengthGroups)
  if(run$conf_l$sunAlt[2] ==2){
    legend(legend=c(paste0(minLength, " cm"), paste0(maxLength, " cm")),col=c("red","blue"),"topright",lty=1,cex =1.5)
  }

  abline(v=0,lty=3)
  abline(v=pi/2,lty=3)
  abline(v=pi,lty=3)
  abline(v=pi*3/2,lty=3)
  abline(v=2*pi,lty=3)
}
