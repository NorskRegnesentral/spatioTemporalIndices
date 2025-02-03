

#' plotResults Plot results
#' @param run fitted object returned by \code{\link{fitModel}}
#' @param what what to plot, options in first element "sunAlt", "depth", "S", "ST" and "SST". If spatial effect is plotted, year and length must be provided
#' @param legend include legend
#' @importFrom graphics abline legend lines
#' @export
plotResults  <- function(run,what=NULL, legend = FALSE){
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
  }

  #TODO: spatial plot

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
