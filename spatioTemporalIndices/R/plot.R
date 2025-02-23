

#' plotResults Plot results
#' @param run fitted object returned by \code{\link{fitModel}}
#' @param what What to plot. Options in first element: "sunAlt", "depth", "ALK", or "space". If "ALK" is selected, year also must be provided, e.g., what = c("ALK",2020). If a spatial plot is selected; year and length or age must also provided. Eg.; what = c("space", 2020,5,"age") or what = c("space", 2020,50,"length").
#' @param year year of interest
#' @param age age of interest for spatial CPUE plot
#' @param length length of interest for spatial CPUE plot
#' @param lon_lat longitude and latitude of interest for age-length-key plot
#' @param xlim optional xlim sent to fields::image.plot
#' @param ylim optional xlim sent to fields::image.plot
#' @param zlim optional zlim sent to fields::image.plot
#' @importFrom graphics abline legend lines
#' @export
plotResults  <- function(run,what=NULL,year = NULL,age = NULL,length = NULL,lon_lat = NULL, xlim = NULL, ylim = NULL, zlim = NULL){
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

    if(!is.null(year)){
      yearInternal = which(run$conf_l$years == year)
    }else{
      stop("Year must be provided for spatial CPUE plot")
    }

    if(!is.null(length)){
      report = run$obj$report()$lengthIndexDetailed
      lengthInternal = max(which(run$conf_l$lengthGroups < length))
      if(is.null(zlim)) zlim = c(0,log(max(report[yearInternal,lengthInternal,])))
    }else if(!is.null(age)){
      report = run$obj$report()$ageIndexDetailed
      ageInternal = age - run$conf_alk$minAge + 2
      if(is.null(zlim)) zlim = c(0,log(max(report[yearInternal,ageInternal,])))
    }else{
      stop("length or age must be provided for spatial CPUE plot")
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
          if(!is.null(length)){
            z[xC,yC] = log(report[yearInternal,lengthInternal,which(run$data$xInt ==xx & run$data$yInt ==yy)])
          }else{
            z[xC,yC] = log(report[yearInternal,ageInternal,which(run$data$xInt ==xx & run$data$yInt ==yy)])
          }
        }
        yC = yC+1
      }
      xC = xC +1
    }

    if(!is.null(length)){
      fields::image.plot(x,y, z,col =  col,
            cex = 1.6,zlim = zlim, xlim = xlim,ylim = ylim,
            xlab = "Eastern direction (km)", ylab = "Northern direction (km)",
            main = paste0("Log CPUE at length ", length, " in year ", year ),
            cex.lab = 1.4,cex.main = 1.5)
    }else{
      fields::image.plot(x,y, z,col =  col,
            cex = 1.6,zlim = zlim, xlim = xlim,ylim = ylim,
            xlab = "Eastern direction (km)", ylab = "Northern direction (km)",
            main = paste0("Log CPUE at age ", age, " in year ",year ),
            cex.lab = 1.4,cex.main = 1.5)
    }

  }else if(what[1] == "ALK"){

    if(is.null(year))stop("Year must be provided for ALK plot")

    lengthInt = seq(min(run$data$length), max(run$data$length), length.out = 200)
    nAges = run$data$ageRange[2]- run$data$ageRange[1] + 1

    col <- grDevices::rainbow(nAges, start = 0, end = 0.95)
    linPredMatrix = matrix(0,length(lengthInt), nAges-1)

    if(!is.null(lon_lat)){
      scale = 1/((4*pi)*exp(run$pl$logKappa_alk))
      sigma = exp(run$pl$logSigma_alk)

      lon_lat_df = data.frame(longitude = lon_lat[1],latitude = lon_lat[2])
      point = sf::st_as_sf(lon_lat_df,coords=c("longitude","latitude"),crs="+proj=longlat")
      pointUTM = sf::st_coordinates(sf::st_transform(point,crs=paste0("+proj=utm +zone=", conf_l$zone," +datum=WGS84 +units=km +no_defs")))
      A_alk = fmesher::fm_basis(attributes(run$data)$meshS, pointUTM)

      for(a in 1:(nAges-1)){
        linPredMatrix[,a] = linPredMatrix[,a] +  (A_alk %*% run$pl$xS_alk[, a]*sigma[1]/sqrt(scale[1])  +
                                                    A_alk %*% run$pl$xST_alk[,which(run$conf_alk$years == year), a])[1,1]*sigma[2]/sqrt(scale[2])
      }
    }

    for(a in 1:(nAges-1)){
      if(run$conf_alk$betaLength==0){
        linPredMatrix[,a] = linPredMatrix[,a] + run$pl$beta0_alk[which(run$conf_alk$years == year), a] + run$pl$betaLength_alk[a]*lengthInt
      }else{
        linPredMatrix[,a] = linPredMatrix[,a] + run$pl$beta0_alk[which(run$conf_alk$years == year), a] + run$pl$betaLength_alk[(which(run$conf_alk$years == year) -1)*(nAges-1) +  a]*lengthInt
      }
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

    main = paste0("Spatially averaged ALK in  year ",year)
    if(!is.null(lon_lat)){
      main = paste0("ALK in  year ",year, " at coordinates (",lon_lat[1],", ",lon_lat[2],")")
    }

    plot(lengthInt,ALK[,1], type = "l", ylab = "Probability", xlab = "Length",
         ylim = c(0,1), main = main,
         col = col[1],lwd = 3, cex.main = 1.3,cex.lab = 1.4, cex.axis =1.2)
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
  }else if(what == "variance"){
    if(run$conf_l$applyALK==1){
      par(mar = c(5, 4, 4, 8) + 0.1)  # Increase right margin for legend
      par(xpd = TRUE)  # Allow legend outside the plot

      sd = data.frame(run$rlSd$logAgeIndex)
      if(length(run$conf_alk$minAge:run$conf_alk$maxAge) != dim(sd)[2]){
        sd = sd[,-1]#Note: Youngest age is not in survey
      }
      sd[sd>100] = NA #Remove infinite variances.
      var = sd^2
      colnames(var) = run$conf_alk$minAge:run$conf_alk$maxAge
      rownames(var) = run$conf_l$years
      line_colors <- 1:ncol(var)
      line_types <- 1:ncol(var)

      matplot(var, xaxt = "n", main = "",ylab = "",xlab = "",
              lwd=3, type = "l", las=1, col = line_colors, lty = line_types)

      axis(side =1, at=1:dim(var)[1], labels=rownames(var))
      mtext(text="Year",cex=1.5,side=1,line=2.1,outer=FALSE)
      mtext(text="Variance",cex=1.5,side=2,line=2.3,outer=FALSE)
      mtext("Variance of log index-at-age", outer=FALSE,  cex=1.6, line=0.5)
      legend("topright", inset = c(-0.3, 0), legend = colnames(var),  col = line_colors, lty = line_types, lwd = 2.8, title = "Age", cex = 1.2)

    }else{

    }
  }else if(what == "correlation"){
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
      ccolors <- c("#A50F15","#DE2D26","#FB6A4A","#FCAE91","#FEE5D9","white","#EFF3FF","#BDD7E7","#6BAED6","#3182BD","#08519C")

      ncols <- ceiling(sqrt(length(covYears)+1))  # Columns based on square root of n
      nrows <- ceiling((length(covYears)+1) / ncols)  # Rows to fit the plots

      mfrow = c(ncols,nrows)
      par(oma=c(1,1,1,1),mar=c(0.5,0.5,0.0,0.0),mfrow = mfrow)
      for(i in 1:length(covYears)){
        xx = stats::cov2cor(covYears[[i]])
        ellipse::plotcorr(xx, col = ccolors[5*xx+6],main =paste0("Year ", run$conf_l$years[i]) , mar = c(0,0,0.7,0.7))
      }
      plot(c(0,2),c(0,1),type = 'n', axes = FALSE,xlab = '', ylab = '', main = '')
      text(x=1.3, y = seq(0,1,l=5), labels = round(seq(-1,1,l=5),2))
      graphics::rasterImage(rev(ccolors), 0, 0, 1,1)
    }else{
      stop("Plot of correlation structures not implemented yet for indices-at-length")
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

#' plotArea
#'
#' Plot the survey area, mesh and integration points
#' @param run Fitted model returned by, \code{\link{fitModel}}
#' @param add Add plot to existing plot.
#' @param cex Sent to plot
#' @param pch Sent to points
#' @param extend Proportion extends xlim and ylim for visual clearity
#' @importFrom graphics abline axis legend lines
#' @export
#plotArea<-function(run, add = FALSE,cex = 0.5, pch = 20, extend = 0.2){
#  utmCRS = paste0("+proj=utm +zone=", run$conf_l$zone," +datum=WGS84 +units=km +no_defs")
#
#  mesh = attributes(run$data)$meshS
#  intPoints = data.frame(UTMX = run$data$xInt, UTMY = run$data$yInt)
#  strata_utm <- sf::st_transform(run$conf_l$strata,utmCRS)
#
#  xlim = range(intPoints[,1]) +c(-diff(range(intPoints[,1]))*extend, diff(range(intPoints[,1]))*extend)
#  ylim = range(intPoints[,2]) +c(-diff(range(intPoints[,2]))*extend, diff(range(intPoints[,2]))*extend)
#
# if(add){
#    points(intPoints,pch = pch,cex = cex)
#  }else{
#    plot(intPoints, xlim = xlim, ylim = ylim ,pch = pch,cex = cex)
#  }
#  plot(strata_utm,col = 0,add = TRUE)
#  plot(mesh,add = TRUE)
#}

