library(spatioTemporalIndices)
library(spatioTemporalALK)
rm(list=ls())

#Read data
dat_l = readRDS("catch_at_length_data_ex_rus.rds")
dat_alk = readRDS("catch_at_age_data_ex_rus.rds")

#Configurations length part
conf_l = defConf(years = 1994:2020, # years to use, use all years with data by default
                 maxLength = 75, # Numeric = use directly; NULL = use input data to determine
                 minLength = 20, # Numeric = use directly; NULL = use input data to determine
                 spatioTemporal = 0,
                 spatial = 1,
                 dLength = 5,
                 reduceLength = 3,
                 stratasystem = list(dsn="shapefiles/Vintertokt_strata_system", layer = "Vintertoktet_nye_strata"),
                 minDepth=100,maxDepth=500,
                 applyALK = 1,
                 cutoff = 100, cbound = 130)

#Define configurations age part
conf_alk = defConf_alk(maxAge = 10,
                       minAge = 3,
                       readability = 1,
                       spatioTemporal = 0,
                       rwBeta0 = 1)


confPred = defConfPred(conf=conf_l,Depth="Data",nInt=3000)

# run model
start_time <- Sys.time()
run = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,ignore.parm.uncertainty = TRUE)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed


#Plot covariate effects
plotResults(what = "sunAlt",run = run)
plotResults(what = "depth",run = run)

#Save indices allong with variance and covariance structures
saveIndex(run, file = "index.txt", folder = "")

#Print parameters
spatioTemporalIndices::partable(run)

#Plot ALK in all years
x11(width = 20, height = 15)
par(oma=c(3,3,2,0.5),mar=c(2,2,2,0.5),mfrow = c(5,6))
for(year in 1994:2020){
  plotALK(run,year = year)
}


#Extract ad-reported values
rl = as.list(run$rep,what = "Est", report = TRUE)
rlSd = as.list(run$rep, "Std", report = TRUE)

#Extract indices
minAge = 3
nYears = length(conf_l$years)
ageLogIndex = rl$logAgeIndex[,-1]
ageLogIndexSd = rlSd$logAgeIndex[,-1]

ageIndexUse = round(as.data.frame(cbind(rep(1,dim(exp(ageLogIndex))[1]),exp(ageLogIndex))),3)
ageLogIndexSdUse = round(as.data.frame(ageLogIndexSd),3)


#Extract yearly covariance matrices
minAge = 3

load("cov_index.Rda")
#Illustrate yearly correlation
library(ellipse)
ccolors <- c("#A50F15","#DE2D26","#FB6A4A","#FCAE91","#FEE5D9","white",
             "#EFF3FF","#BDD7E7","#6BAED6","#3182BD","#08519C")

x11()
par(mfrow = c(5,6))
for(i in 1:length(covYears)){
  xx = cov2cor(covYears[[i]])
  plotcorr(xx, col = ccolors[5*xx+6],main = paste0("Year ", 1993 + i))
}


#Contour plot length#############################
ll = run$obj$report()$lengthIndexDetailed
length = 45
l = which(conf_l$lengthGroups== length)

zlim = c(0,log(max(ll[,l,])))
nYears = length(conf_l$years)

mfrow = c(6,5)
par(oma=c(5,5,4,0.5),mar=c(0.5,0.5,0.0,0.0),mfrow = mfrow)
yaxt = rep(c("s",rep("n",4)),6)
xaxt = rep("n",nYears)
xaxt[1:5] = "s"
xaxt= rev(xaxt)


for(year in 1:nYears){
  x = sort(unique(run$data$xInt))
  y = sort(unique(run$data$yInt))
  z = matrix(NA,nrow = length(x), ncol = length(y))

  xC = 1
  for(xx in x){
    yC = 1
    for(yy in y){
      if(length(which(run$data$xInt ==xx & run$data$yInt ==yy))>0){
        z[xC,yC] = log(ll[year,l,which(run$data$xInt ==xx & run$data$yInt ==yy)])
      }
      yC = yC+1
    }
    xC = xC +1
  }
  breaks = c(-1,seq(0,zlim[2],by  = 1))
  image(x,y, z,col =  colorRampPalette(c("white","yellow", "red"))(length(breaks)-1),
        legend.width = 2, cex = 1.6,zlim = zlim,
        axis.args = list(cex.axis = 1.3),yaxt = yaxt[year],xaxt = xaxt[year],
        main = "")

  title(1993 + year,line = -2,cex.main=2)
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
  points(attributes(run$data)$locObs[which(attributes(run$data)$year == (year +1993)),], cex = 0.1)

  if(year==nYears){
    plot(1,1,axes=F,ylim = c(-99,-98))
    image.plot(x,y, z,col =  colorRampPalette(c("white","yellow", "red"))(length(breaks)-1),
               legend.width = 2, cex = 1.6,zlim = zlim,
               axis.args = list(cex.axis = 1.3),yaxt = yaxt[year],xaxt = xaxt[year],
               main = "" ,legend.cex = 2,legend.only = TRUE,
               smallplot = c(0.1, 0.15, .1, .9))
  }

}

mtext(text="Eastern direction (km)",cex=2,side=1,line=2,outer=TRUE)
mtext(text="Northern direction (km)",cex=2,side=2,line=2.4,outer=TRUE)
mtext(paste0("Predicted log CPUE of length ", length, " cm "  ), outer=TRUE,  cex=2, line=0.5)





#Contour plot ALK
alk = run$obj$report()$ALK_int
length = 45
l = which(conf_l$lengthGroups== length)
age = 5

zlim = c(0,1)
mfrow = c(6,5)
par(oma=c(5,5,4,0.5),mar=c(0.5,0.5,0.0,0.0),mfrow = mfrow)
nYears = length(conf_l$years)

yaxt = rep(c("s",rep("n",4)),6)
xaxt = rep("n",nYears)
xaxt[1:5] = "s"
xaxt= rev(xaxt)

for(year in 1:nYears){
  x = sort(unique(run$data$xInt))
  y = sort(unique(run$data$yInt))
  z = matrix(NA,nrow = length(x), ncol = length(y))

  xC = 1
  for(xx in x){
    yC = 1
    for(yy in y){
      if(length(which(run$data$xInt ==xx & run$data$yInt ==yy))>0){
        #       z[xC,yC] = sum(alk[l,1:(age-1),which(run$data$xInt ==xx & run$data$yInt ==yy),year])
        z[xC,yC] = alk[l,age-1,which(run$data$xInt ==xx & run$data$yInt ==yy),year]
      }
      yC = yC+1
    }
    xC = xC +1
  }
  breaks = c(seq(0,zlim[2],by  = 0.01))
  image(x,y, z,col =  colorRampPalette(c("white","yellow", "red"))(length(breaks)-1),
        legend.width = 2, cex = 1.6,zlim = zlim,
        axis.args = list(cex.axis = 1.3),yaxt = yaxt[year],xaxt = xaxt[year],
        main = "")
  title(1993 + year,line = -2,cex.main=2)


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
  points(attributes(run$data)$locObs[which(attributes(run$data)$year == (year +1993)),], cex = 0.1)

  if(year==27){
    plot(1,1,axes=F,ylim = c(-99,-98))
    image.plot(x,y, z,col =  colorRampPalette(c("white","yellow", "red"))(length(breaks)-1),
               legend.width = 2, cex = 1.6,zlim = zlim,
               axis.args = list(cex.axis = 1.3),yaxt = yaxt[year],xaxt = xaxt[year],
               main = "" ,legend.cex = 2,legend.only = TRUE,
               smallplot = c(0.1, 0.15, .1, .9))
  }

}



mtext(text="Eastern direction (km)",cex=2,side=1,line=2,outer=TRUE)
mtext(text="Northern direction (km)",cex=2,side=2,line=2.4,outer=TRUE)
mtext(paste0("P(age = ", age, "|length = ",length,") in each year "  ), outer=TRUE,  cex=2, line=0.5)







#Contour plot age
aa = run$obj$report()$ageIndexDetailed
for(age in 3:10){
  x11()
  zlim = c(0,log(max(aa[,age,])))

  mfrow = c(6,5)
  par(oma=c(5,5,4,0.5),mar=c(0.5,0.5,0.0,0.0),mfrow = mfrow)

  yaxt = rep(c("s",rep("n",4)),6)
  xaxt = rep("n",nYears)
  xaxt[1:5] = "s"
  xaxt= rev(xaxt)

  for(year in 1:nYears){
    x = sort(unique(run$data$xInt))
    y = sort(unique(run$data$yInt))
    z = matrix(NA,nrow = length(x), ncol = length(y))

    xC = 1
    for(xx in x){
      yC = 1
      for(yy in y){
        if(length(which(run$data$xInt ==xx & run$data$yInt ==yy))>0){
          z[xC,yC] = log(aa[year,age-1,which(run$data$xInt ==xx & run$data$yInt ==yy)])
        }
        yC = yC+1
      }
      xC = xC +1
    }
    breaks = c(-1,seq(0,zlim[2],by  = 1))
    image(x,y, z,col =  colorRampPalette(c("white","yellow", "red"))(length(breaks)-1),
          legend.width = 2, cex = 1.6,zlim = zlim,
          axis.args = list(cex.axis = 1.3),yaxt = yaxt[year],xaxt = xaxt[year],
          main = "")
    title(1993 + year,line = -2,cex.main=2)

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
    points(attributes(run$data)$locObs[which(attributes(run$data)$year == (year +1993)),], cex = 0.1)

    points(1000,8600, cex = rlSd$logAgeIndex[year,age-1]^2 *5, col = "blue")

    if(year==nYears){
      plot(1,1,axes=F,ylim = c(-99,-98))
      image.plot(x,y, z,col =  colorRampPalette(c("white","yellow", "red"))(length(breaks)-1),
                 legend.width = 2, cex = 1.6,zlim = zlim,
                 axis.args = list(cex.axis = 1.3),yaxt = yaxt[year],xaxt = xaxt[year],
                 main = "" ,legend.cex = 2,legend.only = TRUE,
                 smallplot = c(0.1, 0.15, .1, .9))
    }


  }


  mtext(text="Eastern direction (km)",cex=2,side=1,line=2,outer=TRUE)
  mtext(text="Northern direction (km)",cex=2,side=2,line=2.4,outer=TRUE)
  mtext(paste0("Predicted log CPUE of age ", age  ), outer=TRUE,  cex=2, line=0.5)
}

