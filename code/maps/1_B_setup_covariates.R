if(!require(raster)){install.packages('raster');library(raster)}
if(!require(sp)){install.packages('sp');library(sp)}

setwd("~/Dropbox/zika/code")
source("0_numfunctions.R")

### vector occurance probability grid
mskctoff<-c(.681,.397, .681, .397)  ## vector occurance lower bounds
vecf = paste("../data/aegypti.tif",sep="")
vectormx <- as.matrix(raster(vecf),ncol=nx,nrow=ny)
vectormx<- rbind(matrix(NA,120,8640),vectormx)
vectormx<-rbind(vectormx, matrix(NA,720,8640))
vectormx[(vectormx<mskctoff[1] & vectormx>=0) ]<- NA 
min(vectormx[vectormx>0], na.rm=T)

#### Econ values based on GCP (Gross cell product). 
#### Within a country all cells add up to the country's GDPPPP (in Billion USD)
########## Economic grid

allecon<- matrix(as.vector(unlist(read.csv("../data/Gecon_lat_lon.csv", header=T)[,-1])),27445,8)
head(allecon)
allecon[which(is.na(allecon[,3])),3]<-0 
allecon[which(is.na(allecon[,4])),4]<-0 
grdpop<-aggregate(allecon[,3],list(lat=allecon[,1], lon=allecon[,2]), FUN=sum)
grdecon<-aggregate(allecon[,4],list(lat=allecon[,1], lon=allecon[,2]), FUN=sum)

grid.y<- seq(-90,89)
grid.x<- seq(-180,179) 
ecogrid<- matrix(NA,nrow=180,ncol=360)
for (i in 1: (dim(grdecon)[1])) {
  mxloc<- c(which(grid.y == grdecon[i,1]), which(grid.x== grdecon[i,2]))
  ecogrid[mxloc[1], mxloc[2]]<- grdecon[i,3]
}

ecogrid<-ecogrid[180:1,]
ecogrid[which(is.na(ecogrid))]<-0
which(ecogrid==max(ecogrid,na.rm=T), arr.ind=T) ## richest grid

### population
###### note that the econgrid product used 2005 GPW estimates 
popf = paste("../data/glp05ag.bil",sep="")   ### gridded population GPW estimates (2005)
pop2005 <- as.matrix(raster(popf),ncol=nx,nrow=ny)
dim(pop2005)
pop2005<- rbind(matrix(NA,120,8640),pop2005)
pop2005<-rbind(pop2005, matrix(NA,768,8640))
dim(pop2005)

### aggregate population to the resolution of econgrid (1deg) 
pop2005crs<-glbcoarser1deg(pop2005,24,2)

## per capita
econ_pc<- (ecogrid*10^9 / pop2005crs)
econ_pc[which(econ_pc==Inf)]<-0

###### grids that have missing Econgrid values (by serial number)
##### imputed by taking the average of their eight neihgbours (in two iterations)
nills<- c(19532,19356,19357,19536,19891,19892,20072, 20432,20468,20792,20972,21152,20611,20791,20971,20250,20610,20790,20970,20249,20429,20609,20789,20969,20608,
          20788,20607,20606,21351,21352,21531,21532,21711,22765,22766,22767,22945,22946,22947)

## 1st iteration
  econ_pc[nills]
  ecomx<-as.vector(unlist(econ_pc))
  newval <- NULL
  for (bx in nills)
    newval<-c(newval,mean(ecomx[c(bx+1, bx-1,bx+180,bx-180,bx+180+1,bx+180-1,bx-180+1,bx-180-1)],na.rm=T))
  ecomx[nills]<- newval

## 2nd interation
  nills<- nills[which(newval==0)]
  newval <- NULL
  #bx=nills 
  for (bx in nills)
    newval<-c(newval,mean(ecomx[c(bx+1, bx-1,bx+180,bx-180,bx+180+1,bx+180-1,bx-180+1,bx-180-1)],na.rm=T))
  ecomx[nills]<- newval

  length(ecomx)
  ecomx<- matrix(ecomx,180,360)
  
  outdf1 = paste("../generated/economic_cdp_1deg_v3",sep="")  ## grid output
  output1deg(ecomx,outdf1)
  
  ### resample to 2.5 min resolution
  ecogridxmx2<-nearestnb(ecomx)
  outdf1 = paste("../generated/economic_cdp_v3",sep="")
  output2_5(ecogridxmx2,outdf1)    # grid output

  ecopcpp<-ecogridxmx2

  ################
  #### generate the covariates for R0
  
  ###### Americas extents
  newgrid<- as.vector(unlist(read.csv("../data/americas_serial.csv", header=T)[,1]))

  
  grid.x= seq(-180,180,0.041666667)
  grid.y = seq(-90,90,0.041666667)

  allcontent<- cbind(expand.grid((grid.y[4320:1]),grid.x)[,2:1], as.vector(unlist(vectormx)), NA)[newgrid,]
  for (mon in 1:12) {
    climf = paste("../data/tmean",mon,".bil",sep="")
    climdata <- as.matrix(raster(climf),ncol=nx,nrow=ny)
    climdata<-rbind(climdata, matrix(NA,720,8640))
    dim(climdata)
    climdata<- climdata/10
    allcontent<- cbind(allcontent,as.vector(unlist(climdata))[newgrid])
  }

  dim(allcontent)

  ##### small islands missing GDP values  data from CIA https://www.cia.gov/library/publications/the-world-factbook/geos/xx.html
  #### and http://www.indexmundi.com/g/g.aspx?v=67&c=re&l=en
  gdp_pcpppextra<- matrix(as.vector(unlist(read.csv("../data/extra_gdp_values.csv", header=T))),11,3)
  ##gdp_pcpppextra[,1]<- gdp_pcpppextra[,1]* gdp_pcpppextra[,3]/gdp_pcpppextra[,4] 
  contf2 = paste("../data/country_code_2_5min.tif",sep="")
  cntry <- as.matrix(raster(contf2),ncol=nx,nrow=ny)  
  cntry<- rbind(matrix(NA,120,8640),cntry)
  cntry<-rbind(cntry, matrix(NA,720,8640))
  dim(cntry)
  # gd=1

  for(gd in 1:length(gdp_pcpppextra[,1])) {
    cntgrd<- which(cntry==gdp_pcpppextra[gd,2])[-1]
    ecopcpp[cntgrd]<- gdp_pcpppextra[gd,1]
  }
  ecopcpp[which(ecopcpp==Inf)]<-0

  allcontent<- data.frame(cbind(allcontent, as.vector(unlist(ecopcpp))[newgrid]))
  allcontent[which(allcontent[,17]<0),17]<-0

  names(allcontent)<- c('lon', 'lat', 'aegypti','albopictus', 'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','gdp_pcppp2005')
  write.csv(allcontent,"../generated/allcovs_v4_compact.csv", row.names=F, quote=F)
