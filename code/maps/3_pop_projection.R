
setwd("~/Dropbox/zika/code")

#### merge worldpop (South America) and WGP (US) data 
######### population from world pop   
demogf = paste("../data/am2010ppp.tif",sep="")  ##  extentswpop<- c(-118.4497823079,0,0, 32.7115384420) ###
demogdata <- as.matrix(raster(demogf),ncol=nx,nrow=ny)
demogdatacr<-glbcoarser(demogdata,5,2)  ## reduce resolution to 5km

#######
## obtain 2010 GRUMP data
popf = paste("../data/gr10_5k_CLEAN.tif",sep="")
glbpop <- as.matrix(raster(popf),ncol=nx,nrow=ny)
glbpop<- rbind(matrix(NA,120,8640),glbpop)
glbpop<-rbind(glbpop, matrix(NA,720,8640))
sum(glbpop, na.rm=T)
glbpop[which(is.na(glbpop))]<- 0
dim(glbpop)
wpopgrump<- glbpop

### location of the americas (wpop grid)   
wpopgrump[1376:3587,1479:3691] <- demogdatacr
sum(wpopgrump,na.rm=T)
cntf = paste("../generated/country_codes_2_5m.tif",sep="")
cntcode <- as.matrix(raster(cntf),ncol=nx,nrow=ny)
mssin<-c(259,63,251,20)
for (cnt in mssin) 
  wpopgrump[which(cntcode==cnt)]<- glbpop[which(cntcode==cnt)]

write.csv(wpopgrump,"../generated/wpopgpw2010.csv", row.names=F, quote=F)

# project to 2015 
popf = paste("../generated/wpopgpw2010.csv",sep="")
glbpop <- read.csv(popf,header=T)
dim(glbpop)
glbpop[which(is.na(glbpop))]<- 0

sum(glbpop, na.rm=T)
countryf = "../data/country_code_2_5min.tif"
countrymx <- as.matrix(raster(countryf),ncol=nx,nrow=ny)
countrymx<- rbind(matrix(NA,120,8640),countrymx)
countrymx<-rbind(countrymx, matrix(NA,720,8640))
countrymx[countrymx==0]<- 1
dim(countrymx)
table(as.vector(unlist(countrymx)))
reassign<- read.csv("../data/reassigned.csv")    #i=1
for (i in 1:dim(reassign)[1]) 
  countrymx[which(countrymx==reassign[i,1])]<- reassign[i,2]

countrypop<- read.csv("../data/pop_country.csv")
head(countrypop)
countrypop[which(countrypop[,6]==0),6]<-1
cpopgrate<- cbind(countrypop[,1],countrypop[,6]/countrypop[,4])
pgrowth<- matrix(cpopgrate[as.vector(unlist(countrymx)),2],4320,8640)
max(pgrowth, na.rm=T)
pop2015<- glbpop*pgrowth
sum(pop2015, na.rm=T)
sum(glbpop, na.rm=T)
pop2015<- matrix(as.vector(unlist(pop2015)),4320,8640)
write.csv(pop2015,"../generated/wpopgpw2015.csv", row.names=F, quote=F)

