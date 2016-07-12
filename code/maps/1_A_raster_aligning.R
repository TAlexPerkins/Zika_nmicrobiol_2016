setwd("~/Dropbox/zika/code")


### Merge and align raster data 
############### pregnancies and births in Americas
  extents<-c(-125.0,-56.458334181, -28.87500,49.5)  ## includes US
  demogf = paste("../rasters/america2015adjustedBirths_USA.tif",sep="")
  demogata1 <- as.matrix(raster(demogf),ncol=nx,nrow=ny)
  dim(demogata1)
  demogata<- demogata1[-c(1:2624),-c(1:6497,18039:43073)]  ### crop to SA & US
  dim(demogata)
  #### covert to coarser res by summation (2.5 min)
  demogatacr<-glbcoarser(demogata,5,2)
  dim(demogatacr)
  demogatacr<-demogatacr[-2545,]
  demogatacr<- cbind(matrix(NA,2544,1320),demogatacr,matrix(NA,2544,5012))  ## setup a global grid around SA & US
  demogatacr<- rbind(matrix(NA,971,8640),demogatacr,matrix(NA,805,8640))    ## setup a global grid around SA & US
  
  outdf1 = paste("../generated/births_SA_2_5m",sep="")   ### output in .bil format 
  output2_5(demogatacr,outdf1)

###### pregnancies
  demogf = paste("../rasters/america2015pregnancies_USA.tif",sep="")
  demogata1 <- as.matrix(raster(demogf),ncol=nx,nrow=ny)
  demogata<- demogata1[-c(1:2624),-c(1:6497,18039:43073)]
  dim(demogata)

  demogatacr<-glbcoarser(demogata,5,2)
  dim(demogatacr)
  demogatacr<-demogatacr[-2545,]
  demogatacr<- cbind(matrix(NA,2544,1320),demogatacr,matrix(NA,2544,5012))
  demogatacr<- rbind(matrix(NA,971,8640),demogatacr,matrix(NA,805,8640))

  outdf1 = paste("../generated/pregnanciesUS_2_5m",sep="")
  output2_5(demogatacr,outdf1)
  print(outdf1)
  
  ### extent for South America only grids
  extents<-c(-118.3752339997139131,-56.5248320530909467, -28.8371489997140742,32.7215979469088865)
  
  
