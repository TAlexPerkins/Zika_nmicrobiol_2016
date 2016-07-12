library(rgdal)
.libPaths("/afs/crc.nd.edu/user/a/asiraj/raster")
library(raster)
setwd("/afs/crc.nd.edu/user/a/asiraj/zika/code")
source("0_numfunctions.R")

# read population data
pop2015<- read.csv("../generated/wpopgpw2015.csv",header=T)
newgrid<- as.vector(unlist(read.csv("../generated/americas_serial.csv", header=T)[,1]))
#################################
print("reading births raster...")
resf = paste("../generated/birthsUS_2_5m.bil",sep="")
globalbirths <- as.matrix(raster(resf),ncol=nx,nrow=ny)

##### Exclude non-region grids
NAcells<-  as.vector(unlist(read.csv("../data/NA_cells.csv", header=T)))

genid=1
print("reading pregnancies raster...")
resf = paste("../generated/pregnanciesUS_2_5m.bil",sep="")
globalpregs <- as.matrix(raster(resf),ncol=nx,nrow=ny)
remer=c(8,59,92)
n.cores <- 1000
n.lines <- 1000
seqP <- seq(1,n.lines,n.lines/n.cores)
xx <- as.numeric(Sys.getenv("SGE_TASK_ID"))
row.num <- seqP[xx] 
#row.num<- remer[row.num]

for (repl in row.num:row.num) {

  ###### assemble R0 from cluster
  allR0i<- NULL
  print("assembling R0s...")
  for (i in 1:13){
    allR0i<-c(allR0i,as.vector(unlist(read.csv(paste("../generated2/repl_",repl,"_R0x",i,".csv",sep=""),header=T))))
        print(c("repl:",repl,"R0 chunk:",i))
  }

  allR0<- matrix(NA,4320,8640)
  allR0[newgrid]<- allR0i
#  write.csv(allR0,paste("../generated/repl_",repl,"_R0global.csv",sep=""), row.names=F, quote=F)
  #  globalR0<-  matrix(as.vector(unlist(read.csv(paste("../data/R0global.csv",sep=""), header=T))),4320,8640)
  
    ###### assemble AR from cluster
print("assembling ARs...")
  allARi<- NULL
  for (i in 1:13){
    allARi<-c(allARi,as.vector(unlist(read.csv(paste("../gen",genid,"/repl_",repl,"_ARx",i,".csv",sep=""),header=T))))
        print(c("repl:",repl,"AR chunk:",i))
  }

  allAR<- matrix(NA,4320,8640)
  allAR[newgrid]<- allARi

  globalR0 <- allR0
  globalAR <- allAR 
  ########

  ##### map of R0
print("generating R0 maps...") 
  outdf1 = paste("../output2/R0_SA_repl_",repl,sep="")
  output2_5_americas(americas(globalR0,NAcells),outdf1)
  
  ##### map of AR
print("generating AR maps...") 
  outdf1 = paste("../output2/AR_US_SA_repl_",repl,sep="")
  output2_5_americas(americas(globalAR,NAcells),outdf1)

  
  ###### # Infecteds
print("generatng Infecteds map...") 
  outdf1 = paste("../outputbk/cumI_US_SA_repl_",repl,sep="")
  output2_5_americas(americas(globalAR,NAcells)*  americas(pop2015),outdf1)


    
  ###### # infected pregnants
#print("generating infected pregnancies map...") 
#  demogdatacr<-globalpregs[-c(1:971,3516:4320),-c(1:1320,3629:8640)]
#  outdf1 = paste("../outputbk/cumPI_US_SA_repl_",repl,sep="")
#  output2_5_americas((americas(globalAR,NAcells)* demogdatacr),outdf1)
  

  ###### # births by infected
print("generating infected births map...") 
  demogdatacr<-globalbirths[-c(1:971,3516:4320),-c(1:1320,3629:8640)]
  outdf1 = paste("../outputbk/cumBI_US_SA_repl_",repl,sep="")
  output2_5_americas((americas(globalAR,NAcells)* demogdatacr),outdf1)
  
} ### loop through all replicates 
  
