#.libPaths("/afs/crc.nd.edu/user/a/asiraj/R/gdal")
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

print("reading pregnancies raster...")
resf = paste("../generated/pregnanciesUS_2_5m.bil",sep="")
globalpregs <- as.matrix(raster(resf),ncol=nx,nrow=ny)
n.cores <- 1000
n.lines <- 1000
seqP <- seq(1,n.lines,n.lines/n.cores)
xx <- as.numeric(Sys.getenv("SGE_TASK_ID"))
row.num <- seqP[xx] 

for (repl in row.num:row.num) {

  ###### assemble AR from cluster
print("assembling ARs...")
  allARi<- NULL
  for (i in 1:13){
    allARi<-c(allARi,as.vector(unlist(read.csv(paste("../generated/repl_",repl,"_ARx_stat",i,".csv",sep=""),header=T))))
        print(c("repl:",repl,"AR chunk:",i))
  }

  allAR<- matrix(NA,4320,8640)
  allAR[newgrid]<- allARi

  globalAR <- allAR 
  ########
  ##### Americas
  
  ##### map of AR
print("generating AR maps...") 
  outdf1 = paste("../output2/AR_SA_repl_",repl,"_stat",sep="")
  output2_5_americas(americas(globalAR,NAcells),outdf1)
  
  ###### # Infecteds
print("generatng Infecteds map...") 
  outdf1 = paste("../outputbk/cumI_SA_repl_",repl,"_stat",sep="")
  output2_5_americas(americas(globalAR,NAcells)*  americas(pop2015),outdf1)
  
  ###### # infected pregnants
print("generating infected pregnancies map...") 
  demogatacr<-globalpregs[-c(1:971,3516:4320),-c(1:1320,3629:8640)]
  outdf1 = paste("../outputbk/cumPI_SA_repl_",repl,"_stat",sep="")
  output2_5_americas((americas(globalAR,NAcells)* demogatacr),outdf1)
  

  ###### # births by infected
print("generating infected births map...") 
  demogatacr<-globalbirths[-c(1:971,3516:4320),-c(1:1320,3629:8640)]
  outdf1 = paste("../outputbk/cumBI_SA_repl_",repl,"_stat",sep="")
  output2_5_americas((americas(globalAR,NAcells)* demogatacr),outdf1)
  
} ### loop through all replicates 
  
