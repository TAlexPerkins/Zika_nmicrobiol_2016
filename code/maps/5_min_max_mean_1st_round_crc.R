.libPaths("/afs/crc.nd.edu/user/a/asiraj/R/gdal")
library(rgdal)
.libPaths("/afs/crc.nd.edu/user/a/asiraj/raster")
library(raster)
setwd("/afs/crc.nd.edu/user/a/asiraj/zika/code")
source("0_numfunctions.R")

tenth<-seq(0,1000,10)
n.cores <- 100
n.lines <- 100
seqP <- seq(1,n.lines,n.lines/n.cores)
xx <- as.numeric(Sys.getenv("SGE_TASK_ID"))
row.num <- seqP[xx] 

allsum_cumI<- 0
allsum_cumBI<- 0
allsum_R0<- 0
allsum_AR<- 0

allmin_cumI<- 99
allmin_cumBI<- 99
allmin_R0<- 99
allmin_AR<- 99

allmax_cumI<- 0
allmax_cumBI<- 0
allmax_R0<- 0
allmax_AR<- 0

print("sum of ten cumI, PI and BI values...")

for (repl in (tenth[row.num]+1):(tenth[row.num+1])) {

  resf = paste("../output2/R0_SA_repl_",repl,".bil",sep="")
  allsum_R0<- allsum_R0 + as.matrix(raster(resf),ncol=nx,nrow=ny)
  allmin_R0<- pmin(allmin_R0,as.matrix(raster(resf),ncol=nx,nrow=ny))
  allmax_R0<- pmax(allmax_R0,as.matrix(raster(resf),ncol=nx,nrow=ny))

  resf = paste("../output2/AR_SA_repl_",repl,".bil",sep="")
  allsum_AR<- allsum_AR + as.matrix(raster(resf),ncol=nx,nrow=ny)
  allmin_AR<- pmin(allmin_AR,as.matrix(raster(resf),ncol=nx,nrow=ny))
  allmax_AR<- pmax(allmax_AR,as.matrix(raster(resf),ncol=nx,nrow=ny))

  resf = paste("../outputbk/cumI_SA_repl_",repl,".bil",sep="")
  allsum_cumI<- allsum_cumI + as.matrix(raster(resf),ncol=nx,nrow=ny)
  allmin_cumI<- pmin(allmin_cumI,as.matrix(raster(resf),ncol=nx,nrow=ny))
  allmax_cumI<- pmax(allmax_cumI,as.matrix(raster(resf),ncol=nx,nrow=ny))

  resf = paste("../outputbk/cumBI_SA_repl_",repl,".bil",sep="")
  allsum_cumBI<- allsum_cumBI + as.matrix(raster(resf),ncol=nx,nrow=ny)
  allmin_cumBI<- pmin(allmin_cumBI,as.matrix(raster(resf),ncol=nx,nrow=ny))
  allmax_cumBI<- pmax(allmax_cumBI,as.matrix(raster(resf),ncol=nx,nrow=ny))

} ### loop through all replicates 

write.csv(allsum_R0,paste("../output2/R0_SA_total_",row.num,"_1.csv",sep=""), row.names=F, quote=F)
write.csv(allsum_AR,paste("../output2/AR_SA_total_",row.num,"_1.csv",sep=""), row.names=F, quote=F)
write.csv(allsum_cumI,paste("../outputbk/cumI_SA_total_",row.num,"_1.csv",sep=""), row.names=F, quote=F)
write.csv(allsum_cumBI,paste("../outputbk/cumBI_SA_total_",row.num,"_1.csv",sep=""), row.names=F, quote=F)

write.csv(matrix(as.vector(unlist(allmin_R0)),2544,2308),paste("../output2/R0_SA_min_",row.num,"_1.csv",sep=""), row.names=F, quote=F)
write.csv(matrix(as.vector(unlist(allmin_AR)),2544,2308),paste("../output2/AR_SA_min_",row.num,"_1.csv",sep=""), row.names=F, quote=F)
write.csv(matrix(as.vector(unlist(allmin_cumI)),2544,2308),paste("../outputbk/cumI_SA_min_",row.num,"_1.csv",sep=""), row.names=F, quote=F)
write.csv(matrix(as.vector(unlist(allmin_cumBI)),2544,2308),paste("../outputbk/cumBI_SA_min_",row.num,"_1.csv",sep=""), row.names=F, quote=F)

write.csv(matrix(as.vector(unlist(allmax_R0)),2544,2308),paste("../output2/R0_SA_max_",row.num,"_1.csv",sep=""), row.names=F, quote=F)
write.csv(matrix(as.vector(unlist(allmax_AR)),2544,2308),paste("../output2/AR_SA_max_",row.num,"_1.csv",sep=""), row.names=F, quote=F)
write.csv(matrix(as.vector(unlist(allmax_cumI)),2544,2308),paste("../outputbk/cumI_SA_max_",row.num,"_1.csv",sep=""), row.names=F, quote=F)
write.csv(matrix(as.vector(unlist(allmax_cumBI)),2544,2308),paste("../outputbk/cumBI_SA_max_",row.num,"_1.csv",sep=""), row.names=F, quote=F)

