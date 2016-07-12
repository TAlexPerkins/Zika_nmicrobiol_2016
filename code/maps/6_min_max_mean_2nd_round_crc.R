#.libPaths("/afs/crc.nd.edu/user/a/asiraj/R/gdal")
library(rgdal)
.libPaths("/afs/crc.nd.edu/user/a/asiraj/raster")
library(raster)
setwd("/afs/crc.nd.edu/user/a/asiraj/zika/code")
source("0_numfunctions.R")

tenth<-seq(0,100,10)
n.cores <- 10
n.lines <- 10
seqP <- seq(1,n.lines,n.lines/n.cores)
xx <- as.numeric(Sys.getenv("SGE_TASK_ID"))
row.num <- seqP[xx] 

allsum_cumI<- 0
allsum_cumBI<- 0
allsum_R0<-0
allsum_AR<-0

allmin_cumI<- 99
allmin_cumBI<- 99
allmin_R0<-99
allmin_AR<-99

allmax_cumI<- 0
allmax_cumBI<- 0
allmax_R0<-0
allmax_AR<-0

print("sum of ten cumI, PI and BI values...")

for (repl in (tenth[row.num]+1):(tenth[row.num+1])) {
  ###### assemble R0 from cluster
  
  resf = paste("../output2/R0_SA_total_",repl,"_1.csv",sep="")
  allsum_R0<- allsum_R0 + read.csv(resf,header=T)

  resf = paste("../output2/AR_SA_total_",repl,"_1.csv",sep="")
  allsum_AR<- allsum_AR + read.csv(resf,header=T)

  resf = paste("../outputbk/cumI_SA_total_",repl,"_1.csv",sep="")
  allsum_cumI<- allsum_cumI + read.csv(resf,header=T)

  resf = paste("../outputbk/cumBI_SA_total_",repl,"_1.csv",sep="")
  allsum_cumBI<- allsum_cumBI + read.csv(resf,header=T)
  
  resf = paste("../output2/R0_SA_min_",repl,"_1.csv",sep="")
  allmin_R0<- pmin(matrix(allmin_R0,2544,2308) ,as.numeric(unlist(read.csv(resf,header=T))))

  resf = paste("../output2/AR_SA_min_",repl,"_1.csv",sep="")
  allmin_AR<- pmin(matrix(allmin_AR,2544,2308) ,as.numeric(unlist(read.csv(resf,header=T))))

  resf = paste("../outputbk/cumI_SA_min_",repl,"_1.csv",sep="")
  allmin_cumI<- pmin(matrix(allmin_cumI,2544,2308) ,as.numeric(unlist(read.csv(resf,header=T))))

  resf = paste("../outputbk/cumBI_SA_min_",repl,"_1.csv",sep="")
  allmin_cumBI<- pmin(matrix(allmin_cumBI,2544,2308),as.numeric(unlist(read.csv(resf,header=T))))


  resf = paste("../output2/R0_SA_max_",repl,"_1.csv",sep="")
  allmax_R0<- pmax(matrix(allmax_R0,2544,2308),as.numeric(unlist(read.csv(resf,header=T))))

  resf = paste("../output2/AR_SA_max_",repl,"_1.csv",sep="")
  allmax_AR<- pmax(matrix(allmax_AR,2544,2308),as.numeric(unlist(read.csv(resf,header=T))))


  resf = paste("../outputbk/cumI_SA_max_",repl,"_1.csv",sep="")
  allmax_cumI<- pmax(matrix(allmax_cumI,2544,2308) ,as.numeric(unlist(read.csv(resf,header=T))))

  resf = paste("../outputbk/cumBI_SA_max_",repl,"_1.csv",sep="")
  allmax_cumBI<- pmax(matrix(allmax_cumBI,2544,2308) ,as.numeric(unlist(read.csv(resf,header=T))))


} ### loop through all cumulative replicates 

write.csv(allsum_R0,paste("../output2/R0_SA_total_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(allsum_AR,paste("../output2/AR_SA_total_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(allsum_cumI,paste("../outputbk/cumI_SA_total_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(allsum_cumBI,paste("../outputbk/cumBI_SA_total_",row.num,"_0.csv",sep=""), row.names=F, quote=F)

write.csv(matrix(allmin_R0,2544,2308),paste("../output2/R0_SA_min_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(matrix(allmin_AR,2544,2308),paste("../output2/AR_SA_min_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(matrix(allmin_cumI,2544,2308),paste("../outputbk/cumI_SA_min_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(matrix(allmin_cumBI,2544,2308),paste("../outputbk/cumBI_SA_min_",row.num,"_0.csv",sep=""), row.names=F, quote=F)

write.csv(matrix(allmax_R0,2544,2308),paste("../output2/R0_SA_max_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(matrix(allmax_AR,2544,2308),paste("../output2/AR_SA_max_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(matrix(allmax_cumI,2544,2308),paste("../outputbk/cumI_SA_max_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(matrix(allmax_cumBI,2544,2308),paste("../outputbk/cumBI_SA_max_",row.num,"_0.csv",sep=""), row.names=F, quote=F)


