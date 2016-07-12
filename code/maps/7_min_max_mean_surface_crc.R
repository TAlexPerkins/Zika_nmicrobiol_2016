bPaths("/afs/crc.nd.edu/user/a/asiraj/R/gdal")
library(rgdal)
.libPaths("/afs/crc.nd.edu/user/a/asiraj/raster")
library(raster)
setwd("/afs/crc.nd.edu/user/a/asiraj/zika/code")
source("0_numfunctions.R")

n.cores <- 1
n.lines <- 1
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

for (repl in 1:10) {
  ###### assemble R0 from cluster
  
  resf = paste("../output2/R0_SA_total_",repl,"_0.csv",sep="")
  allsum_R0<- allsum_R0 + read.csv(resf,header=T)

  resf = paste("../output2/AR_SA_total_",repl,"_0.csv",sep="")
  allsum_AR<- allsum_AR + read.csv(resf,header=T)

  resf = paste("../outputbk/cumI_SA_total_",repl,"_0.csv",sep="")
  allsum_cumI<- allsum_cumI + read.csv(resf,header=T)

  resf = paste("../outputbk/cumBI_SA_total_",repl,"_0.csv",sep="")
  allsum_cumBI<- allsum_cumBI + read.csv(resf,header=T)

  resf = paste("../output2/R0_SA_min_",repl,"_0.csv",sep="")
  allmin_R0<- pmin(matrix(allmin_R0,2544,2308),as.numeric(unlist(read.csv(resf,header=T))))

  resf = paste("../output2/AR_SA_min_",repl,"_0.csv",sep="")
  allmin_AR<- pmin(matrix(allmin_AR,2544,2308),as.numeric(unlist(read.csv(resf,header=T))))

  resf = paste("../outputbk/cumI_SA_min_",repl,"_0.csv",sep="")
  allmin_cumI<- pmin(matrix(allmin_cumI,2544,2308),as.numeric(unlist(read.csv(resf,header=T))))

  resf = paste("../outputbk/cumBI_SA_min_",repl,"_0.csv",sep="")
  allmin_cumBI<- pmin(matrix(allmin_cumBI,2544,2308),as.numeric(unlist(read.csv(resf,header=T))))

  resf = paste("../output2/R0_SA_max_",repl,"_0.csv",sep="")
  allmax_R0<- pmax(matrix(allmax_R0,2544,2308),as.numeric(unlist(read.csv(resf,header=T))))

  resf = paste("../output2/AR_SA_max_",repl,"_0.csv",sep="")
  allmax_AR<- pmax(matrix(allmax_AR,2544,2308),as.numeric(unlist(read.csv(resf,header=T))))

  resf = paste("../outputbk/cumI_SA_max_",repl,"_0.csv",sep="")
  allmax_cumI<- pmax(matrix(allmax_cumI,2544,2308),as.numeric(unlist(read.csv(resf,header=T))))

  resf = paste("../outputbk/cumBI_SA_max_",repl,"_0.csv",sep="")
  allmax_cumBI<- pmax(matrix(allmax_cumBI,2544,2308),as.numeric(unlist(read.csv(resf,header=T))))
  
} ### loop through all replicates 

print(dim(allsum_cumBI))


  ####### # R0
print("generatng mean R0 map...")
  outdf1 = paste("../output2/R0_SA_mean_1k",sep="")
  output2_5_americas(allsum_R0/1000,outdf1)

  ######## # Attack rate
print("generatng mean AR map...")
  outdf1 = paste("../output2/AR_SA_mean_1k",sep="")
  output2_5_americas(allsum_AR/1000,outdf1)

  ###### # Infecteds
print("generatng mean Infecteds map...") 
  outdf1 = paste("../outputbk/cumI_SA_mean_1k",sep="")
  output2_5_americas(allsum_cumI/1000,outdf1)
  
  
  ###### # births by infected
print("generating mean Infected births map...") 
  outdf1 = paste("../outputbk/cumBI_SA_mean_1k",sep="")
  output2_5_americas(allsum_cumBI/1000,outdf1)


  ####### # R0
print("generatng min R0 map...")
  outdf1 = paste("../output2/R0_SA_min_1k",sep="")
  output2_5_americas( matrix(as.vector(unlist(allmin_R0)),2544,2308),outdf1)

  ######## # Attack rate
print("generatng min AR map...")
  outdf1 = paste("../output2/AR_SA_min_1k",sep="")
  output2_5_americas( matrix(as.vector(unlist(allmin_AR)),2544,2308),outdf1)

  ###### # Infecteds
print("generatng min Infecteds map...") 
  outdf1 = paste("../outputbk/cumI_SA_min_1k",sep="")
  output2_5_americas( matrix(as.vector(unlist(allmin_cumI)),2544,2308),outdf1)

  ###### # births by infected
print("generating min Infected births map...") 
  outdf1 = paste("../outputbk/cumBI_SA_min_1k",sep="")
  output2_5_americas( matrix(as.vector(unlist(allmin_cumBI)),2544,2308),outdf1)

 ####### # R0
print("generatng max R0 map...")
  outdf1 = paste("../output2/R0_SA_max_1k",sep="")
  output2_5_americas( matrix(as.vector(unlist(allmax_R0)),2544,2308),outdf1)

  ######## # Attack rate
print("generatng max AR map...")
  outdf1 = paste("../output2/AR_SA_max_1k",sep="")
  output2_5_americas(matrix(as.vector(unlist(allmax_AR)),2544,2308),outdf1)

  ###### # Infecteds
print("generatng max Infecteds map...") 
  outdf1 = paste("../outputbk/cumI_SA_max_1k",sep="")
  output2_5_americas(matrix(as.vector(unlist(allmax_cumI)),2544,2308),outdf1)

  
  ###### # births by infected
print("generating max Infected births map...") 
  outdf1 = paste("../outputbk/cumBI_SA_max_1k",sep="")
  output2_5_americas(matrix(as.vector(unlist(allmax_cumBI)),2544,2308),outdf1)



write.csv(allsum_R0/1000,paste("../output2/R0_SA_mean_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(allsum_AR/1000,paste("../output2/AR_SA_mean_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(allsum_cumI/1000,paste("../outputbk/cumI_SA_mean_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(allsum_cumBI/1000,paste("../outputbk/cumBI_SA_mean_",row.num,"_0.csv",sep=""), row.names=F, quote=F)

write.csv(allmin_R0,paste("../output2/R0_SA_min_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(allmin_AR,paste("../output2/AR_SA_min_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(allmin_cumI,paste("../outputbk/cumI_SA_min_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(allmin_cumBI,paste("../outputbk/cumBI_SA_min_",row.num,"_0.csv",sep=""), row.names=F, quote=F)

write.csv(allmax_R0,paste("../output2/R0_SA_max_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(allmax_AR,paste("../output2/AR_SA_max_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(allmax_cumI,paste("../outputbk/cumI_SA_max_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(allmax_cumBI,paste("../outputbk/cumBI_SA_max_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
