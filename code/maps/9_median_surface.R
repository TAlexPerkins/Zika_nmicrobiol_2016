.libPaths("/afs/crc.nd.edu/user/a/asiraj/R/gdal")
library(rgdal)
.libPaths("/afs/crc.nd.edu/user/a/asiraj/raster")
library(raster)
setwd("/afs/crc.nd.edu/user/a/asiraj/zika/code")
source("0_numfunctions.R")

tenth<-seq(0,1000,10)
n.cores <- 1
n.lines <- 1

print("median of ten cumI, PI and BI values...")

allmed_R0<- NULL
allmed_AR<- NULL
allmed_cumI<- NULL
allmed_cumBI<- NULL

for (repl in (1:100)) {

	allmed_R0<- c(allmed_R0,as.vector(unlist(read.csv(paste("../output2/R0_US_SA_med_",repl,"_1.csv",sep=""), header=T))))
	allmed_AR<- c(allmed_AR,as.vector(unlist(read.csv(paste("../output2/AR_US_SA_med_",repl,"_1.csv",sep=""), header=T))))
	allmed_cumI<- c(allmed_cumI,as.vector(unlist(read.csv(paste("../outputbk/cumI_US_SA_med_",repl,"_1.csv",sep=""), header=T))))
	allmed_cumBI<- c(allmed_cumBI,as.vector(unlist(read.csv(paste("../outputbk/cumBI_US_SA_med_",repl,"_1.csv",sep=""), header=T))))
}
length(allmed_R0)

allmed_R0 <- matrix(allmed_R0[-c(5871553:5872552)], 2544,2308)
allmed_AR <- matrix(allmed_AR[-c(5871553:5872552)], 2544,2308)
allmed_cumI<- matrix(allmed_cumI[-c(5871553:5872552)], 2544,2308)
allmed_cumBI<- matrix(allmed_cumBI[-c(5871553:5872552)], 2544,2308)

write.csv(allmed_R0,paste("../output2/R0_US_SA_med_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(allmed_AR,paste("../output2/AR_US_SA_med_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(allmed_cumI,paste("../outputbk/cumI_US_SA_med_",row.num,"_0.csv",sep=""), row.names=F, quote=F)
write.csv(allmed_cumBI,paste("../outputbk/cumBI_US_SA_med_",row.num,"_0.csv",sep=""), row.names=F, quote=F)

  ####### # R0
print("generatng mean R0 map...")
  outdf1 = paste("../output2/R0_US_SA_median_1k",sep="")
  output2_5_americas(allmed_R0,outdf1)

  ######## # Attack rate
print("generatng mean AR map...")
  outdf1 = paste("../output2/AR_US_SA_median_1k",sep="")
  output2_5_americas(allmed_AR,outdf1)

  ###### # Infecteds
print("generatng mean Infecteds map...") 
  outdf1 = paste("../outputbk/cumI_US_SA_median_1k",sep="")
  output2_5_americas(allmed_cumI,outdf1)
  
  ###### # births by infected
print("generating mean Infected births map...") 
  outdf1 = paste("../outputbk/cumBI_US_SA_median_1k",sep="")
  output2_5_americas(allmed_cumBI,outdf1)



