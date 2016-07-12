library(rgdal)
if(!require(scam)) {install.packages('scam'); library(scam)}

setwd("/afs/crc.nd.edu/user/a/asiraj/zika/code")
source("0_numfunctions.R")

newgrid<- as.vector(unlist(read.csv("../generated/americas_serial.csv", header=T)[,1]))
cluster.setup<-read.csv("../data/clustersetup3.csv", header=T)

###########
n.cores <- 13000
n.lines <- 13000
seqP <- seq(1,n.lines,n.lines/n.cores)
xx <- as.numeric(Sys.getenv("SGE_TASK_ID"))
row.num <- seqP[xx]

#registerDoParallel(cores=4)
repl<- cluster.setup[row.num,3]	
print(c("loading functions:", repl))
load("../data/functions_stat.RData")
	rep1 = rep.master[repl,1]
	rep2 = rep.master[repl,2]
#	print(R0)

print(c("replicate",repl))

skiper<-seq(0,1496521,115117)
sid<- cluster.setup[row.num,2]
print(c("reading chunk",sid))
x<- read.csv("../generated/allcovs_v4_compact.csv",skip=skiper[sid], nrows=115117, header=T)
x<- matrix(as.vector(unlist(x)),115117,17)
x[,3]<- as.vector(unlist(read.csv(paste("../generated/aegypti_",rep2,".csv",sep=""), header=T)))[newgrid][(skiper[sid]+1):(skiper[sid+1])]
x<- data.frame(cbind(x[,3],apply(x[,5:16],1,mean, na.rm=T),x[,17]))
names(x)<- c('aegypti','avg_temp','gdp_pcppp2005')
x$factor.adj =1
ik<-which(is.na(x[,3]) | x[,3]==0) # | x[,1]<.681)
x[ik,3]<-NA

print("running AR...")
AR.predicted<-as.vector(pnorm(predict(lm.rand[[rep1]],x)))

print("saving AR...")
write.csv(AR.predicted,paste("../generated/repl_",repl,"_ARx_stat",sid,".csv",sep=""), row.names=F, quote=F)


