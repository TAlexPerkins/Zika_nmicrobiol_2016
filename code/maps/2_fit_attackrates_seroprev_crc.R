if(!require(scam)) {install.packages('scam'); library(scam)}
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
#	load(paste("../generated/functions_R0_AR_1k_prev.RData",sep""))
	load(paste("../generated/functions_R0_AR_1k.RData", sep=""))
	rep1 = rep.master[repl,1]
	rep2 = rep.master[repl,2]
	rep3 = rep.master[repl,3]
#	print(R0)

	##function to obtain predicted adjustment factors for random scam

	scam.pred = function(econ){
	  predict(scam.est.list[[repl]],newdata=data.frame(econ = econ))
	}

	# define functions to calculate R0 and AR
	mort = mort.fun[[rep2]]
	eip = eip.fun[[rep3]]
	R0 = function(data){
	  as.numeric(
	    a ^ 2 * b * c.r * -log(1 - data$aegypti) *
	      exp(scam.pred(log(data$gdp_pcppp2005))) *
	      sapply(1:nrow(data),function(row){mean(sort(sapply(data[
	        row,c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')],
	        function(T){exp(-1/mort(T)*eip(T))*mort(T)}),decreasing=TRUE)[1:6])}))
	}
	R0.vec = seq(0,1000,.01)
	AR.vec = numeric(length(R0.vec))
	for(ii in 1:length(AR.vec)){
	  AR.vec[ii] = 1 - optimize(f=function(S){(S-exp(-R0.vec[ii]*(1-S)))^2},interval=c(0,1))$minimum
	}
	AR.fun = approxfun(R0.vec,AR.vec)
	AR = function(R0.in,h){
	  sapply(R0.in,function(RR)ifelse(RR<=1000,AR.fun(RR^h),1))
	}

	print(c("replicate",repl))

        skiper<-seq(0,1496521,115117)
        sid<- cluster.setup[row.num,2]
        print(c("reading chunk",sid))

        x<- read.csv("../generated/allcovs_v4_compact.csv",skip=skiper[sid], nrows=115117, header=T)
        x<- matrix(as.vector(unlist(x)),115117,17)
        x[,3]<- as.vector(unlist(read.csv(paste("../generated/aegypti_",rep1,".csv",sep=""), header=T)))[newgrid][(skiper[sid]+1):(skiper[sid+1])]
        x<- data.frame(x)
        names(x)<- c('log','lat','aegypti','albopictus','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','gdp_pcppp2005')

        x$factor.adj =1
        ik<-which(is.na(x[,17]) | x[,17]==0) # | x[,3]<.681)
        x[ik,17]<-100
        print(c(a,b,c.r))
        print(c("chunk:",sid, "running R0..."))
        R0.predicted<- R0(x)
        R0.predicted[ik]<-NA
        print("running AR...")
        AR.predicted = AR(R0.predicted,h.list[[repl]]) 

        print("saving R0...")
        write.csv(R0.predicted,paste("../generated/repl_",repl,"_R0x",sid,".csv",sep=""), row.names=F, quote=F)
        print("saving AR...")
        write.csv(AR.predicted,paste("../generated/repl_",repl,"_ARx",sid,".csv",sep=""), row.names=F, quote=F)

