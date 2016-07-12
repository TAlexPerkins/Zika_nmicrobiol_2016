# number of replicates
reps = 1000


# get random samples from confidence interval for EIP-temperature relationship
library(pomp)
beta.0 = c(6,10)
beta.T = c(-.29,-.12)
scalars = sobol(vars = list(beta.0 = beta.0, beta.T = beta.T),1e5)
eip = function(T,beta.0,beta.T){exp(beta.0 + beta.T * T)}
eip.vec = sapply(1:nrow(scalars),function(ii)eip(30,scalars[ii,1],scalars[ii,2]))
eip.sd = optimize(
  f = function(sd){(qnorm(0.05,6.1,sd) - 3.4) ^ 2 + (qnorm(0.95,6.1,sd) - 9.9) ^ 2},
  interval = c(0,10)
)$minimum
beta.samples = matrix(0,reps,2)
eip.fun = list()
for(rr in 1:reps){
  repeat{
    eip.draw = rnorm(1,6.1,eip.sd)
    if(eip.draw > 0){
      break
    }    
  }
  which.range = which(eip.vec > (eip.draw - 0.05) & eip.vec < (eip.draw + 0.05))
  beta.samples[rr,] = colMeans(scalars[which.range,])
  eip.fun[[rr]] =  approxfun(seq(0,50,.1), sapply(seq(0,50,.1), function(T){exp(beta.samples[rr,1] + beta.samples[rr,2] * T)}))
}


# get random samples from confidence interval of mortality-temperature relationship
load('../../data/algam_85re.Rdata')
lifespan <- function(Temperature,se.mult){
  dd <- seq(0,120,length.out=(120*24+2))
  nwdd <- data.frame(Days=dd,Temperature=rep(Temperature, (120*24+2)), Study_number=5, Feed_B=2, Feed_S=1)
  nwdd <- cbind(nwdd,logDay=log(nwdd$Days+1), logTemp=log(nwdd$Temperature+1))  ## +1 avoids log(0) 
  prediction <- as.vector(unlist(predict(algam,newdata = (nwdd), se.fit = TRUE, type = "response")$fit))
  prediction.se <- as.vector(unlist(predict(algam,newdata = (nwdd), se.fit = TRUE, type = "response")$se.fit))
  prediction <- prediction + se.mult * prediction.se
  prediction <- prediction[-1]
  prediction[1:24] <- prediction[1:24]/prediction[1]
  prediction[which(prediction>1)] <- 1
  prediction[which(prediction<=0.001)] <- 0
  diffDeath <- (prediction[1: (length(prediction)-1)] - prediction[2:length(prediction)])
  diffDeath <- diffDeath/sum(diffDeath)
  return(pmax(0, sum(dd[2:length(prediction)]*diffDeath)))
}
fielddata <- rbind(c(20, 34, 1 - 0.91))
tgrd <- 0.1
lifespan.range <- sapply(seq(fielddata[1],fielddata[2],tgrd), function (tr) {lifespan(tr,0)})
fieldcorxn <- fielddata[3] - 1/mean(lifespan.range)
mort.fun = list()
for(rr in 1 : reps){
  normdev = rnorm(1,0,1)
  mort.fun[[rr]] = approxfun(seq(0,50,tgrd), sapply(seq(0,50,tgrd), function(TT) 1/(1/lifespan(TT,normdev)+fieldcorxn)))
}


# save random draws from EIP and mortality relationships with temperature
save(reps,eip.fun,mort.fun,file='../../generated/params.RData')
