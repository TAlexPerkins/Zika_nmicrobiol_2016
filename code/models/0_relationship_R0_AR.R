# load libraries and seroprevalence data and metadata
library(scam)
x = read.csv('../../data/seroprev_metadata.csv')
x$factor.adj = 1

# relationship between temperature and mortality
load('../../data/algam_85re.Rdata')
lifespan = function(Temperature,se.mult){
  dd = seq(0,120,length.out=(120*24+2))
  nwdd = data.frame(Days=dd,Temperature=rep(Temperature, (120*24+2)), Study_number=5, Feed_B=2, Feed_S=1)
  nwdd = cbind(nwdd,logDay=log(nwdd$Days+1), logTemp=log(nwdd$Temperature+1))  ## +1 avoids log(0) 
  prediction = as.vector(unlist(predict(algam,newdata = (nwdd), se.fit = TRUE, type = "response")$fit))
  prediction.se = as.vector(unlist(predict(algam,newdata = (nwdd), se.fit = TRUE, type = "response")$se.fit))
  prediction = prediction + se.mult * prediction.se
  prediction = prediction[-1]
  prediction[1:24] = prediction[1:24]/prediction[1]
  prediction[which(prediction>1)] = 1
  prediction[which(prediction<=0.001)] = 0
  diffDeath = (prediction[1: (length(prediction)-1)] - prediction[2:length(prediction)])
  diffDeath = diffDeath/sum(diffDeath)
  return(pmax(0, sum(dd[2:length(prediction)]*diffDeath)))
}
fielddata = rbind(c(20, 34, 1 - 0.91))
tgrd = 0.1
lifespan.range = sapply(seq(fielddata[1],fielddata[2],tgrd), function (tr) {lifespan(tr,0)})
fieldcorxn = fielddata[3] - 1/mean(lifespan.range)
mort.fun = approxfun(seq(0,50,tgrd), sapply(seq(0,50,tgrd), function(TT) 1/(1/lifespan(TT,0)+fieldcorxn)))

# relationship between temperature and extrinsic incubation period
eip = function(T){exp(8 - .2 * T)}

# constant parameters
b = 0.4
c.r = 3.5
a = 1 / 1.5

# R0 as a function of mosquito-human ratio and temperature
R0 = function(m,T){
  g = 1 / mort.fun(T)
  m * a ^ 2 * b * c.r * exp(-g * eip(T)) / g
}

# attack rate as a function of R0
R0.vec = seq(0,1000,.01)
AR.vec = numeric(length(R0.vec))
for(ii in 1:length(AR.vec)){
  AR.vec[ii] = 1 - optimize(f=function(S){(S-exp(-R0.vec[ii]*(1-S)))^2},interval=c(0,1))$minimum
}
AR.fun.vec = approxfun(R0.vec,AR.vec)
AR.fun = function(R0,h){
  if(R0 < 1){
    return(0)
  }
  if(R0 >= 1){
    return(AR.fun.vec(R0 ^ h))
  }
}

# R0 for a given location
R0.base = function(row){
  vecs = c(0,0)
  vecs[1]=-log(1-max(x$aegypti[row],x$albopictus[row]))
  m.fun = max(vecs)
  mean(sort(sapply(x[
    row,c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')],
    function(TT) R0(m.fun,TT)),decreasing=TRUE)[1:6])
}

# attack rate for a given location
AR.avg = function(row,h){
  R0.fun = x$factor.adj[row] * R0.base(row)
  ifelse(R0.fun<1000,AR.fun(R0.fun,h),1)
}

# R0 needed to make attack rate consistent with seroprevalence data
R0.needed = function(row,h){
  optim(par=1,fn=function(par){abs(AR.fun(par,h)-x$seroprev[row])})$par
}

# factor change in R0 necessary to make each attack rate == seroprevalence data
h.vec = seq(.01,1,.01)
factor.needed = matrix(0,nrow(x),length(h.vec))
for(hh in 1:length(h.vec)){
  factor.needed[,hh] =
    sapply(1:nrow(x),function(rr)R0.needed(rr,h.vec[hh])) /
    sapply(1:nrow(x),function(rr)R0.base(rr))
}


# plot match between attack rates and seroprevalence data for given R0
pdf(file='../../outputs/attackrate_vs_R0.pdf',width=6.5,height=2.25)
  layout(matrix(1:3,1,3))
  par(oma=c(2.75,4,0,0),mar=c(1,.5,.5,.5))
  
  h = 1
  df = data.frame(econ=log(x$gdp_pcppp2005),factor=log(factor.needed[,ncol(factor.needed)]))
  x$factor.adj = 1
  R0.predicted = sapply(1:nrow(x),function(rr) x$factor.adj[rr] * R0.base(rr))
  plot(seq(0,25,.01), sapply(seq(0,25,.01),function(RR)AR.fun(RR,h)),
       type='l',xlim=c(0,7),ylim=c(0,1),xlab='',ylab='',main='',las=1)
  segments(
    R0.predicted,sapply(1:nrow(x),function(rr)AR.avg(rr,h)),
    R0.predicted,x$seroprev,col=rgb(0,0,0,.5))
  points(R0.predicted,x$seroprev,pch=19,cex=.5)
  mtext('a',side=3,at=.25,line=-1.5)
  mtext('Epidemic attack rate',side=2,line=2.75)
  
  scam.est = scam(factor~s(econ,bs='mdcx'),data=df)
  x$factor.adj = exp(scam.est$fitted.values)
  R0.predicted = sapply(1:nrow(x),function(rr) x$factor.adj[rr] * R0.base(rr))
  plot(seq(0,25,.01), sapply(seq(0,25,.01),function(RR)AR.fun(RR,h)),
       type='l',xlim=c(0,7),ylim=c(0,1),xlab='',ylab='',main='',las=1,yaxt='n')
  segments(
    R0.predicted,sapply(1:nrow(x),function(rr)AR.avg(rr,h)),
    R0.predicted,x$seroprev,col=rgb(0,0,0,.5))
  points(R0.predicted,x$seroprev,pch=19,cex=.5)
  mtext('b',side=3,at=.25,line=-1.5)
  mtext(expression('Basic reproduction number '*R[0]),1,line=2.75)

  res.vec = numeric(length=length(h.vec))
  for(hh in 1:length(h.vec)){
    h = h.vec[hh]
    df = data.frame(econ=log(x$gdp_pcppp2005),factor=log(factor.needed[,hh]))
    x$factor.adj = 1
    scam.est = scam(factor~s(econ,bs='mdcx'),data=df)
    x$factor.adj = exp(scam.est$fitted.values)
    res.vec[hh] = sum((sapply(1:nrow(x),function(rr)AR.avg(rr,h)) - x$seroprev)^2)
  }
  h = h.vec[which.min(res.vec)]
  df = data.frame(econ=log(x$gdp_pcppp2005),factor=log(factor.needed[,which.min(res.vec)]))
  x$factor.adj = 1
  scam.est = scam(factor~s(econ,bs='mdcx'),data=df)
  x$factor.adj = exp(scam.est$fitted.values)
  
  plot(seq(0,25,.01),sapply(seq(0,25,.01),function(RR)AR.fun(RR,h)),
       type='l',xlim=c(0,18),ylim=c(0,1),xlab='',ylab='',main='',las=1,yaxt='n')
  R0.predicted = sapply(1:nrow(x),function(rr) x$factor.adj[rr] * R0.base(rr))
  segments(
    R0.predicted,sapply(1:nrow(x),function(rr)AR.avg(rr,h)),
    R0.predicted,x$seroprev,col=rgb(0,0,0,.5))
  points(R0.predicted,x$seroprev,pch=19,cex=.5)
  mtext('c',side=3,at=.25,line=-1.5)
dev.off()


# gather up functions and other stuff to save for showing relationships in Figure 4
R0 = function(data){
  as.numeric(a ^ 2 * b * c.r *
               -log(1-data$aegypti) *
               exp(predict(scam.est,list(econ=log(data$gdp_pcppp2005)))) *
               sapply(1:nrow(data),function(row){mean(sort(sapply(data[
                 row,c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')],
                 function(T){exp(-1/mort.fun(T)*eip(T))*mort.fun(T)}),decreasing=TRUE)[1:6])}))
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

save(R0,AR,AR.fun,AR.fun,scam.est,mort.fun,eip,a,b,c.r,h,file='../../generated/functions_R0_AR_mean.RData')
