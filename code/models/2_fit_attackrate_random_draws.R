# load functions with replicates of EIP and mortality functions
load('../../generated/params.RData')

# constant parameters
b = 0.4
c.r = 3.5
a = 1 / 1.5

# R0 as a function of mosquito-human ratio and temperature
R0 = function(m,T,rep2,rep3){
  g = 1 / mort.fun[[rep2]](T)
  e = eip.fun[[rep3]](T)
  m * a ^ 2 * b * c.r * exp(-g * e) / g
}

# load replicate aegypti occurrence probabilities
aegypti.maps.pts.data = as.matrix(read.csv('../../data/aegypti_pt_vals_all.csv')[,-(1:4)])

# seroprevalence data and metadata
x.data = read.csv('../../data/seroprev_metadata.csv')
x.data$factor.adj = 1

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
R0.base = function(row,rep1,rep2,rep3){
  vecs = c(0,0)
  vecs[1]=-log(1-pmax(aegypti.maps.pts.data[row,rep1],x.data$albopictus))
  m.fun = max(vecs)
  mean(sort(sapply(x[
    row,c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')],
    function(TT) R0(m.fun,TT,rep2,rep3)),decreasing=TRUE)[1:6])
}

# attack rate for a given location
AR.avg = function(row,h,rep1,rep2,rep3){
  R0.fun = x$factor.adj[row] * R0.base(row,rep1,rep2,rep3)
  ifelse(R0.fun<1000,AR.fun(R0.fun,h),1)
}

# R0 needed to make attack rate consistent with seroprevalence data
R0.needed = function(row,h){
  optim(par=1,fn=function(par){abs(AR.fun(par,h)-x$seroprev[row])})$par
}


# fit a SCAM for the multiplication factor - economic index relationship
# for each of 1,000 random draws of EIP, mortality, and Aedes occurrence
library(scam)
scam.est.list = list()
h.list = list()
h.list[[reps+1]] = 1
rep.master = cbind(
  sample(ncol(aegypti.maps.pts.data),reps,replace=T),
  1:reps,
  1:reps)
for(rr in 1 : reps){
  repeat{
    if(!is.null(h.list[[rr]])){
      break
    }
    
    # select random combinations of random draws from each of the EIP, occurrence, and mortality distributions
    rep.master[rr,] = c(
      sample(ncol(aegypti.maps.pts.data),1,replace=T),
      sample(reps,1,replace=T),
      sample(reps,1,replace=T))
    rep1 = rep.master[rr,1]
    rep2 = rep.master[rr,2]
    rep3 = rep.master[rr,3]
    
    # get serological data and Aedes aegypti occurrence probabilities by site
    sites = 1:nrow(x.data)
    x = x.data[sites,]
    aegypti.maps.pts = aegypti.maps.pts.data[sites,rep1]
    
    # factor change in R0 necessary to make each attack rate == seroprevalence data for each h value
    h.vec = seq(.01,1,.01)
    factor.needed = matrix(0,nrow(x),length(h.vec))
    for(hh in 1:length(h.vec)){
      factor.needed[,hh] =
        sapply(1:nrow(x),function(rr)R0.needed(rr,h.vec[hh])) /
        sapply(1:nrow(x),function(rr)R0.base(rr,rep1,rep2,rep3))
    }
  
    # fit SCAMs for each h value in the range from 0.01 to 1 and record attack rate residuals
    res.vec = numeric(length=length(h.vec))
    for(hh in 1:length(h.vec)){
      h = h.vec[hh]
      df = data.frame(econ=log(x$gdp_pcppp2005),factor=log(factor.needed[,hh]))
      x$factor.adj = 1
      if(sum(is.infinite(df$factor))){
        next
      }
      scam.est = scam(factor~s(econ,bs='mdcx'),data=df)
      x$factor.adj = exp(scam.est$fitted.values)
      res.vec[hh] = sum((sapply(1:nrow(x),function(rr)AR.avg(rr,h,rep1,rep2,rep3)) - x$seroprev)^2)
    }
    
    # select h value that minimizes residuals and re-fit associated SCAM
    h = h.vec[which.min(res.vec)]
    df = data.frame(econ=log(x$gdp_pcppp2005),factor=log(factor.needed[,which.min(res.vec)]))
    x$factor.adj = 1
    if(sum(is.infinite(df$factor))){
      next
    }
    scam.est = scam(factor~s(econ,bs='mdcx'),data=df)
    x$factor.adj = exp(scam.est$fitted.values)
    
    # store best-fit SCAM model and associated h value
    scam.est.list[[rr]] = scam.est
    h.list[[rr]] = h
  }
}


# save data
save(reps,rep.master,scam.est.list,mort.fun,eip.fun,a,b,c.r,h.list,file='../../generated/functions_R0_AR_random_draws.RData')
