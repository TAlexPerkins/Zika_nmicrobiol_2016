# load functions with replicates of EIP and mortality functions
load('../../generated/params.RData')
load('../../generated/functions_R0_AR_random_draws.RData')

# constant parameters
b = 0.4
c.r = 3.5
a = 1 / 1.5

# R0 as a function of mosquito-human ratio and temperature
R0 = function(m,T,rep2,rep3){
  g = mort.fun[[rep2]](T)
  e = eip.fun[[rep3]](T)
  m * a ^ 2 * b * c.r * exp(-g * e) / g
}

# load replicate aegypti occurrence probabilities
aegypti.maps.pts.data = as.matrix(read.csv('../../data/aegypti_pt_vals_all.csv')[,-(1:4)])

# load replicates of infections and such
proj = read.csv('../../outputs/country_output_model.csv')
proj.inf = proj[,1+(1:1000)]
proj.inf = matrix(as.numeric(unlist(proj.inf)),nrow(proj.inf),ncol(proj.inf))
proj.births = proj[,1+1000+(1:1000)]
proj.births = matrix(as.numeric(unlist(proj.births)),nrow(proj.births),ncol(proj.births))
countries = as.character(proj$NAME)


# make Figure 2
pdf('../../outputs/posteriors_mechanistic.pdf',width=6.5,5)

  layout(t(matrix(1:15,5,3)))
  par(oma=c(0,2,1.5,0.5))
  par(mar=c(4.5,3,0,0.5))
  
  # occurrence probability
  occur.hist = matrix(0,2,length(seq(0,1,.01))-1)
  occur.mids = occur.hist
  for(ii in 1:2){
    occur.hist[ii,] =
      hist(aegypti.maps.pts.data[ii,],breaks=seq(0,1,.01),plot=F)$density /
      ncol(aegypti.maps.pts.data)
    occur.mids[ii,] =
      hist(aegypti.maps.pts.data[ii,],breaks=seq(0,1,.01),plot=F)$mids
  }
  plot(-100,-100,type='l',xlim=c(.7,1),ylim=c(0,.16),xlab='',ylab='')
  polygon(c(occur.mids[1,],0),c(occur.hist[1,],0),col=rgb(0,0,1,.3))
  polygon(c(occur.mids[2,],0),c(occur.hist[2,],0),col=rgb(0,1,0,.3))
  mtext('Occurrence prob.',1,cex=.7,line=2)
  mtext('Density',2,cex=.7,line=1.9)
  mtext('Model parameters',2,cex=.7,line=3.4)
  mtext('a',3,at=.7,cex=.7)
  
  # lifespan
  temp = seq(15,35,0.01)
  plot(temp,mort.fun[[1]](temp),type='l',col=rgb(0,0,0,.2),ylim=c(4,12),xlab='',ylab='')
  for(ii in 2:length(mort.fun)){
    lines(temp,mort.fun[[ii]](temp),col=rgb(0,0,0,.2))
  }
  mtext('Temperature',1,cex=.7,line=2)
  mtext('Mean lifespan',2,cex=.7,line=1.9)
  mtext('b',3,at=15,cex=.7)
  
  # EIP
  temp = seq(15,35,0.01)
  plot(temp,eip.fun[[1]](temp),type='l',col=rgb(0,0,0,.05),xlab='',ylab='')
  for(ii in 2:length(eip.fun)){
    lines(temp,eip.fun[[ii]](temp),col=rgb(0,0,0,.05))
  }
  mtext('Temperature',1,cex=.7,line=2)
  mtext('EIP',2,cex=.7,line=1.9)
  mtext('c',3,at=15,cex=.7)
  
  # economics
  library(scam)
  econ = seq(6.5,10.5,0.1)
  plot(
    econ,predict(scam.est.list[[1]],newdata=data.frame(econ=econ)),
    type='l',col=rgb(0,0,0,.05),ylim=c(-1,6),xlab='',ylab='')
  for(ii in 2:length(scam.est.list)){
    lines(econ,predict(scam.est.list[[ii]],newdata=data.frame(econ=econ)),col=rgb(0,0,0,.1))
  }
  mtext('ln(econ. index)',1,cex=.7,line=2)
  mtext('ln(m multiplier)',2,cex=.7,line=2)
  mtext('d',3,at=6.5,cex=.7)
  
  # alpha (or h here in code)
  R0.vec = seq(0,1000,.01)
  AR.vec = numeric(length(R0.vec))
  for(ii in 1:length(AR.vec)){
    AR.vec[ii] = 1 - optimize(f=function(S){(S-exp(-R0.vec[ii]*(1-S)))^2},interval=c(0,1))$minimum
  }
  AR.fun = approxfun(R0.vec,AR.vec)
  AR = function(R0.in,h){
    sapply(R0.in,function(RR)ifelse(RR<=1000,AR.fun(RR^h),1))
  }
  R0.in = seq(0,15,.1)
  plot(R0.in,AR(R0.in,h.list[[1]]),col=rgb(0,0,0,.1),type='l',ylim=c(0,.6),xlab='',ylab='')
  for(ii in 2:(length(h.list)-1)){
    lines(R0.in,AR(R0.in,jitter(h.list[[ii]],factor=2)),col=rgb(0,0,0,.1))
  }
  mtext(expression(R[0]),1,cex=.7,line=2)
  mtext('Epidemic attack rate',2,cex=.7,line=1.9)
  mtext('e',3,at=0,cex=.7)
  
  par(mar=c(5,0,0,.5))
  
  # infections by country
  # all
  plot(
    -100,-100,type='l',xlim=range(hist(colSums(proj.inf[,]/1e6),breaks=50,plot=F)$mids),
    ylim=c(0,max(hist(colSums(proj.inf[,]/1e6),breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
  polygon(
    c(hist(colSums(proj.inf[,]/1e6),breaks=50,plot=F)$mids,rev(hist(colSums(proj.inf[,]/1e6),breaks=50,plot=F)$mids)),
    c(hist(colSums(proj.inf[,]/1e6),breaks=50,plot=F)$density,rep(0,length(hist(colSums(proj.inf[,]/1e6),breaks=50,plot=F)$mids))),
    col=rgb(1,0,0,.5))
  axis(1,at=seq(85,115,10),labels=c('','','',''))
  text(seq(85,115,10),par("usr")[3]-.008, labels = c('85.0','95.0','105.0','115.0'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
  mtext('Normalized density',2,cex=.7,line=.4)
  mtext('f',3,at=min(hist(colSums(proj.inf[,]/1e6),breaks=50,plot=F)$mids),cex=.7)
  mtext('Total',3,cex=.7)
  
  # Brazil
  plot(
    -100,-100,type='l',xlim=range(hist(proj.inf[11,]/1e6,breaks=50,plot=F)$mids),
    ylim=c(0,max(hist(proj.inf[11,]/1e6,breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
  polygon(
    c(hist(proj.inf[11,]/1e6,breaks=50,plot=F)$mids,rev(hist(proj.inf[11,]/1e6,breaks=50,plot=F)$mids)),
    c(hist(proj.inf[11,]/1e6,breaks=50,plot=F)$density,rep(0,length(hist(proj.inf[11,]/1e6,breaks=50,plot=F)$mids))),
    col=rgb(1,0,0,.5))
  axis(1,at=seq(35,45,5),labels=c('','',''))
  text(seq(35,45,5),par("usr")[3]-.015, labels = c('35.0','40.0','45.0'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
  mtext('g',3,at=min(hist(proj.inf[11,]/1e6,breaks=50,plot=F)$mids),cex=.7)
  mtext('Brazil',3,cex=.7)
  
  # Venezuela
  plot(
    -100,-100,type='l',xlim=range(hist(proj.inf[49,]/1e6,breaks=50,plot=F)$mids),
    ylim=c(0,max(hist(proj.inf[49,]/1e6,breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
  polygon(
    c(hist(proj.inf[49,]/1e6,breaks=50,plot=F)$mids,rev(hist(proj.inf[49,]/1e6,breaks=50,plot=F)$mids)),
    c(hist(proj.inf[49,]/1e6,breaks=50,plot=F)$density,rep(0,length(hist(proj.inf[49,]/1e6,breaks=50,plot=F)$mids))),
    col=rgb(1,0,0,.5))
  axis(1,at=seq(7,8,.5),labels=c('','',''))
  text(seq(7,8,.5),par("usr")[3]-.11, labels = c('7.0','7.5','8.0'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
  mtext('Infections among all people (millions)',1,cex=.7,line=2.25)
  mtext('h',3,at=min(hist(proj.inf[49,]/1e6,breaks=50,plot=F)$mids),cex=.7)
  mtext('Venezuela',3,cex=.7)

  # Colombia
  plot(
    -100,-100,type='l',xlim=range(hist(proj.inf[16,]/1e6,breaks=50,plot=F)$mids),
    ylim=c(0,max(hist(proj.inf[16,]/1e6,breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
  polygon(
    c(hist(proj.inf[16,]/1e6,breaks=50,plot=F)$mids,rev(hist(proj.inf[16,]/1e6,breaks=50,plot=F)$mids)),
    c(hist(proj.inf[16,]/1e6,breaks=50,plot=F)$density,rep(0,length(hist(proj.inf[16,]/1e6,breaks=50,plot=F)$mids))),
    col=rgb(1,0,0,.5))
  axis(1,at=seq(6.5,7.5,.5),labels=c('','',''))
  text(seq(6.5,7.5,.5),par("usr")[3]-.15, labels = c('6.5','7.0','7.5'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
  mtext('i',3,at=min(hist(proj.inf[16,]/1e6,breaks=50,plot=F)$mids),cex=.7)
  mtext('Colombia',3,cex=.7)
  
  # Puerto Rico
  plot(
    -100,-100,type='l',xlim=range(hist(proj.inf[39,]/1e6,breaks=50,plot=F)$mids),
    ylim=c(0,max(hist(proj.inf[39,]/1e6,breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
  polygon(
    c(hist(proj.inf[39,]/1e6,breaks=50,plot=F)$mids,rev(hist(proj.inf[39,]/1e6,breaks=50,plot=F)$mids)),
    c(hist(proj.inf[39,]/1e6,breaks=50,plot=F)$density,rep(0,length(hist(proj.inf[39,]/1e6,breaks=50,plot=F)$mids))),
    col=rgb(1,0,0,.5))
  axis(1,at=seq(.9,1.1,.1),labels=c('','',''))
  text(seq(.9,1.1,.1),par("usr")[3]-1.1, labels = c('0.9','1.0','1.1'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
  mtext('j',3,at=min(hist(proj.inf[39,]/1e6,breaks=50,plot=F)$mids),cex=.7)
  mtext('Puerto Rico',3,cex=.7)
  
  # infections among births by country
  # all
  plot(
    -100,-100,type='l',xlim=range(hist(colSums(proj.births[,]/1e6),breaks=50,plot=F)$mids),
    ylim=c(0,max(hist(colSums(proj.births[,]/1e6),breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
  polygon(
    c(hist(colSums(proj.births[,]/1e6),breaks=50,plot=F)$mids,rev(hist(colSums(proj.births[,]/1e6),breaks=50,plot=F)$mids)),
    c(hist(colSums(proj.births[,]/1e6),breaks=50,plot=F)$density,rep(0,length(hist(colSums(proj.births[,]/1e6),breaks=50,plot=F)$mids))),
    col=rgb(0,0,0,.5))
  axis(1,at=seq(1.6,2,.2),labels=c('','',''))
  text(seq(1.6,2,.2),par("usr")[3]-.5, labels = c('1.6','1.8','2.0'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
  mtext('Normalized density',2,cex=.7,line=.5)
  mtext('k',3,at=min(hist(colSums(proj.births[,]/1e6),breaks=50,plot=F)$mids),cex=.7)
  mtext('Total',3,cex=.7)
  
  # Brazil
  plot(
    -100,-100,type='l',xlim=range(hist(proj.births[11,]/1e6,breaks=50,plot=F)$mids),
    ylim=c(0,max(hist(proj.births[11,]/1e6,breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
  polygon(
    c(hist(proj.births[11,]/1e6,breaks=50,plot=F)$mids,rev(hist(proj.births[11,]/1e6,breaks=50,plot=F)$mids)),
    c(hist(proj.births[11,]/1e6,breaks=50,plot=F)$density,rep(0,length(hist(proj.births[11,]/1e6,breaks=50,plot=F)$mids))),
    col=rgb(0,0,0,.5))
  axis(1,at=seq(.5,.7,.1),labels=c('','',''))
  text(seq(.5,.7,.1),par("usr")[3]-1.5, labels = c('0.5','0.6','0.7'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
  mtext('l',3,at=min(hist(proj.births[11,]/1e6,breaks=50,plot=F)$mids),cex=.7)
  mtext('Brazil',3,cex=.7)
  
  # Venezuela
  plot(
    -100,-100,type='l',xlim=range(hist(proj.births[49,]/1e6,breaks=50,plot=F)$mids),
    ylim=c(0,max(hist(proj.births[49,]/1e6,breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
  polygon(
    c(hist(proj.births[49,]/1e6,breaks=50,plot=F)$mids,rev(hist(proj.births[49,]/1e6,breaks=50,plot=F)$mids)),
    c(hist(proj.births[49,]/1e6,breaks=50,plot=F)$density,rep(0,length(hist(proj.births[49,]/1e6,breaks=50,plot=F)$mids))),
    col=rgb(0,0,0,.5))
  axis(1,at=seq(.14,.16,.01),labels=c('','',''))
  text(seq(.14,.16,.01),par("usr")[3]-4.5, labels = c('0.14','0.15','0.16'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
  mtext('Infections among childbearing women (millions)',1,cex=.7,line=2.2)
  mtext('m',3,at=min(hist(proj.births[49,]/1e6,breaks=50,plot=F)$mids),cex=.7)
  mtext('Venezuela',3,cex=.7)

  # Colombia
  plot(
    -100,-100,type='l',xlim=range(hist(proj.births[16,]/1e6,breaks=50,plot=F)$mids),
    ylim=c(0,max(hist(proj.births[16,]/1e6,breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
  polygon(
    c(hist(proj.births[16,]/1e6,breaks=50,plot=F)$mids,rev(hist(proj.births[16,]/1e6,breaks=50,plot=F)$mids)),
    c(hist(proj.births[16,]/1e6,breaks=50,plot=F)$density,rep(0,length(hist(proj.births[16,]/1e6,breaks=50,plot=F)$mids))),
    col=rgb(0,0,0,.5))
  axis(1,at=seq(.12,.15,.01),labels=c('','','',''))
  text(seq(.12,.15,.01),par("usr")[3]-10, labels = c('0.12','0.13','0.14','0.15'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
  mtext('n',3,at=min(hist(proj.births[16,]/1e6,breaks=50,plot=F)$mids),cex=.7)
  mtext('Colombia',3,cex=.7)
  
  # Puerto Rico
  plot(
    -100,-100,type='l',xlim=range(hist(proj.births[39,]/1e6,breaks=50,plot=F)$mids),
    ylim=c(0,max(hist(proj.births[39,]/1e6,breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
  polygon(
    c(hist(proj.births[39,]/1e6,breaks=50,plot=F)$mids,rev(hist(proj.births[39,]/1e6,breaks=50,plot=F)$mids)),
    c(hist(proj.births[39,]/1e6,breaks=50,plot=F)$density,rep(0,length(hist(proj.births[39,]/1e6,breaks=50,plot=F)$mids))),
    col=rgb(0,0,0,.5))
  axis(1,at=seq(.0105,.0135,.001),labels=c('','','',''))
  text(seq(.0105,.0135,.001),par("usr")[3]-70, labels = c('0.0105','0.0115','0.0125','0.0135'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
  mtext('o',3,at=min(hist(proj.births[39,]/1e6,breaks=50,plot=F)$mids),cex=.7)
  mtext('Puerto Rico',3,cex=.7)
  
dev.off()
