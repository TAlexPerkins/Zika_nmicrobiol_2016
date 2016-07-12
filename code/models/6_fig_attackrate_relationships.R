# load and process data
library(scam)
load('../../generated/functions_R0_AR_mean.RData')
load('../../generated/functions_stat_mean.RData')
x = read.csv('../../data/seroprev_metadata.csv')
x$aegypti = pmax(x$aegypti,x$albopictus)
x$temp = rowMeans(as.matrix(x[,8:19]))
x$avg_temp = rowMeans(as.matrix(x[,8:19]))
x$seroprev[x$seroprev==0] = 1e-3
x$seroprobit = qnorm(x$seroprev)


# calculate projected attack rates under mechanistic and statistical models
R0 = function(data){
  as.numeric(a ^ 2 * b * c.r *
               -log(1 - data$aegypti) *
               exp(predict(scam.est,list(econ=data$gdp_pcppp2005))) *
               exp(-1/mort.fun(data$temp)*eip(data$temp))*mort.fun(data$temp))
}
df = expand.grid(
  gdp_pcppp2005 = seq(6.5,10.5,.5),
  aegypti = seq(.35,.95,length.out=3),
  temp = seq(18,38,.1)
)
df$R0 = R0(df)
df$AR.model = AR.fun(df$R0^h)
df$avg_temp = df$temp
df$AR.stat.full = pnorm(predict(lm.est,newdata=df))
df$AR.stat.step = pnorm(predict(lm.step,newdata=df))


# make Figure 4
pdf('../../outputs/AR_relationships.pdf',width=6.5,height=3.5)
layout(t(matrix(1:9,3,3)))
par(oma=c(4,2,2,0),mar=c(0,2,.25,.25))

cols = log(unique(df$gdp_pcppp2005))
cols = (cols - min(cols)) / (max(cols) - min(cols))
cols = rgb(1-cols,cols,0,.5)

d = df[df$aegypti==unique(df$aegypti)[1],]
plot(
  unique(d$temp),d$AR.model[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[1]],
  type='l',las=1,col=cols[1],lwd=2,ylim=c(0,1),xaxt='n')
for(ii in 2:length(unique(d$gdp_pcppp2005))){
  lines(unique(d$temp),d$AR.model[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[ii]],col=cols[ii],lwd=2)
}
mtext('a',3,line=-1.5,at=18.5)
mtext('Pr(occurrence) = 0.35',3,line=.5)
cols.legend = seq(6.5,10.5,1)
cols.legend = (cols.legend - min(cols.legend)) / (max(cols.legend) - min(cols.legend))
cols.legend = rgb(1-cols.legend,cols.legend,0,.5)

d = df[df$aegypti==unique(df$aegypti)[2],]
plot(
  unique(d$temp),d$AR.model[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[1]],
  type='l',las=1,col=cols[1],lwd=2,ylim=c(0,1),xaxt='n')
for(ii in 2:length(unique(d$gdp_pcppp2005))){
  lines(unique(d$temp),d$AR.model[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[ii]],col=cols[ii],lwd=2)
}
mtext('b',3,line=-1.5,at=18.5)
mtext('Pr(occurrence) = 0.65',3,line=.5)

d = df[df$aegypti==unique(df$aegypti)[3],]
plot(
  unique(d$temp),d$AR.model[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[1]],
  type='l',las=1,col=cols[1],lwd=2,ylim=c(0,1),xaxt='n')
for(ii in 2:length(unique(d$gdp_pcppp2005))){
  lines(unique(d$temp),d$AR.model[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[ii]],col=cols[ii],lwd=2)
}
mtext('c',3,line=-1.5,at=18.5)
mtext('Pr(occurrence) = 0.95',3,line=.5)


d = df[df$aegypti==unique(df$aegypti)[1],]
plot(
  unique(d$temp),d$AR.stat.full[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[1]],
  type='l',las=1,col=cols[1],lwd=2,ylim=c(0,1),xaxt='n')
for(ii in 2:length(unique(d$gdp_pcppp2005))){
  lines(unique(d$temp),d$AR.stat.full[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[ii]],col=cols[ii],lwd=2)
}
mtext('Attack rate',2,line=2.5)
mtext('d',3,line=-1.5,at=18.5)
cols.legend = seq(6.5,10.5,1)
cols.legend = (cols.legend - min(cols.legend)) / (max(cols.legend) - min(cols.legend))
cols.legend = rgb(1-cols.legend,cols.legend,0,.5)
legend(
  x = 30.7, y = 18.5, legend = seq(6.5,10.5,1), title = 'ln(GCP)',
  lwd = 2, lty = 1, bty = 'n', col = cols.legend
)

d = df[df$aegypti==unique(df$aegypti)[2],]
plot(
  unique(d$temp),d$AR.stat.full[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[1]],
  type='l',las=1,col=cols[1],lwd=2,ylim=c(0,1),xaxt='n')
for(ii in 2:length(unique(d$gdp_pcppp2005))){
  lines(unique(d$temp),d$AR.stat.full[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[ii]],col=cols[ii],lwd=2)
}
mtext('e',3,line=-1.5,at=18.5)

d = df[df$aegypti==unique(df$aegypti)[3],]
plot(
  unique(d$temp),d$AR.stat.full[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[1]],
  type='l',las=1,col=cols[1],lwd=2,ylim=c(0,1),xaxt='n')
for(ii in 2:length(unique(d$gdp_pcppp2005))){
  lines(unique(d$temp),d$AR.stat.full[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[ii]],col=cols[ii],lwd=2)
}
mtext('f',3,line=-1.5,at=18.5)


d = df[df$aegypti==unique(df$aegypti)[1],]
plot(
  unique(d$temp),d$AR.stat.step[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[1]],
  type='l',las=1,col=cols[1],lwd=2,ylim=c(0,1))
for(ii in 2:length(unique(d$gdp_pcppp2005))){
  lines(unique(d$temp),d$AR.stat.step[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[ii]],col=cols[ii],lwd=2)
}
mtext('g',3,line=-1.5,at=18.5)
cols.legend = seq(6.5,10.5,1)
cols.legend = (cols.legend - min(cols.legend)) / (max(cols.legend) - min(cols.legend))
cols.legend = rgb(1-cols.legend,cols.legend,0,.5)
legend(
  x = 30.7, y = 18.5, legend = seq(6.5,10.5,1), title = 'ln(GCP)',
  lwd = 2, lty = 1, bty = 'n', col = cols.legend
)

d = df[df$aegypti==unique(df$aegypti)[2],]
plot(
  unique(d$temp),d$AR.stat.step[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[1]],
  type='l',las=1,col=cols[1],lwd=2,ylim=c(0,1))
for(ii in 2:length(unique(d$gdp_pcppp2005))){
  lines(unique(d$temp),d$AR.stat.step[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[ii]],col=cols[ii],lwd=2)
}
mtext(expression('Temperature (' * degree * 'C)'),1,line=2.75)
mtext('h',3,line=-1.5,at=18.5)

d = df[df$aegypti==unique(df$aegypti)[3],]
plot(
  unique(d$temp),d$AR.stat.step[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[1]],
  type='l',las=1,col=cols[1],lwd=2,ylim=c(0,1))
for(ii in 2:length(unique(d$gdp_pcppp2005))){
  lines(unique(d$temp),d$AR.stat.step[d$gdp_pcppp2005==unique(d$gdp_pcppp2005)[ii]],col=cols[ii],lwd=2)
}
mtext('i',3,line=-1.5,at=18.5)

dev.off()
