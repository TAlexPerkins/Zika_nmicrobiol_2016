# load replicates of numbers infected under the statistical model
proj = as.matrix(read.csv('../../outputs/country_output_stat.csv')[,-1])
proj.inf = proj[,substr(colnames(proj),1,1)=='I']
proj.births = proj[,substr(colnames(proj),1,2)=='BI']
countries = as.character(read.csv('../../outputs/country_output_stat.csv')[,1])


# make Figure 3
pdf('../../outputs/posteriors_stat.pdf',width=6.5,3)

layout(t(matrix(1:10,5,2)))
par(oma=c(0,2,1.5,.5))
par(mar=c(4.5,0,0,.5))

# infections by country
# all
plot(
  -100,-100,type='l',xlim=range(hist(colSums(proj.inf[,]/1e6),breaks=50,plot=F)$mids),
  ylim=c(0,max(hist(colSums(proj.inf[,]/1e6),breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
polygon(
  c(hist(colSums(proj.inf[,]/1e6),breaks=50,plot=F)$mids,rev(hist(colSums(proj.inf[,]/1e6),breaks=50,plot=F)$mids)),
  c(hist(colSums(proj.inf[,]/1e6),breaks=50,plot=F)$density,rep(0,length(hist(colSums(proj.inf[,]/1e6),breaks=50,plot=F)$mids))),
  col=rgb(1,0,0,.5))
axis(1,at=seq(0,600,by=100),labels=c('','','','','','',''))
text(seq(0,600,by=100), par("usr")[3]-.007, labels = c('0','100','200','300','400','500','600'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
mtext('Normalized density',2,cex=.7,line=.4)
mtext('a',3,at=min(hist(colSums(proj.inf[,]/1e6),breaks=50,plot=F)$mids),cex=.7)
mtext('Total',3,cex=.7)

# Brazil
plot(
  -100,-100,type='l',xlim=range(hist(proj.inf[11,]/1e6,breaks=50,plot=F)$mids),
  ylim=c(0,max(hist(proj.inf[11,]/1e6,breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
polygon(
  c(hist(proj.inf[11,]/1e6,breaks=50,plot=F)$mids,rev(hist(proj.inf[11,]/1e6,breaks=50,plot=F)$mids)),
  c(hist(proj.inf[11,]/1e6,breaks=50,plot=F)$density,rep(0,length(hist(proj.inf[11,]/1e6,breaks=50,plot=F)$mids))),
  col=rgb(1,0,0,.5))
axis(1,at=seq(0,max(proj.inf[11,]),by=30e6)/1e6,labels=c('','','','','','',''))
text(seq(0,max(proj.inf[11,]),by=30e6)/1e6, par("usr")[3]-.015, labels = c('0','30','60','90','120','150','180'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
mtext('b',3,at=min(hist(proj.inf[11,]/1e6,breaks=50,plot=F)$mids),cex=.7)
mtext('Brazil',3,cex=.7)

# Venezuela
plot(
  -100,-100,type='l',xlim=range(hist(proj.inf[49,]/1e6,breaks=50,plot=F)$mids),
  ylim=c(0,max(hist(proj.inf[49,]/1e6,breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
polygon(
  c(hist(proj.inf[49,]/1e6,breaks=50,plot=F)$mids,rev(hist(proj.inf[49,]/1e6,breaks=50,plot=F)$mids)),
  c(hist(proj.inf[49,]/1e6,breaks=50,plot=F)$density,rep(0,length(hist(proj.inf[49,]/1e6,breaks=50,plot=F)$mids))),
  col=rgb(1,0,0,.5))
axis(1,at=seq(0,max(proj.inf[49,]),by=5e6)/1e6,labels=c('','','','','',''))
text(seq(0,max(proj.inf[49,]),by=5e6)/1e6, par("usr")[3]-.15, labels = c('0','5','10','15','20','25'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
mtext('Infections among all people (millions)',1,cex=.7,line=2)
mtext('c',3,at=min(hist(proj.inf[49,]/1e6,breaks=50,plot=F)$mids),cex=.7)
mtext('Venezuela',3,cex=.7)

# Colombia
plot(
  -100,-100,type='l',xlim=range(hist(proj.inf[16,]/1e6,breaks=50,plot=F)$mids),
  ylim=c(0,max(hist(proj.inf[16,]/1e6,breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
polygon(
  c(hist(proj.inf[16,]/1e6,breaks=50,plot=F)$mids,rev(hist(proj.inf[16,]/1e6,breaks=50,plot=F)$mids)),
  c(hist(proj.inf[16,]/1e6,breaks=50,plot=F)$density,rep(0,length(hist(proj.inf[16,]/1e6,breaks=50,plot=F)$mids))),
  col=rgb(1,0,0,.5))
axis(1,at=seq(0,max(proj.inf[16,]),by=10e6)/1e6,labels=c('','','','',''))
text(seq(0,max(proj.inf[16,]),by=10e6)/1e6, par("usr")[3]-.07, labels = c('0','10','20','30','40'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
mtext('d',3,at=min(hist(proj.inf[16,]/1e6,breaks=50,plot=F)$mids),cex=.7)
mtext('Colombia',3,cex=.7)

# Puerto Rico
plot(
  -100,-100,type='l',xlim=range(hist(proj.inf[39,]/1e6,breaks=50,plot=F)$mids),
  ylim=c(0,max(hist(proj.inf[39,]/1e6,breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
polygon(
  c(hist(proj.inf[39,]/1e6,breaks=50,plot=F)$mids,rev(hist(proj.inf[39,]/1e6,breaks=50,plot=F)$mids)),
  c(hist(proj.inf[39,]/1e6,breaks=50,plot=F)$density,rep(0,length(hist(proj.inf[39,]/1e6,breaks=50,plot=F)$mids))),
  col=rgb(1,0,0,.5))
axis(1,at=seq(0,max(proj.inf[39,]),by=.5e6)/1e6,labels=c('','','','','','',''))
text(seq(0,max(proj.inf[39,]),by=.5e6)/1e6, par("usr")[3]-1.5, labels = c('0','0.5','1.0','1.5','2.0','2.5','3.0'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
mtext('e',3,at=min(hist(proj.inf[39,]/1e6,breaks=50,plot=F)$mids),cex=.7)
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
axis(1,at=seq(0,max(colSums(proj.births[,])),by=2000000)/1e6,labels=c('','','','','',''))
text(seq(0,max(colSums(proj.births[,])),by=2000000)/1e6, par("usr")[3]-.33, labels = c('0','2','4','6','8','10'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
mtext('Normalized density',2,cex=.7,line=.4)
mtext('f',3,at=min(hist(colSums(proj.births[,]/1e6),breaks=50,plot=F)$mids),cex=.7)
mtext('Total',3,cex=.7)

# Brazil
plot(
  -100,-100,type='l',xlim=range(hist(proj.births[11,]/1e6,breaks=50,plot=F)$mids),
  ylim=c(0,max(hist(proj.births[11,]/1e6,breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
polygon(
  c(hist(proj.births[11,]/1e6,breaks=50,plot=F)$mids,rev(hist(proj.births[11,]/1e6,breaks=50,plot=F)$mids)),
  c(hist(proj.births[11,]/1e6,breaks=50,plot=F)$density,rep(0,length(hist(proj.births[11,]/1e6,breaks=50,plot=F)$mids))),
  col=rgb(0,0,0,.5))
axis(1,at=seq(0,max(proj.births[11,]),by=500000)/1e6,labels=c('','','','','',''))
text(seq(0,max(proj.births[11,]),by=500000)/1e6, par("usr")[3]-1.5, labels = c('0','0.5','1.0','1.5','2.0','2.5'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
mtext('g',3,at=min(hist(proj.births[11,]/1e6,breaks=50,plot=F)$mids),cex=.7)
mtext('Brazil',3,cex=.7)

# Venezuela
plot(
  -100,-100,type='l',xlim=range(hist(proj.births[49,]/1e6,breaks=50,plot=F)$mids),
  ylim=c(0,max(hist(proj.births[49,]/1e6,breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
polygon(
  c(hist(proj.births[49,]/1e6,breaks=50,plot=F)$mids,rev(hist(proj.births[49,]/1e6,breaks=50,plot=F)$mids)),
  c(hist(proj.births[49,]/1e6,breaks=50,plot=F)$density,rep(0,length(hist(proj.births[49,]/1e6,breaks=50,plot=F)$mids))),
  col=rgb(0,0,0,.5))
axis(1,at=seq(0,max(proj.births[49,]),by=100000)/1e6,labels=c('','','','','',''))
text(seq(0,max(proj.births[49,]),by=100000)/1e6, par("usr")[3]-7, labels = c('0','0.1','0.2','0.3','0.4','0.5'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
mtext('Infections among childbearing women (millions)',1,cex=.7,line=2)
mtext('h',3,at=min(hist(proj.births[49,]/1e6,breaks=50,plot=F)$mids),cex=.7)
mtext('Venezuela',3,cex=.7)

# Colombia
plot(
  -100,-100,type='l',xlim=range(hist(proj.births[16,]/1e6,breaks=50,plot=F)$mids),
  ylim=c(0,max(hist(proj.births[16,]/1e6,breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
polygon(
  c(hist(proj.births[16,]/1e6,breaks=50,plot=F)$mids,rev(hist(proj.births[16,]/1e6,breaks=50,plot=F)$mids)),
  c(hist(proj.births[16,]/1e6,breaks=50,plot=F)$density,rep(0,length(hist(proj.births[16,]/1e6,breaks=50,plot=F)$mids))),
  col=rgb(0,0,0,.5))
axis(1,at=seq(0,max(proj.births[16,]),by=200000)/1e6,labels=c('','','','',''))
text(seq(0,max(proj.births[16,]),by=200000)/1e6, par("usr")[3]-3.4, labels = c('0','0.2','0.4','0.6','0.8'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
mtext('i',3,at=min(hist(proj.births[16,]/1e6,breaks=50,plot=F)$mids),cex=.7)
mtext('Colombia',3,cex=.7)

# Puerto Rico
plot(
  -100,-100,type='l',xlim=range(hist(proj.births[39,]/1e6,breaks=50,plot=F)$mids),
  ylim=c(0,max(hist(proj.births[39,]/1e6,breaks=50,plot=F)$density)),xlab='',ylab='',xaxt='n',yaxt='n')
polygon(
  c(hist(proj.births[39,]/1e6,breaks=50,plot=F)$mids,rev(hist(proj.births[39,]/1e6,breaks=50,plot=F)$mids)),
  c(hist(proj.births[39,]/1e6,breaks=50,plot=F)$density,rep(0,length(hist(proj.births[39,]/1e6,breaks=50,plot=F)$mids))),
  col=rgb(0,0,0,.5))
axis(1,at=seq(0,max(proj.births[39,]),by=10000)/1e6,labels=c('','','','',''))
text(seq(0,max(proj.births[39,]),by=10000)/1e6, par("usr")[3]-80, labels = c('0','0.01','0.02','0.03','0.04'), srt = 45, pos = 1, xpd = TRUE,cex=.9)
mtext('j',3,at=min(hist(proj.births[39,]/1e6,breaks=50,plot=F)$mids),cex=.7)
mtext('Puerto Rico',3,cex=.7)

dev.off()
