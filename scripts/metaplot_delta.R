args = commandArgs(trailingOnly=TRUE)
samplelist=read.table("samples_list.txt",sep='\t')


metafile1 = read.table('metaplot_delta_CG.txt',sep='\t')
metafile2 = read.table('metaplot_delta_CHG.txt',sep='\t')
metafile3 = read.table('metaplot_delta_CHH.txt',sep='\t')

# y_cg=max(apply(as.matrix(metafile1[,7:86]),2,mean,na.rm=T))*100
# y_chg=max(apply(as.matrix(metafile2[,7:86]),2,mean,na.rm=T))*100
# y_chh=max(apply(as.matrix(metafile3[,7:86]),2,mean,na.rm=T))*100

#mean each coulmn
##CG
av=apply(as.matrix(metafile1[,7:86]),2,mean,na.rm=T)
outfile=paste("metaplot_delta_CG.pdf",sep="")
pdf(outfile,height=5,width=8)
par(mar=c(2.5,4.5,2.5,2.5))
hy_cg=ifelse(max(av)>=0, max(av)*100, max(av)*100)
ly_cg=ifelse(min(av)>=0, min(av)*100, min(av)*100)

if(hy_cg<2 && ly_cg>-2){
  buf_cg=0
  hy_cg=2
  ly_cg=-2
}else if(hy_cg>=0 && ly_cg<0){
  buf_cg=ceiling((abs(ly_cg)+hy_cg)/10)
}else if((hy_cg>=0 && ly_cg>=0)){
  buf_cg=ceiling((hy_cg-ly_cg)/10)
}else if((hy_cg<0 && ly_cg>=0)){
  buf_cg=ceiling((abs(hy_cg)+ly_cg)/10)
}else if((hy_cg<0 && ly_cg<0)){
  buf_cg=ceiling((abs(ly_cg)-abs(hy_cg))/10)
}

plot(NA,ylab=expression(paste(Delta,' CG methylation levels (%)')),xlab=NA,type='l',cex.main = 1,cex.axis=1.1,cex.lab=1.3,xlim=c(0,ncol(metafile1)),ylim=c(ly_cg-buf_cg,hy_cg+buf_cg),xaxt="n",cex=1.3,las=1)
lines(av*100,xlab=NA,type='l',cex.main = 1,xaxt="n",cex=1.5,col= '#054C7F',lwd = 2.5)
lines(seq(-20,length(av)),vector("numeric", length(seq(-20,length(av)))),type='h',lwd=0.5,col='black') #draw zero line

xx<-c(10,20,40,60,70)
#xx<-c(20,40,60)
ind<-c("-2kb"," ","Genebody"," ","2kb")

abline(v=20,lty=2,col='grey70')
abline(v=60,lty=2,col='grey70')
axis(1,at=xx,tick = FALSE, labels=ind)

legend("topright", legend = paste(args[1],args[2],sep=' - '), col = '#054C7F', lwd = 2,cex = 1,seg.len= 1.1,box.lty=0)
dev.off()

####CHG
av1=apply(as.matrix(metafile2[,7:86]),2,mean,na.rm=T)


outfile=paste("metaplot_delta_CHG.pdf",sep="")
pdf(outfile,height=5,width=8)
par(mar=c(2.5,4.5,2.5,2.5))
#plot(NA,ylab='CHG methylation levels(%)',xlab=NA,type='l',cex.main = 1,xlim=c(0,ncol(metafile2)),ylim=c(0,y_chg),xaxt="n",cex=1.5)
hy_chg=ifelse(max(av1)>=0, max(av1)*100, max(av1)*100)
ly_chg=ifelse(min(av1)>=0, min(av1)*100, min(av1)*100)


if(hy_chg<2 && ly_chg>-2){
  buf_chg=0
  hy_chg=2
  ly_chg=-2
}else if(hy_chg>=0 && ly_chg<0){
  buf_chg=ceiling((abs(ly_chg)+hy_chg)/10)
}else if((hy_chg>=0 && ly_chg>=0)){
  buf_chg=ceiling((hy_chg-ly_chg)/10)
}else if((hy_chg<0 && ly_chg>=0)){
  buf_chg=ceiling((abs(hy_chg)+ly_chg)/10)
}else if((hy_chg<0 && ly_chg<0)){
  buf_chg=ceiling((abs(ly_chg)-abs(hy_chg))/10)
}

plot(NA,ylab=expression(paste(Delta,' CHG methylation levels (%)')),xlab=NA,type='l',cex.main = 1,cex.axis=1.1,cex.lab=1.3,cex.lab=1.5,xlim=c(0,ncol(metafile2)),ylim=c(ly_chg-buf_chg,hy_chg+buf_chg),xaxt="n",cex=1.5,las=1)
lines(av1*100,xlab=NA,type='l',cex.main = 1,xaxt="n",cex=1.5,col= '#054C7F',lwd = 2.5)
lines(seq(-20,length(av1)),vector("numeric", length(seq(-20,length(av1)))),type='h',lwd=0.5,col='black') #draw zero line

xx<-c(10,20,40,60,70)
#xx<-c(20,40,60)
ind<-c("-2kb"," ","Genebody"," ","2kb")

abline(v=20,lty=2,col='grey70')
abline(v=60,lty=2,col='grey70')
axis(1,at=xx,tick = FALSE, labels=ind)

legend("topright", legend = paste(args[1],args[2],sep=' - '), col = '#054C7F', lwd = 2,cex = 1,seg.len= 1.1,box.lty=0)
dev.off()


##CHH
av2=apply(as.matrix(metafile3[,7:86]),2,mean,na.rm=T)

outfile=paste("metaplot_delta_CHH.pdf",sep="")
pdf(outfile,height=5,width=8)
par(mar=c(2.5,4.5,2.5,2.5))
hy_chh=ifelse(max(av2)>=0, max(av2)*100, max(av2)*100)
ly_chh=ifelse(min(av2)>=0, min(av2)*100, min(av2)*100)


if(hy_chh<2 && ly_chh>-2){
  buf_chh=0
  hy_chh=2
  ly_chh=-2
}else if(hy_chh>=0 && ly_chh<0){
  buf_chh=ceiling((abs(ly_chh)+hy_chh)/10)
}else if((hy_chh>=0 && ly_chh>=0)){
  buf_chh=ceiling((hy_chh-ly_chh)/10)
}else if((hy_chh<0 && ly_chh>=0)){
  buf_chh=ceiling((abs(hy_chh)+ly_chh)/10)
}else if((hy_chh<0 && ly_chh<0)){
  buf_chh=ceiling((abs(ly_chh)-abs(hy_chh))/10)
}

plot(NA,ylab=expression(paste(Delta,' CHH methylation levels (%)')),xlab=NA,type='l',cex.main = 1,cex.axis=1.1,cex.lab=1.3,cex.lab=1.5,xlim=c(0,ncol(metafile3)),ylim=c(ly_chh-buf_chh,hy_chh+buf_chh),xaxt="n",cex=1.5,las=1)
lines(av2*100,xlab=NA,type='l',cex.main = 1,xaxt="n",cex=1.5,col= '#054C7F',lwd = 2.5)
lines(seq(-20,length(av2)),vector("numeric", length(seq(-20,length(av2)))),type='h',lwd=0.5,col='black') #draw zero line

xx<-c(10,20,40,60,70)
#xx<-c(20,40,60)
ind<-c("-2kb"," ","Genebody"," ","2kb")

abline(v=20,lty=2,col='grey70')
abline(v=60,lty=2,col='grey70')
axis(1,at=xx,tick = FALSE, labels=ind)
legend("topright", legend = paste(args[1],args[2],sep=' - '), col = '#054C7F', lwd = 2,cex = 1,seg.len= 1.1,box.lty=0)
dev.off()
