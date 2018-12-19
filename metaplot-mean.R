args = commandArgs(trailingOnly=TRUE)
samplelist=read.table("samples_list.txt",sep='\t')

#n_rain=nrow(samplelist)
#cl <- rainbow(n_rain)
metafile1=read.table(paste(samplelist[1,1],"_CG.matrix.gz",sep=''),skip=1)
metafile2=read.table(paste(samplelist[1,1],"_CHG.matrix.gz",sep=''),skip=1)
metafile3=read.table(paste(samplelist[1,1],"_CHH.matrix.gz",sep=''),skip=1)
y_cg=max(apply(as.matrix(metafile1[,7:86]),2,mean,na.rm=T))*100
y_chg=max(apply(as.matrix(metafile2[,7:86]),2,mean,na.rm=T))*100
y_chh=max(apply(as.matrix(metafile3[,7:86]),2,mean,na.rm=T))*100

#---

outfile=paste("./metaplot_mean_CG.png",sep="")
bitmap(outfile,type="png16m",height=3.5,width=3.5,res=450)
#plot(NA,ylab='CG methylation levels(%)',xlab=NA,type='l',cex.main = 1,xlim=c(0,ncol(metafile1)),ylim=c(0,y_cg),xaxt="n",cex=1.5)

#palette()[1] "black"   "red"     "green3"  "blue"    "cyan"    "magenta" "yellow" "gray"
si=0
for (i in 1:nrow(samplelist)){
  name=samplelist[i,3]
  if(i>1)
    if(name!=samplelist[i-1,3])
      si=i
}
av=0
av1=0
av2=0
for (i in 1:nrow(samplelist)){
  name=samplelist[i,1]
  file=paste(name,"_CG.matrix.gz",sep='')
  metafile=read.table(file,sep='\t',skip=1)
  if(i<si){
    av1=av1+apply(as.matrix(metafile[,7:86]),2,mean,na.rm=T)
  }else{
    av2=av2+apply(as.matrix(metafile[,7:86]),2,mean,na.rm=T)
  }
}
av1=av1/(si-1)
av2=av2/((length(samplelist)+1)-(si-1))
if(args[2]==samplelist[si-1,3] && args[2]!=args[3]){
  av=av1-av2
}else if (args[2]==samplelist[si,3] && args[2]!=args[3]) {
  av=av2-av1
}
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

plot(NA,ylab=expression(paste(Delta,' CG methylation levels (%)')),xlab=NA,type='l',cex.main = 1,xlim=c(0,ncol(metafile1)),ylim=c(ly_cg-buf_cg,hy_cg+buf_cg),xaxt="n",cex=1.5,las=1)
lines(av*100,xlab=NA,type='l',cex.main = 1,xaxt="n",cex=1.5,col=palette()[i])
lines(seq(-20,length(av)),vector("numeric", length(seq(-20,length(av)))),type='h',lwd=0.5,col='black') #draw zero line

xx<-c(20,40,60)
ind<-c("TSS",args[1],"TES")

abline(v=20,lty=2,col='grey70')
abline(v=60,lty=2,col='grey70')
axis(1,at=xx,labels=ind)

legend("topright", legend = paste(args[2],args[3],sep=' - '), col = t(palette()[i]), lwd = 1,cex = 0.5)
dev.off()

#---

#---chg
outfile=paste("./metaplot_mean_CHG.png",sep="")
bitmap(outfile,type="png16m",height=3.5,width=3.5,res=450)
#plot(NA,ylab='CHG methylation levels(%)',xlab=NA,type='l',cex.main = 1,xlim=c(0,ncol(metafile2)),ylim=c(0,y_chg),xaxt="n",cex=1.5)

#palette()[1] "black"   "red"     "green3"  "blue"    "cyan"    "magenta" "yellow" "gray"
si=0
for (i in 1:nrow(samplelist)){
  name=samplelist[i,3]
  if(i>1)
    if(name!=samplelist[i-1,3])
      si=i
}
av=0
av1=0
av2=0
for (i in 1:nrow(samplelist)){
  name=samplelist[i,1]
  file=paste(name,"_CHG.matrix.gz",sep='')
  metafile=read.table(file,sep='\t',skip=1)
  if(i<si){
    av1=av1+apply(as.matrix(metafile[,7:86]),2,mean,na.rm=T)
  }else{
    av2=av2+apply(as.matrix(metafile[,7:86]),2,mean,na.rm=T)
  }
}
av1=av1/(si-1)
av2=av2/((length(samplelist)+1)-(si-1))
if(args[2]==samplelist[si-1,3] && args[2]!=args[3]){
  av=av1-av2
}else if (args[2]==samplelist[si,3] && args[2]!=args[3]) {
  av=av2-av1
}
hy_chg=ifelse(max(av)>=0, max(av)*100, max(av)*100)
ly_chg=ifelse(min(av)>=0, min(av)*100, min(av)*100)

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

plot(NA,ylab=expression(paste(Delta,' CHG methylation levels (%)')),xlab=NA,type='l',cex.main = 1,xlim=c(0,ncol(metafile1)),ylim=c(ly_chg-buf_chg,hy_chg+buf_chg),xaxt="n",cex=1.5,las=1)
lines(av*100,xlab=NA,type='l',cex.main = 1,xaxt="n",cex=1.5,col=palette()[i])
lines(seq(-20,length(av)),vector("numeric", length(seq(-20,length(av)))),type='h',lwd=0.5,col='black') #draw zero line

xx<-c(20,40,60)
ind<-c("TSS",args[1],"TES")

abline(v=20,lty=2,col='grey70')
abline(v=60,lty=2,col='grey70')
axis(1,at=xx,labels=ind)

legend("topright", legend = paste(args[2],args[3],sep=' - '), col = t(palette()[i]), lwd = 1,cex = 0.5)
dev.off()

#---

#CHH
outfile=paste("./metaplot_mean_CHH.png",sep="")
bitmap(outfile,type="png16m",height=3.5,width=3.5,res=450)
#plot(NA,ylab='CHH methylation levels(%)',xlab=NA,type='l',cex.main = 1,xlim=c(0,ncol(metafile3)),ylim=c(0,y_chh),xaxt="n",cex=1.5)

#palette()[1] "black"   "red"     "green3"  "blue"    "cyan"    "magenta" "yellow" "gray"
#palette()[1] "black"   "red"     "green3"  "blue"    "cyan"    "magenta" "yellow" "gray"
si=0
for (i in 1:nrow(samplelist)){
  name=samplelist[i,3]
  if(i>1)
    if(name!=samplelist[i-1,3])
      si=i
}
av=0
av1=0
av2=0
for (i in 1:nrow(samplelist)){
  name=samplelist[i,1]
  file=paste(name,"_CHH.matrix.gz",sep='')
  metafile=read.table(file,sep='\t',skip=1)
  if(i<si){
    av1=av1+apply(as.matrix(metafile[,7:86]),2,mean,na.rm=T)
  }else{
    av2=av2+apply(as.matrix(metafile[,7:86]),2,mean,na.rm=T)
  }
}
av1=av1/(si-1)
av2=av2/((length(samplelist)+1)-(si-1))
if(args[2]==samplelist[si-1,3] && args[2]!=args[3]){
  av=av1-av2
}else if (args[2]==samplelist[si,3] && args[2]!=args[3]) {
  av=av2-av1
}
hy_chh=ifelse(max(av)>=0, max(av)*100, max(av)*100)
ly_chh=ifelse(min(av)>=0, min(av)*100, min(av)*100)

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

plot(NA,ylab=expression(paste(Delta,' CHH methylation levels (%)')),xlab=NA,type='l',cex.main = 1,xlim=c(0,ncol(metafile1)),ylim=c(ly_chh-buf_chh,hy_chh+buf_chh),xaxt="n",cex=1.5,las=1)
lines(av*100,xlab=NA,type='l',cex.main = 1,xaxt="n",cex=1.5,col=palette()[i])
lines(seq(-20,length(av)),vector("numeric", length(seq(-20,length(av)))),type='h',lwd=0.5,col='black') #draw zero line

xx<-c(20,40,60)
ind<-c("TSS",args[1],"TES")

abline(v=20,lty=2,col='grey70')
abline(v=60,lty=2,col='grey70')
axis(1,at=xx,labels=ind)

legend("topright", legend = paste(args[2],args[3],sep=' - '), col = t(palette()[i]), lwd = 1,cex = 0.5)
dev.off()
