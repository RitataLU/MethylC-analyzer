#metaplot
args = commandArgs(trailingOnly=TRUE)

samplelist=read.table("samples_list.txt",sep='\t')
colors = c('#7360A5',  '#5DBD82')


metafile1=read.table(paste(samplelist[1,1],"_CG.matrix.gz",sep=''),skip=1)
metafile2=read.table(paste(samplelist[1,1],"_CHG.matrix.gz",sep=''),skip=1)
metafile3=read.table(paste(samplelist[1,1],"_CHH.matrix.gz",sep=''),skip=1)
y_cg=max(apply(as.matrix(metafile1[,7:86]),2,mean,na.rm=T),na.rm=T)*1.2*100
y_chg=max(apply(as.matrix(metafile2[,7:86]),2,mean,na.rm=T),na.rm=T)*1.2*100
y_chh=max(apply(as.matrix(metafile3[,7:86]),2,mean,na.rm=T),na.rm=T)*1.2*100

#new add (for ylim scale)
for (i in 1:nrow(samplelist)){
    if(max(apply(as.matrix(read.table(paste(samplelist[i,1],"_CG.matrix.gz",sep=''),skip=1)[,7:86]),2,mean,na.rm=T),na.rm=T)*1.2*100 > y_cg)
    {
        y_cg = max(apply(as.matrix(read.table(paste(samplelist[i,1],"_CG.matrix.gz",sep=''),skip=1)[,7:86]),2,mean,na.rm=T),na.rm=T)*1.2*100
       
    }

    if(max(apply(as.matrix(read.table(paste(samplelist[i,1],"_CHG.matrix.gz",sep=''),skip=1)[,7:86]),2,mean,na.rm=T),na.rm=T)*1.2*100 > y_chg)
    {
        y_chg = max(apply(as.matrix(read.table(paste(samplelist[i,1],"_CHG.matrix.gz",sep=''),skip=1)[,7:86]),2,mean,na.rm=T),na.rm=T)*1.2*100
    }

    if(max(apply(as.matrix(read.table(paste(samplelist[i,1],"_CHH.matrix.gz",sep=''),skip=1)[,7:86]),2,mean,na.rm=T),na.rm=T)*1.2*100 > y_chh)
    {
        y_chh = max(apply(as.matrix(read.table(paste(samplelist[i,1],"_CHH.matrix.gz",sep=''),skip=1)[,7:86]),2,mean,na.rm=T),na.rm=T)*1.2*100
    }
}

if(y_cg<4){
  md_y_cg=4
}else{
  md_y_cg=y_cg
}

outfile=paste("./metaplot_CG.pdf",sep="")
pdf(outfile,,height=5,width=8)
#bitmap(outfile,type="png16m",height=3.5,width=3.8,res=450)
par(mar=c(2.5,4.5,2.5,2.5))
plot(NA,ylab='CG methylation levels (%)',xlab=NA,type='l',cex.main = 1,cex.axis=1.5,cex.lab=1.5,cex.lab=1.5,xlim=c(0,ncol(metafile1)),ylim=c(0,md_y_cg),xaxt="n",cex=1.5,las=1)

#palette()[1] "black"   "red"     "green3"  "blue"    "cyan"    "magenta" "yellow" "gray"
for (i in 1:nrow(samplelist)){
    name=samplelist[i,1]
    file=paste(name,"_CG.matrix.gz",sep='')
    group=samplelist[i,3]
    if (group == unique(samplelist[c('V3')])[[1]][1]) {
    clr = colors[1]
    } else {
    clr = colors[2]
    }
    metafile=read.table(file,sep='\t',skip=1)
    av1=apply(as.matrix(metafile[,7:86]),2,mean,na.rm=T)
    lines(av1*100,xlab=NA,type='l',cex.main = 1,xaxt="n",cex=1.5,col=clr)
}

xx<-c(10,20,40,60,70)
#xx<-c(20,40,60)
ind<-c("-2kb","TSS","Genebody","TES","2kb")

abline(v=20,lty=2,col='grey70')
abline(v=60,lty=2,col='grey70')
axis(1,at=xx,tick = FALSE, labels=ind)
legend('topright', legend = t(unique(samplelist[c('V3')])[[1]]), col = t(colors), lwd = 1.2,cex = 1,seg.len= 1.1,box.lty=0)
dev.off()

#---CHG
if(y_chg<4){
  md_y_chg=4
}else{
  md_y_chg=y_chg
}
outfile=paste("./metaplot_CHG.pdf",sep="")
pdf(outfile,,height=5,width=8)
#bitmap(outfile,type="png16m",height=3.5,width=3.8,res=450)
par(mar=c(2.5,4.5,2.5,2.5))
plot(NA,ylab='CHG methylation levels (%)',xlab=NA,type='l',cex.main = 1,cex.axis=1.5,cex.lab=1.5,cex.lab=1.5,xlim=c(0,ncol(metafile2)),ylim=c(0,md_y_chg),xaxt="n",cex=1.5,las=1)

#palette()[1] "black"   "red"     "green3"  "blue"    "cyan"    "magenta" "yellow" "gray"
for (i in 1:nrow(samplelist)){
    name=samplelist[i,1]
    file=paste(name,"_CHG.matrix.gz",sep='')
    group=samplelist[i,3]
    if (group == unique(samplelist[c('V3')])[[1]][1]) {
    clr = colors[1]
    } else {
    clr = colors[2]
    }

    metafile=read.table(file,sep='\t',skip=1)
    av1=apply(as.matrix(metafile[,7:86]),2,mean,na.rm=T)
    lines(av1*100,xlab=NA,type='l',cex.main = 1,xaxt="n",cex=1.5,col=clr)
}

xx<-c(10,20,40,60,70)
#xx<-c(20,40,60)
ind<-c("-2kb","TES","Genebody","TES","2kb")

abline(v=20,lty=2,col='grey70')
abline(v=60,lty=2,col='grey70')
axis(1,at=xx,tick = FALSE, labels=ind)
dev.off()

#---CHH
if(y_chh<4){
  md_y_chh=4
}else{
  md_y_chh=y_chh
}
outfile=paste("./metaplot_CHH.pdf",sep="")
pdf(outfile,,height=5,width=8)
#bitmap(outfile,type="png16m",height=3.5,width=3.8,res=450)
par(mar=c(2.5,4.5,2.5,2.5))
plot(NA,ylab='CHH methylation levels (%)',xlab=NA,type='l',cex.main = 1,cex.axis=1.5,cex.lab=1.5,cex.lab=1.5,xlim=c(0,ncol(metafile3)),ylim=c(0,md_y_chh),xaxt="n",cex=1.5,las=1)

#palette()[1] "black"   "red"     "green3"  "blue"    "cyan"    "magenta" "yellow" "gray"
for (i in 1:nrow(samplelist)){
    name=samplelist[i,1]
    file=paste(name,"_CHH.matrix.gz",sep='')
    group=samplelist[i,3]
    if (group == unique(samplelist[c('V3')])[[1]][1]) {
    clr = colors[1]
    } else {
    clr = colors[2]
    }

    metafile=read.table(file,sep='\t',skip=1)
    av1=apply(as.matrix(metafile[,7:86]),2,mean,na.rm=T)
    lines(av1*100,xlab=NA,type='l',cex.main = 1,xaxt="n",cex=1.5,col=clr)
}

xx<-c(10,20,40,60,70)
#xx<-c(20,40,60)
ind<-c("-2kb","TSS","Genebody","TES","2kb")

abline(v=20,lty=2,col='grey70')
abline(v=60,lty=2,col='grey70')
axis(1,at=xx,tick = FALSE, labels=ind)
legend('topright', legend = t(unique(samplelist[c('V3')])[[1]]), col = t(colors), lwd = 1.2,cex = 1,seg.len= 1.1,box.lty=0)
dev.off()

