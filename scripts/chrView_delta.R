#chrView-mean plot

chrvwlist=read.table("chrView_delta_list.txt",sep='\t', stringsAsFactors=FALSE)


viewfile1=read.table('chrView_delta.txt',sep='\t',header=T)
ym_cg=max(viewfile1$deltaCG,na.rm=TRUE)*100
ym_chg=max(viewfile1$deltaCHG,na.rm=TRUE)*100
ym_chh=max(viewfile1$deltaCHH,na.rm=TRUE)*100

ys_cg=min(viewfile1$deltaCG,na.rm=TRUE)*100
ys_chg=min(viewfile1$deltaCHG,na.rm=TRUE)*100
ys_chh=min(viewfile1$deltaCHH,na.rm=TRUE)*100

buf_cg=0
buf_chg=0
buf_chh=0

if(ym_cg<2 && ys_cg>-2){
  buf_cg=0
  ym_cg=2
  ys_cg=-2
}else if(ym_cg>=0 && ys_cg<0){
  buf_cg=ceiling((abs(ys_cg)+ym_cg)/10)
}else if((ym_cg>=0 && ys_cg>=0)){
  buf_cg=ceiling((ym_cg-ys_cg)/10)
}else if((ym_cg<0 && ys_cg>=0)){
  buf_cg=ceiling((abs(ym_cg)+ys_cg)/10)
}else if((ym_cg<0 && ys_cg<0)){
  buf_cg=ceiling((abs(ys_cg)-abs(ym_cg))/10)
}

#CG
outfile=paste("./chrView_delta_CG.pdf",sep="")
pdf(outfile,height= 4,width=20)
par(mar=c(5.1,5.1,2.1,2.1))
par(las=2)
plot(NA,ylab=expression(paste(Delta,' CG methylation levels (%)')),xlab="",col='white',pch=19,type='p',cex.main = 1,cex.axis=1.5,cex.lab=1.5,cex.lab=1.5,frame = FALSE,ylim=c(ys_cg-buf_cg,ym_cg+buf_cg),xlim=c(0,nrow(viewfile1)),xaxt="n",cex=1.5,las=1)
#axis(2,at=seq(ys_cg,ym_cg))
#axis(2,at=seq(ys_cg,ym_cg))

name=chrvwlist[1,1]
#file=chrvwlist[i,2]
#viewfile=read.table(file,sep='\t',header=T)

chr=unique(viewfile1$chromosome)

lines(seq(-20,length(viewfile1$Position)),vector("numeric", length(seq(-20,length(viewfile1$Position)))),type='h',lwd=0.5,col='black') #draw zero line
for (j in chr){
  sub=viewfile1[viewfile1$chromosome==j,]
  sub$Position = row.names(sub)
  lines(sub$Position,sub$deltaCG*100,lwd=1.5,col= '#054C7F')
        #last row
  abline(v=as.numeric(sub[nrow(sub),1])+0.5,col="grey",lwd=0.5)
  #label
  labsite=sub[as.integer(nrow(sub)/2),1]
  lab=sub$chromosome[1]
  sta_p=sub[1,1]
  end_p=sub[(nrow(sub)),1]

axis(1,at=labsite ,label=lab,tick=F,cex.axis=1.5)
axis(1,at=sta_p ,label=FALSE)
axis(1,at=end_p ,label=FALSE)

xx<-c(sta_p,end_p)
axis(1,at=xx,labels=F)
}

legend("topright",legend = name, col =  '#054C7F', lwd = 1.5,cex = 1,seg.len= 1.1,box.lty=0)
#legend("topright", legend = t(chrvwlist[1]), col = t(cl), lwd = 1,cex = 1.1)
dev.off()

outfile=paste("./chrView_delta_CHG.pdf",sep="")
pdf(outfile,height= 4,width=20)
par(mar=c(5.1,5.1,2.1,2.1))
par(las=2)
plot(NA,ylab=expression(paste(Delta,' CHG methylation levels (%)')),xlab="",col='white',pch=19,type='p',cex.main = 1,cex.axis=1.5,cex.lab=1.5,cex.lab=1.5,frame = FALSE,ylim=c(ys_chg-buf_chg,ym_chg+buf_chg),xlim=c(0,nrow(viewfile1)),xaxt="n",cex=1.5,las=1)
#axis(2,at=seq(ys_chg,ym_chg))
#axis(2,at=seq(ys_chg,ym_chg))

name=chrvwlist[1,1]
#file=chrvwlist[i,2]
#viewfile=read.table(file,sep='\t',header=T)
chr=unique(viewfile1$chromosome)
lines(seq(-20,length(viewfile1$Position)),vector("numeric", length(seq(-20,length(viewfile1$Position)))),type='h',lwd=0.5,col='black') #draw zero line
for (j in chr){
  sub=viewfile1[viewfile1$chromosome==j,]
  sub$Position = row.names(sub)
  lines(sub$Position,sub$deltaCHG*100,lwd=1.5,col= '#054C7F')
        #last row
  abline(v=as.numeric(sub[nrow(sub),1])+0.5,col="grey",lwd=0.5)
  #label
  labsite=sub[as.integer(nrow(sub)/2),1]
  lab=sub$chromosome[1]
  sta_p=sub[1,1]
  end_p=sub[(nrow(sub)),1]

axis(1,at=labsite ,label=lab,tick=F,cex.axis=1.5)
axis(1,at=sta_p ,label=FALSE)
axis(1,at=end_p ,label=FALSE)

xx<-c(sta_p,end_p)
axis(1,at=xx,labels=F)
}

legend("topright",legend = name, col =  '#054C7F', lwd = 1.5,cex = 1,seg.len= 1.1,box.lty=0)
#legend("topright", legend = t(chrvwlist[1]), col = t(cl), lwd = 1,cex = 1.1)
dev.off()

outfile=paste("./chrView_delta_CHH.pdf",sep="")
pdf(outfile,height= 4,width=20)
par(mar=c(5.1,5.1,2.1,2.1))
par(las=2)
plot(NA,ylab=expression(paste(Delta,' CHH methylation levels (%)')),xlab="",col='white',pch=19,type='p',cex.main = 1,cex.axis=1.5,cex.lab=1.5,cex.lab=1.5,frame = FALSE,ylim=c(ys_chh-buf_chh,ym_chh+buf_chh),xlim=c(0,nrow(viewfile1)),xaxt="n",cex=1.5,las=1)
#axis(2,at=seq(ys_chh,ym_chh))
#axis(2,at=seq(ys_chh,ym_chh))

name=chrvwlist[1,1]
#file=chrvwlist[i,2]
#viewfile=read.table(file,sep='\t',header=T)

chr=unique(viewfile1$chromosome)

lines(seq(-20,length(viewfile1$Position)),vector("numeric", length(seq(-20,length(viewfile1$Position)))),type='h',lwd=0.5,col='black') #draw zero line
for (j in chr){
  sub=viewfile1[viewfile1$chromosome==j,]
  sub$Position = row.names(sub)
  lines(sub$Position,sub$deltaCHH*100,lwd=1.5,col= '#054C7F')
        #last row
  abline(v=as.numeric(sub[nrow(sub),1])+0.5,col="grey",lwd=0.5)
  #label
  labsite=sub[as.integer(nrow(sub)/2),1]
  lab=sub$chromosome[1]
  sta_p=sub[1,1]
  end_p=sub[(nrow(sub)),1]

axis(1,at=labsite ,label=lab,tick=F,cex.axis=1.5)
axis(1,at=sta_p ,label=FALSE)
axis(1,at=end_p ,label=FALSE)

xx<-c(sta_p,end_p)
axis(1,at=xx,labels=F)
}

legend("topright",legend = name, col =  '#054C7F', lwd = 1.5,cex = 1,seg.len= 1.1,box.lty=0)
#legend("topright", legend = t(chrvwlist[1]), col = t(cl), lwd = 1,cex = 1.1)
dev.off()
