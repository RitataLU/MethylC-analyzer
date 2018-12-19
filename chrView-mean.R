#chrView-mean plot

chrvwlist=read.table("chrView_mean_list.txt",sep='\t', stringsAsFactors=FALSE)

n_rain=nrow(chrvwlist)
cl <- rainbow(n_rain)
viewfile1=read.table(chrvwlist[1,2],sep='\t',header=T)
ym_cg=max(viewfile1$meanCG,na.rm=TRUE)*100
ym_chg=max(viewfile1$meanCHG,na.rm=TRUE)*100
ym_chh=max(viewfile1$meanCHH,na.rm=TRUE)*100

ys_cg=min(viewfile1$meanCG,na.rm=TRUE)*100
ys_chg=min(viewfile1$meanCHG,na.rm=TRUE)*100
ys_chh=min(viewfile1$meanCHH,na.rm=TRUE)*100

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
outfile=paste("./chrView_mean_CG.png",sep="")
bitmap(outfile,type="png16m",height=2,width=8,res=400)

plot(NA,ylab=expression(paste(Delta,' CG methylation levels (%)')),xlab="Chromosome",col='white',pch=19,type='p',cex.main = 1,frame = FALSE,ylim=c(ys_cg-buf_cg,ym_cg+buf_cg),xlim=c(0,nrow(viewfile1)),xaxt="n",cex=1.5,las=1)
#axis(2,at=seq(ys_cg,ym_cg))
for (i in 1:nrow(chrvwlist)){
	name=chrvwlist[i,1]
	file=chrvwlist[i,2]
	viewfile=read.table(file,sep='\t',header=T)
	
	chr=unique(viewfile$chromosome)
	
	lines(seq(-20,length(viewfile$Position)),vector("numeric", length(seq(-20,length(viewfile$Position)))),type='h',lwd=0.5,col='black') #draw zero line
	for (j in chr){
		sub=viewfile[viewfile$chromosome==j,]
		lines(sub$Position,sub$meanCG*100,lwd=0.5,col= cl[i])
		#last row
		abline(v=sub[nrow(sub),1]+0.5,col="grey",lwd=0.5)
		#label
		labsite=sub[as.integer(nrow(sub)/2),1]
		lab=sub$chromosome[1]
		sta_p=sub[1,1]
		end_p=sub[(nrow(sub)),1]

axis(1,at=labsite ,label=lab,tick=F)
axis(1,at=sta_p ,label=FALSE)
axis(1,at=end_p ,label=FALSE)

xx<-c(sta_p,end_p)
axis(1,at=xx,labels=F)
}
}
legend("topright", legend = t(chrvwlist[1]), col = t(cl), lwd = 1,cex = 0.5)
dev.off()

if(ym_chg<2 && ys_chg>-2){
  buf_chg=0
  ym_chg=2
  ys_chg=-2
}else if(ym_chg>=0 && ys_chg<0){
  buf_chg=ceiling((abs(ys_chg)+ym_chg)/10)
}else if((ym_chg>=0 && ys_chg>=0)){
  buf_chg=ceiling((ym_chg-ys_chg)/10)
}else if((ym_chg<0 && ys_chg>=0)){
  buf_chg=ceiling((abs(ym_chg)+ys_chg)/10)
}else if((ym_chg<0 && ys_chg<0)){
  buf_chg=ceiling((abs(ys_chg)-abs(ym_chg))/10)
}

#CHG
outfile=paste("./chrView_mean_CHG.png",sep="")
bitmap(outfile,type="png16m",height=2,width=8,res=400)

plot(NA,ylab=expression(paste(Delta,' CHG methylation levels (%)')),xlab="Chromosome",col='white',pch=19,type='p',cex.main = 1,frame = FALSE,ylim=c(ys_chg-buf_chg,ym_chg+buf_chg),xlim=c(0,nrow(viewfile1)),xaxt="n",cex=1.5,las=1)

for (i in 1:nrow(chrvwlist)){
	name=chrvwlist[i,1]
	file=chrvwlist[i,2]
	viewfile=read.table(file,sep='\t',header=T)
	
	chr=unique(viewfile$chromosome)
	lines(seq(-20,length(viewfile$Position)),vector("numeric", length(seq(-20,length(viewfile$Position)))),type='h',lwd=0.5,col='black') #draw zero line
	for (j in chr){
		sub=viewfile[viewfile$chromosome==j,]
		lines(sub$Position,sub$meanCHG*100,lwd=0.5,col= cl[i])
		
		#last row
		abline(v=sub[nrow(sub),1]+0.5,col="grey",lwd=0.5)
		#label
		labsite=sub[as.integer(nrow(sub)/2),1]
		lab=sub$chromosome[1]
		sta_p=sub[1,1]
		end_p=sub[(nrow(sub)),1]

axis(1,at=labsite ,label=lab,tick=F)
axis(1,at=sta_p ,label=FALSE)
axis(1,at=end_p ,label=FALSE)

xx<-c(sta_p,end_p)
axis(1,at=xx,labels=F)
}
}
legend("topright", legend = t(chrvwlist[1]), col = t(cl), lwd = 1,cex = 0.5)
dev.off()


if(ym_chh<2 && ys_chh>-2){
  buf_chh=0
  ym_chh=2
  ys_chh=-2
}else if(ym_chh>=0 && ys_chh<0){
  buf_chh=ceiling((abs(ys_chh)+ym_chh)/10)
}else if((ym_chh>=0 && ys_chh>=0)){
  buf_chh=ceiling((ym_chh-ys_chh)/10)
}else if((ym_chh<0 && ys_chh>=0)){
  buf_chh=ceiling((abs(ym_chh)+ys_chh)/10)
}else if((ym_chh<0 && ys_chh<0)){
  buf_chh=ceiling((abs(ys_chh)-abs(ym_chh))/10)
}

#CHH
outfile=paste("./chrView_mean_CHH.png",sep="")
bitmap(outfile,type="png16m",height=2,width=8,res=400)

plot(NA,ylab=expression(paste(Delta,' CHH methylation levels (%)')),xlab="Chromosome",col='white',pch=19,type='p',cex.main = 1,frame = FALSE,ylim=c(ys_chh-buf_chh,ym_chh+buf_chh),xlim=c(0,nrow(viewfile1)),xaxt="n",cex=1.5,las=1)

for (i in 1:nrow(chrvwlist)){
	name=chrvwlist[i,1]
	file=chrvwlist[i,2]
	viewfile=read.table(file,sep='\t',header=T)
	
	chr=unique(viewfile$chromosome)
	
	lines(seq(-20,length(viewfile$Position)),vector("numeric", length(seq(-20,length(viewfile$Position)))),type='h',lwd=0.5,col='black') #draw zero line
	for (j in chr){
		sub=viewfile[viewfile$chromosome==j,]
		lines(sub$Position,sub$meanCHH*100,lwd=0.5,col= cl[i])
		#last row
		abline(v=sub[nrow(sub),1]+0.5,col="grey",lwd=0.5)
		#label
		labsite=sub[as.integer(nrow(sub)/2),1]
		lab=sub$chromosome[1]
		sta_p=sub[1,1]
		end_p=sub[(nrow(sub)),1]

axis(1,at=labsite ,label=lab,tick=F)
axis(1,at=sta_p ,label=FALSE)
axis(1,at=end_p ,label=FALSE)

xx<-c(sta_p,end_p)
axis(1,at=xx,labels=F)
}
}
legend("topright", legend = t(chrvwlist[1]), col = t(cl), lwd = 1,cex = 0.5)
dev.off()


