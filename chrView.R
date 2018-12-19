
#chrView plot

chrvwlist=read.table("chrView_list.txt",sep='\t', stringsAsFactors=FALSE)

n_rain=nrow(chrvwlist)
cl <- rainbow(n_rain)
viewfile1=read.table(chrvwlist[1,2],sep='\t',header=T)
y_cg=max(viewfile1$meanCG,na.rm=TRUE)*1.5*100
y_chg=max(viewfile1$meanCHG,na.rm=TRUE)*1.5*100
y_chh=max(viewfile1$meanCHH,na.rm=TRUE)*1.5*100

#CG
outfile=paste("./chrView_CG.png",sep="")
bitmap(outfile,type="png16m",height=2,width=8,res=400)

plot(NA,ylab='CG methylation levels (%)',xlab="Chromosome",col='white',pch=19,type='p',cex.main = 1,frame = FALSE,ylim=c(0,100),xlim=c(0,nrow(viewfile1)),xaxt="n",cex=1.5,las=1)

for (i in 1:nrow(chrvwlist)){
	name=chrvwlist[i,1]
	file=chrvwlist[i,2]
	viewfile=read.table(file,sep='\t',header=T)
	
	chr=unique(viewfile$chromosome)
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


#CHG
outfile=paste("./chrView_CHG.png",sep="")
bitmap(outfile,type="png16m",height=2,width=8,res=400)

if(y_chg<4){
  md_y_chg=4
}else{
  md_y_chg=y_chg
}
plot(NA,ylab='CHG methylation levels (%)',xlab="Chromosome",col='white',pch=19,type='p',cex.main = 1,frame = FALSE,ylim=c(0,md_y_chg),xlim=c(0,nrow(viewfile1)),xaxt="n",cex=1.5,las=1)

for (i in 1:nrow(chrvwlist)){
	name=chrvwlist[i,1]
	file=chrvwlist[i,2]
	viewfile=read.table(file,sep='\t',header=T)
	
	chr=unique(viewfile$chromosome)
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


#CHH
outfile=paste("./chrView_CHH.png",sep="")
bitmap(outfile,type="png16m",height=2,width=8,res=400)

if(y_chh<4){
  md_y_chh=4
}else{
  md_y_chh=y_chh
}
plot(NA,ylab='CHH methylation levels (%)',xlab="Chromosome",col='white',pch=19,type='p',cex.main = 1,frame = FALSE,ylim=c(0,md_y_chh),xlim=c(0,nrow(viewfile1)),xaxt="n",cex=1.5,las=1)

for (i in 1:nrow(chrvwlist)){
	name=chrvwlist[i,1]
	file=chrvwlist[i,2]
	viewfile=read.table(file,sep='\t',header=T)
	
	chr=unique(viewfile$chromosome)
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

