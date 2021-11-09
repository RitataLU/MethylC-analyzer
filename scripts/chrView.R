
chrvwlist=read.table("chrView_list.txt", sep = '\t',stringsAsFactors=FALSE)

colors = c('#7360A5',  '#5DBD82')
#select colors

viewfile1=read.table(chrvwlist[1,2],sep='\t',header=T)

y_cg=max(viewfile1$meanCG,na.rm=TRUE)*1.2*100
y_chg=max(viewfile1$meanCHG,na.rm=TRUE)*1.2*100
y_chh=max(viewfile1$meanCHH,na.rm=TRUE)*1.2*100

#CG plot
outfile=paste("chrView_CG.pdf",sep="")
pdf(outfile,,height=4,width=20)
#bitmap(outfile,type="png16m",height=2.5,width=10,res=400)
par(mar=c(5.1,5.2,2.1,2.1))
par(las=2)
# par(xaxt="n")
plot(NA,ylab='CG methylation levels (%)',xlab="",col='white',pch=19,type='p',
    cex.main = 1,cex.axis=1.5,cex.lab=1.5,cex.sub=1.5,frame = FALSE,ylim=c(0,100),
    xlim=c(0,nrow(viewfile1)),xaxt="n",cex=1.5,las=1)

for (i in 1:nrow(chrvwlist)){
  name=chrvwlist[i,1]
  file=chrvwlist[i,2]
  group=chrvwlist[i,3]

  if (group == unique(chrvwlist[c('V3')])[[1]][1]) {
  clr = colors[1]
  } else {
  clr = colors[2]
  }

  viewfile=read.table(file, sep='\t',header=T)
  viewfile$Position = row.names(viewfile)
  
  chr=unique(viewfile$chromosome)
  for (j in chr){
    sub=viewfile[viewfile$chromosome==j,]
    #sub$Position = row.names(sub)
    
    lines(sub$Position,sub$meanCG*100,lwd=1.2,col= clr)
    #lines(rollapply(zoo(sub$meanCG, order.by =sub$Position ),FUN=mean, width=3,by = 1, align = "center"), col = clr,lwd = 2)


    #last row
    abline(v=as.numeric(sub[nrow(sub),1])+0.5,col="grey",lwd=0.5)
    #label
    labsite=sub[as.integer(nrow(sub)/2),1]
    lab=sub$chromosome[1]
    sta_p=sub[1,1]
    end_p=sub[(nrow(sub)),1]

    axis(1,at=labsite,labels = lab, tick=F,srt = 90, cex.axis = 1.5)
    #axis(1, at=seq(1, 10, by=1), labels = FALSE)
    axis(1,at=sta_p ,label=FALSE)
    axis(1,at=end_p ,label=FALSE)

    xx<-c(sta_p,end_p)
    axis(1,at=xx,labels=F)
  }
}
legend(x = "topright",legend = t(unique(chrvwlist[c('V3')])[[1]]), col = t(colors), lwd = 2,cex = 1.2, seg.len= 1.5, box.lty=0)
dev.off()

outfile=paste("chrView_CHG.pdf",sep="")
pdf(outfile,,height=4,width=20)
#bitmap(outfile,type="png16m",height=2.5,width=10,res=400)
par(mar=c(5.1,5.2,2.1,2.1))
par(las=2)
if(y_chg<4){
md_y_chg=4
}else{ 
md_y_chg=y_chg
}
plot(NA,ylab='CHG methylation levels (%)',xlab="",col='white',pch=19,type='p',
  cex.main = 1,cex.axis=1.5,cex.lab=1.5,cex.sub=1.5,frame = FALSE,ylim=c(0,md_y_chg+5),
  xlim=c(0,nrow(viewfile1)),xaxt="n",cex=1.5,las=1)

for (i in 1:nrow(chrvwlist)){
  name=chrvwlist[i,1]
  file=chrvwlist[i,2]
  group=chrvwlist[i,3]

  if (group == unique(chrvwlist[c('V3')])[[1]][1]) {
  clr = colors[1]
  } else {
  clr = colors[2]
  }

  viewfile=read.table( file,sep='\t',header=T)
  
  chr=unique(viewfile$chromosome)
  for (j in chr){
    sub=viewfile[viewfile$chromosome==j,]
    sub$Position = row.names(sub)
    lines(sub$Position,sub$meanCHG*100,lwd=1.2,col= clr)
    #last row
    abline(v=as.numeric(sub[nrow(sub),1])+0.5,col="grey",lwd=0.5)
    #label
    labsite=sub[as.integer(nrow(sub)/2),1]
    lab=sub$chromosome[1]
    sta_p=sub[1,1]
    end_p=sub[(nrow(sub)),1]

    axis(1,at=labsite,label=lab,tick=F,srt= 90,adj=0.6,cex.axis=1.5)
    axis(1,at=sta_p ,label=FALSE)
    axis(1,at=end_p ,label=FALSE)

    xx<-c(sta_p,end_p)
    axis(1,at=xx,labels=F)
  }
}
legend(x = "topright",legend = t(unique(chrvwlist[c('V3')])[[1]]), col = t(colors), lwd = 2,cex = 1.2, seg.len= 1.5, box.lty=0)
dev.off()


#CHH plot
outfile=paste("chrView_CHH.pdf",sep="")
pdf(outfile,,height=4,width=20)
#bitmap(outfile,type="png16m",height=2.5,width=10,res=400)
par(mar=c(5.1,5.2,2.1,2.1))
par(las=2)
if(y_chh<4){
md_y_chh=4
}else{ 
md_y_chh=y_chh
}
plot(NA,ylab='CHH methylation levels (%)',xlab="",col='white',pch=19,type='p',
  cex.main = 1,cex.axis=1.5,cex.lab=1.5,cex.sub=1.5,frame = FALSE,ylim=c(0,md_y_chh+5),
  xlim=c(0,nrow(viewfile1)),xaxt="n",cex=1.5,las=1)

for (i in 1:nrow(chrvwlist)){
  name=chrvwlist[i,1]
  file=chrvwlist[i,2]
  group=chrvwlist[i,3]

  if (group == unique(chrvwlist[c('V3')])[[1]][1]) {
  clr = colors[1]
  } else {
  clr = colors[2]
  }

  viewfile=read.table(file,sep='\t',header=T)
  
  chr=unique(viewfile$chromosome)
  for (j in chr){
    sub=viewfile[viewfile$chromosome==j,]
    sub$Position = row.names(sub)
    lines(sub$Position,sub$meanCHH*100,lwd=1.2,col= clr)
    #last row
    abline(v=as.numeric(sub[nrow(sub),1])+0.5,col="grey",lwd=0.5)
    #label
    labsite=sub[as.integer(nrow(sub)/2),1]
    lab=sub$chromosome[1]
    sta_p=sub[1,1]
    end_p=sub[(nrow(sub)),1]

    axis(1,at=labsite,label=lab,tick=F,srt= 90,adj=0.6,cex.axis=1.5)
    axis(1,at=sta_p ,label=FALSE)
    axis(1,at=end_p ,label=FALSE)

    xx<-c(sta_p,end_p)
    axis(1,at=xx,labels=F)
  }
}
legend(x = "topright",legend = t(unique(chrvwlist[c('V3')])[[1]]), col = t(colors), lwd = 2,cex = 1.2, seg.len= 1.5, box.lty=0)
dev.off()



