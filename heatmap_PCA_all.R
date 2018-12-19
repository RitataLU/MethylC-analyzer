library(graphics)
library(gplots)

args <- commandArgs(trailingOnly = TRUE)

#file
input1=args[1]

dat=read.table(input1 ,header = T)

a=strsplit(input1,'_')
b=strsplit(a[[1]][2],'.txt')

end=ncol(dat)
ma=apply(dat[4:end],1,max)
mi=apply(dat[4:end],1,min)
cutoff=args[2]
dat2=dat[(ma-mi) > cutoff,]
dim(dat2)
#15253     6


###
mycol <- colorpanel(n=40,low="white",high="red4")
outfile=paste("heatmap_",b,"_",cutoff,".png",sep = "",collapse='')
bitmap(outfile,type="png16m",height=5,width=3,res=400)
heatmap.2(as.matrix(dat2[4:end]),scale='none',dendrogram = c("both"),Colv=T,Rowv=T,key=T,trace="none",tracecol="azure",col=mycol ,labRow=NA,xlab=NA,margins=c(10,4),cex.lab=0.1,main= paste(b,cutoff,sep="_"))

dev.off()

#PCA
library("ggplot2")

pca <- prcomp(t(dat2[4:end]))
summary(pca)
p <- ggplot(data=as.data.frame(pca$x), aes(PC1,PC2, label="PCA_",b,"_",cutoff))
p +geom_hline(yintercept=0, colour="gray25")+geom_vline(xintercept=0, colour="gray65")+theme(legend.position = "none")+geom_point(size=3)+geom_text(aes(label=rownames(pca$x)),hjust=0.5, vjust=2, size=4)+xlab(paste0("PC1=",round(summary(pca)$importance[2,1],3)*100,"%"))+ylab(paste0("PC2=",round(summary(pca)$importance[2,2],3)*100,"%"))
ggsave(file=paste("PCA_",b,"_",cutoff,".jpeg",sep = "",collapse=''),dpi=300)
dev.off()


