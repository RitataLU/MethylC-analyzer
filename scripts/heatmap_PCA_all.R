#install.packages("viridis")
library(viridis)
library(graphics)
library(gplots)
library(ggplot2)
library(ComplexHeatmap)

args <- commandArgs(trailingOnly = TRUE)

#file
input1= args[1]

dat=read.table(input1 ,header = T)

a=strsplit(input1,'_')
b=strsplit(a[[1]][2],'.txt')

end=ncol(dat)
ma=apply(dat[4:end],1,max)
mi=apply(dat[4:end],1,min)
cutoff=args[2]
dat2=dat[(ma-mi) > cutoff,]

dim(dat2)


###
h = Heatmap(as.matrix(dat2[4:end]), clustering_method_rows  = 'ward.D2',clustering_method_columns = 'ward.D2',
               column_title_gp = gpar(fontsize = 10), name = "methylation",row_dend_width = unit(2, "cm"), column_dend_height = unit(4, "cm"),
               col=viridis(5),na_col = "black",show_row_names = FALSE)
  
pdf(paste("Heatmap_",b,"_",cutoff,".pdf",sep = ""),width=6,height=8,useDingbats=FALSE)
draw(h)
dev.off()

##PCA
pca <- prcomp( t(dat2[4:end]))
summary(pca)

theme<-theme(strip.background=element_blank(),
             plot.title = element_text(hjust = 0.5,face = "bold", size =15),
             axis.text.x=element_text(face = "bold",colour="black",size=18),
             axis.text.y=element_text(face = "bold",colour="black",size=18),
             axis.title.x = element_text( colour="black",size=20, face="bold"),
             axis.title.y = element_text( colour="black",size=20, face="bold"),
             axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line") ,
             legend.key = element_rect(fill = "white", colour = "black"),
             legend.text = element_text(size = 15),
             legend.title = element_text(face = "bold"))

p <- ggplot(data=as.data.frame(pca$x),
                  aes(x=PC1, y=PC2))+
  geom_point(alpha = 0.8,size=2)+  theme_minimal() +
  xlab(paste0("PC1=",round(summary(pca)$importance[2,1],3)*100,"%"))+
  ylab(paste0("PC2=",round(summary(pca)$importance[2,2],3)*100,"%"))

p2 <- p+ geom_hline(yintercept=0, colour="gray55")+
geom_vline(xintercept=0, colour="gray55")+ theme+
geom_text(aes(label=rownames(pca$x)),vjust= 1.5,hjust="inward", size=4)


ggsave(plot=p2  ,height=8,width=10,dpi=300,
         filename=paste("PCA_",b,"_",cutoff,".pdf",sep=''), useDingbats=FALSE)

