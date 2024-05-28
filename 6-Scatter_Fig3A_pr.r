load("../mat.rda")

sf=apply(mat,2,sum)
sf=sf/min(sf)
mat=sweep(mat,2,sf,FUN="/")

peak<-read.table(file="../all.peak.bed",header=F)[,c(1:3)]
mat<-cbind(peak,mat)
names(mat)<-c("chr","start","end","case1","case2","control1","control2","control3")
mat$loci<-paste(mat$chr,mat$start,mat$end,sep="_")

DESeq<-read.table(file="../DESeq.xls",header=T)

sig<-DESeq[DESeq$fdr.self<2,c(1:3,13)]
names(sig)<-c("chr","start","end","fdr")
sig$loci<-paste(sig$chr,sig$start,sig$end,sep="_")

data<-merge(mat,sig,by="loci")
data<-data[,c(5:9,13)]

data$case<-(data[,1]+data[,2])/2+0.1
data$control<-(data[,3]+data[,4]+data[,5])/3+0.1

data$color<-"grey"
data[data$fdr<0.05,]$color<-"black"

library(ggplot2)
library(dplyr)

p<-ggplot(data, aes(x=log(case,10), y=log(control,10),color=color))+
   geom_point()+
   theme_classic()+
   scale_colour_manual(values = c("black","grey"))+
   theme( axis.title.x = element_text(size = 30,color="black"),
            axis.text.y = element_text(size = 30,color="black"),
            axis.text.x = element_text(size = 30,color="black"),
            axis.title.y = element_text(size = 30,color="black"))+
  geom_abline(slope=1,intercept=0)+
  xlab("control")+
  ylab("case")+
  theme(legend.position="none")

pdf(file="Scatter_plot_log10.pdf",6,6)
p
dev.off()
