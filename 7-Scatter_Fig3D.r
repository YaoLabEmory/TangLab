# Significantly upregulated ISGs or Non-ISGs in Bru-seq. Upregulated genes that are also detected by DRIP-seq are highlighted. 

.libPaths("/home/yli99/R/x86_64-pc-linux-gnu-library/4.0")

PRV<-read.table(file="PRV.DESeq.uniqgene.xls",header=F)[,c(13,14)]
MR<-read.table(file="MR.DESeq.uniqgene.xls",header=F)[,c(13,14)]

names(PRV)<-c("FDR_PRV","Gene")
names(MR)<-c("FDR_MR","Gene")

data<-merge(PRV,MR,by="Gene")

ISG<-read.table(file="../CommonUpISG",header=F)
NonISG<-read.table(file="../CommonUpNonISG",header=F)

names(ISG)<-"Gene"
names(NonISG)<-"Gene"

ISGFDR<-merge(data,ISG,by="Gene")
NonISGFDR<-merge(data,NonISG,by="Gene")

row.names(ISGFDR)<-ISGFDR$Gene
row.names(NonISGFDR)<-NonISGFDR$Gene

ISGFDR<-ISGFDR[,-1]
NonISGFDR<-NonISGFDR[,-1]

ISGFDR<-(-log(ISGFDR,10))
NonISGFDR<-(-log(NonISGFDR,10))

NonISGFDR$FDR_MR<-(-NonISGFDR$FDR_MR)

data<-rbind(ISGFDR,NonISGFDR)
data$Gene<-row.names(data)

FourUp<-read.table(file="/projects/compbio/users/yli/DNAseq/DRIP_TangLab/Mapping/mr_New/Scatter/FourUp",header=F)
names(FourUp)<-"Gene"
FourUp$color<-"red"

data<-merge(data,FourUp,by="Gene",all.x=T)
data <- replace(data, is.na(data), "grey")

data$size=2
data[data$color=="red",]$size=3

data[data$color=="red" & data$FDR_MR<0,]$color="blue"

library(ggplot2)
library(ggallin)
#p<-ggplot(data, aes(x=FDR_MR, y=FDR_PRV,color=color))+
p<-ggplot()+
#   geom_point(size=data$size)+
   theme_classic()+
#   scale_colour_manual(values = c("royalblue2","grey","darkorange"))+
   theme( axis.title.x = element_text(size = 30,color="black"),
            axis.text.y = element_text(size = 30,color="black"),
            axis.text.x = element_text(size = 30,color="black"),
            axis.title.y = element_text(size = 30,color="black"),
            legend.position = "none")+
#  geom_vline(xintercept=0,linetype="dashed")+
  xlab("mr_fdr")+
  ylab("prv_fdr")+
#  theme(legend.position="none")+
  ylim(0,250)+
  scale_x_continuous(trans = pseudolog10_trans,limits = c(-150, 150))
p<-p+
geom_point(data = data[data$color == "grey", ], aes(x = FDR_MR, y = FDR_PRV, color = color, size = size)) +
geom_point(data = data[data$color == "blue", ], aes(x = FDR_MR, y = FDR_PRV, color = color, size = size),shape = 15) +
geom_point(data = data[data$color == "red", ], aes(x = FDR_MR, y = FDR_PRV, color = color, size = size), shape = 15) +
scale_color_manual(values = c("grey" = "grey", "black" = "black", "red" = "darkorange", "blue" = "royalblue2")) +
  scale_size_identity()  +
geom_vline(xintercept=0,linetype="dashed")
pdf(file="ISGNonISGBru.pdf",6,6)
p
dev.off()
