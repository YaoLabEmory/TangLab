# Genome-wide significantly upregulated and downregulated sites detected by DRIP-seq in IFN-infected cells (FDR < 0.05).
# DRIP-seq: Heatmap of comparison between IFN vs Mock. 

.libPaths("/mnt/icebreaker/data2/home/yli/R/x86_64-redhat-linux-gnu-library/3.6")
load("../mat.rda")

sf=apply(mat,2,sum)
sf=sf/min(sf)
mat=sweep(mat,2,sf,FUN="/")

peak<-read.table(file="../all.peak.bed",header=F)[,c(1:3)]
mat<-cbind(peak,mat)
names(mat)<-c("chr","start","end","case1","case2","case3","control1","control2","control3")
mat$loci<-paste(mat$chr,mat$start,mat$end,sep="_")

DESeq<-read.table(file="../DESeq.xls",header=T)
sig<-DESeq[DESeq$fdr.self<0.05,c(1:3)]
names(sig)<-c("chr","start","end")
sig$loci<-paste(sig$chr,sig$start,sig$end,sep="_")

data<-merge(mat,sig,by="loci")
data<-data[,c(5:10)]

library(pheatmap)

library(RColorBrewer)
color = colorRampPalette(rev(brewer.pal(n = 7, name ="PiYG")))(100)

data$FC<-((data[,1]+data[,2]+data[,3])/3+0.1)/((data[,6]+data[,4]+data[,5])/3+0.1)
data<-data[order(data$FC),c(4:6,1:3)]
pdf(file="sig_heatmap.pdf",4,6)
pheatmap(data,cluster_row=F,cluster_col=F,scale="row",show_rownames=F,col=color)


dev.off()
