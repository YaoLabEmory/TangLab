# Genome-wide significantly upregulated and downregulated R-loop regions detected by DRIP-seq in ZIKVMR-infected U-251 MG/IFNAR1KO cells (FDR < 0.05).

.libPaths("/home/yli99/R/x86_64-pc-linux-gnu-library/4.0")
load("../mat.rda")

sf=apply(mat,2,sum)
sf=sf/min(sf)
mat=sweep(mat,2,sf,FUN="/")

peak<-read.table(file="../all.peak.bed",header=F)[,c(1:3)]
mat<-cbind(peak,mat)
names(mat)<-c("chr","start","end","case1","case2","control1")
mat$loci<-paste(mat$chr,mat$start,mat$end,sep="_")

DESeq<-read.table(file="../DESeq.xls",header=T)
sig<-DESeq[DESeq$fdr.self<0.05,c(1:3)]
names(sig)<-c("chr","start","end")
sig$loci<-paste(sig$chr,sig$start,sig$end,sep="_")

data<-merge(mat,sig,by="loci")
data<-data[,c(7,5,6)]

library(pheatmap)

library(RColorBrewer)
color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdGy")))(100)

data$FC<-((data[,2]+data[,3])/2+0.1)/(data[,1]+0.1)
data<-data[order(data$FC),c(1:3)]
pdf(file="sig_heatmap_KO.pdf",4,6)
pheatmap(data,cluster_row=F,cluster_col=F,scale="row",show_rownames=F,col=color)


dev.off()
