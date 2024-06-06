# Heatmap demonstrate the normalized Bru-seq reads in 5 R-loop processing genes in ZIKVMR and ZIKVPR-infected cells.

.libPaths("/mnt/icebreaker/data2/home/yli/R/x86_64-redhat-linux-gnu-library/3.6")

load("ZMR_vs_M_New_mat.rda")
mr<-mat
colnames(mr)<-c("MR1","MR2","Mock1","Mock2")
load("ZPRV_vs_M_New_mat.rda")
prv<-mat
colnames(prv)<-c("PR1","PR2","Mock1","Mock2")
mat<-cbind(mr,prv)[,c(1,2,5,6,3,4)]

sf=apply(mat,2,sum)
sf=sf/min(sf)
mat=sweep(mat,2,sf,FUN="/")

peak<-read.table(file="hg19.gene.bed",header=F)[,4]
mat<-cbind(peak,mat)
names(mat)<-c("gene","MR1","MR2","PR1","PR2","Mock1","Mock2")
write.table(mat,file="gene.xlsx",row.names=F,col.names=T,quote=F)
#mat$loci<-paste(mat$chr,mat$start,mat$end,sep="_")

#DESeq<-read.table(file="../DESeq.xls",header=T)
#sig<-DESeq[DESeq$fdr.self<0.05,c(1:3)]
sig<-read.table(file="genelist1.txt",header=F)
names(sig)<-"gene"
#names(sig)<-c("chr","start","end")
#sig$loci<-paste(sig$chr,sig$start,sig$end,sep="_")
data<-merge(mat,sig,by="gene")
rownames(data)<-data$gene
data<-data[,c(2:7)]

library(pheatmap)

library(RColorBrewer)
color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)

data$FC<-((data[,1]+data[,2]+data[,3]+data[,4])/4+0.1)/((data[,5]+data[,6])/2+0.1)
data<-data[order(data$FC),1:6]
pdf(file="sig_heatmap1.pdf",4,3)
pheatmap(data,cluster_row=F,cluster_col=F,scale="row",show_rownames=T,col=color)


dev.off()
