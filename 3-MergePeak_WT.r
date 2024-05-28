#!/usr/bin/env Rscript
.libPaths("/home/yli99/R/x86_64-pc-linux-gnu-library/4.0")
library(Rsamtools)
library(DESeq2)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

# Step1: Process and obtain data
setwd("/projects/compbio/users/ywa3222/DNA_Tang_050924/Output/")
case1="U251_WT_ZV1/U251_WT_ZV1_peaks.narrowPeak"
case2="U251_WT_ZV2/U251_WT_ZV2_peaks.narrowPeak"
control1="U251_WT_Mock/U251_WT_Mock_peaks.narrowPeak"

case1.gr=import(case1)
case2.gr=import(case2)
control1.gr=import(control1)
case.peak.gr<-Reduce(union, list(
  case1.gr,
  case2.gr))

control.peak.gr=Reduce(union, list(
  control1.gr))

case.peak<-as.data.frame(case.peak.gr)
control.peak<-as.data.frame(control.peak.gr)
write.table(case.peak,file="WT/case.peak.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(control.peak,file="WT/control.peak.bed",row.names=F,col.names=F,sep="\t",quote=F)

all.peak.gr=Reduce(union, list(
  case1.gr,
  case2.gr,
  control1.gr))

all.peak<-as.data.frame(all.peak.gr)
write.table(all.peak,file="WT/all.peak.bed",row.names=F,col.names=F,sep="\t",quote=F)

case1="U251_WT_ZV1/U251_WT_ZV1.sorted.bed"
case2="U251_WT_ZV2/U251_WT_ZV2.sorted.bed"
control1="U251_WT_Mock/U251_WT_Mock.sorted.bed"
case1.gr=import(case1)
case2.gr=import(case2)
control1.gr=import(control1)

case1c=countOverlaps(all.peak.gr, case1.gr)
case2c=countOverlaps(all.peak.gr, case2.gr)
control1c=countOverlaps(all.peak.gr, control1.gr)

mat=data.frame(case1c, case2c, control1c)

sf=apply(mat,2,sum)
sf=sf/min(sf)
matnorm=sweep(mat,2,sf,FUN="/")

save(mat,file="WT/mat.rda")

caseread=apply(matnorm[,1:2],1,mean)
controlread=matnorm[,3]

MAT<-mat

mat<-mat[(mat$control1c>20) |
               (mat$case1c>20 & mat$case2c>20),]

condition=c(rep("treated",2),rep("control",1))
colData=data.frame(condition=condition)
dds=DESeqDataSetFromMatrix(countData =mat,colData = colData,design= ~condition)
dds=DESeq(dds)
result=results(dds)
result<-as.data.frame(result)

id=as.numeric(rownames(result))
resultmat=MAT[id,]
caseread=caseread[id]
controlread=controlread[id]
fc=(caseread+0.001)/(controlread+0.001)
all.peak.gr=as.data.frame(all.peak.gr)
all.peak.gr=all.peak.gr[id,1:3]
allmat=data.frame(all.peak.gr=all.peak.gr,controlread=controlread,caseread=caseread,fc=fc,result=result)

allmat<-allmat[order(allmat$result.pval),]
allmat$fdr.self<-p.adjust(allmat$result.pval,method="fdr")

write.table(allmat,file="WT/DESeq.xls",sep="\t",quote=FALSE,row.names=FALSE)
