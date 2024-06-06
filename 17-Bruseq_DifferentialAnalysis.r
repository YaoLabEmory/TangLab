#Bru-seq reads were counted in gene body and stored in a matrix called "mat.rda". Reads were normalized to total mapped reads and call diff by DESeq2

#!/apps/R-4.0.2/bin/R

library(Rsamtools)
library(DESeq2)
library(rtracklayer)

bins<-import("hg19.gene.bed")
all.peak.gr<-bins

load("mat.rda")
MAT<-mat
sf<-c(47409981,23932360,164079119,54083552)
sf=sf/min(sf)
matnorm=sweep(mat,2,sf,FUN="/")

caseread=apply(matnorm[,1:2],1,mean)
controlread=apply(matnorm[,3:4],1,mean)

mat<-mat[(mat$control1c>20 & mat$control2c>20) |
               (mat$case1c>20 & mat$case2c>20),]

condition=c(rep("treated",2),rep("control",2))
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
allmatFDR<-allmat[allmat$fdr.self<0.1,]
allmatFDR<-allmatFDR[!is.na(allmatFDR$fdr.self),]

write.table(allmat,file="DESeq.xls",sep="\t",quote=FALSE,row.names=FALSE)
write.table(allmatFDR,file="DESeq.FDR0.1.xls",sep="\t",quote=FALSE,row.names=FALSE)