#Count Bruseq reads in gene body for MR and Mock

#!/apps/R-4.0.2/bin/R

library(Rsamtools)
library(DESeq2)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

bins<-import("hg19.gene.bed")

case1="../sortedbed/ZMR1BRU.sorted.bed"
case2="../sortedbed/ZMR2BRU.sorted.bed"

control1="../sortedbed/M1BRU.sorted.bed"
control2="../sortedbed/M2BRU.sorted.bed"

case1.gr=import(case1)
case2.gr=import(case2)
control1.gr=import(control1)
control2.gr=import(control2)

case1c=countOverlaps(bins, case1.gr)
case2c=countOverlaps(bins, case2.gr)
control1c=countOverlaps(bins, control1.gr)
control2c=countOverlaps(bins, control2.gr)

mat=data.frame(case1c, case2c, control1c,control2c)

save(mat,file="mat.rda")
