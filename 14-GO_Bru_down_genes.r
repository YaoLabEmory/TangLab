.libPaths("/home/yli99/R/x86_64-pc-linux-gnu-library/4.0")

library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)

GO=read.table("GO.xls",header=T,sep="\t")
GO<-GO[,c(1,3,5,6,8)]
GO<-GO %>%separate(GO.biological.process.complete, c("name", "ID"), "\\(")
colnames(GO) = c("GoTerm","GOID","genes","dir","FoldEnrichment","FDR")
GO<-GO[GO$dir=="+",]
GO<-GO[GO$FDR<0.05,]
GO$FoldEnrichment<-as.numeric(GO$FoldEnrichment)
GO<-GO[order(GO$FoldEnrichment, decreasing = TRUE),]
#GO$log10FDR = -log10(as.numeric(as.character(GO$FDR)))

GO=head(GO,n=5)
GO$GoTerm <- factor(GO$GoTerm, levels = GO$GoTerm[order(GO$FoldEnrichment)])
GO$enrichFactor=GO$genes

p = ggplot(GO,aes(FoldEnrichment,GoTerm))+
           geom_point(aes(size=genes,color=FoldEnrichment))+
           scale_colour_gradient(low="blue",high="blue")+
           labs(color="FoldEnrichment",size="gene number",x="FoldEnrichment",y="",title=""
           )
p=p+expand_limits(x=c(0,5))
p=p+theme_bw()
p=p+theme(
        axis.title.x=element_text(size=15,color="black"),
        axis.title.y=element_text(size=15,color="black"),
                axis.text.x=element_text(size=12,color="black"),
        axis.text.y=element_text(size=12,color="black"),
        legend.text = element_text(size = rel(1.2)),
                legend.title = element_text(size = rel(1.2)),
                )
pdf(file="GO.pdf",height=4,width=9)
p
dev.off()
