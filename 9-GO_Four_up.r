.libPaths("/mnt/icebreaker/data2/home/yli/R/x86_64-redhat-linux-gnu-library/3.6")

library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)

GO=read.table("GO.xls",header=T,sep="\t")
GO<-GO[,c(1,3,5,8)]
GO<-GO %>%separate(GO.biological.process.complete, c("name", "ID"), "\\(")
colnames(GO) = c("GoTerm","GOID","genes","dir","FDR")
GO<-GO[GO$dir=="+",]
GO$log10FDR = -log10(as.numeric(as.character(GO$FDR)))

GO=head(GO,n=5)
GO<-GO[order(GO$log10FDR),]
GO$GoTerm <- factor(GO$GoTerm, levels = GO$GoTerm[order(GO$log10FDR)])
GO$enrichFactor=GO$genes

p = ggplot(GO,aes(log10FDR,GoTerm))+
           geom_point(aes(size=genes,color=log10FDR))+
           scale_colour_gradient(low="red",high="red")+
           labs(color="-log10FDR",size="gene number",x="-log10FDR",y="",title=""
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
pdf(file="GO.pdf",height=4,width=7)
p
dev.off()
