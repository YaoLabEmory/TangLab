# DRIP-seq logFC for Bru-seq upregulated ISGs and Non-ISGs. **** indicates P ＜ 0.0001 (Unpaired t test).

.libPaths("/home/yli99/R/x86_64-pc-linux-gnu-library/4.0")
mr<-read.table(file="/projects/compbio/users/yli/DNAseq/DRIP_TangLab/Mapping/mr_New/DESeq.gene.color.xls",header=T)
names(mr)<-c(names(mr)[-c(14,15)],"gene","color")

prv<-read.table(file="/projects/compbio/users/yli/DNAseq/DRIP_TangLab/Mapping/prv_New/DESeq.gene.color.xls",header=T)
names(prv)<-c(names(prv)[-c(14,15)],"gene","color")

ISG<-read.table(file="../CommonUpISG",header=F)[,1]
NonISG<-read.table(file="../CommonUpNonISG",header=F)[,1]

mrISG<-mr[mr$gene %in% ISG,]
mrNonISG<-mr[mr$gene %in% NonISG,]
prvISG<-prv[prv$gene %in% ISG,]
prvNonISG<-prv[prv$gene %in% NonISG,]

mrISGlogFC<-log(mrISG$fc,2)
mrNonISGlogFC<-log(mrNonISG$fc,2)
prvISGlogFC<-log(prvISG$fc,2)
prvNonISGlogFC<-log(prvNonISG$fc,2)

write.table(mrISGlogFC,file="mrISGlogFC.txt",quote=F,sep="\t",row.name=F, col.name=F)
write.table(mrNonISGlogFC,file="mrNonISGlogFC.txt",quote=F,sep="\t",row.name=F, col.name=F)
write.table(prvISGlogFC,file="prvISGlogFC.txt",quote=F,sep="\t",row.name=F, col.name=F)
write.table(prvNonISGlogFC,file="prvNonISGlogFC.txt",quote=F,sep="\t",row.name=F, col.name=F)

t.test(mrISGlogFC,mrNonISGlogFC)
t.test(prvISGlogFC,prvNonISGlogFC)

data<-as.data.frame(stack(list(mrISG=mrISGlogFC, mrNonISG=mrNonISGlogFC, prvISG=prvISGlogFC, prvNonISG=prvNonISGlogFC)))
#data<-as.data.frame(stack(list(ISG=c(mrISGlogFC,prvISGlogFC), NonISG=c(mrNonISGlogFC,prvNonISGlogFC))))
names(data)<-c("logFC","type")

library(ggplot2)

p<-ggplot(data, aes(x=type, y=logFC, color=type,fill=type)) +
   geom_boxplot(outlier.shape = NA)+
   scale_colour_manual(values = c("darkorange","royalblue2","darkorange","royalblue2"))+
   scale_fill_manual(values = c("darkorange","royalblue2","darkorange","royalblue2"))
   #geom_jitter(color="black", size=0.4, alpha=0.9)

p<-p+theme()+
            xlab("")+
            ylab("logFC")+
            theme_classic()+
            theme( axis.title.x = element_text(size = 0),
            axis.text.y = element_text(size = 20,colour="black"),
            axis.text.x = element_text(size = 20,colour="black"),
            axis.title.y = element_text(size = 20,colour="black"))+
            theme(axis.text.x = element_text(angle = 90,hjust =1,vjust=0.5))+
            theme(legend.text=element_text(size=20,colour="black"))+
            theme(legend.title=element_text(size=0))+
            theme(legend.position='none')+
            theme(axis.line = element_line(colour = 'black', size = 1))

pdf(file="boxplot.pdf",4,5)
p
dev.off()
