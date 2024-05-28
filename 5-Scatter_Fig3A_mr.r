.libPaths("/home/yli99/R/x86_64-pc-linux-gnu-library/4.0")
data<-read.table(file="../DESeq.gene.color.xls",header=T)
names(data)<-c(names(data)[-c(14,15)],"gene","color")

library(ggplot2)
library(dplyr)

p <- ggplot(data, aes(x = log(controlread, 10), y = log(caseread, 10))) +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Control R-loop") +
  ylab("Case R-loop") +
  theme(
    axis.title.x = element_text(size = 30, color = "black"),
    axis.text.y = element_text(size = 30, color = "black"),
    axis.text.x = element_text(size = 30, color = "black"),
    axis.title.y = element_text(size = 30, color = "black"),
    legend.position = "none"
  )

p <- p + geom_point(data = subset(data, color == "grey"), aes(color = color),
                    color = "grey")

p <- p + geom_point(data = subset(data, color == "blue"), aes(color = color),
                    color = "blue")

p <- p + geom_point(data = subset(data, color == "red"), aes(color = color),
                    color = "red")


pdf(file="Scatter_plot_log10.pdf",6,6)
p
dev.off()
