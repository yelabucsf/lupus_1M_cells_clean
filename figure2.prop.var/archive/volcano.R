library(data.table)
library(ggplot2)

args=commandArgs(TRUE)
ct=args[1]

res=fread(paste('v5..', ct , '.deseq.txt', sep=''))

ifn=read.table('/ye/yelabstore/10x.lupus/batch4.analysis/ifn.lupus.crow.etal.txt', header=F)$V1


res$ifn=res$V1 %in% ifn
res$sig=res$qval < 0.05

print(head(res))

lim=max(abs(res$log2FoldChange)) + 0.2
min_lim=-lim

pdf(paste(ct, '.volcano.pdf', sep=''), useDingbats=F)
ggplot(res, aes(x=res$log2FoldChange, y=pvalue, color=ifn, shape=sig))+ theme_bw() + xlim(-5, 5) + geom_point()
dev.off()
