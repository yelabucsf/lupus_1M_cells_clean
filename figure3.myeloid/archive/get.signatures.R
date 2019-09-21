library(data.table)

args=commandArgs(TRUE)
pbmc=fread(args[1])
out=args[2]
#samples=read.table('/ye/yelabstore2/10x.lupus/disease/v2/v2.interaction.colData.txt', row.names=NULL)
load('/ye/yelabstore2/10x.lupus/disease/demux.v2/norep.noadjrace.colData.rda')
samples.use=colData

#given a file with genes and clusters that they belong to, output a signature for each cell type for all clusters
genes=read.table('/ye/yelabstore2/10x.lupus/disease/correlate.props/gene.e.clusters.txt', header=T, sep='\t')
clusters=unique(genes$X8.Clusters)

#load gene expression
#samples.use=samples[!duplicated(samples$row.names), ]


all(colnames(pbmc)[-1] ==samples.use$row.names) #make sure everything is in the same order

for(c in clusters){
    print(c)
    ifn=genes[which(genes$X8.Clusters==c),]$id
#    print(ifn)
    pbmc.use=t(t(pbmc[,-1])/colSums(pbmc[, -1]))*median(colSums(pbmc[, -1]))
    pbmc.norm=log2(pbmc.use + 1)

    expr.ifn=pbmc.norm[match(ifn, pbmc$gene), ]
    expr.adj=t(apply(expr.ifn, 1, function(x){return(lm(x~as.factor(samples.use$batch))$residuals)}))
    pcs=prcomp(t(expr.adj))

    #print proportion of variance explained by first
    print('Proportion of Variance explained by first PC')
    eigs <- pcs$sdev^2
    print(eigs[1] / sum(eigs))

    pcs=pcs$x
    ifn.sig=data.frame(ifn=pcs[,1])
    rownames(ifn.sig)=rownames(pcs)
    ifn.binary=data.frame(ifn=ifn.sig$ifn > 0)
    ifn.binary$ifn=as.factor(ifn.binary$ifn)
    levels(ifn.binary$ifn)=c(0,1)
    rownames(ifn.binary)=rownames(pcs)

    write.table(ifn.sig, paste(out, 'cluster.', c, '.signature.txt', sep=''), quote=F, col.names=T, row.names=T)
#    write.table(ifn.binary, '/ye/yelabstore2/10x.lupus/disease/v2/v2.ifn.binary.txt', quote=F, col.names=T, row.names=T)

}



