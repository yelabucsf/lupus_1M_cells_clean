pca_cell_type_plot <- function(cell_type, covars, ifn, clinical) {
    healthy <- fread(paste("v2.healthy.",cell_type,".expr.matrix.eqtl.txt",sep=""))
    sle <- fread(paste("v2.",cell_type,".expr.matrix.eqtl.txt",sep=""))

    gene.names <- healthy$gene;

    all <- cbind(healthy[,-1], sle[,-1]);
    ##browser();
    covars <- covars[na.omit(match(colnames(all),covars$ind)),];
    ifn_ifn <- ifn$ifn[na.omit(match(colnames(all),ifn$V1))]

    all.log <- log2(all);
    ##sle.log <- log2(sle[,-1]);
    if(length(which(is.na(rowSums(all.log)))) > 0) {
        all.log <- all.log[-which(is.na(rowSums(all.log)))];
    }

    all.log.norm <- apply(all.log, 2, function(x) {(x-mean(x))})
    all.log.std <- t(apply(all.log.norm, 1, function(x) {(x-mean(x))/sd(x)}))
    all.log.std.res <- t(lm(t(all.log.std)~covars$well*covars$disease)$residual);

    all.prcomp <- prcomp(all.log.std.res);

    plot.df <- data.frame(all.prcomp$rotation[,c(1,2)], covars);

    ##   browser();

    ggplot <- ggplot(aes(PC1,PC2,color=disease),data=plot.df)+geom_point()
    ggsave(ggplot, file=paste("prcomp.",cell_type,".disease.pdf",sep=""))

    ggplot <- ggplot(aes(PC1,PC2,color=well),data=plot.df)+geom_point()
    ggsave(ggplot, file=paste("prcomp.",cell_type,".well.pdf",sep=""))

    ggplot <- ggplot(aes(PC1,PC2,color=pop),data=plot.df)+geom_point()
    ggsave(ggplot, file=paste("prcomp.",cell_type,".pop.pdf",sep=""))

    ggplot <- ggplot(aes(PC1,PC2,color=batch),data=plot.df)+geom_point()
    ggsave(ggplot, file=paste("prcomp.",cell_type,".batch.pdf",sep=""))


    plot.df2 <- data.frame(all.prcomp$rotation[,c(1,2)], ifn=ifn_ifn);

    ggplot <- ggplot(aes(PC1,PC2,color=ifn),data=plot.df2)+geom_point()
    ggsave(ggplot, file=paste("prcomp.",cell_type,".ifn.pdf",sep=""))

    ##all.prcomp <- all.prcomp[!is.na(match(colnames(all),clinical$genotypeid)),]
    ##clinical <- clinical[na.omit(match(colnames(all),clinical$genotypeid))]

    plot.df3 <- data.frame(all.prcomp$rotation[!is.na(match(colnames(all),clinical$genotypeid)),c(1,2)], 
    clinical[na.omit(match(colnames(all),clinical$genotypeid))],
    ifn=ifn$ifn[na.omit(match(clinical$genotypeid,ifn$V1))]);

    ggplot <- ggplot(aes(PC1,PC2,color=sledaiscore),data=plot.df3)+geom_point()
    ggsave(ggplot, file=paste("prcomp.",cell_type,".sledaiscore.pdf",sep=""))

    ggplot <- ggplot(aes(PC1,PC2,color=sleactivity),data=plot.df3)+geom_point()
    ggsave(ggplot, file=paste("prcomp.",cell_type,".sleactivity.pdf",sep=""))

    return(list(plot.df, plot.df2, plot.df3))
}

