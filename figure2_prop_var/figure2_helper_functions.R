library(nnls);
qualitative_prediction <- function(trait = trait, indices = NULL, method="svmLinear") {
    trait_index = match(trait, colnames(joined_pivot))
    
    if(is.null(indices)) {
        indices <- 1:nrow(joined_pivot);
    }

    performance <- NULL;
    
    ##print(indices);

    df <- data.frame(as.matrix(comp_lineage), trait=factor(unlist(joined_pivot[,..trait_index])))
    model <- train( trait ~ ., df[indices,], method=method, trControl = trainControl(method="cv", number =10))
    ##cat("comp_lineage: ", mean(model$results$Accuracy), "\n")
    performance <- rbind(performance, c(type="comp_lineage", performance=mean(model$results$Accuracy)));

    df <- data.frame(as.matrix(comp), trait=factor(unlist(joined_pivot[,..trait_index])))
    model <- train( trait ~ ., df[indices,], method=method, trControl = trainControl(method="cv", number =10))
    ##cat("comp: ", mean(model$results$Accuracy), "\n")
    performance <- rbind(performance, c(type="comp", performance=mean(model$results$Accuracy)));

    df <- data.frame(t(pbmc_expr_std[match(ifn_genes_use,expr_gene_names),]), trait=factor(unlist(joined_pivot[,..trait_index])))
    model <- train( trait ~ ., df[indices,], method=method, trControl = trainControl(method="cv", number =10))
    ##cat("pbmc IFN: ", mean(model$results$Accuracy), "\n")
    performance <- rbind(performance, c(type="PBMC IFN", performance=mean(model$results$Accuracy)));

    df <- data.frame(t(all_expr_cts_var), trait=factor(unlist(joined_pivot[,..trait_index])))
    model <- train( trait ~ ., df[indices,], method=method, trControl = trainControl(method="cv", number =10))
    ##cat("ct VAR: ", mean(model$results$Accuracy), "\n")
    performance <- rbind(performance, c(type="CT VAR", performance=mean(model$results$Accuracy)));

    df <- data.frame(t(all_includect_comp_expr_var), trait=factor(unlist(joined_pivot[,..trait_index])))
    model <- train( trait ~ ., df[indices,], method=method, trControl = trainControl(method="cv", number =10))
    ##cat("ct VAR + comp: ", mean(model$results$Accuracy), "\n")
    performance <- rbind(performance, c(type="CT VAR COMP", performance=mean(model$results$Accuracy)));
       
    ##t(performance)
    return(performance)
}

quantitative_prediction <- function(trait = trait, indices = NULL, method="svmLinear") {
    trait_index = match(trait, colnames(joined_pivot))
    
    if(is.null(indices)) {
        indices <- 1:nrow(joined_pivot);
    }

    performance <- NULL;
    
    df <- data.frame(as.matrix(comp_lineage), trait=unlist(joined_pivot[,..trait_index]))
    model <- train( trait ~ ., df[indices,], method=method, trControl = trainControl(method="cv", number =10))
    ##cat("comp_lineage: ", mean(model$results$Rsquared), "\n")
    performance <- rbind(performance, c(type="comp_lineage", performance=mean(model$results$Rsquared)));

    df <- data.frame(as.matrix(comp), trait=unlist(joined_pivot[,..trait_index]))
    model <- train( trait ~ ., df[indices,], method=method, trControl = trainControl(method="cv", number =10))
    ##cat("comp: ", mean(model$results$Rsquared), "\n")
    performance <- rbind(performance, c(type="comp", performance=mean(model$results$Rsquared)));

    df <- data.frame(t(pbmc_expr_std[match(ifn_genes_use,expr_gene_names),]), trait=unlist(joined_pivot[,..trait_index]))
    model <- train( trait ~ ., df[indices,], method=method, trControl = trainControl(method="cv", number =10))
    ##cat("pbmc IFN: ", mean(model$results$Rsquared), "\n")
    performance <- rbind(performance, c(type="PBMC IFN", performance=mean(model$results$Rsquared)));

    df <- data.frame(t(all_expr_cts_var), trait=unlist(joined_pivot[,..trait_index]))
    model <- train( trait ~ ., df[indices,], method=method, trControl = trainControl(method="cv", number =10))
    ##cat("ct VAR: ", mean(model$results$Rsquared), "\n")
    performance <- rbind(performance, c(type="CT VAR", performance=mean(model$results$Rsquared)));

    df <- data.frame(t(all_includect_comp_expr_var), trait=unlist(joined_pivot[,..trait_index]))
    ##print(dim(df));
    model <- train( trait ~ ., df[indices,], method=method, trControl = trainControl(method="cv", number =10))
    ##cat("ct VAR + comp: ", mean(model$results$Rsquared), "\n")
    performance <- rbind(performance, c(type="CT VAR COMP", performance=mean(model$results$Rsquared)));
    
    ##t(performance);
    return(performance)
}

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

standardize <- function(x) {
    a <- (x-mean(x))/sd(x);
    a[is.na(a)] <- 0;
    return(a)
}

test_single_gene <- function(expr, props, exprs, i) {
    df = NULL;
    cors = NULL;
    
    model_prop <- lm(expr[i,]~props-1)
    predicted_prop <- model_prop$fitted
    df <- data.frame(observed=expr[i,],predicted=predicted_prop,type="prop");
    cors <- cor(expr[i,], predicted_prop);
    print(summary(model_prop)$r.squared)

    model_expr <- lm(expr[i,]~exprs[i,,]-1)
    predicted_expr <- model_expr$fitted
    df <- rbind(df, data.frame(observed=expr[i,],predicted=predicted_expr,type="expr"))
    cors <- c(cors,cor(expr[i,], predicted_expr))
    print(summary(model_expr)$r.squared)
    
    predicted <- 10^(exprs[i,,1])*props[,1]+
                 10^(exprs[i,,2])*props[,2]+
                 10^(exprs[i,,3])*props[,3]+
                 10^(exprs[i,,4])*props[,4]+
                 10^(exprs[i,,5])*props[,5]+
                 10^(exprs[i,,6])*props[,6]+
                 10^(exprs[i,,7])*props[,7]
    
    predicted_log <- log10(predicted);
    
    df <- rbind(df, data.frame(observed=expr[i,],predicted=predicted_log,type="factor"))
    cors <- c(cors,cor(expr[i,], predicted_log))
    
    print(ggplot(aes(observed, predicted), data=df)+geom_point()+facet_wrap(~type, scales="fixed")+theme_bw()+geom_abline(slope=1))
    names(cors) <- c("prop","expr","factor")
    print(cors);
    ##return(test);
}

plot_single_gene_prop <- function(expr, props, i) {
    df = NULL;
    cors = NULL;
    
    df <- reshape(data.frame(expr=expr[i,], as.data.frame(props)), 
        varying=c("cM","ncM","cDC","pDC","Th","Tc","NK","B"),
        times=c("cM","ncM","cDC","pDC","Th","Tc","NK","B"),
        timevar="ct_name",
        v.names="ct_perc",
        direction='long')

    df$ct_name = factor(df$ct_name, levels=c("cM","ncM","cDC","pDC","Th","Tc","NK","B"))
    
    print(ggplot(aes(ct_perc, expr, color=pop), data=df)+geom_point()
    +facet_wrap(~ct_name, scales="free_x",nrow=1)
    +theme_bw()+geom_smooth(method = "lm", aes(fill=pop))+scale_colour_manual(values=c("#D82B29","#1B78B6","#F57F20"),aesthetics = c("colour", "fill")))
}


plot_single_gene_expr <- function(expr, exprs_list, i) {
      df = NULL;
    cors = NULL;
    
    
            exprs <- rbind(data.frame(exprs_list[[1]][i,,],pop="broad_white"),
                           data.frame(exprs_list[[2]][i,,],pop="sle_asian"),
                            data.frame(exprs_list[[3]][i,,],pop="sle_white"))
        
    df <- reshape(data.frame(expr=expr[i,], as.data.frame(exprs)), 
        varying=c("cM","ncM","cDC","pDC","Th","Tc","NK","B"),
        times=c("cM","ncM","cDC","pDC","Th","Tc","NK","B"),
        timevar="ct_name",
        v.names="ct_expr",
        direction='long')

    df$ct_name = factor(df$ct_name, levels=c("cM","ncM","cDC","pDC","Th","Tc","NK","B"))
    
    print(ggplot(aes(ct_expr, expr, color=pop), data=df)+geom_point()
    +facet_wrap(~ct_name, scales="free_x",nrow=1)
    +theme_bw()+geom_smooth(method = "lm", aes(fill=pop))+scale_colour_manual(values=c("#D82B29","#1B78B6","#F57F20"),aesthetics = c("colour", "fill")))
    ##return(test);
        
}

## account for covariance
##var(a+b) = var(a)+var(b)+2*cov(a,b)

model_prop_old <- function(expr, props) {
    output_prop <- NULL;
    for(i in 1:nrow(expr)) {
        model <- lm(expr[i,]~props-1)
        r2 = summary(model)$coef[,1]^2*apply(props, 2,var,na.rm=T)/var(expr[i,])*100;
        output_prop <- rbind(output_prop, 
                               c(r2,sum(r2),summary(model)$adj.r.squared*100))
    }
    colnames(output_prop) <- c(colnames(props), "total r2","adj r2");
    rownames(output_prop) <- dimnames(expr)[[1]]
    return(output_prop)
}

model_prop <- function(expr, props) {
    output_prop <- NULL;
    cov_mat <- cov(props);

    for(i in 1:nrow(expr)) {
        ##cat("Processing: ",i,"\n")
        model <- lm(expr[i,]~props-1)
        coef = summary(model)$coef[,1];

        vars = NULL;
        
        #j is each cell type
        for(j in 1:ncol(props)) {
            var_j = sum(cov_mat[j,]*coef[j]*coef);
            vars <- c(vars, var_j);
        }
        
        r2s = vars/sum(vars);
        
        ##print(r2);
        r2 = sum(vars)/var(expr[i,]);
        r2_adj = r2-(1-r2)*(ncol(props))/(ncol(expr)-ncol(props)-1);
        
        output_prop <- rbind(output_prop, 
                               c(r2s*100,
                                 r2_adj*100,
                                 summary(model)$adj.r.squared*100))
        
        ##cat(dim(output_prop),"\n")
    }
    ##print(c(colnames(props), "total r2","adj r2"))
    colnames(output_prop) <- c(colnames(props), "total r2","adj r2");
    rownames(output_prop) <- dimnames(expr)[[1]]
    return(output_prop)
}

model_expr_old <- function(expr, exprs) {
    output_expr <- NULL;
    for(i in 1:nrow(expr)) {
        model <- lm(expr[i,]~exprs[i,,]-1);
    
        r2 = summary(model)$coef[,1]^2*apply(exprs[i,,], 2,var,na.rm=T)/var(expr[i,])*100;
    
        output_expr <- rbind(output_expr, 
                       c(r2,sum(r2),summary(model)$adj.r.squared*100))
    }
    colnames(output_expr) <- c(dimnames(exprs)[[3]], "total r2","adj r2");
    rownames(output_expr) <- dimnames(exprs)[[1]]
    return(output_expr)
}


model_expr <- function(expr, exprs) {
    output_expr <- NULL;
    
    for(i in 1:nrow(expr)) {
        ##cat("Processing: ",i,"\n")
        cov_mat <- cov(exprs[i,,]);
        model <- lm(expr[i,]~exprs[i,,]-1);
        coef <- summary(model)$coef[,1];
        
        vars <- NULL;
        
        #j is each cell type
        for(j in 1:dim(exprs)[[3]]) {
            var_j = sum(cov_mat[j,]*coef[j]*coef);
            vars <- c(vars, var_j);
        }
        
        ##vars <- unlist(vars);
        ##print(vars);
        r2s = vars/sum(vars);
        
        ##print(r2);
        r2 = sum(vars)/var(expr[i,]);
        r2_adj = r2-(1-r2)*dim(exprs)[[3]]/(ncol(expr)-dim(exprs)[[3]]-1);
        
        output_expr <- rbind(output_expr, 
                               c(r2s*100,
                                 r2_adj*100,summary(model)$adj.r.squared*100))
        ##cat(dim(output_expr),"\n")
        
    }
    colnames(output_expr) <- c(dimnames(exprs)[[3]], "total r2","adj r2");
    rownames(output_expr) <- dimnames(expr)[[1]]
    return(output_expr)
}


## modeling proportion and expr together
model_prop_expr <- function(expr, props, exprs) {
    output_prop <- NULL;
    output_expr <- NULL;
    
    for(i in 1:nrow(expr)) {
        
        cov_mat <- cov(cbind(props, exprs[i,,]));
        model <- lm(expr[i,]~props+exprs[i,,]-1);
        coef <- summary(model)$coef[,1];
        
        vars <- NULL;
        
        #j is each cell type
        for(j in 1:(ncol(props)+dim(exprs)[[3]])) {
            var_j = sum(cov_mat[j,]*coef[j]*coef);
            vars <- c(vars, var_j);
        }
        
        ##vars <- unlist(vars);
        ##print(vars);
        r2s = vars/sum(vars);
        
        ##print(r2);
        
        prop_indices <- 1:ncol(props);
        expr_indices <- (ncol(props)+1):(ncol(props)+dim(exprs)[[3]]);
        
        r2_prop = sum(vars[prop_indices])/var(expr[i,]);
        r2_prop_adj = r2_prop;
        ##r2_prop_adj = 1-((1-r2_prop)*(ncol(expr)-1))/(ncol(expr)-ncol(props)-1);             
        ##r2_prop_adj = r2_prop-(1-r2_prop)*(dim(exprs)[[3]]+ncol(props))/(ncol(expr)-(dim(exprs)[[3]]+ncol(props))-1);             
        
        r2_expr = sum(vars[expr_indices])/var(expr[i,]);
        r2_expr_adj = r2_expr;
        ##r2_expr_adj = 1-((1-r2_expr)*(ncol(expr)-1))/(ncol(expr)-dim(exprs)[[3]]-1);             
        ##r2_expr_adj = r2_expr-(1-r2_expr)*(dim(exprs)[[3]]+ncol(props))/(ncol(expr)-(dim(exprs)[[3]]+ncol(props))-1);
        
        r2 = sum(vars)/var(expr[i,]);
        ##r2_adj = 1-((1-r2)*(ncol(expr)-1))/(ncol(expr)-dim(exprs)[[3]]+ncol(props)-1); 
        r2_adj = r2-(1-r2)*(dim(exprs)[[3]]+ncol(props))/(ncol(expr)-(dim(exprs)[[3]]+ncol(props))-1);
        
        output_prop <- rbind(output_prop, 
                               c(r2s[prop_indices]*100,
                                 r2_prop_adj*100,r2_adj*100,summary(model)$adj.r.squared*100))
        
        output_expr <- rbind(output_expr, 
                               c(r2s[expr_indices]*100,
                                 r2_expr_adj*100, r2_adj*100,summary(model)$adj.r.squared*100))
        
#         if(i == 1 || i == 2) {
#             cat("Processing: ",i,"\n")
#             cat(prop_indices,"\n")
#             cat(expr_indices,"\n")
#             cat(r2_prop,"\n");
#             cat(r2_prop_adj,"\n");
#             cat(r2_expr,"\n");
#             cat(r2_expr_adj,"\n")
#             cat(output_prop,"\n")
#             cat(output_expr,"\n")
#             cat(r2s*100,"\n")
#             cat(r2s[prop_indices],"\n");     
#             cat(r2s[expr_indices],"\n");    
#             cat(vars,"\n")

#         }

        ##cat(dim(output_expr),"\n")
        
    }
    colnames(output_prop) <- c(colnames(props), "adj r2","total prop_expr r2", "adj prop_expr r2");
    rownames(output_prop) <- dimnames(expr)[[1]];
    
    colnames(output_expr) <- c(dimnames(exprs)[[3]], "adj r2", "total prop_expr r2","adj prop_expr r2");
    rownames(output_expr) <- dimnames(expr)[[1]];
    
    return(list(output_prop, output_expr))
}



model_prop_expr_noah <- function(expr, props, exprs) {    
    output_prop <- NULL;
    output_expr <- NULL;
    output_covar_expr <- NULL;
    output_var2_expr <- NULL;
    output_expr_pbmc_count <- NULL;
    output_covar_expr_pbmc_count <- NULL;
    output_var2_expr_pbmc_count <- NULL;

    cov_prop_mat <- cov(props);
    expr_res <- expr;

    for(i in 1:nrow(expr)) {
        ## 1. calculate contribution from cell type together
        cov_expr_pbmc_count_mat <- cov(exprs[i,,]);
        
        vars_expr_pbmc_count <- NULL;
        covars_expr_pbmc_count <- NULL;
        vars2_expr_pbmc_count <- NULL; ## vars2 is adding up vars and covars and compartmentalizing to specific cell types
        
        for(j in 1:dim(exprs)[[3]]) {
            var_j_expr_pbmc_count = diag(cov_expr_pbmc_count_mat)[[j]]; ##sum(cov_expr_mat[j,]);
            covar_j_expr_pbmc_count = sum(cov_expr_pbmc_count_mat[j,])-diag(cov_expr_pbmc_count_mat)[[j]];
            vars_expr_pbmc_count <- c(vars_expr_pbmc_count, var_j_expr_pbmc_count);
            covars_expr_pbmc_count <- c(covars_expr_pbmc_count, covar_j_expr_pbmc_count)
            vars2_expr_pbmc_count <- c(vars2_expr_pbmc_count, sum(cov_expr_pbmc_count_mat[j,]))
        }   
  
        r2s_expr_pbmc_count = vars_expr_pbmc_count/var(expr[i,]);
        r2_expr_pbmc_count = sum(vars_expr_pbmc_count)/var(expr[i,]);
        r2_expr_pbmc_count_adj = r2_expr_pbmc_count-(1-r2_expr_pbmc_count)*dim(exprs)[[3]]/(ncol(expr)-dim(exprs)[[3]]-1);
    
        r2s_covar_expr_pbmc_count = covars_expr_pbmc_count/var(expr[i,]);
        r2_covar_expr_pbmc_count = sum(covars_expr_pbmc_count)/var(expr[i,]);
        r2_covar_expr_pbmc_count_adj = r2_covar_expr_pbmc_count-(1-r2_covar_expr_pbmc_count)*dim(exprs)[[3]]/(ncol(expr)-dim(exprs)[[3]]-1);
        
        r2s_var2_expr_pbmc_count = vars2_expr_pbmc_count/var(expr[i,]);
        r2_var2_expr_pbmc_count = sum(vars2_expr_pbmc_count)/var(expr[i,]);
        r2_var2_expr_pbmc_count_adj = r2_var2_expr_pbmc_count-(1-r2_var2_expr_pbmc_count)*dim(exprs)[[3]]/(ncol(expr)-dim(exprs)[[3]]-1);
    
        output_expr_pbmc_count <- rbind(output_expr_pbmc_count, 
                               c(r2s_expr_pbmc_count*100,
                                 r2_expr_pbmc_count*100,
                                 r2_expr_pbmc_count_adj*100))

        output_covar_expr_pbmc_count <- rbind(output_covar_expr_pbmc_count, 
                               c(r2s_covar_expr_pbmc_count*100,
                                 r2_covar_expr_pbmc_count*100,
                                 r2_covar_expr_pbmc_count_adj*100))
    
        output_var2_expr_pbmc_count <- rbind(output_var2_expr_pbmc_count, 
                               c(r2s_var2_expr_pbmc_count*100,
                                 r2_var2_expr_pbmc_count*100,
                                 r2_var2_expr_pbmc_count_adj*100))
    
        ## 2. calculate contribution from proportion only        
        model <- lm(expr[i,]~props-1)
        coef = summary(model)$coef[,1];
        expr_res[i,] <- residuals(model);

        vars_prop = NULL;
        
        for(j in 1:ncol(props)) {
            var_j_prop = sum(cov_prop_mat[j,]*coef[j]*coef);
            vars_prop <- c(vars_prop, var_j_prop);
        }
        
        r2s_prop = vars_prop/var(expr[i,]);
        r2_prop = sum(vars_prop)/var(expr[i,]);
        r2_prop_adj = r2_prop-(1-r2_prop)*(ncol(props))/(ncol(expr)-ncol(props)-1);
        
        ## 3. calculate contribution from expression after residualizing for proportion
        cov_expr_mat <- cov(exprs[i,,]-props*matrix(coef,nrow(props),length(coef),byrow=T));
        
        vars_expr <- NULL;
        covars_expr <- NULL;
        vars2_expr <- NULL; ## vars2 is adding up vars and covars and compartmentalizing to specific cell types
        
        for(j in 1:dim(exprs)[[3]]) {
            var_j_expr = diag(cov_expr_mat)[[j]]; ##sum(cov_expr_mat[j,]);
            covar_j_expr = sum(cov_expr_mat[j,])-diag(cov_expr_mat)[[j]];
            vars_expr <- c(vars_expr, var_j_expr);
            covars_expr <- c(covars_expr, covar_j_expr)
            vars2_expr <- c(vars2_expr, sum(cov_expr_mat[j,]))
        }   
  
        r2s_expr = vars_expr/var(expr[i,]);
        r2_expr = sum(vars_expr)/var(expr[i,]);
        r2_expr_adj = r2_expr-(1-r2_expr)*dim(exprs)[[3]]/(ncol(expr)-dim(exprs)[[3]]-1);
    
        r2s_covar_expr = covars_expr/var(expr[i,]);
        r2_covar_expr = sum(covars_expr)/var(expr[i,]);
        r2_covar_expr_adj = r2_covar_expr-(1-r2_covar_expr)*dim(exprs)[[3]]/(ncol(expr)-dim(exprs)[[3]]-1);
        
        r2s_var2_expr = vars2_expr/var(expr[i,]);
        r2_var2_expr = sum(vars2_expr)/var(expr[i,]);
        r2_var2_expr_adj = r2_var2_expr-(1-r2_var2_expr)*dim(exprs)[[3]]/(ncol(expr)-dim(exprs)[[3]]-1);

        output_prop <- rbind(output_prop, 
                            c(r2s_prop*100,
                            r2_prop*100,
                            r2_prop_adj*100,
                            (r2_var2_expr+r2_prop)*100,
                            (r2_var2_expr_adj+r2_prop_adj)*100))

        output_expr <- rbind(output_expr, 
                               c(r2s_expr*100,
                                 r2_expr*100,
                                 r2_expr_adj*100,
                                 (r2_expr+r2_prop)*100,
                                 (r2_expr_adj+r2_prop_adj)*100))##, sum(vars_expr), var(expr[i,])))
        
        output_covar_expr <- rbind(output_covar_expr, 
                               c(r2s_covar_expr*100,
                                 r2_covar_expr*100,
                                 r2_covar_expr_adj*100,
                                 (r2_expr+r2_covar_expr+r2_prop)*100,
                                 (r2_expr_adj+r2_covar_expr_adj+r2_prop_adj)*100))##,sum(covars_expr), var(expr[i,])))

        output_var2_expr <- rbind(output_var2_expr, 
                               c(r2s_var2_expr*100,
                                 r2_var2_expr*100,
                                 r2_var2_expr_adj*100,
                                 (r2_var2_expr+r2_prop)*100,
                                 (r2_var2_expr_adj+r2_prop_adj)*100))##,sum(covars_expr), var(expr[i,])))
        
     }

    # colnames(output_expr_pbmc_count) <- c(dimnames(exprs)[[3]], "total r2", "adj r2");
    # rownames(output_expr_pbmc_count) <- dimnames(expr)[[1]]

    # colnames(output_covar_expr_pbmc_count) <- c(dimnames(exprs)[[3]], "total r2", "adj r2");
    # rownames(output_covar_expr_pbmc_count) <- dimnames(expr)[[1]]

    # colnames(output_var2_expr_pbmc_count) <- c(dimnames(exprs)[[3]], "total r2", "adj r2");
    # rownames(output_var2_expr_pbmc_count) <- dimnames(expr)[[1]]

    # colnames(output_prop) <- c(colnames(props), "total r2","adj r2","adj prop_expr r2");
    # rownames(output_prop) <- dimnames(expr)[[1]]

    # colnames(output_expr) <- c(dimnames(exprs)[[3]], "total r2", "adj r2", "adj prop_expr r2");
    # rownames(output_expr) <- dimnames(expr)[[1]]
    
    # colnames(output_covar_expr) <- c(dimnames(exprs)[[3]], "total r2", "adj r2", "adj prop_expr r2");
    # rownames(output_covar_expr) <- dimnames(expr)[[1]]

    # colnames(output_var2_expr) <- c(dimnames(exprs)[[3]], "total r2", "adj r2", "adj prop_expr r2");
    # rownames(output_var2_expr) <- dimnames(expr)[[1]]

    colnames(output_expr_pbmc_count) <- c(dimnames(exprs)[[3]], "total r2", "adj r2");
    rownames(output_expr_pbmc_count) <- dimnames(expr)[[1]]

    colnames(output_covar_expr_pbmc_count) <- c(dimnames(exprs)[[3]], "total r2", "adj r2");
    rownames(output_covar_expr_pbmc_count) <- dimnames(expr)[[1]]

    colnames(output_var2_expr_pbmc_count) <- c(dimnames(exprs)[[3]], "total r2", "adj r2");
    rownames(output_var2_expr_pbmc_count) <- dimnames(expr)[[1]]

    colnames(output_prop) <- c(colnames(props), "total r2", "adj r2", "total prop_expr r2", "adj prop_expr r2");
    rownames(output_prop) <- dimnames(expr)[[1]]

    colnames(output_expr) <- c(dimnames(exprs)[[3]], "total r2", "adj r2", "total prop_expr r2", "adj prop_expr r2");
    rownames(output_expr) <- dimnames(expr)[[1]]
    
    colnames(output_covar_expr) <- c(dimnames(exprs)[[3]], "total r2", "adj r2", "total prop_expr r2", "adj prop_expr r2");
    rownames(output_covar_expr) <- dimnames(expr)[[1]]

    colnames(output_var2_expr) <- c(dimnames(exprs)[[3]], "total r2", "adj r2", "total prop_expr r2", "adj prop_expr r2");
    rownames(output_var2_expr) <- dimnames(expr)[[1]]

    
    return(list(output_prop, output_expr, output_covar_expr, output_var2_expr, output_expr_pbmc_count, output_covar_expr_pbmc_count, output_var2_expr_pbmc_count))
}

## modeling proportion and expr together
model_prop_expr_nnls <- function(expr, props, exprs) {
    output_prop <- NULL;
    output_expr <- NULL;
    
    for(i in 1:nrow(expr)) {
        
        cov_mat <- cov(cbind(props, exprs[i,,]));
        model <- nnls(cbind(props,exprs[i,,]), expr[i,]);
        ##model <- lm(expr[i,]~props+exprs[i,,]-1);
        coef <- model$x;
        
        vars <- NULL;
        
        #j is each cell type
        for(j in 1:(ncol(props)+dim(exprs)[[3]])) {
            var_j = sum(cov_mat[j,]*coef[j]*coef);
            vars <- c(vars, var_j);
        }
        
        ##vars <- unlist(vars);
        ##print(vars);
        r2s = vars/sum(vars);
        
        ##print(r2);
        
        prop_indices <- 1:ncol(props);
        expr_indices <- (ncol(props)+1):(ncol(props)+dim(exprs)[[3]]);
        
        r2_prop = sum(vars[prop_indices])/var(expr[i,]);
        r2_prop_adj = r2_prop;
        ##r2_prop_adj = 1-((1-r2_prop)*(ncol(expr)-1))/(ncol(expr)-ncol(props)-1);             
        ##r2_prop_adj = r2_prop-(1-r2_prop)*(dim(exprs)[[3]]+ncol(props))/(ncol(expr)-(dim(exprs)[[3]]+ncol(props))-1);             
        
        r2_expr = sum(vars[expr_indices])/var(expr[i,]);
        r2_expr_adj = r2_expr;
        ##r2_expr_adj = 1-((1-r2_expr)*(ncol(expr)-1))/(ncol(expr)-dim(exprs)[[3]]-1);             
        ##r2_expr_adj = r2_expr-(1-r2_expr)*(dim(exprs)[[3]]+ncol(props))/(ncol(expr)-(dim(exprs)[[3]]+ncol(props))-1);
        
        r2 = sum(vars)/var(expr[i,]);
        ##r2_adj = 1-((1-r2)*(ncol(expr)-1))/(ncol(expr)-dim(exprs)[[3]]+ncol(props)-1); 
        r2_adj = r2-(1-r2)*(dim(exprs)[[3]]+ncol(props))/(ncol(expr)-(dim(exprs)[[3]]+ncol(props))-1);
        
        output_prop <- rbind(output_prop, 
                               c(r2s[prop_indices]*100,
                                 r2_prop_adj*100,r2_adj*100))
        
        output_expr <- rbind(output_expr, 
                               c(r2s[expr_indices]*100,
                                 r2_expr_adj*100, r2_adj*100))
        
#         if(i == 1 || i == 2) {
#             cat("Processing: ",i,"\n")
#             cat(prop_indices,"\n")
#             cat(expr_indices,"\n")
#             cat(r2_prop,"\n");
#             cat(r2_prop_adj,"\n");
#             cat(r2_expr,"\n");
#             cat(r2_expr_adj,"\n")
#             cat(output_prop,"\n")
#             cat(output_expr,"\n")
#             cat(r2s*100,"\n")
#             cat(r2s[prop_indices],"\n");     
#             cat(r2s[expr_indices],"\n");    
#             cat(vars,"\n")

#         }

        ##cat(dim(output_expr),"\n")
        
    }
    colnames(output_prop) <- c(colnames(props), "adj r2","total prop_expr r2");
    rownames(output_prop) <- dimnames(expr)[[1]];
    
    colnames(output_expr) <- c(dimnames(exprs)[[3]], "adj r2", "total prop_expr r2");
    rownames(output_expr) <- dimnames(expr)[[1]];
    
    return(list(output_prop, output_expr))
}

plot_bar_chart <- function(output_prop, output_expr, pheatmap_out, cuttree_out, normalize_r2 = TRUE) {
    colors = c("#E58606", "#5D69B1", "#52BCA3", "#99C945", "#CC61B0", "#24796C", "#DAA51B", "#764E9F")
    
#     colors = c('#1f77b4','#ff7f0e',
#  '#279e68',
# ## '#d62728',
#  '#aa40fc',
#  '#8c564b',
#  '#e377c2',
# ## '#b5bd61',
#  '#17becf',
#  '#aec7e8',
#  '#ffbb78',
#  '#98df8a',
#  '#ff9896',
#  '#c5b0d5',
#  '#c49c94',
#  '#f7b6d2',
#  '#dbdb8d',
#  '#9edae5',
#  '#ad494a',
#  '#8c6d31')
    
    options(repr.plot.width = 6, repr.plot.height = 15)
    
    output_prop <- output_prop[match(cuttree_out$gene, rownames(output_prop)),]
    
    if(normalize_r2) {
        prop_r2 = output_prop[,"total r2"]/100;
    } else {
        prop_r2 = 1;
    }
    
    
#     out <- rbind(data.frame(gene_names=rownames(output_prop), prop=output_prop[,"cM"]*prop_r2, cell="cM"),
#                  data.frame(gene_names=rownames(output_prop), prop=output_prop[,"Th"]*prop_r2, cell="Th"),
#                  data.frame(gene_names=rownames(output_prop), prop=output_prop[,"B"]*prop_r2, cell="B"),
#                  data.frame(gene_names=rownames(output_prop), prop=output_prop[,"NK"]*prop_r2, cell="NK"),
#                  data.frame(gene_names=rownames(output_prop), prop=output_prop[,"ncM"]*prop_r2, cell="ncM"),
#                  data.frame(gene_names=rownames(output_prop), prop=output_prop[,"Tc"]*prop_r2, cell="Tc"),
#                  data.frame(gene_names=rownames(output_prop), prop=output_prop[,"cDC"]*prop_r2, cell="cDC"),
#                   data.frame(gene_names=rownames(output_prop), prop=output_prop[,"pDC"]*prop_r2, cell="pDC")
#                 )
    
       out <- rbind(data.frame(gene_names=rownames(output_prop), prop=output_prop[,"cM"], prop_r2 = prop_r2, cell="cM"),
                 data.frame(gene_names=rownames(output_prop), prop=output_prop[,"Th"], prop_r2 = prop_r2, cell="Th"),
                 data.frame(gene_names=rownames(output_prop), prop=output_prop[,"B"], prop_r2 = prop_r2, cell="B"),
                 data.frame(gene_names=rownames(output_prop), prop=output_prop[,"NK"], prop_r2 = prop_r2, cell="NK"),
                 data.frame(gene_names=rownames(output_prop), prop=output_prop[,"ncM"], prop_r2 = prop_r2, cell="ncM"),
                 data.frame(gene_names=rownames(output_prop), prop=output_prop[,"Tc"], prop_r2 = prop_r2, cell="Tc"),
                 data.frame(gene_names=rownames(output_prop), prop=output_prop[,"cDC"], prop_r2 = prop_r2, cell="cDC"),
                  data.frame(gene_names=rownames(output_prop), prop=output_prop[,"pDC"], prop_r2 = prop_r2, cell="pDC")
                )

    ##matched <- match(out$gene_names, rev(pheatmap_out$tree_row$labels[pheatmap_out$tree_row$order]))
    ##matched <- match(out$gene_names, cuttree_out$gene)
    ##out <- out[which(!is.na(matched)),]
    ##out_ordered <- out[order(na.omit(matched)),]
#     out_ordered <- cbind(out_ordered, cluster=cuttree_out[match(out_ordered$gene_names,cuttree_out[,"gene"]),"cluster"])

    ##out_ordered <- out[match(cuttree_out$gene, out$gene_names),]
    out_ordered <- cbind(out, cluster=cuttree_out$cluster);
    
    out_ordered$gene_names <- factor(out_ordered$gene_names,levels=c(unique(as.character(out_ordered$gene_names))))
    out_ordered$cell = factor(out_ordered$cell, levels=c('cM', 'ncM', 'cDC', 'pDC', 'Th', 'Tc', 'NK','B'));
    
    
    prop_plot <- ggplot(out_ordered)+geom_bar(aes(x=factor(gene_names,levels=rev(levels(gene_names))),y=prop,fill=cell),stat="identity")+
    ##geom_point(aes(x=gene_names, y=as.numeric(prop_r2)*100, group=1), stat='summary', fun.y=sum)+
    ##stat_summary(fun.y=sum, geom="line")
    geom_line(aes(x=factor(gene_names,levels=rev(levels(gene_names))), y=as.numeric(prop_r2)*100, group=1),stat="identity")+
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.spacing.y=unit(0.3, "lines"),
          ##aspect.ratio = 2,
          plot.margin = unit(c(1,-2,1,0), "mm"),
          legend.position="none")+
    facet_grid(cluster~.,scales = "free_y",space="free_y",switch="y")+
    ##facet_grid(~cluster,scales = "free_x",space="free_x")+
    scale_y_reverse(limits=c(120,-10),breaks=c(0, 25, 50, 75, 100))+
    coord_flip()+
    scale_fill_manual(values=colors)
    
#     prop_plot <- ggplot(out_ordered)
#     +geom_bar(aes(x=as.numeric(gene_names),y=prop,fill=cell),stat="identity",width=1)
#     +geom_line(aes(x=as.numeric(gene_names), y=as.numeric(prop_r2)*100),stat="identity")+theme(axis.title.x = element_blank(), 
#                                                                                                axis.title.y = element_blank(), 
#                                                                                                axis.text.y = element_blank(), 
#                                                                                                axis.ticks.y = element_blank(), 
#                                                                                                plot.margin = unit(c(1,-2,1,0), "mm"),
#                                                                                                legend.position="none")
#     +facet_grid(cluster~.,scales = "free_y",space="free_y",switch="y")
#     +scale_y_reverse(limits=c(120,-10),breaks=c(0, 25, 50, 75, 100))
#     +coord_flip()
#     +scale_fill_manual(values=colors)

            ##panel.spacing.y=unit(0, "lines"), 
        ##panel.spacing.x=unit(0, "lines"),

#     prop_plot <- ggplot(aes(x=gene_names,y=prop,fill=cell),data=out_ordered)+geom_bar(stat="identity")+theme(axis.title.x = element_blank(), 
#         axis.title.y = element_blank(), 
#         axis.text.y = element_blank(), 
#         ##axis.text.y = element_text(size=5),
#         axis.ticks.y = element_blank(), 
#         plot.margin = unit(c(1,-1,1,0), "mm"),
#         panel.spacing.y=unit(0.35, "lines"),                                                                                                           
#         legend.position="none")+facet_grid(cluster~.,scales = "free_y",space="free_y",switch="y")+scale_y_reverse(limits=c(100,-10),breaks=c(0, 25, 50, 75, 100))+coord_flip()+
#         scale_fill_manual(values=colors)
    
#     prop_plot <- ggplot(aes(x=gene_names,y=prop,fill=cell),data=out_ordered)+geom_bar(stat="identity")+theme(axis.title.x = element_blank(), 
#         axis.title.y = element_blank(), 
#         axis.text.y = element_blank(), 
#         ##axis.text.y = element_text(size=5),
#         axis.ticks.y = element_blank(), 
#         plot.margin = unit(c(1,-1,1,0), "mm"),
#         legend.position="none")+facet_grid(cluster~.,scales = "free_y",space="free_y",switch="y")+scale_y_reverse()+coord_flip()

    options(repr.plot.width = 6, repr.plot.height = 15)
    
    output_expr <- output_expr[match(cuttree_out$gene, rownames(output_expr)),]
    
    if(normalize_r2) {
        expr_r2 = output_expr[,"total r2"]/100;
        ##expr_r2 = (100-output_prop[,"adj r2"])/100;
    } else {
        expr_r2 = 1;
    }
    
#     out_expr <- rbind(data.frame(gene_names=rownames(output_expr), prop=output_expr[,"cM"]*expr_r2, cell="cM"),
#              data.frame(gene_names=rownames(output_expr), prop=output_expr[,"Th"]*expr_r2, cell="Th"),
#              data.frame(gene_names=rownames(output_expr), prop=output_expr[,"B"]*expr_r2, cell="B"),
#              data.frame(gene_names=rownames(output_expr), prop=output_expr[,"NK"]*expr_r2, cell="NK"),
#              data.frame(gene_names=rownames(output_expr), prop=output_expr[,"ncM"]*expr_r2, cell="ncM"),
#              data.frame(gene_names=rownames(output_expr), prop=output_expr[,"Tc"]*expr_r2, cell="Tc"),
#              data.frame(gene_names=rownames(output_expr), prop=output_expr[,"cDC"]*expr_r2, cell="cDC"),
#             data.frame(gene_names=rownames(output_expr), prop=output_expr[,"pDC"]*expr_r2, cell="pDC")
#                      )
   
        out_expr <- rbind(data.frame(gene_names=rownames(output_expr), prop=output_expr[,"cM"], expr_r2=expr_r2, cell="cM"),
             data.frame(gene_names=rownames(output_expr), prop=output_expr[,"Th"], expr_r2=expr_r2, cell="Th"),
             data.frame(gene_names=rownames(output_expr), prop=output_expr[,"B"], expr_r2=expr_r2, cell="B"),
             data.frame(gene_names=rownames(output_expr), prop=output_expr[,"NK"], expr_r2=expr_r2, cell="NK"),
             data.frame(gene_names=rownames(output_expr), prop=output_expr[,"ncM"], expr_r2=expr_r2, cell="ncM"),
             data.frame(gene_names=rownames(output_expr), prop=output_expr[,"Tc"], expr_r2=expr_r2, cell="Tc"),
             data.frame(gene_names=rownames(output_expr), prop=output_expr[,"cDC"], expr_r2=expr_r2, cell="cDC"),
            data.frame(gene_names=rownames(output_expr), prop=output_expr[,"pDC"], expr_r2=expr_r2, cell="pDC")
                     )

            ##cuttree_ordered <- cuttree_out[order(cuttree_out[,"cluster"]),]
#     ##matched_expr <- match(out_expr$gene_names, rev(pheatmap_out$tree_row$labels[pheatmap_out$tree_row$order]))
#     matched_expr <- match(out_expr$gene_names, cuttree_out$gene)
#     out_expr <- out_expr[which(!is.na(matched)),]
#     out_expr_ordered <- out_expr[order(na.omit(matched)),]
#     out_expr_ordered <- cbind(out_expr_ordered, cluster=cuttree_out[match(out_expr_ordered$gene_names,cuttree_out[,"gene"]),"cluster"])

    ##out_expr_ordered <- out_expr[match(cuttree_out$gene, out_expr$gene_names),]
    out_expr_ordered <- cbind(out_expr, cluster=cuttree_out$cluster);

    ##out_expr_ordered$cluster <- factor(out_expr_ordered$cluster, levels=c('dc_specific','tc_specific','th_specific','all_ifn','myeloid_ifn','ncM_ifn','cM_DC_ifn','cM_ifn'))

    out_expr_ordered$gene_names <- factor(out_expr_ordered$gene_names,levels=c(unique(as.character(out_expr_ordered$gene_names))))
    out_expr_ordered$cell = factor(out_expr_ordered$cell, levels=c('cM', 'ncM', 'cDC', 'pDC', 'Th', 'Tc', 'NK', 'B'));

    expr_plot <- ggplot(out_expr_ordered)+geom_bar(aes(x=factor(gene_names,levels=rev(levels(gene_names))),y=prop,fill=cell),stat="identity")+##,width=1)+
    ##geom_point(aes(x=gene_names, y=as.numeric(expr_r2)*100,group=1), stat='summary', fun.y=sum)+
    ##stat_summary(fun.y=sum, geom="line")
    geom_line(aes(x=factor(gene_names,levels=rev(levels(gene_names))), y=as.numeric(expr_r2)*100, group=1),stat="identity")+
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.spacing.y=unit(0.3, "lines"),
          ##aspect.ratio = 2,
          plot.margin = unit(c(1,-2,1,0), "mm"),
          legend.position="none")+
    facet_grid(cluster~.,scales = "free_y",space="free_y")+
    ##facet_grid(~cluster,scales = "free_x",space="free_x")+
    scale_y_continuous(limits=c(-10, 120),breaks=c(0, 25, 50, 75, 100))+
    coord_flip()+
    scale_fill_manual(values=colors)
    
#     expr_plot <- ggplot(out_expr_ordered)+geom_bar(stat="identity",aes(x=as.numeric(gene_names),y=prop,fill=cell),width=1)+geom_line(aes(x=as.numeric(gene_names), y=expr_r2*100), stat="identity")+theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
#         axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#         plot.margin = unit(c(1,0,1,-2), "mm"),legend.position="none")
#     +facet_grid(cluster~.,scales = "free_y",space="free_y")+coord_flip()+scale_fill_manual(values=colors)+scale_y_continuous(limits=c(-10,120),breaks=c(0, 25, 50, 75, 100))

            #panel.spacing.y=unit(0, "lines"),
        #panel.spacing.x=unit(0, "lines"),

    
    ##expr_plot <- expr_plot + ggplot(aes(x=gene_names, y=expr_r2), data=out_expr)+geom_point();

#     expr_plot <- ggplot(aes(x=gene_names,y=prop,fill=cell),data=out_expr_ordered)+geom_bar(stat="identity")+theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
#         axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#         plot.margin = unit(c(1,0,1,-1), "mm"),legend.position="none")+facet_grid(cluster~.,scales = "free_y",space="free_y")+coord_flip()

    gene_name <- ggplot(aes(x=gene_names,y=prop,fill=cell),data=out_expr_ordered)+geom_bar(stat="identity")+ylim(-10,120)+facet_grid(cluster~.,scales = "free_y",space="free_y")+coord_flip()+theme(panel.spacing.y=unit(0.3, "lines"));
    
    library(gridExtra)
    gg1 <- ggplot_gtable(ggplot_build(prop_plot))
    gg2 <- ggplot_gtable(ggplot_build(expr_plot))

    ##plot_out = grid.arrange(gg1, gg2,ncol=2,widths=c(4/9,4/9))
    plot_out = arrangeGrob(gg1, gg2,ncol=2,widths=c(4/9,4/9))
    return(plot_out);
}

##ggplot(aes(x=gene_names,y=prop,fill=cell),data=out_ordered)+geom_bar(stat="identity")+facet_grid(~cluster,scales = "free_x", space = "free_x")
