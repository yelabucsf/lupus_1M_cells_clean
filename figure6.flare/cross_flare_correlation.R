library(pheatmap);

lymphoid <- fread("~/Downloads/Lymphocyte_CrossX_Flare_Correlations.csv")
lymphoid_small <- lymphoid[,-1]
rownames(lymphoid_small) <- lymphoid$V1
pheatmap(lymphoid_small)


myeloid <- fread("~/Downloads/Myeloid_CrossX_Flare_Correlations.csv")
myeloid_small <- myeloid[,-1]
rownames(myeloid_small) <- myeloid$V1
pheatmap(myeloid_small)

myeloid_clusters <- c("cM", "iM", "ncM", "actM", "compM", "Mac1", "Mac2", "cDC1", "cDC2","pDC", "2","0","10", "9", "22","13","16","21");
cross_myeloid_indices <- match(myeloid_clusters, rownames(myeloid_small))[1:10];
flare_myeloid_indices <- match(myeloid_clusters, rownames(myeloid_small))[11:18];
myeloid_names <- c("Cross_cM", "Cross_iM", "Cross_ncM", "Cross_actM", "Cross_compM", "Cross_Mac1", "Cross_Mac2", "Cross_cDC1", "Cross_cDC2", "Cross_pDC",
                   "Flare_cM", "Flare_iM", "Flare_ncM", "Flare_actM", "Flare_compM", "Flare_Mac", "Flare_cDC", "Flare_pDC");
myeloid_small <- myeloid_small[cross_myeloid_indices,..flare_myeloid_indices];
rownames(myeloid_small) <- myeloid_names[1:10];
colnames(myeloid_small) <- myeloid_names[11:18];

pheatmap(myeloid_small, cluster_rows=F, cluster_cols=F, scale="column")




lymphoid_clusters <- c("B atypical", "B memory", "B naive", "MK contam", "NK1", "NK2", "NK3", "Plasmablast", "Progenitor", "Prolif T", "T4 IFN","T4 memory","T4 naive", "T4 regulatory", "T8 cyto1","T8 cyto2","T8 cyto3","T8 memory", "T8 naive",
                       "19", "12", "4", "20", "6", "17", "27", "29", "23", "3", "5", "14", "1", "11", "7");
cross_lymphoid_indices <- match(lymphoid_clusters, rownames(lymphoid_small))[1:10];
flare_lymphoid_indices <- match(lymphoid_clusters, rownames(lymphoid_small))[11:18];
lymphoid_names <- c("Cross_cM", "Cross_iM", "Cross_ncM", "Cross_actM", "Cross_compM", "Cross_Mac1", "Cross_Mac2", "Cross_cDC1", "Cross_cDC2", "Cross_pDC",
                   "Flare_cM", "Flare_iM", "Flare_ncM", "Flare_actM", "Flare_compM", "Flare_Mac", "Flare_cDC", "Flare_pDC");
lymphoid_small <- lymphoid_small[cross_lymphoid_indices,..flare_lymphoid_indices];
rownames(lymphoid_small) <- lymphoid_names[1:10];
colnames(lymphoid_small) <- lymphoid_names[11:18];

pheatmap(lymphoid_small, cluster_rows=F, cluster_cols=F, scale="column")
