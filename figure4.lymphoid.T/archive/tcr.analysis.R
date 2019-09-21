library(data.table)
library(ggplot2)

th <- fread("Th_giniIndex_hiseq.csv")
tc <- fread("Tc_giniIndex_hiseq.csv")

df <- rbind(cbind(th, type="T4"), cbind(tc, type="T8"))
df$disease = "SLE";
df$disease[grep("IGTB",df$sample)] = "Ctrl";

plot = ggplot(aes(disease, gini, color=type),data=df)+geom_boxplot()+theme_bw()
ggsave(plot, file="figures.V6/gini.png", height=3, width=3)


th_clones <- fread("Th_sharedclones_hiseq_CDR3.csv")
tc_clones <- fread("Tc_sharedclones_hiseq_CDR3.csv")

tcr <- fread("TCR_paired_final_hiseq.csv")


lineage <- fread("TCR_gini_lymphlineage.csv")
lineage$disease = "SLE";
lineage$disease[grep("IGTB",lineage$Sample_ID)] = "Ctrl";

plot = ggplot(aes(disease, coefficient, color=ct_cov),data=lineage)+geom_boxplot()+theme_bw()
ggsave(plot, file="figures.V6/gini_lineage.png", height=3, width=9)

high_res <- fread("TCR_gini_highreslymph.csv")
##high_res <- high_res[which(high_res$cell_count>),]
high_res$disease = "SLE";
high_res$disease[grep("IGTB",high_res$Sample_ID)] = "Ctrl";

plot = ggplot(aes(disease, coefficient, color=lymp_ct_cov),data=high_res)+geom_boxplot()+theme_bw()
ggsave(plot, file="figures.V6/gini_high_res.png", height=3, width=9)
