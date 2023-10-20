#plot heatmap showing (# gene-metab pairs implicated in both tissues)/(# gene-metab pairs tested in both tissues) for all tissue pairs

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)
#library(gplots)


pairs_tested_file <- "results/ukbb_gtex_multi_intact/n_pairs_tested_pairwise.txt"

implicated_pairs_dir <- "results/ukbb_gtex_multi_intact/gppc_rst/sig_rst/"

outfile <- "results/ukbb_gtex_multi_intact/prop_pairwise_overlap_over_union_heatmap_clustered.pdf"


pairs_tested <- fread(pairs_tested_file,h=T)

pairs_implicated <- fread(cmd = paste0("cat ",implicated_pairs_dir,"*_sig_rst.txt"),h=F)

pairs_implicated <- pairs_implicated %>% dplyr::select(V1,V19,V20)

colnames(pairs_implicated) <- c("gene","tissue","cid")

pairs_implicated$cid_gene <- paste0(pairs_implicated$cid,"_",pairs_implicated$gene)

pairs_tested$n_overlap_implicated <- NA

for (i in 1:nrow(pairs_tested)){

	implicated_tis1 <- pairs_implicated %>%
		filter(tissue == pairs_tested$Tissue1[i]) 

	implicated_tis2 <- pairs_implicated %>%
		filter(tissue == pairs_tested$Tissue2[i])

	n_overlap <- length(intersect(implicated_tis1$cid_gene,implicated_tis2$cid_gene))

	n_union <- length(union(implicated_tis1$cid_gene,implicated_tis2$cid_gene))

	pairs_tested$n_overlap_implicated[i] <- n_overlap

	pairs_tested$n_union[i] <- n_union
}

pairs_tested$prop_implicated_both_tissues <- pairs_tested$n_overlap_implicated/pairs_tested$n_overlap

#summary(pairs_tested$prop_implicated_both_tissues)



#Try #(gene-metab pairs implicated in both tissues)/#(union of gene-metab pairs implicated in both tissues)

pairs_tested$union_implicated_tissues <- pairs_tested$n_overlap_implicated/pairs_tested$n_union

#summary(pairs_tested$union_implicated_tissues)


#base r function for clustering

plot_mat <- pairs_tested %>%
	dplyr::select(Tissue1,Tissue2,union_implicated_tissues) %>%
	pivot_wider(names_from = Tissue2, values_from = union_implicated_tissues) %>%
	dplyr::select(-Tissue1)

plot_mat <- as.matrix(plot_mat)
rownames(plot_mat) <- colnames(plot_mat)


pdf(outfile)
print(pheatmap(plot_mat,
         fontsize_row=5,fontsize_col=5))
dev.off()
