#Compute the number of gene-metabolite pairs tested in pairwise tissues  


library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)

dir <- "results/ukbb_gtex_multi_intact/gppc_rst/"

tissues <- fread("data/gtex_v8_tissue.txt",h=F)$V1

outfile <- "results/ukbb_gtex_multi_intact/n_pairs_tested_pairwise.txt"

pair_df <- data.frame(tissue1 = rep(tissues,1,each=49),tissue2 = rep(tissues,49))

outdf <- NULL

for (i in 1:nrow(pair_df)){

        tis1 <- fread(paste0(dir,pair_df$tissue1[i],"/C100000007_",pair_df$tissue1[i],"_gppc_rst.txt"),h=T)

	tis2 <- fread(paste0(dir,pair_df$tissue2[i],"/C100000007_",pair_df$tissue2[i],"_gppc_rst.txt"),h=T)

	n_overlap <- length(intersect(tis1$ensid_protein_id_chr,tis2$ensid_protein_id_chr))

        tmp <- data.frame(Tissue1 = pair_df$tissue1[i],Tissue2 = pair_df$tissue2[i],n_overlap = n_overlap*1408)

        outdf <- rbind.data.frame(outdf,tmp)
}

write.table(outdf,file=outfile,quote=F,sep = '\t',row.names = F, col.names = T)
