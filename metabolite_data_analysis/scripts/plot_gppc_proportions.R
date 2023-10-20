#Visualize proprotions of multi-INTACT genes implicated by marginal analyses

library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)


tissue_file <- "data/gtex_v8_tissue.txt"

out_file <- "results/ukbb_gtex_multi_intact/gppc_proportion_plot.pdf"

tissues <- fread(tissue_file,h=F)$V1

df <- NULL

for (i in 1:length(tissues)){

        tmp <- fread(paste0("results/ukbb_gtex_multi_intact/gppc_rst/gppc_summary/", tissues[i],"_gppc_summary.txt"))

        tmp <- tmp %>% dplyr::select(V1,V2,V3,V6,V7,V8,V9)

        colnames(tmp) <- c("CID","Tissue","n_sig_scan","n_scan_only","n_expr_only","n_protein_only","n_expr_and_protein")

        df <- rbind.data.frame(df,tmp)
}

#set tissue order
tissue_order <- df %>%
        group_by(Tissue) %>%
        summarise(n_pairs_scan = sum(n_sig_scan)) %>%
        arrange(desc(n_pairs_scan)) %>%
        dplyr::select(Tissue) %>% unlist() %>% as.vector()


pdf(out_file)
df %>%
        group_by(Tissue) %>%
        summarise(n_pairs_scan_only = sum(n_scan_only),
                  n_pairs_expr_only = sum(n_expr_only),
                  n_pairs_protein_only = sum(n_protein_only),
		  n_pairs_expr_protein = sum(n_expr_and_protein)) %>%
        dplyr::rename("Multi_INTACT_scan_only"="n_pairs_scan_only",
                      "INTACT_expression"="n_pairs_expr_only",
                      "INTACT_protein"="n_pairs_protein_only",
		      "INTACT_expression_and_protein" = "n_pairs_expr_protein") %>%
        data.frame() %>%
	pivot_longer(cols = Multi_INTACT_scan_only:INTACT_expression_and_protein,
			names_to = "Type",
			values_to = "N_pairs") %>%
        mutate(Type = case_when(Type == "Multi_INTACT_scan_only" ~ "Multi-INTACT only",
                                  Type == "INTACT_expression" ~ "INTACT (Expression)",
                                  Type == "INTACT_protein" ~ "INTACT (Protein)",
				  Type == "INTACT_expression_and_protein" ~ "INTACT (Expression & Protein)")) %>%
	mutate(Type = factor(Type)) %>%
	mutate(Type = fct_relevel(Type, c("Multi-INTACT only","INTACT (Expression)","INTACT (Protein)","INTACT (Expression & Protein)"))) %>%
	mutate(Tissue = factor(Tissue)) %>%
        mutate(Tissue = fct_relevel(Tissue, rev(tissue_order))) %>%
        ggplot(aes(x=Tissue,y=N_pairs,fill = Type)) +
        geom_bar(stat = 'identity') +
	coord_flip() +
#        geom_text(aes(label=n_pairs),hjust = -0.15,size = 2.75) +
        ylab("Number of significant gene-metabolite pairs implicated by Multi-INTACT") +
#        scale_y_continuous(limits = c(0,450)) +
        theme_bw() +
        theme(text = element_text(size = 10, face = 'bold'),
				  legend.title=element_blank())
dev.off()


