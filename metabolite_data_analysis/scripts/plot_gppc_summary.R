#Visualize GPPC summary results by tissue

library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)


tissue_file <- "data/gtex_v8_tissue.txt"

out_file <- "results/ukbb_gtex_multi_intact/summary_gppc.pdf"

tissues <- fread(tissue_file,h=F)$V1

df <- NULL

for (i in 1:length(tissues)){

	tmp <- fread(paste0("results/ukbb_gtex_multi_intact/gppc_rst/gppc_summary/", tissues[i],"_gppc_summary.txt"))

	tmp <- tmp %>% dplyr::select(V1,V2,V3,V4,V5)
	
	colnames(tmp) <- c("CID","Tissue","n_sig_scan","n_sig_intact_expr","n_sig_intact_protein")

	df <- rbind.data.frame(df,tmp)
}


#set tissue order
tissue_order <- df %>%
        group_by(Tissue) %>%
        summarise(n_pairs_scan = sum(n_sig_scan),
                  n_pairs_intact_expr = sum(n_sig_intact_expr),
                  n_pairs_intact_protein = sum(n_sig_intact_protein)) %>%
        dplyr::rename("Multi_INTACT_scan"="n_pairs_scan",
                      "INTACT_expression"="n_pairs_intact_expr",
                      "INTACT_protein"="n_pairs_intact_protein") %>%
        data.frame() %>%
        arrange(desc(Multi_INTACT_scan)) %>%
	dplyr::select(Tissue) %>% unlist() %>% as.vector()


pdf(out_file)
df %>%
	group_by(Tissue) %>%
	summarise(n_pairs_scan = sum(n_sig_scan),
	          n_pairs_intact_expr = sum(n_sig_intact_expr),
		  n_pairs_intact_protein = sum(n_sig_intact_protein)) %>%
	dplyr::rename("Multi_INTACT"="n_pairs_scan",
		      "INTACT_Expression"="n_pairs_intact_expr",
		      "INTACT_Protein"="n_pairs_intact_protein") %>%
	data.frame() %>%
	arrange(desc(Multi_INTACT)) %>%
	pivot_longer(cols = c("Multi_INTACT","INTACT_Expression","INTACT_Protein"),
			names_to = "Method",
			values_to = "n_pairs") %>%
	mutate(Method = case_when(Method == "Multi_INTACT" ~ "Multi-INTACT",
				  Method == "INTACT_Expression" ~ "INTACT (Expression)",
				  Method == "INTACT_Protein" ~ "INTACT (Protein)")) %>%
        mutate(Tissue = factor(Tissue)) %>%
	mutate(Tissue = fct_relevel(Tissue, rev(tissue_order))) %>%
	ggplot(aes(x=Tissue,y=n_pairs)) + 
	geom_bar(stat = 'identity',fill='dark blue') +
	geom_text(aes(label=n_pairs),hjust = -0.15,size = 2.25) +
	ylab("Number of significant gene-metabolite pairs") + 
	scale_y_continuous(limits = c(0,450)) +
	facet_wrap(~Method) + 
	coord_flip() + 
	theme_bw() +
	theme(text = element_text(size = 10, face = 'bold'))
dev.off()
