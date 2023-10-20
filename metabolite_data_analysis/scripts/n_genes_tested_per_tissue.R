#summarise the number of genes tested in Multi-INTACT analysis by tissue

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)

dir <- "results/ukbb_gtex_multi_intact/gppc_rst/"

tissues <- fread("data/gtex_v8_tissue.txt",h=F)$V1

outfile <- "results/ukbb_gtex_multi_intact/n_tested_m_intact.pdf"

plotdf <- NULL

for (i in 1:length(tissues)){

	rst <- fread(paste0(dir,tissues[i],"/C100000007_",tissues[i],"_gppc_rst.txt"),h=T)

	n_tested <- nrow(rst)

	tmp <- data.frame(Tissue = tissues[i],n_tested = n_tested)

	plotdf <- rbind.data.frame(plotdf,tmp)
}


tissue_order <- plotdf %>%
	arrange(desc(n_tested)) %>%
	dplyr::select(Tissue) %>% unlist() %>% as.vector()

pdf(outfile)
plotdf %>%
	mutate(Tissue = factor(Tissue)) %>%
	mutate(Tissue = fct_relevel(Tissue,tissue_order)) %>%
	ggplot(aes(x = Tissue,y=n_tested)) +
	geom_col() +
	xlab("Tissue") +
	ylab("Number of genes tested by Multi-INTACT") +
	theme_bw() + 
	theme(text = element_text(size = 10, face = "bold"),
	      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()
