#plot histogram showing: x-axis: #tissues, y-axis: number of metabolite-PCG pairs identified

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)


implicated_pairs_dir <- "results/ukbb_gtex_multi_intact/gppc_rst/sig_rst/"

outfile <- "results/ukbb_gtex_multi_intact/tissue_concordance_histogram.pdf"

pairs_implicated <- fread(cmd = paste0("cat ",implicated_pairs_dir,"*_sig_rst.txt"),h=F)

pairs_implicated <- pairs_implicated %>% dplyr::select(V1,V19,V20)

colnames(pairs_implicated) <- c("gene","tissue","cid")

pairs_implicated$protein_id <- str_split_fixed(pairs_implicated$gene,":",n = 3)[,2]
pairs_implicated$symbol <- gsub("\\_.*", "", pairs_implicated$protein_id)
pairs_implicated$cid_gene <- paste0(pairs_implicated$cid,"_",pairs_implicated$symbol)


pdf(outfile)
pairs_implicated %>%
	group_by(cid_gene) %>%
	summarise(n = n()) %>%
	ggplot(aes(x = n)) + 
	geom_histogram(stat="count") +
	xlab("Number of tissues in which the PCG-metabolite pair is identified") +
	ylab("Number of PCG-metabolite pairs") +
	theme_bw() +
	theme(text = element_text(size = 10,face= "bold"))
dev.off()


#which tissue has the most metabolite-PCG pairs?

pairs_implicated %>%
	group_by(tissue) %>%
	summarise(n = n()) %>%
	arrange(desc(n))

#What percent of pairs implicated in a single tissue?

single_tiss_pairs <- pairs_implicated %>%
        group_by(cid_gene) %>%
        summarise(n = n()) %>%
	group_by(n) %>%
	summarise(n_pairs = n()) %>%
	data.frame() %>% 
	filter(n == 1) %>%
	dplyr::select(n_pairs) %>%
	unlist() %>% as.numeric()

print(single_tiss_pairs/length(unique(pairs_implicated$cid_gene)))

#What percent of pairs implicated in 2-10 tissues?

few_tiss_pairs <- pairs_implicated %>%
        group_by(cid_gene) %>%
        summarise(n = n()) %>%
        group_by(n) %>%
        summarise(n_pairs = n()) %>%
        data.frame() %>%
        filter(n > 1 & n <= 10)

print(sum(few_tiss_pairs$n_pairs)/length(unique(pairs_implicated$cid_gene)))

#What percent of pairs implicated 40+ tissues?

many_tiss_pairs <- pairs_implicated %>%
        group_by(cid_gene) %>%
        summarise(n = n()) %>%
        group_by(n) %>%
        summarise(n_pairs = n()) %>%
        data.frame() %>%
        filter(n > 40)

print(sum(many_tiss_pairs$n_pairs)/length(unique(pairs_implicated$cid_gene)))

