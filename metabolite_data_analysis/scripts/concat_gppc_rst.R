#Concatenate scan summaries into tissue-specific files
library(data.table)

tissue_file <- "data/gtex_v8_tissue.txt"

out_dir <- "results/ukbb_gtex_multi_intact/gppc_rst/gppc_summary"

tissues <- fread(tissue_file,h=F)$V1

for (i in 1:length(tissues)){

	tissue_dir <- paste0('results/ukbb_gtex_multi_intact/gppc_rst/',tissues[i],'/')

	system(paste0('cat ',tissue_dir,"*_gppc_summary.txt > ",out_dir,"/",tissues[i],"_gppc_summary.txt"))
}

