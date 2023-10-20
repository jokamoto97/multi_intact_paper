args  <- commandArgs(trailingOnly=T);


suppressPackageStartupMessages(require(data.table));
suppressPackageStartupMessages(require(stringr));
suppressPackageStartupMessages(require(dplyr));
source('scripts/multi_intact.R')


metab_ids <- read.table('data/metab_names.txt',h=F)

metabs <- metab_ids$V1

#cid <- metabs[position]
cid <- args[1]

tissue_file <- "data/gtex_v8_tissue.txt"

tissues <- fread(tissue_file,h=F)$V1

for (i in 1:length(tissues)){

if(!dir.exists("results/ukbb_gtex_multi_intact/gppc_rst/sig_rst")){
        dir.create("results/ukbb_gtex_multi_intact/gppc_rst/sig_rst");
};


infile <- paste0("results/ukbb_gtex_multi_intact/gppc_rst/", tissues[i],"/",cid,"_",tissues[i],"_gppc_rst.txt")

outfile <- paste0("results/ukbb_gtex_multi_intact/gppc_rst/sig_rst/",cid,"_",tissues[i],"_sig_rst.txt")

current_files <- paste0("results/ukbb_gtex_multi_intact/gppc_rst/sig_rst/",list.files("results/ukbb_gtex_multi_intact/gppc_rst/sig_rst"))

if(outfile %in% current_files){
	next
}

dat <- fread(infile,h=T)

out <- dat %>% filter(FDR05_sig_scan == TRUE) %>%
	mutate(Tissue = tissues[i]) %>%
	mutate(CID = cid)


if(nrow(out) > 0){

	write.table(out,file = outfile,quote = F, sep ='\t',row.names = F, col.names = F)
}

}


