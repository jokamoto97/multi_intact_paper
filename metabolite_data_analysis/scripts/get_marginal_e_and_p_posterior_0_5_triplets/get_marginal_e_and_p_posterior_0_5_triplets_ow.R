args  <- commandArgs(trailingOnly=T);


suppressPackageStartupMessages(require(data.table));
suppressPackageStartupMessages(require(stringr));
suppressPackageStartupMessages(require(dplyr));
suppressPackageStartupMessages(require(tidyr));

source('scripts/multi_intact.R')


metab_ids <- read.table('data/metab_names.txt',h=F)

metabs <- metab_ids$V1

#cid <- metabs[position]
cid <- args[1]

tissue_file <- "data/gtex_v8_tissue.txt"

tissues <- fread(tissue_file,h=F)$V1


for (i in 1:length(tissues)){


if(!dir.exists(paste0("results/ukbb_gtex_multi_intact/em_algorithm_rst/high_marginal_e_and_p_posterior_triplets"))){
        dir.create(paste0("results/ukbb_gtex_multi_intact/em_algorithm_rst/high_marginal_e_and_p_posterior_triplets"));
};

infile <- paste0("results/ukbb_gtex_multi_intact/em_algorithm_rst/", tissues[i],"/",cid,"_",tissues[i],"_posterior_rst.txt")

outfile <- paste0("results/ukbb_gtex_multi_intact/em_algorithm_rst/high_marginal_e_and_p_posterior_triplets/",cid,"_",tissues[i],"_high_marginal_e_and_p_posterior_triplets.txt")

current_files <- paste0("results/ukbb_gtex_multi_intact/em_algorithm_rst/high_marginal_e_and_p_posterior_triplets/",list.files("results/ukbb_gtex_multi_intact/em_algorithm_rst/high_marginal_e_and_p_posterior_triplets/"))

if(outfile %in% current_files){
       next
}

dat <- fread(infile,h=T)

out <- dat %>% 
	mutate(marg_posterior_E = posterior_E + posterior_EP,
	       marg_posterior_P = posterior_P + posterior_EP) %>%
	filter(marg_posterior_E >= 0.5 & marg_posterior_P >= 0.5)

if (nrow(out) > 0){

        out$Tissue <- tissues[i]

        out$cid <- cid

        write.table(out,file = outfile,sep = '\t',quote = F, row.names = F, col.names = F)

}
if (nrow(out) == 0){
        next
}

}

