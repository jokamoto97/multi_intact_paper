args  <- commandArgs(trailingOnly=T);


suppressPackageStartupMessages(require(data.table));
suppressPackageStartupMessages(require(stringr));
suppressPackageStartupMessages(require(dplyr));
suppressPackageStartupMessages(require(tidyr));

source('scripts/multi_intact.R')


if(!dir.exists(paste0("results/ukbb_gtex_multi_intact/twas_pwas_sign_concordance"))){
        dir.create(paste0("results/ukbb_gtex_multi_intact/twas_pwas_sign_concordance"));
};


metab_ids <- read.table('data/metab_names.txt',h=F)

metabs <- metab_ids$V1

#cid <- metabs[position]
cid <- args[1]

tissue_file <- "data/gtex_v8_tissue.txt"

tissues <- fread(tissue_file,h=F)$V1

outfile <- paste0("results/ukbb_gtex_multi_intact/twas_pwas_sign_concordance/",cid,"_twas_pwas_sign_log10p_spearman_cor.txt")

current_files <- paste0("results/ukbb_gtex_multi_intact/twas_pwas_sign_concordance/",list.files("results/ukbb_gtex_multi_intact/twas_pwas_sign_concordance/"))


if(!(outfile %in% current_files)){

outdf <- NULL

for (i in 1:length(tissues)){


infile <- paste0("data/multi_intact_input/", tissues[i],"/",cid,"_",tissues[i],"_multi_intact_input.txt")


current_files <- paste0("results/ukbb_gtex_multi_intact/twas_pwas_sign_concordance/",list.files("results/ukbb_gtex_multi_intact/twas_pwas_sign_concordance/"))

if(outfile %in% current_files){
       next
}

dat <- fread(infile,h=T)

tmp <- dat %>%
	mutate(twas_p = 2*pnorm(-abs(PTWAS_Z)),
               pwas_p = 2*pnorm(-abs(pwas_z))) %>%
        mutate(log10_twasp = -log10(twas_p),
               log10_pwasp = -log10(pwas_p)) %>%
        mutate(sign_log10twasp = sign(PTWAS_Z) * log10_twasp,
               sign_log10pwasp = sign(pwas_z) * log10_pwasp) %>%
        summarise(Estimate = cor.test(sign_log10twasp,sign_log10pwasp,method = "spearman")$estimate,
                  pval = cor.test(sign_log10twasp,sign_log10pwasp,method = "spearman")$p.value,
		  statistic = cor.test(sign_log10twasp,sign_log10pwasp,method = "spearman")$statistic) %>%
	mutate(Tissue = tissues[i],CID = cid) %>%
	dplyr::select(CID,Tissue,Estimate,pval,statistic) 

outdf <- rbind.data.frame(outdf,tmp)
}
write.table(outdf,file = outfile,row.names=F, col.names = F, quote = F, sep = '\t')
}
