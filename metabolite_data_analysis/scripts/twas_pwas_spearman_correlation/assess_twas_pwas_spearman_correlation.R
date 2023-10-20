#Analyze metabolite-tissue signed log pval correlation coefficients

suppressPackageStartupMessages(require(data.table));
suppressPackageStartupMessages(require(stringr));
suppressPackageStartupMessages(require(dplyr));
suppressPackageStartupMessages(require(tidyr));
suppressPackageStartupMessages(require(tidyverse));
suppressPackageStartupMessages(require(ggplot2));

dir <- "results/ukbb_gtex_multi_intact/twas_pwas_sign_concordance/"

dat <- fread(cmd = paste0("cat ",dir,"*_twas_pwas_sign_log10p_spearman_cor.txt"),h=F)

colnames(dat) <- c("CID","Tissue","correlation","pval","statistic")


#Set tissue order
tissue_order <- dat %>%
	group_by(Tissue) %>%
	summarise(medians = median(correlation)) %>%
	arrange(desc(medians)) %>%
	dplyr::select(Tissue) %>% unlist() %>% as.vector()


#Mean and SE across all metab-tissue pairs

print(paste0("Mean correlation = ",mean(dat$correlation)))
print(paste0("SE of the mean = ",sd(dat$correlation)/sqrt(nrow(dat))))
print(paste0("Variance of the correlation = ",sd(dat$correlation)^2))


