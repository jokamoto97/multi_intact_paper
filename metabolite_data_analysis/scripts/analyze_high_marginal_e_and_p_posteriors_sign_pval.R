#Compute correlation of signed -log10 TWAS and PWAS pvalues

suppressPackageStartupMessages(require(data.table));
suppressPackageStartupMessages(require(dplyr));
suppressPackageStartupMessages(require(ggplot2));
suppressPackageStartupMessages(require(tidyr));
suppressPackageStartupMessages(require(tidyverse));



dir <- "results/ukbb_gtex_multi_intact/em_algorithm_rst/high_marginal_e_and_p_posterior_triplets/"

#outfile1 <- "results/ukbb_gtex_multi_intact/sign_pval_cor.pdf"

#Metabolite annotations

#annot <- fread("data/metsim_merge_metabolite_annotation_kegg_hmdb.txt",h=T) %>% filter(SUPER_PATHWAY == "Lipid")


dat <- fread(cmd = paste0("cat ",dir,"*_high_marginal_e_and_p_posterior_triplets.txt"),h=F)
colnames(dat) <- c("ensid_protein_id_chr",
                   "PTWAS_Z",
                   "pwas_z",
                   "GLCP_expr",
                   "GLCP_protein",
                   "chisq_stat",
                   "chisq_df",
                   "fp_coloc",
                   "posterior_EP",
                   "posterior_P",
                   "posterior_E",
                   "posterior_0",
                   "marg_posterior_E",
                   "marg_posterior_P",
                   "Tissue",
                   "CID")


#Marginal across tissues

dat %>%
        mutate(twas_p = 2*pnorm(-abs(PTWAS_Z)),
               pwas_p = 2*pnorm(-abs(pwas_z))) %>%
        mutate(log10_twasp = -log10(twas_p),
               log10_pwasp = -log10(pwas_p)) %>%
        mutate(sign_log10twasp = sign(PTWAS_Z) * log10_twasp,
               sign_log10pwasp = sign(pwas_z) * log10_pwasp) %>%
        summarise(estimate = cor.test(sign_log10twasp,sign_log10pwasp,method = "spearman")$estimate,
                  pval = cor.test(sign_log10twasp,sign_log10pwasp,method = "spearman")$p.value) 

	




