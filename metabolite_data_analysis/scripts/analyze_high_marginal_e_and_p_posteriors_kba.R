#Summarize gene-metabolite-tissue triplets with marginal posteriors E and P >= 0.5

suppressPackageStartupMessages(require(data.table));
suppressPackageStartupMessages(require(dplyr));
suppressPackageStartupMessages(require(ggplot2));
suppressPackageStartupMessages(require(tidyr));
suppressPackageStartupMessages(require(tidyverse));



dir <- "results/ukbb_gtex_multi_intact/em_algorithm_rst/high_marginal_e_and_p_posterior_triplets/"

#outfile1 <- "results/ukbb_gtex_multi_intact/high_marginal_e_and_p_posterior_triplets_kba.txt"

outfile2 <- "results/ukbb_gtex_multi_intact/marginal_ep_posterior_triplets_0_5_sign_concordance_kba.txt"



m_intact_rst <- fread(cmd = paste0("cat ",dir,"*_high_marginal_e_and_p_posterior_triplets.txt"),h=F)
colnames(m_intact_rst) <- c("ensid_protein_id_chr",
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




m_intact_rst$ENSID <- str_split_fixed(m_intact_rst$ensid_protein_id_chr,":",n = 3)[,1]
m_intact_rst$protein_id <- str_split_fixed(m_intact_rst$ensid_protein_id_chr,":",n = 3)[,2]
m_intact_rst$chrom <- str_split_fixed(m_intact_rst$ensid_protein_id_chr,":",n = 3)[,3]
m_intact_rst$symbol <- gsub("\\_.*", "", m_intact_rst$protein_id)
m_intact_rst$cid_gene <- paste0(m_intact_rst$CID,"_",m_intact_rst$symbol)

#Read in KBA data

eric_file <- 'data/tab_kba_genes_4_dap_signal_clusters.txt'
eric_dat <- fread(eric_file) %>%
        #filter(CID %in% intact_rst$metab_id) %>% 
        dplyr::select(CHROM,START,STOP,CID,eric_gene_group)

tmp <- rbindlist(
  lapply(strsplit(eric_dat$eric_gene_group, "\\W"), function(x) data.table(t(x))),
  fill = TRUE
)

eric_dat <- cbind(eric_dat,tmp) %>%
        pivot_longer(cols = paste0("V",seq(1,9)),names_to = 'tmp',values_to = 'gene') %>%
        filter(is.na(gene) == F & gene != 'unknown') %>%
        mutate(cid_gene = paste0(CID,"_",gene))


dat <- merge(m_intact_rst,eric_dat,by = 'cid_gene') %>%
        dplyr::select(cid_gene,CHROM,PTWAS_Z,pwas_z,Tissue) %>%
        unique()

#write.table(dat,file = outfile1, quote = F, sep = '\t',row.names = F, col.names = T)


#Assess sign concordance

dat %>%
        mutate(sign_match_flag = case_when(sign(PTWAS_Z) == sign(pwas_z) &
                                           PTWAS_Z != 0 & pwas_z != 0 ~ TRUE,
                                           TRUE ~ FALSE)) %>%
        mutate(sign_opposite_flag = case_when(sign(PTWAS_Z) != sign(pwas_z) &
                                              PTWAS_Z != 0 & pwas_z != 0 ~ TRUE,
                                                TRUE ~ FALSE)) %>%
        summarise(sign_matches = sum(sign_match_flag),
                  opposite_signs = sum(sign_opposite_flag))


#set tissue order
tissue_order <- dat %>%
        group_by(Tissue) %>%
        summarise(n_pairs_scan = n()) %>%
        arrange(desc(n_pairs_scan)) %>%
        dplyr::select(Tissue) %>% unlist() %>% as.vector()



#save plot data
dat %>%
	mutate(sign_match_flag = case_when(sign(PTWAS_Z) == sign(pwas_z) &
                                           PTWAS_Z != 0 & pwas_z != 0 ~ TRUE,
                                           TRUE ~ FALSE)) %>%
        mutate(sign_opposite_flag = case_when(sign(PTWAS_Z) != sign(pwas_z) &
                                              PTWAS_Z != 0 & pwas_z != 0 ~ TRUE,
                                                TRUE ~ FALSE)) %>%
        group_by(Tissue) %>%
        summarise(n_pairs_conc = sum(sign_match_flag),
                  n_pairs_disc = sum(sign_opposite_flag)) %>%
        dplyr::rename("Concordant_marginal_z_signs" = "n_pairs_conc",
                      "Discordant_marginal_z_signs"="n_pairs_disc") %>%
        data.frame() %>%
        pivot_longer(cols = Concordant_marginal_z_signs:Discordant_marginal_z_signs,
                        names_to = "Type",
                        values_to = "N_pairs") %>%
	write.table(file = outfile2, quote = F, sep = '\t',row.names = F, col.names = T)


#Test whether proportion of sign concordance in liver is significantly different from the rest of the tissue "population"

liver_df <- dat %>%
        group_by(Tissue) %>%
        mutate(sign_match_flag = case_when(sign(PTWAS_Z) == sign(pwas_z) &
                                           PTWAS_Z != 0 & pwas_z != 0 ~ TRUE,
                                           TRUE ~ FALSE)) %>%
        mutate(sign_opposite_flag = case_when(sign(PTWAS_Z) != sign(pwas_z) &
                                              PTWAS_Z != 0 & pwas_z != 0 ~ TRUE,
                                                TRUE ~ FALSE)) %>%
        summarise(n_pairs_conc = sum(sign_match_flag),
                  n_pairs_disc = sum(sign_opposite_flag)) %>%
        dplyr::rename("Concordant_marginal_z_signs" = "n_pairs_conc",
                      "Discordant_marginal_z_signs"="n_pairs_disc") %>%
        data.frame() %>%
       filter(Tissue == "Liver")

non_liver_df <- dat %>%
        group_by(Tissue) %>%
        mutate(sign_match_flag = case_when(sign(PTWAS_Z) == sign(pwas_z) &
                                           PTWAS_Z != 0 & pwas_z != 0 ~ TRUE,
                                           TRUE ~ FALSE)) %>%
        mutate(sign_opposite_flag = case_when(sign(PTWAS_Z) != sign(pwas_z) &
                                              PTWAS_Z != 0 & pwas_z != 0 ~ TRUE,
                                                TRUE ~ FALSE)) %>%
        summarise(n_pairs_conc = sum(sign_match_flag),
                  n_pairs_disc = sum(sign_opposite_flag)) %>%
        dplyr::rename("Concordant_marginal_z_signs" = "n_pairs_conc",
                      "Discordant_marginal_z_signs"="n_pairs_disc") %>%
        data.frame() %>%
       filter(Tissue != "Liver") %>%
        summarise(total_n_conc = sum(Concordant_marginal_z_signs),
                        total_n_disc = sum(Discordant_marginal_z_signs))



prop.test(c(liver_df$Concordant_marginal_z_signs,non_liver_df$total_n_conc),
                c(liver_df$Concordant_marginal_z_signs + liver_df$Discordant_marginal_z_signs,
                        non_liver_df$total_n_conc + non_liver_df$total_n_disc))



