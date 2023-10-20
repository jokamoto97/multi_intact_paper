#Plot agreement between gene-tissue-metabolite triplets that are identified by Multi-INTACT or KBA.

suppressPackageStartupMessages(require(data.table));
suppressPackageStartupMessages(require(stringr));
suppressPackageStartupMessages(require(dplyr));
suppressPackageStartupMessages(require(tidyr));
suppressPackageStartupMessages(require(ggplot2));
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(eulerr))



dir <- "results/ukbb_gtex_multi_intact/gppc_rst/sig_rst/"


#Limit to genes tested by both TWAS in at least 1 tissue and PWAS

intersect_symbols_file <- "data/genes_tested_in_pwas_and_twas.txt"

intersect_symbols <- read.table(intersect_symbols_file,h=F,sep='\t')$V1

#Tested metabolites

metab_ids <- read.table('data/metab_names.txt',h=F)$V1

#Multi-INTACT results

m_intact_rst <- fread(cmd = paste0("cat ",dir,"*_sig_rst.txt"),h=F)
colnames(m_intact_rst) <- c("ensid_protein_id_chr",
                           "PTWAS_Z",
                           "pwas_z",
                           "GLCP_expr",
                           "GLCP_protein",
                           "chisq_stat",
                           "chisq_df",
                           "scan_posterior",
                           "intact_expr_posterior",
                           "intact_protein_posterior",
                           "FDR05_sig_scan",
                           "FDR05_sig_intact_expr",
                           "FDR05_sig_intact_protein",
                           "scan_only",
                           "expr_only",
                           "protein_only",
                           "expr_and_protein",
                           "consistent_marginal_z_sign",
                           "Tissue",
                           "CID")

m_intact_rst$ENSID <- str_split_fixed(m_intact_rst$ensid_protein_id_chr,":",n = 3)[,1]
m_intact_rst$protein_id <- str_split_fixed(m_intact_rst$ensid_protein_id_chr,":",n = 3)[,2]
m_intact_rst$chrom <- str_split_fixed(m_intact_rst$ensid_protein_id_chr,":",n = 3)[,3]
m_intact_rst$symbol <- gsub("\\_.*", "", m_intact_rst$protein_id)
m_intact_rst$cid_gene <- paste0(m_intact_rst$CID,"_",m_intact_rst$symbol)
m_intact_rst <- m_intact_rst %>% filter(symbol %in% intersect_symbols) #%>%
#	filter(sign(PTWAS_Z) == sign(pwas_z))

#KBA data

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
        mutate(cid_gene = paste0(CID,"_",gene)) %>%
	filter(gene %in% intersect_symbols) %>%
	filter(CID %in% metab_ids)





eric_pairs <- unique(eric_dat$cid_gene)

m_intact_pairs<- unique(m_intact_rst$cid_gene)

eric_only <- length(setdiff(eric_pairs,m_intact_pairs))

m_intact_only <- length(setdiff(m_intact_pairs,eric_pairs))

eric_m_intact <- length(intersect(eric_pairs, m_intact_pairs))


pdf('results/venn_diagrams/m_intact_kba_gene_metabolite_venn_diagram.pdf')
s2 <- c("Knowledge\nBased" = eric_only,
        "Multi-INTACT" = m_intact_only,
        "Knowledge\nBased&Multi-INTACT" = eric_m_intact)
plot(euler(s2),quantities = list(type = "counts",cex=2),fills = list(fill = c("grey", "purple"), alpha = 0.5),
     labels = list(col = "black",cex = 1.25))
dev.off()


#intersect_eric_m_intact <- intersect(eric_pairs, m_intact_pairs)
