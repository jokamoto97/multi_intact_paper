#Analyze INTACT results

library(data.table)
library(ggplot2)
library(INTACT)
library(dplyr)
library(qvalue)
library(stringr)
library(tidyr)
#library(venneuler)
library(eulerr)



fdr_rst2 <- function (posterior, alpha = 0.05)
{
    gene_num <- seq(1, length(posterior))
    lfdr_sort = sort(1 - posterior)
    FDR = cumsum(lfdr_sort)/(1:length(lfdr_sort))
    thresh = 1 - lfdr_sort[max(which(FDR <= alpha))]
    rej_gene = as.numeric(gene_num[which(posterior >= thresh)])
    out_tmp <- rep(FALSE, length(posterior))
    out_tmp[rej_gene] <- TRUE
    out <- data.frame(posterior = posterior, sig = out_tmp)
    return(out)
}


#Limit to genes tested by both TWAS in at least 1 tissue and PWAS

intersect_symbols_file <- "data/genes_tested_in_pwas_and_twas.txt"

intersect_symbols <- read.table(intersect_symbols_file,h=F,sep='\t')$V1

#Tested metabolites

metab_ids <- fread('data/metab_names.txt',h=F)


#Protein INTACT results

intact_rst_file <- 'data/intact_sig_genes_prelim.txt'

intact_rst <- fread(intact_rst_file) 

intact_rst$gene_name <- gsub("\\_.*", "", intact_rst$protein_id)

intact_rst$cid_gene <- paste0(intact_rst$metab_id,"_",intact_rst$gene_name)

intact_rst <- intact_rst %>% 
	filter(gene_name %in% intersect_symbols & metab_id %in% metab_ids$V1)

#KBA

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
	filter(gene %in% intersect_symbols & CID %in% metab_ids$V1) 

merged <- merge(eric_dat,intact_rst, by = 'cid_gene') %>%
	dplyr::select(cid_gene,CHROM,START,STOP,protein_cis_region) %>%
	unique() 


eqtl_intact_file <- 'data/tab.intact.dap.cid.pcg.sig.gz'

eqtl_intact <- fread(eqtl_intact_file) %>%
	filter(SYMBOL %in% intersect_symbols & CID %in% metab_ids$V1)

eqtl_intact$cid_gene <- paste0(eqtl_intact$CID,'_',eqtl_intact$SYMBOL)

#Make Venn Diagram

pqtl_pairs <- unique(intact_rst$cid_gene)

eqtl_pairs <- unique(eqtl_intact$cid_gene)

kba_pairs <- unique(eric_dat$cid_gene)

pqtl_eqtl <- intersect(pqtl_pairs, eqtl_pairs)

pqtl_kba <- intersect(pqtl_pairs, kba_pairs)

eqtl_kba <- intersect(eqtl_pairs, kba_pairs)

##Diagram entries

pqtl_only <- length(setdiff(pqtl_pairs, union(eqtl_pairs,kba_pairs)))

eqtl_only <- length(setdiff(eqtl_pairs, union(pqtl_pairs,kba_pairs)))

kba_only <- length(setdiff(kba_pairs, union(pqtl_pairs,eqtl_pairs)))

pqtl_eqtl_kba <- length(intersect(pqtl_eqtl, kba_pairs))

pqtl_and_eqtl <- length(setdiff(pqtl_eqtl,intersect(pqtl_eqtl,kba_pairs)))

pqtl_and_kba <- length(setdiff(pqtl_kba,intersect(pqtl_eqtl,kba_pairs)))

eqtl_and_kba <- length(setdiff(eqtl_kba,intersect(pqtl_eqtl,kba_pairs)))

length(union(kba_pairs,union(pqtl_pairs,eqtl_pairs)))

sum(c(pqtl_only,eqtl_only,kba_only,pqtl_and_eqtl,pqtl_and_kba,eqtl_and_kba,pqtl_eqtl_kba))



#KBA - expression

pdf('results/venn_diagrams/sig_genes_venn_diagram_kba_expr.pdf')
s2 <- c("Knowledge\nBased" = length(setdiff(kba_pairs,eqtl_pairs)),
	"INTACT\n(Expression)" = length(setdiff(eqtl_pairs,kba_pairs)),
        "INTACT\n(Expression)&Knowledge\nBased" = length(intersect(kba_pairs,eqtl_pairs)))

plot(euler(s2),quantities = list(type = "counts",cex = 2),fills = list(fill = rev(c("steelblue4","grey")), alpha = 0.5),
     labels = list(col = "black",cex = 1.25))
dev.off()


#KBA - protein

pdf('results/venn_diagrams/sig_genes_venn_diagram_kba_protein.pdf')
s2 <- c("Knowledge\nBased" = length(setdiff(kba_pairs,pqtl_pairs)),
	"INTACT\n(Protein)" = length(setdiff(pqtl_pairs,kba_pairs)),
        "INTACT\n(Protein)&Knowledge\nBased" = length(intersect(kba_pairs,pqtl_pairs)))

plot(euler(s2),quantities = list(type = "counts",cex = 2),fills = list(fill = rev(c("red","grey")), alpha = 0.5),
     labels = list(col = "black",cex = 1.25))
dev.off()


#union of pqtl_kba and eqtl_kba

#pqtl_kba_union_eqtl_kba <- union(pqtl_kba,eqtl_kba)
