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

if(!dir.exists(paste0("results/ukbb_gtex_multi_intact/gppc_rst/", tissues[i]))){
        dir.create(paste0("results/ukbb_gtex_multi_intact/gppc_rst/", tissues[i]));
};

infile <- paste0("data/multi_intact_input/",tissues[i],"/",cid,"_",tissues[i],"_multi_intact_input.txt")

out_file1 <- paste0("results/ukbb_gtex_multi_intact/gppc_rst/", tissues[i],"/",cid,"_",tissues[i],"_gppc_rst.txt")

out_file2 <- paste0("results/ukbb_gtex_multi_intact/gppc_rst/", tissues[i],"/",cid,"_",tissues[i],"_gppc_summary.txt")

current_files <- paste0("results/ukbb_gtex_multi_intact/gppc_rst/", tissues[i],"/",list.files(paste0("results/ukbb_gtex_multi_intact/gppc_rst/", tissues[i])))

#if(out_file1 %in% current_files & out_file2 %in% current_files){
#       next
#}

dat <- fread(infile,h=T)

##Run Multi-INTACT scan

coloc_dat <- dat %>% dplyr::select(ensid_protein_id_chr,GLCP_expr,GLCP_protein)

rst <- multi_intact_scan(coloc_rst = coloc_dat,
                         chisq_vec = dat$chisq_stat,
                         chisq_vec_dof = dat$chisq_df)


rst$FDR05_sig <- fdr_rst2(posterior = rst$posterior)$sig

rst <- rst %>% dplyr::select(ensid_protein_id_chr,posterior,FDR05_sig) %>%
        dplyr::rename("scan_posterior" = "posterior",
                      "FDR05_sig_scan" = "FDR05_sig")

dat <- merge(dat,rst,by = "ensid_protein_id_chr")

#Run marginal INTACT analyses

dat$intact_expr_posterior <- intact(GLCP_vec = dat$GLCP_expr,z_vec = dat$PTWAS_Z)

dat$FDR05_sig_intact_expr <- fdr_rst2(posterior = dat$intact_expr_posterior)$sig

dat$intact_protein_posterior <- intact(GLCP_vec = dat$GLCP_protein,z_vec = dat$pwas_z)

dat$FDR05_sig_intact_protein <- fdr_rst2(posterior = dat$intact_protein_posterior)$sig

out <- dat %>% arrange(desc(scan_posterior)) %>%
        mutate(scan_only = case_when(FDR05_sig_scan == TRUE &
                                     FDR05_sig_intact_expr == FALSE &
                                     FDR05_sig_intact_protein == FALSE ~ TRUE,
                             TRUE ~ FALSE)) %>%
        mutate(expr_only = case_when(FDR05_sig_scan == TRUE &
                                     FDR05_sig_intact_expr == TRUE &
                                     FDR05_sig_intact_protein == FALSE ~ TRUE,
                             TRUE ~ FALSE)) %>%
        mutate(protein_only = case_when(FDR05_sig_scan == TRUE &
                                     FDR05_sig_intact_expr == FALSE &
                                     FDR05_sig_intact_protein == TRUE ~ TRUE,
                             TRUE ~ FALSE)) %>%
        mutate(expr_and_protein = case_when(FDR05_sig_scan == TRUE &
                                     FDR05_sig_intact_expr == TRUE &
                                     FDR05_sig_intact_protein == TRUE ~ TRUE,
                             TRUE ~ FALSE)) %>%
        mutate(consistent_marginal_z_sign = case_when(FDR05_sig_scan == TRUE & (sign(PTWAS_Z) == sign(pwas_z)) ~ TRUE,
                                                     FDR05_sig_scan == TRUE & !(sign(PTWAS_Z) == sign(pwas_z)) ~ FALSE,
                                                     FDR05_sig_scan == FALSE ~ NA)) %>%
        dplyr::select(ensid_protein_id_chr,PTWAS_Z,pwas_z,GLCP_expr,GLCP_protein,chisq_stat,
        chisq_df,scan_posterior,intact_expr_posterior,intact_protein_posterior,FDR05_sig_scan,FDR05_sig_intact_expr,FDR05_sig_intact_protein,scan_only,expr_only,protein_only,expr_and_protein,consistent_marginal_z_sign)

#Summarize results

summary_rst <- out %>%
        summarise(n_sig_scan = sum(FDR05_sig_scan),
                  n_sig_intact_expr = sum(FDR05_sig_intact_expr),
                  n_sig_intact_protein = sum(FDR05_sig_intact_protein),
                  n_scan_only = sum(scan_only),
                  n_expr_only = sum(expr_only),
                  n_protein_only = sum(protein_only),
                  n_expr_and_protein = sum(expr_and_protein),
                  n_consistent_marginal_z_sign = sum(consistent_marginal_z_sign,na.rm=T)) %>%
        mutate(n_inconsistent_marginal_z_sign = n_sig_scan - n_consistent_marginal_z_sign)

summary_rst <- summary_rst %>% mutate(Tissue = tissues[i],CID = cid) %>%
        dplyr::select(CID,Tissue,n_sig_scan,n_sig_intact_expr,n_sig_intact_protein,n_scan_only,n_expr_only,n_protein_only,n_expr_and_protein,n_consistent_marginal_z_sign,n_inconsistent_marginal_z_sign)


write.table(out,file=out_file1,quote = F, sep = '\t',row.names = F, col.names = T)

write.table(summary_rst,file=out_file2,quote = F, sep = '\t',row.names = F, col.names = F)

}

