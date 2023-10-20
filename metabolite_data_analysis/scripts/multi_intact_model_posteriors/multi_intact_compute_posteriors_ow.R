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

tissue_priors_file <- "data/tissue_priors.txt"

tissue_priors <- fread(tissue_priors_file,h=F)

colnames(tissue_priors) <- c("Tissue","pi_0","pi_E","pi_P","pi_EP","converged")

for (i in 1:length(tissues)){

if(!dir.exists(paste0("results/ukbb_gtex_multi_intact/em_algorithm_rst/", tissues[i]))){
        dir.create(paste0("results/ukbb_gtex_multi_intact/em_algorithm_rst/", tissues[i]));
};

infile <- paste0("data/multi_intact_input/",tissues[i],"/",cid,"_",tissues[i],"_multi_intact_input.txt")

out_file <- paste0("results/ukbb_gtex_multi_intact/em_algorithm_rst/", tissues[i],"/",cid,"_",tissues[i],"_posterior_rst.txt")

current_files <- paste0("results/ukbb_gtex_multi_intact/em_algorithm_rst/", tissues[i],"/",list.files(paste0("results/ukbb_gtex_multi_intact/em_algorithm_rst/", tissues[i])))

if(out_file %in% current_files){
       next
}

dat <- fread(infile,h=T)
dat$fp_coloc <- linear(pmax(dat$GLCP_protein,dat$GLCP_expr))
pi0 <- 1 - .pi1_fun_multi(chisq_vec = dat$chisq_stat,df = dat$chisq_df)

bf_df <- dat %>%
        mutate(xwas_z = qnorm(pchisq(chisq_stat,df = chisq_df,lower.tail = F,log.p = T) - log(2),lower.tail = F,log.p = T)) %>%
        mutate(BF_EP_ss = wakefield_bf_z_ln(xwas_z)) %>%
        mutate(BF_E_ss = wakefield_bf_z_ln(PTWAS_Z)) %>%
        mutate(BF_P_ss = wakefield_bf_z_ln(pwas_z)) %>%
        mutate(BF_0 = log(1)) %>%
        pivot_longer(cols = c(BF_0,BF_E_ss,BF_P_ss,BF_EP_ss),
                     names_to = "BF_Type",
                     values_to = "BF")



em_rst <- tissue_priors[tissue_priors$Tissue == tissues[i]]

#Use estimates if EM converged; else, use (1/3,1/3,1/3) for alternative scenarios

if (em_rst$converged == TRUE){

        pi_hats <- as.vector(unlist(em_rst[1,2:5]))

}
if (em_rst$converged == FALSE){

        pi_0 <- as.numeric(em_rst[1,2])

        pi_hats <- c(pi_0,rep(1-pi_0,3)/3)
}



dat %>%
        mutate(posterior_EP = class_posteriors(w = pi_hats, bf = bf_df$BF, fp_coloc = fp_coloc)[,4]) %>%
        mutate(posterior_P = class_posteriors(w = pi_hats, bf = bf_df$BF, fp_coloc = fp_coloc)[,3]) %>%
        mutate(posterior_E = class_posteriors(w = pi_hats, bf = bf_df$BF, fp_coloc = fp_coloc)[,2]) %>%
        mutate(posterior_0 = class_posteriors(w = pi_hats, bf = bf_df$BF, fp_coloc = fp_coloc)[,1]) %>%
        write.table(file = out_file,quote = F, sep = '\t', row.names = F, col.names =T)

}

