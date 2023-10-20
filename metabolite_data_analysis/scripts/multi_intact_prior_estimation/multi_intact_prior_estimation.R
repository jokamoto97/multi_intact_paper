#!/usr/bin/env Rscript
pos=as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))


suppressPackageStartupMessages(require(data.table));
suppressPackageStartupMessages(require(stringr));
suppressPackageStartupMessages(require(dplyr));
suppressPackageStartupMessages(require(tidyr));

source('scripts/multi_intact.R')

tissue_file <- "data/gtex_v8_tissue.txt"

tissues <- fread(tissue_file,h=F)$V1

tissue_dir <- paste0("data/multi_intact_input/",tissues[pos],"/")


if(!dir.exists("results/ukbb_gtex_multi_intact/em_algorithm_rst/prior_estimates")){
        dir.create("results/ukbb_gtex_multi_intact/em_algorithm_rst/prior_estimates");
};

out_file <- paste0("results/ukbb_gtex_multi_intact/em_algorithm_rst/prior_estimates/",tissues[pos],"_em_algorithm_priors.txt")


current_files <- paste0("results/ukbb_gtex_multi_intact/em_algorithm_rst/prior_estimates/",list.files(paste0("results/ukbb_gtex_multi_intact/em_algorithm_rst/prior_estimates/")))

if(!(out_file %in% current_files)){
	
	dat <- fread(cmd = paste0("awk FNR!=1 ",tissue_dir,"*_multi_intact_input.txt"),h=F)
	colnames(dat) <- c("ensid_protein_id_chr",
				  "PTWAS_Z",
				  "pwas_z",
				  "GLCP_expr",
				  "GLCP_protein",
				  "chisq_stat",
				  "chisq_df")


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

	pi_start = c(pi0,rep(1-pi0,3)/3)


	pi_hat_coloc_pi0 <- SQUAREM::squarem(par = pi_start,
        	                   fixptfn = bf.em_coloc_pi0,
                	           objfn = bf.loglik_coloc,
                        	   control = list(tol = 1.e-08,
                                	          minimize=FALSE,
                                        	  maxiter=50),
                           	bf = bf_df$BF,
                           	fp_coloc = dat$fp_coloc)

	pi_hats <- pi_hat_coloc_pi0$par
	
	out <- data.frame(Tissue = tissues[pos],
			  pi_0 = pi_hats[1],
			  pi_E = pi_hats[2],
			  pi_P = pi_hats[3],
			  pi_EP = pi_hats[4],
			  pi_hat_coloc_pi0$convergence)

	write.table(out,file = out_file,quote = F, sep = '\t',row.names = F, col.names = F)
}

