#Compare Multi-INTACT EM algorithm posteriors to INTACT posterior as a E-->Y or P-->Y effect classifier

library(dplyr)
library(ggplot2)
library(tidyr)
library(SQUAREM)
library(pROC)

source('scripts/multi_intact.R')

pi1_fun <- function(z_vec,lambda = 0.5){

  p_vec <- 2*pnorm(abs(z_vec),lower.tail = FALSE)

  p_vec <- p_vec[which(p_vec != 1)]

  pi0 <- length(which(p_vec > lambda))/(length(p_vec)*(1-lambda))

  pi0_max <-  0.99

  pi1 <- 1- min(pi0_max,pi0)

  return(pi1)
}

simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels, scores = scores[order(scores,decreasing = T)])
}



#Function to compute marginal posteriors and return ROC function input given a simulated data set


marginal_posteriors <- function(file_num){

	df <- read.table(filenames[file_num],h=T,sep = '\t')

	df$fp_coloc <- linear(pmax(df$glcp_protein,df$glcp_expr))

	pi0 <- 1 - pi1_fun(z_vec = qnorm(pchisq(df$chisq,df = 2,lower.tail = F)/2))

	bf_df <- df %>%
        	mutate(xwas_z = qnorm(pchisq(chisq,df = 2,lower.tail = F)/2,lower.tail = F)) %>%
        	mutate(BF_EP_ss = wakefield_bf_z_ln(xwas_z)) %>%
        	mutate(BF_E_ss = wakefield_bf_z_ln(z_ptwas_expr)) %>%
        	mutate(BF_P_ss = wakefield_bf_z_ln(z_ptwas_protein)) %>%
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
        	                   fp_coloc = df$fp_coloc)

	pi_hats_qvalue <- pi_hat_coloc_pi0$par


	out_df <- df %>%
        	mutate(xwas_z = qnorm(pchisq(chisq,df = 2,lower.tail = F)/2,lower.tail = F)) %>%
        	mutate(BF_EP_ss = wakefield_bf_z_ln(xwas_z)) %>%
        	mutate(BF_E_ss = wakefield_bf_z_ln(z_ptwas_expr)) %>%
       		mutate(BF_P_ss = wakefield_bf_z_ln(z_ptwas_protein)) %>%
        	mutate(BF_0 = log(1)) %>%
        	mutate(posterior_EP = class_posteriors(w = pi_hats_qvalue, bf = bf_df$BF, fp_coloc = fp_coloc)[,4]) %>%
        	mutate(posterior_P = class_posteriors(w = pi_hats_qvalue, bf = bf_df$BF, fp_coloc = fp_coloc)[,3]) %>%
        	mutate(posterior_E = class_posteriors(w = pi_hats_qvalue, bf = bf_df$BF, fp_coloc = fp_coloc)[,2]) %>%
        	mutate(posterior_0 = class_posteriors(w = pi_hats_qvalue, bf = bf_df$BF, fp_coloc = fp_coloc)[,1]) %>%
		mutate(intact_posterior_E = intact(GLCP_vec = glcp_expr,z_vec = z_ptwas_expr)) %>%
		mutate(intact_posterior_P = intact(GLCP_vec = glcp_protein,z_vec = z_ptwas_protein)) %>%
		mutate(marginal_posterior_E = posterior_EP + posterior_E) %>%
		mutate(marginal_posterior_P = posterior_EP + posterior_P) %>%
		mutate(P_truth = case_when(delta != 0 ~ T,
				   TRUE ~ F)) %>%
		mutate(E_truth = case_when(gamma != 0 ~ T,
				   TRUE ~ F)) %>%
		mutate(abs_z_expr = abs(z_ptwas_expr),
		       abs_z_protein = abs(z_ptwas_protein)) %>%
		dplyr::select(gene,P_truth,E_truth,marginal_posterior_P,marginal_posterior_E,intact_posterior_P,intact_posterior_E,glcp_expr,glcp_protein,abs_z_expr,abs_z_protein)
		
	return(out_df)
}

filenames <-c(paste0("sim_data/pi1_5_pi2_25_pi3_25/sim.phi_0_6.",seq(1,33),".summary"),
              paste0("sim_data/pi1_25_pi2_5_pi3_25/sim.phi_0_6.",seq(1,33),".summary"),
              paste0("sim_data/pi1_25_pi2_25_pi3_5/sim.phi_0_6.",seq(1,34),".summary"))



outfile <- "sim_rst/roc_m_intact_vary_priors.txt"

#Compute marginal posteriors for all simulated data

roc_df <- NULL

for (i in 1:length(filenames)){

	tmp <- marginal_posteriors(i)

	roc_df <- rbind.data.frame(roc_df,tmp)

}

write.table(roc_df,file = outfile,row.names = F, col.names = T, sep = '\t',quote = F)


