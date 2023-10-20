#Bayes factor EM algorithm to estimate priors

library(dplyr)
library(ggplot2)
library(tidyr)
library(SQUAREM)

source('scripts/multi_intact.R')



pi1_fun <- function(z_vec,lambda = 0.5){

  p_vec <- 2*pnorm(abs(z_vec),lower.tail = FALSE)

  p_vec <- p_vec[which(p_vec != 1)]

  pi0 <- length(which(p_vec > lambda))/(length(p_vec)*(1-lambda))

  pi0_max <-  0.99

  pi1 <- 1- min(pi0_max,pi0)

  return(pi1)
}


e_filename = "sim_data/pi1_5_pi2_25_pi3_25/sim.summary.v4.rst"

p_filename = "sim_data/pi1_25_pi2_5_pi3_25/sim.summary.v4.rst"

ep_filename = "sim_data/pi1_25_pi2_25_pi3_5/sim.summary.v4.rst"




df_e <- read.table(e_filename,header=T,sep = '\t')

df_p <- read.table(p_filename,header=T,sep = '\t')

df_ep <- read.table(ep_filename,header=T,sep = '\t')


df_e$fp_coloc <- linear(pmax(df_e$glcp_protein,df_e$glcp_expr))

df_p$fp_coloc <- linear(pmax(df_p$glcp_protein,df_p$glcp_expr))

df_ep$fp_coloc <- linear(pmax(df_ep$glcp_protein,df_ep$glcp_expr))

pi0_e <- 1 - pi1_fun(z_vec = qnorm(pchisq(df_e$chisq,df = 2,lower.tail = F)/2))

pi0_p <- 1 - pi1_fun(z_vec = qnorm(pchisq(df_p$chisq,df = 2,lower.tail = F)/2))

pi0_ep <- 1 - pi1_fun(z_vec = qnorm(pchisq(df_ep$chisq,df = 2,lower.tail = F)/2))



bf_df_e <- df_e %>%
        mutate(xwas_z = qnorm(pchisq(chisq,df = 2,lower.tail = F)/2,lower.tail = F)) %>%
        mutate(BF_EP_ss = wakefield_bf_z_ln(xwas_z)) %>%
        mutate(BF_E_ss = wakefield_bf_z_ln(z_ptwas_expr)) %>%
        mutate(BF_P_ss = wakefield_bf_z_ln(z_ptwas_protein)) %>%
	mutate(BF_0 = log(1)) %>%
	pivot_longer(cols = c(BF_0,BF_E_ss,BF_P_ss,BF_EP_ss),
                     names_to = "BF_Type",
                     values_to = "BF")

pi_start_e = c(pi0_e,rep(1-pi0_e,3)/3)

pi_hat_coloc_pi0_e <- SQUAREM::squarem(par = pi_start_e,
                           fixptfn = bf.em_coloc_pi0,
                           objfn = bf.loglik_coloc,
                           control = list(tol = 1.e-08,
                                          minimize=FALSE,
                                          maxiter=50),
                           bf = bf_df_e$BF,
                           fp_coloc = df_e$fp_coloc)



pi_hats_e <- pi_hat_coloc_pi0_e$par

######

bf_df_p <- df_p %>%
        mutate(xwas_z = qnorm(pchisq(chisq,df = 2,lower.tail = F)/2,lower.tail = F)) %>%
        mutate(BF_EP_ss = wakefield_bf_z_ln(xwas_z)) %>%
        mutate(BF_E_ss = wakefield_bf_z_ln(z_ptwas_expr)) %>%
        mutate(BF_P_ss = wakefield_bf_z_ln(z_ptwas_protein)) %>%
        mutate(BF_0 = log(1)) %>%
        pivot_longer(cols = c(BF_0,BF_E_ss,BF_P_ss,BF_EP_ss),
                     names_to = "BF_Type",
                     values_to = "BF")

pi_start_p = c(pi0_p,rep(1-pi0_p,3)/3)

pi_hat_coloc_pi0_p <- SQUAREM::squarem(par = pi_start_p,
                           fixptfn = bf.em_coloc_pi0,
                           objfn = bf.loglik_coloc,
                           control = list(tol = 1.e-08,
                                          minimize=FALSE,
                                          maxiter=50),
                           bf = bf_df_p$BF,
                           fp_coloc = df_p$fp_coloc)



pi_hats_p <- pi_hat_coloc_pi0_p$par


#######

bf_df_ep <- df_ep %>%
        mutate(xwas_z = qnorm(pchisq(chisq,df = 2,lower.tail = F)/2,lower.tail = F)) %>%
        mutate(BF_EP_ss = wakefield_bf_z_ln(xwas_z)) %>%
        mutate(BF_E_ss = wakefield_bf_z_ln(z_ptwas_expr)) %>%
        mutate(BF_P_ss = wakefield_bf_z_ln(z_ptwas_protein)) %>%
        mutate(BF_0 = log(1)) %>%
        pivot_longer(cols = c(BF_0,BF_E_ss,BF_P_ss,BF_EP_ss),
                     names_to = "BF_Type",
                     values_to = "BF")

pi_start_ep = c(pi0_ep,rep(1-pi0_ep,3)/3)

pi_hat_coloc_pi0_ep <- SQUAREM::squarem(par = pi_start_ep,
                           fixptfn = bf.em_coloc_pi0,
                           objfn = bf.loglik_coloc,
                           control = list(tol = 1.e-08,
                                          minimize=FALSE,
                                          maxiter=50),
                           bf = bf_df_ep$BF,
                           fp_coloc = df_ep$fp_coloc)



pi_hats_ep <- pi_hat_coloc_pi0_ep$par




#Posterior distributions


posterior_plots <- function(df,long_df,pis,prop_prefix){
	rst_df <- df %>%
        mutate(xwas_z = qnorm(pchisq(chisq,df = 2,lower.tail = F)/2,lower.tail = F)) %>%
        mutate(BF_EP_ss = wakefield_bf_z_ln(xwas_z)) %>%
        mutate(BF_E_ss = wakefield_bf_z_ln(z_ptwas_expr)) %>%
        mutate(BF_P_ss = wakefield_bf_z_ln(z_ptwas_protein)) %>%
        mutate(BF_0 = log(1))

	rst_df$posterior_EP <- class_posteriors(w = pis, bf = long_df$BF, fp_coloc = df$fp_coloc)[,4]
	rst_df$posterior_P <- class_posteriors(w = pis, bf = long_df$BF, fp_coloc = df$fp_coloc)[,3]
	rst_df$posterior_E <- class_posteriors(w = pis, bf = long_df$BF, fp_coloc = df$fp_coloc)[,2]
	rst_df$posterior_0 <- class_posteriors(w = pis, bf = long_df$BF, fp_coloc = df$fp_coloc)[,1]
	grDevices::cairo_pdf(paste0('sim_rst/',prop_prefix,'_posterior_distributions.pdf'))
	print(rst_df %>%
		filter(gamma != 0 | delta != 0) %>%
		mutate(class_num = case_when(gamma != 0 & delta != 0 ~ 'class1',
					     gamma == 0 & delta != 0 ~ 'class2',
					     gamma != 0 & delta == 0 ~ 'class3'))%>%
		pivot_longer(cols = c(posterior_0,
                              posterior_E,
                              posterior_P,
                              posterior_EP),
                     names_to = "Posterior_Type",
                     values_to = "Posterior") %>%
        filter(Posterior_Type != "posterior_0") %>%
        mutate(Posterior_Type = factor(Posterior_Type,
                                levels = c("posterior_E","posterior_P","posterior_EP"),
                                labels = c("Expression","Protein","Expression +\nProtein"))) %>%
        mutate(class_num = factor(class_num,
                                  levels = c("class1","class2","class3"),
                                  labels = c("(E,P) \u2192 Y",
                                             "P \u2192 Y only",
                                             "E \u2192 Y only"))) %>%
        ggplot(aes(x = Posterior_Type,y = Posterior,fill = Posterior_Type)) +
        geom_violin(trim = T,scale = 'width') +
        #geom_boxplot() +
        stat_summary(fun=mean, geom="crossbar",width = 0.7,colour = "red") +
        #coord_cartesian(ylim = c(0,1)) +
        facet_grid(~class_num) +
        xlab("Posterior Type") +
        ylab(expression("Posterior"))+
	scale_x_discrete(labels = c('Expression' = expression(paste(italic("M"["E"])," (Expression)")),
                                      'Protein'   = expression(paste(italic("M"["P"])," (Protein)")),
                                      "Expression +\nProtein" = expression(paste(italic("M"["E+P"])," (Expression + Protein)")))) +
        theme_bw() +
        theme(text = element_text(size = 13,face="bold"),axis.text.x = element_text(angle=30, vjust=.8, hjust=0.8),legend.position="none",aspect.ratio = 1))
	dev.off()
}


posterior_plots(df = df_e,long_df = bf_df_e,pis = pi_hats_e,prop_prefix = "pi1_5_pi2_25_pi3_25")
posterior_plots(df = df_p,long_df = bf_df_p,pis = pi_hats_p,prop_prefix = "pi1_25_pi2_5_pi3_25")
posterior_plots(df = df_ep,long_df = bf_df_ep,pis = pi_hats_ep,prop_prefix = "pi1_25_pi2_25_pi3_5")

