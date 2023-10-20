#Generate posterior model probability distribution visusalization for simulated data

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

#function to compute model posteriors for causal genes

model_posteriors <- function(file_num){

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
                filter(gamma != 0 | delta != 0) %>%
		mutate(class_num = case_when(gamma != 0 & delta != 0 ~ "class1",
					     gamma == 0 & delta != 0 ~ "class2",
					     gamma != 0 & delta == 0 ~ "class3")) %>%
                dplyr::select(gene,posterior_P,posterior_E,posterior_EP,class_num)
        return(out_df)
}

filenames <-c(paste0("sim_data/pi1_5_pi2_25_pi3_25/sim.phi_0_6.",seq(1,33),".summary"),
              paste0("sim_data/pi1_25_pi2_5_pi3_25/sim.phi_0_6.",seq(1,33),".summary"),
              paste0("sim_data/pi1_25_pi2_25_pi3_5/sim.phi_0_6.",seq(1,34),".summary"))


rst_df <- NULL

for (i in 1:length(filenames)){

        tmp <- model_posteriors(i)

        rst_df <- rbind.data.frame(rst_df,tmp)

}



out_file <- "sim_rst/posterior_ss_violin_qval_vary_priors.pdf"


grDevices::cairo_pdf(out_file)
rst_df %>%
	pivot_longer(cols = c(posterior_E,
                              posterior_P,
                              posterior_EP),
                     names_to = "Posterior_Type",
                     values_to = "Posterior") %>%
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
	#coord_flip() +
        facet_grid(~class_num) +
        xlab("Posterior Type") +
        ylab(expression("Posterior"))+
	scale_x_discrete(labels = c('Expression' = expression(paste(italic("M"["E"])," (Expression)")),
                                      'Protein'   = expression(paste(italic("M"["P"])," (Protein)")),
				      "Expression +\nProtein" = expression(paste(italic("M"["E+P"])," (Expression + Protein)")))) +
        #ggtitle(~ paste("Posterior distributions with ", pi [0], " estimated in q-value procedure")) +
        theme_bw() +
        theme(text = element_text(size = 13,face="bold"),axis.text.x = element_text(angle=45, vjust=1, hjust=1),legend.position="none",aspect.ratio = 1)
dev.off()



