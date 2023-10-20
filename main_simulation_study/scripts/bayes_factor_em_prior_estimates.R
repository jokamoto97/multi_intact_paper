#Compute mean (SE) prior estimates for each model across the simulated data

library(dplyr)
library(ggplot2)
library(tidyr)
library(SQUAREM)
#library(pROC)

source('scripts/multi_intact.R')

pi1_fun <- function(z_vec,lambda = 0.5){

  p_vec <- 2*pnorm(abs(z_vec),lower.tail = FALSE)

  p_vec <- p_vec[which(p_vec != 1)]

  pi0 <- length(which(p_vec > lambda))/(length(p_vec)*(1-lambda))

  pi0_max <-  0.99

  pi1 <- 1- min(pi0_max,pi0)

  return(pi1)
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

        return(pi_hats_qvalue)
}

filenames <-c(paste0("sim_data/pi1_5_pi2_25_pi3_25/sim.phi_0_6.",seq(1,33),".summary"),
              paste0("sim_data/pi1_25_pi2_5_pi3_25/sim.phi_0_6.",seq(1,33),".summary"),
              paste0("sim_data/pi1_25_pi2_25_pi3_5/sim.phi_0_6.",seq(1,34),".summary"))



rst_df <- NULL



for (i in 1:length(filenames)){

        tmp <- marginal_posteriors(i)

	tmp_df <- data.frame("pi0" = tmp[1],"pi1" = tmp[2],"pi2" = tmp[3],"pi3" = tmp[4])

	if (i >= 1 & i < 34 ){tmp_df$most_common_alt_model <- "E"}

	if (i >= 34 & i < 67){tmp_df$most_common_alt_model <- "P"}

	if (i >= 67 & i < 101 ){tmp_df$most_common_alt_model <- "E+P"}

        rst_df <- rbind.data.frame(rst_df,tmp_df)
}



#plot each panel separately

#without legends

plot_fun_no_legend <- function(most_common_mechanism,df){
        df <- df %>% filter(most_common_alt_model == most_common_mechanism)
grDevices::cairo_pdf(paste0("sim_rst/m_intact_prior_estimates_boxplots_vary_priors_",most_common_mechanism,"_no_legend.pdf"))
print(df %>%
        dplyr::select(pi1,pi2,pi3,most_common_alt_model) %>%
        mutate(scaling_factor = pi1 + pi2 + pi3) %>%
        mutate(pi1_scaled = pi1/scaling_factor,
               pi2_scaled = pi2/scaling_factor,
               pi3_scaled = pi3/scaling_factor) %>%
        dplyr::select(pi1_scaled,pi2_scaled,pi3_scaled,most_common_alt_model) %>%
        mutate(true_e_prop = case_when(most_common_alt_model == "E" ~ 0.5,
                                       TRUE ~ 0.247)) %>%
        mutate(true_p_prop = case_when(most_common_alt_model == "P" ~ 0.5,
                                       TRUE ~ 0.25)) %>%
        mutate(true_ep_prop = case_when(most_common_alt_model == "E+P" ~ 0.5,
                                       TRUE ~ 0.253)) %>%
        pivot_longer(cols = c(pi1_scaled,
                              pi2_scaled,
                              pi3_scaled),
                     names_to = "param_type",
                     values_to = "estimate") %>%
        mutate(Posterior_Type = factor(param_type,
                                levels = c("pi1","pi2","pi3"),
                                labels = c("pi1","pi2","pi3"))) %>%
        ggplot(aes(x = param_type,y = estimate,fill = param_type)) +
        geom_boxplot() +
        xlab("") +
        ylab(expression("Estimate"))+
        geom_hline(aes(yintercept = true_e_prop,linetype = "E \u2192 Y only"),col="red",alpha = 0.9) +
        geom_hline(aes(yintercept = true_p_prop,linetype = "P \u2192 Y only"),col="dark green",alpha = 0.9) +
        geom_hline(aes(yintercept = true_ep_prop,linetype = "(E,P) \u2192 Y"),col="blue",alpha = 0.9) +
        scale_linetype_manual(name = "True Mechanism Proportion", values = c(3,2,1),
                      guide = guide_legend(override.aes = list(color = c("blue","red","dark green")))) +
        scale_x_discrete(labels = c('pi1_scaled' = expression(hat(italic(f))[italic(M[E])]),
                                      'pi2_scaled' = expression(hat(italic(f))[italic(M[P])]),
                                      "pi3_scaled" = expression(hat(italic(f))[italic(M[E+P])]))) +
        guides(fill = "none")+
        theme_bw() +
        theme(text = element_text(size = 25,face="bold"),,legend.position="none",aspect.ratio = 1))
dev.off()
}


plot_fun_no_legend(most_common_mechanism = "E",df = rst_df)
plot_fun_no_legend(most_common_mechanism = "P",df = rst_df)
plot_fun_no_legend(most_common_mechanism = "E+P",df = rst_df)



