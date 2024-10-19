#Generate grid plot of PVE for each class 

library(dplyr)
library(tidyr)
library(stringr)
library(tidyr)
library(ggplot2)


filenames <-c(paste0("sim_data/pi1_5_pi2_25_pi3_25/sim.phi_0_6.",seq(1,33),".summary"),
              paste0("sim_data/pi1_25_pi2_5_pi3_25/sim.phi_0_6.",seq(1,33),".summary"),
              paste0("sim_data/pi1_25_pi2_25_pi3_5/sim.phi_0_6.",seq(1,34),".summary"))



output_file <- "sim_rst/sim_pve_vary_priors.pdf"


df <- data.frame()

for (i in 1:length(filenames)){

	tmp <- read.table(filenames[i],h=T,sep='\t')

        df <- rbind.data.frame(df,tmp)
}


#Data frame with mean PVEs for each class

means_df <- df %>%
        pivot_longer(cols = expr_r2:ct_r2_med_protein,
                     names_to = "pve_type",
                     values_to = "pve") %>%
        mutate(pve_type = factor(pve_type,
                                 levels = c("expr_r2",
                                            "protein_r2",
                                            "ct_r2",
                                            "ct_r2_med_expr",
                                            "ct_r2_med_protein"),
                                 labels = c("Expression",
                                            "Protein",
                                            "Complex Trait, per Gene",
                                            "Expression-Mediated Complex Trait, per Gene",
					    "Protein-Mediated Complex Trait, per Gene"))) %>%   
	group_by(pve_type) %>%
	summarise(mean_pve = format(round(mean(pve),digits = 3),nsmall = 3))

means_df$labs <- paste0("Mean = ",means_df$mean_pve)
means_df$mean_pve <- as.numeric(means_df$mean_pve)


grDevices::cairo_pdf(output_file)
df %>%
	pivot_longer(cols = expr_r2:ct_r2_med_protein,
		     names_to = "pve_type",
		     values_to = "pve") %>%
	mutate(pve_type = factor(pve_type,
				 levels = c("expr_r2",
					    "protein_r2",
					    "ct_r2",
					    "ct_r2_med_expr",
					    "ct_r2_med_protein"),
				 labels = c("Expression",
                                            "Protein",
                                            "Complex Trait, per Gene",
                                            "Expression-Mediated Complex Trait, per Gene",
                                            "Protein-Mediated Complex Trait, per Gene"))) %>%
	ggplot(aes(x = pve)) + 
	geom_histogram(bins = 100) +
	xlim(c(0,1)) +
	geom_text(data = means_df,mapping = aes(x = mean_pve + 0.01, y = 10000, 
						label = labs,col = "red"),
                          hjust   = -0.1,
                          vjust   = -1) +
	geom_vline(data = means_df,aes(xintercept = mean_pve),col="red",alpha = 0.2,linetype = "dashed") +
	facet_wrap(~pve_type,ncol = 1) + 
	xlab(expression(paste(R^2))) +
	ylab("Number of Genes Across Simulated Data Sets") +
	theme_bw() + 
	theme(text = element_text(size = 10, face = "bold"),legend.position = "none")
dev.off()


