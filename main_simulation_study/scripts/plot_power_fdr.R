#Generate grid plot of power & FDR, marginal over dags 

library(dplyr)
library(ggh4x)
library(tidyr)
library(stringr)
library(ggpattern)

infile <- "sim_rst/power_fdr_vary_priors.txt"
output_file <- "sim_rst/power_fdr_vary_priors.pdf"

df <- read.table(infile,h=T,sep = '\t')

df$qtl <- c(rep(c("Expression","Protein"),2),"Protein and Expression","Expression","Protein")
df$method <- sub(" .*", "", df$method)
df$method[df$method == "PTWAS"] <- "TWAS"
df$method[df$method == "GLCP"] <- "Colocalization"


#Change protein-TWAS to PWAS

df$method[7] <- "PWAS"


plotdf <- df %>%
	mutate(Width = c(1,1,1,1,0.35,1,1)) %>%
        pivot_longer(cols = c(FDR,Power),
                     names_to = "Type",
                     values_to = "value") %>%
        mutate(method = factor(method,
                               levels = rev(c(
                                          "Multi-INTACT",
                                          "INTACT",
                                          "Colocalization",
					  "TWAS",
					  "PWAS"
                                          )))) %>%
        mutate(num_molecular_pheno = case_when(qtl == "Protein and Expression" ~ "Multi-Molecular Phenotype Method",                                             qtl != "Protein and Expression" ~ "Single Molecular Phenotype Methods")) %>%
        mutate(num_molecular_pheno = factor(num_molecular_pheno, levels = rev(c("Single Molecular Phenotype Methods", "Multi-Molecular Phenotype Method")))) %>%
        mutate(qtl = factor(qtl, levels = rev(c("Expression","Protein", "Protein and Expression"))))


#Asterisk vector for methods with exceeded FDR

annot_vec <- plotdf %>%
	mutate(Type = factor(Type, levels = c("Power","FDR"))) %>%
	pivot_wider(names_from = "Type",values_from = "value") %>%
	mutate(annot = case_when(FDR > 0.05 ~ "*",
				 TRUE ~ " ")) %>%
	pivot_longer(cols = c(FDR,Power),
                     names_to = "Type",
                     values_to = "value") %>%
	mutate(annot = case_when(Type == "Power" ~ annot,
				 Type == "FDR" ~ " ")) %>%
	dplyr::select(annot) %>% unlist() %>% as.vector()


plotdf <- plotdf %>%  
	mutate(SE = case_when(Type == "FDR" ~ FDR_SE,
                              TRUE ~ Power_SE)) %>%
	mutate(lower = value - SE,
	       upper = value + SE)

#Multi_INTACT power

m_intact_power <- as.numeric(plotdf$value[plotdf$method == "Multi-INTACT" & plotdf$Type == "Power"])


#plot rst

grDevices::cairo_pdf(output_file)
        print(
        plotdf %>%
		mutate(Type = factor(Type, levels = c("Power","FDR"))) %>%
		mutate(annot = annot_vec) %>%
                ggplot(aes(x = method,y = value, fill = Type,pattern = Type)) +
#                geom_col(position = position_dodge2(preserve = "single"),width = plotdf$Width) +
		geom_bar_pattern(position = position_dodge2(preserve = "single"),stat="identity",
				 width = plotdf$Width,
		                 pattern_fill = "black",
		      		 pattern_angle = 45,
		 		 pattern_density = 0.1,
			       	 pattern_spacing = 0.025,
				 pattern_key_scale_factor = 0.6)+
		geom_errorbar(aes(ymin = lower,
                                ymax = upper),position=position_dodge2(preserve = "single"),
			      width=plotdf$Width) +
	        scale_pattern_manual(values = c(FDR = "stripe", Power = "none")) +
                scale_fill_manual(values=c("#56B4E9","red")) +
#                facet_nested(num_molecular_pheno + qtl~ class_num,scales = "free_y") +
		facet_nested(~ num_molecular_pheno + qtl,scales='free_x') +
#		force_panelsizes(rows = unit(4, "cm"),
#                   cols = unit(c(5, 5, 5), "cm"),
#                   TRUE) +
		geom_text(aes(label = annot),col="blue",hjust = 3.5,vjust = -0.5) +
                geom_hline(aes(yintercept = 0.05,linetype = "FDR = 0.05"),col = 'red',alpha= 0.5) +
                geom_hline(aes(yintercept = m_intact_power,linetype = "Multi-INTACT Power"),
                        col = "blue",alpha = 0.4) +
                scale_linetype_manual(name = "", values = c(3, 2),
                              guide = guide_legend(override.aes = list(color = c("red", "blue")))) +
                ylab('') +
                xlab("Method") +
#                coord_flip() +
                theme_bw() +
                theme(text = element_text(size = 10,face="bold"),
                axis.text.x = element_text(angle = 45, hjust=1),
		legend.position="none",
		legend.title=element_blank(),aspect.ratio = 1)) 
dev.off()

