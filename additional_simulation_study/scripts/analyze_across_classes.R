#Generate grid plot of power & FDR for each dag and class 

library(dplyr)
library(ggh4x)
library(tidyr)
library(stringr)


class1_filenames = c("sim_rst/e_to_y_and_p_to_y_effects/p_to_e_effects/power_fdr.rst",
		     "sim_rst/e_to_y_and_p_to_y_effects/no_effects_between_e_and_p/power_fdr.rst",
		     "sim_rst/e_to_y_and_p_to_y_effects/e_to_p_effects/power_fdr.rst")

class2_filenames = c("sim_rst/p_to_y_effects_only/p_to_e_effects/power_fdr.rst",
		     "sim_rst/p_to_y_effects_only/no_effects_between_e_and_p/power_fdr.rst",
		     "sim_rst/p_to_y_effects_only/e_to_p_effects/power_fdr.rst")

class3_filenames = c("sim_rst/e_to_y_effects_only/p_to_e_effects/power_fdr.rst",
		    "sim_rst/e_to_y_effects_only/no_effects_between_e_and_p/power_fdr.rst",
		    "sim_rst/e_to_y_effects_only/e_to_p_effects/power_fdr.rst")





output_file_class1 <- "sim_rst/e_to_y_and_p_to_y_effects/rst_barplot.pdf"
output_file_class2 <- "sim_rst/p_to_y_effects_only/rst_barplot.pdf"
output_file_class3 <- "sim_rst/e_to_y_effects_only/rst_barplot.pdf"
output_file_table <- "sim_rst/power_fdr_table.txt"

output_file_class_vec <- c(output_file_class1,output_file_class2,output_file_class3)

df <- data.frame()

for (i in 1:length(class1_filenames)){
        tmp1 <- read.table(class1_filenames[i],header=T,sep = '\t')
        tmp2 <- read.table(class2_filenames[i],header=T,sep = '\t')
	tmp3 <- read.table(class3_filenames[i],header=T,sep = '\t')
	
	tmp1$class_num <- "class1"
	tmp2$class_num <- "class2"
	tmp3$class_num <- "class3"
	
	df <- rbind.data.frame(df,tmp1,tmp2,tmp3)
}

#Write to results table

df %>% filter(method != "Post hoc") %>%
	write.table(file = output_file_table,quote = F, sep = '\t',row.names = F, col.names = T)



###


df$qtl <- rep(c(rep(c("Protein","Expression"),5),rep("Protein and Expression",3)),9)
df$method[df$method == "Post hoc"] <- "Post_hoc"
df$method <- sub(" .*", "", df$method)


#Create dataframe to add line corresponding to power of multi-INTACT
power_df <- df %>%
	filter(method == "Multi-INTACT") %>%
	select(Power,dag,class_num,qtl) 

power_df <- rbind.data.frame(power_df,
			     power_df,
			     power_df) %>%
	mutate(qtl = rep(c("Protein","Expression","Protein and Expression"),each = 9)) %>%
	mutate(dag = factor(dag, levels = c("dag1",
                                            "dag2",
                                            "dag3"),labels = c(expression("P \u2192 E"),
                                                                expression("Pleiotropy"),
                                                                expression("E \u2192 P")))) %>%
        mutate(class_num = factor(class_num, levels = c("class1",
                                                "class2",
                                                "class3"), labels = c(("Class 1 \n E -> Y & P -> Y"),
                                                                        ("Class 2 \n P -> Y"),
                                                                        ("Class 3 \n E -> Y")))) %>%
       mutate(num_molecular_pheno = case_when(qtl == "Protein and Expression" ~ "Multi-Molecular Phenotype Method",
					     qtl != "Protein and Expression" ~ "Single Molecular Phenotype Methods"))	



plotdf <- df %>%
	mutate(Width = rep(c(1,1,1,1,1,1,1,1,1,1,rep(0.3,3)),9)) %>%
	mutate(method = str_replace(method, "GLCP", "Colocalization")) %>%
        select(-total_rej) %>%
        pivot_longer(cols = c(FDR,Power),
                     names_to = "Type",
                     values_to = "value") %>%
        mutate(dag = factor(dag, levels = c("dag1",
                                            "dag2",
                                            "dag3"),labels = c(expression("P \u2192 E"),
                                                                expression("Pleiotropy"),
                                                                expression("E \u2192 P")))) %>%
        mutate(class_num = factor(class_num, levels = c("class1",
                                                "class2",
                                                "class3"), labels = c(("Class 1 \n E -> Y & P -> Y"),
                                                                        ("Class 2 \n P -> Y"),
                                                                        ("Class 3 \n E -> Y")))) %>%
        mutate(method = factor(method,
                               levels = rev(c("Chi-square",
                                          "FOCUS",
                                          "PTWAS",
                                          "Post_hoc",
                                          "SMR",
                                          "Multi-INTACT",
                                          "INTACT",
                                          "Colocalization"
                                          )))) %>%
        mutate(num_molecular_pheno = case_when(qtl == "Protein and Expression" ~ "Multi-Molecular Phenotype Method",                                             qtl != "Protein and Expression" ~ "Single Molecular Phenotype Methods")) %>%
        mutate(num_molecular_pheno = factor(num_molecular_pheno, levels = rev(c("Single Molecular Phenotype Methods", "Multi-Molecular Phenotype Method")))) %>%
        mutate(qtl = factor(qtl, levels = rev(c("Expression","Protein", "Protein and Expression")))) %>%
	filter(!(method %in% c("Post_hoc","Chi-square")))







###Break up the plot by classes

class_plot_fun <- function(dat,class_int){
	if(class_int == 1){
		lab <- "Class 1 \n E -> Y & P -> Y"
	}
	if(class_int == 2){
		lab <- "Class 2 \n P -> Y"
        }
	if(class_int == 3){
		lab <- "Class 3 \n E -> Y"
        }
	dat <- dat %>%
		filter(class_num == lab) 
	grDevices::cairo_pdf(output_file_class_vec[class_int])
	print(
	dat %>%
		filter(method != "Post_hoc") %>%
		ggplot(aes(x = method,y = value, fill = Type)) +
		geom_col(position = position_dodge2(preserve = "single"),width = dat$Width) +
#        	geom_bar(position='dodge', stat = "identity") +
        	scale_fill_manual(values=c("red","#56B4E9")) +
        	facet_nested(num_molecular_pheno + qtl~ dag,scales = "free_y") +
        	geom_hline(aes(yintercept = 0.05,linetype = "FDR = 0.05"),col = 'red',alpha= 0.3) +
        	geom_hline(data = power_df %>% filter(class_num == lab), aes(yintercept = Power,linetype = "Multi-INTACT Power"),
                   	col = "blue",alpha = 0.2) +
        	scale_linetype_manual(name = "", values = c(2, 2),
                              guide = guide_legend(override.aes = list(color = c("red", "blue")))) +
        	ylab('') +
        	xlab("Method") +
		coord_flip() +
        	theme_bw() +
        	theme(text = element_text(size = 10,face="bold"),
              	axis.text.x = element_text(angle = 45, hjust=1),
		legend.title=element_blank()))
	dev.off()
}




class_plot_fun(dat = plotdf,class_int = 1)
class_plot_fun(dat = plotdf,class_int = 2)
class_plot_fun(dat = plotdf,class_int = 3)




