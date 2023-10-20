#check equivalence of X2 statistic using DAP-G notation vs MultiXcan notation
library(dplyr)
library(ggplot2)
library(aod)
library(tidyr)




df <- read.table("sim_data/sum_stats_approaches/compare_x2_estimators.tab",sep = '\t',header=T)


pdf('sim_rst/compare_x2_estimators.pdf')
df %>%
	rename("Approach_1" = "X2_dapg",
	       "Approach_2" = "X2_multi") %>%
	pivot_longer(cols = c("Approach_1","Approach_2"),
		     names_to = "summary_method",
		     values_to = "X2") %>%
	ggplot(aes(x = X2_indiv, y = X2)) +
	geom_point() + 
	ylab(expression(paste(chi^2, " Test Statistic Using Summary Data Method",sep=""))) +
	xlab(expression(paste(chi^2, " Test Statistic Using Individual-level Data",sep=""))) +	
	facet_wrap(~summary_method) +	
	geom_abline(intercept = 0, slope = 1, col = "red") + 
	theme_bw()+
	theme(text = element_text(size = 10,face = 'bold'),aspect.ratio = 1)
dev.off()


