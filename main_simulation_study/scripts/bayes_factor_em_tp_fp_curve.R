#Plot TP vs FP for Multi-INTACT, INTACT, TWAS/PWAS, and colocalization, up to a specified FDR threshold for Multi-INTACT

library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)

d <- read.table("sim_rst/roc_m_intact_vary_priors.txt", head=T)

outfile <- "sim_rst/sim_tp_fp_curve.pdf"

count_by_desc<-function(x, value, truth){

        tp = length(which(value>=x&truth==1))
        fp = length(which(value>=x&truth==0))
        return(c(tp,fp))

}

# Expression data

mpe <- unique(round(sort(d$marginal_posterior_E,decreasing=T),6))
mpe_rst <- t(sapply(mpe, function(x) count_by_desc(x,d$marginal_posterior_E,d$E_truth)))
twase <- unique(round(sort(d$abs_z_expr,decreasing=T),3))
twase_rst <- t(sapply(twase,function(x) count_by_desc(x, d$abs_z_expr, d$E_truth)))
intacte <- unique(round(sort(d$intact_posterior_E,decreasing=T),6))
intacte_rst <- t(sapply(intacte, function(x) count_by_desc(x, d$intact_posterior_E,d$E_truth)))
enloce <-  unique(round(sort(d$glcp_expr,decreasing=T),6))
enloce_rst <- t(sapply(enloce, function(x) count_by_desc(x, d$glcp_expr, d$E_truth)))


# Protein data

mpp <- unique(round(sort(d$marginal_posterior_P,decreasing=T),6))
mpp_rst <- t(sapply(mpp, function(x) count_by_desc(x,d$marginal_posterior_P,d$P_truth)))
twasp <- unique(round(sort(d$abs_z_protein,decreasing=T),3))
twasp_rst <- t(sapply(twasp,function(x) count_by_desc(x, d$abs_z_protein, d$P_truth)))
intactp <- unique(round(sort(d$intact_posterior_P,decreasing=T),6))
intactp_rst <- t(sapply(intactp, function(x) count_by_desc(x, d$intact_posterior_P,d$P_truth)))
enlocp <- unique(round(sort(d$glcp_protein,decreasing=T),6))
enlocp_rst <- t(sapply(enlocp, function(x) count_by_desc(x, d$glcp_protein, d$P_truth)))


fdr_thresh <- 0.20


### Get cutoffs 

cutoffe <- mpe_rst[min(which(mpe_rst[,2]/mpe_rst[,1] > fdr_thresh/(1-fdr_thresh))),]
cutoffe <- round(cutoffe/1000)*1000


cutoffp <- mpp_rst[min(which(mpp_rst[,2]/mpp_rst[,1] > fdr_thresh/(1-fdr_thresh))),]
cutoffp <- round(cutoffp/1000)*1000


### Cutoff for fdr 0.05

mpe_cutoff05 <- mpe_rst[min(which(mpe_rst[,2]/mpe_rst[,1] > 0.05/(1-0.05))),]
twase_cutoff05 <- twase_rst[min(which(twase_rst[,2]/twase_rst[,1] > 0.05/(1-0.05))),]
intacte_cutoff05 <- intacte_rst[min(which(intacte_rst[,2]/intacte_rst[,1] > 0.05/(1-0.05))),]
enloce_cutoff05 <- enloce_rst[min(which(enloce_rst[,2]/enloce_rst[,1] > 0.05/(1-0.05))),]


mpp_cutoff05 <- mpp_rst[min(which(mpp_rst[,2]/mpp_rst[,1] > 0.05/(1-0.05))),]
twasp_cutoff05 <- twasp_rst[min(which(twasp_rst[,2]/twasp_rst[,1] > 0.05/(1-0.05))),]
intactp_cutoff05 <- intactp_rst[min(which(intactp_rst[,2]/intactp_rst[,1] > 0.05/(1-0.05))),]
enlocp_cutoff05 <- enlocp_rst[min(which(enlocp_rst[,2]/enlocp_rst[,1] > 0.05/(1-0.05))),]




#Make plot data frame

mpe_df <- data.frame(mpe_rst) %>% mutate(type = "m_intact")
intacte_df <- data.frame(intacte_rst) %>% mutate(type = "intact")
twase_df <- data.frame(twase_rst) %>% mutate(type = "twas")
enloce_df <- data.frame(enloce_rst) %>% mutate(type = "enloc")

expr_df <- rbind.data.frame(mpe_df,intacte_df,twase_df,enloce_df) %>% mutate(qtl = "Identifying Expression Effects")


mpp_df <- data.frame(mpp_rst) %>% mutate(type = "m_intact")
intactp_df <- data.frame(intactp_rst) %>% mutate(type = "intact")
twasp_df <- data.frame(twasp_rst) %>% mutate(type = "twas")
enlocp_df <- data.frame(enlocp_rst) %>% mutate(type = "enloc")

protein_df <- rbind.data.frame(mpp_df,intactp_df,twasp_df,enlocp_df) %>% mutate(qtl = "Identifying Protein Effects")


plotdf <- rbind.data.frame(expr_df,protein_df)
plotdf <- plotdf %>%
        mutate(type = case_when(type == "m_intact" ~ "Multi-INTACT",
                                type == "intact" ~ "INTACT",
                                type == "twas" ~ "TWAS",
                                TRUE ~ "Colocalization")) %>%
        mutate(type = factor(type)) %>%
        mutate(type = fct_relevel(type,rev(c("TWAS","Colocalization","INTACT","Multi-INTACT")))) %>%
        dplyr::rename("Method" = "type") 


text_df <- data.frame(label = "*",
		      X1 = c(mpe_cutoff05[1],
    			     intacte_cutoff05[1],
			     twase_cutoff05[1],
			     enloce_cutoff05[1],
			     mpp_cutoff05[1],
                             intactp_cutoff05[1],
                             twasp_cutoff05[1],
                             enlocp_cutoff05[1]),
		      X2 = c(mpe_cutoff05[2],
                             intacte_cutoff05[2],
                             twase_cutoff05[2],
                             enloce_cutoff05[2],
                             mpp_cutoff05[2],
                             intactp_cutoff05[2],
                             twasp_cutoff05[2],
                             enlocp_cutoff05[2]),
		      qtl = rep(c("Identifying Expression Effects","Identifying Protein Effects"),each = 4),
		      Method = rep(c("Multi-INTACT","INTACT","TWAS","Colocalization"),2))

grDevices::cairo_pdf(outfile)
print(
      plotdf %>%
	ggplot(aes(x = X2, y = X1, group = Method)) + 
	geom_line(aes(color = Method,linetype=Method)) + 
	facet_wrap(~qtl) +
	xlab("False Discoveries") +
	ylab("True Discoveries") +
	scale_linetype_manual(name = "Method", values = c(1,2,3,4))+
	scale_color_manual(values = c("purple","dark green","blue","red")) + 
	coord_cartesian(xlim = c(0,min(cutoffe[2],cutoffp[2])),ylim = c(0,min(cutoffe[1],cutoffp[1]))) +
	geom_text(data = text_df,aes(x = X2,y=X1,label = label),size=5,vjust   = 0.8) +
	theme_bw() + 
	theme(text = element_text(size = 10,face = "bold"),aspect.ratio = 1)
)
dev.off()

