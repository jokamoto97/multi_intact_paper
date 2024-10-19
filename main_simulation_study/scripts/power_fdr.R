library(qvalue)
library(INTACT)
library(dplyr)
library(ggplot2)
source("scripts/multi_intact.R")


assess<-function(rej,cgene){

    fp_set = setdiff(rej, cgene)
    
    tp_set = intersect(rej, cgene)

    FDR = length(fp_set)/length(rej)

    power = length(tp_set)/length(cgene)

    return(data.frame(total_rej = length(rej),
                      FDR = round(FDR,digits = 3),
                      Power = round(power,digits = 3)))

}



power_fdr <- function(dat){

	cgenes <- dat$gene[dat$gamma != 0 | dat$delta != 0]

	#PTWAS pvals and FDR correction
	
	p_ptwas_expr <- 2*pnorm(-abs(dat$z_ptwas_expr))
	p_ptwas_protein <- 2*pnorm(-abs(dat$z_ptwas_protein))

	ptwas_qrst_expr <- qvalue(p_ptwas_expr,fdr.level=0.05)
	ptwas_qrst_protein <- qvalue(p_ptwas_protein,fdr.level=0.05)
	
	rej_ptwas_expr <- dat$gene[which(ptwas_qrst_expr$sig)]
	rej_ptwas_protein <- dat$gene[which(ptwas_qrst_protein$sig)]
	
	#Colocalization FDR correction

	rej_glcp_protein <- dat$gene[fdr_rst2(dat$glcp_protein)$sig]
	rej_glcp_expr <- dat$gene[fdr_rst2(dat$glcp_expr)$sig]

	#INTACT FDR correction
	
	rej_intact_protein <- dat$gene[fdr_rst2(intact(z_vec = dat$z_ptwas_protein,
								 prior = linear, 
								 GLCP_vec = dat$glcp_protein))$sig]
	rej_intact_expr <- dat$gene[fdr_rst2(intact(z_vec = dat$z_ptwas_expr,
                                                                 prior = linear,
                                                                 GLCP_vec = dat$glcp_expr))$sig]

	#Multi-INTACT FDR correction

	coloc_rst <- dat[,c("gene","glcp_expr","glcp_protein")]
	multi_intact_rst <- multi_intact_scan(chisq_vec_dof = 2,
                                      chisq_vec = dat$chisq,
                                      coloc_rst = coloc_rst)
	multi_rst <- multi_intact_rst[order(match(multi_intact_rst$gene,dat$gene)),]
	m_intact_pprobs <- multi_rst$posterior
	rej_m_intact <- multi_rst$gene[fdr_rst2(m_intact_pprobs)$sig]

	
	#FDR/Power rst

	ptwas_protein_rst <- assess(rej_ptwas_protein,cgene=cgenes)
	ptwas_expr_rst <- assess(rej_ptwas_expr,cgene=cgenes)

	glcp_protein_rst <- assess(rej_glcp_protein,cgene=cgenes)
	glcp_expr_rst <- assess(rej_glcp_expr,cgene=cgenes)

	intact_protein_rst <- assess(rej_intact_protein,cgene=cgenes)
	intact_expr_rst <- assess(rej_intact_expr,cgene=cgenes)

	m_intact_rst <- assess(rej_m_intact,cgene=cgenes)


	out <- rbind.data.frame(ptwas_protein_rst,
			       ptwas_expr_rst,
			       glcp_protein_rst,
			       glcp_expr_rst,
			       intact_protein_rst,
			       intact_expr_rst,
			       m_intact_rst)
	out$method <- c("PTWAS protein",
               "PTWAS expression",
               "GLCP protein",
               "GLCP expression",
               "INTACT protein",
               "INTACT expression",
               "Multi-INTACT")
	
	#out <- out[,c(4,5,1,2,3)]
	
	return(out)
}

filenames <-c(paste0("sim_data/pi1_5_pi2_25_pi3_25/sim.phi_0_6.",seq(1,33),".summary"),
              paste0("sim_data/pi1_25_pi2_5_pi3_25/sim.phi_0_6.",seq(1,33),".summary"),
              paste0("sim_data/pi1_25_pi2_25_pi3_5/sim.phi_0_6.",seq(1,34),".summary"))


output_file <- "sim_rst/power_fdr_vary_priors.txt"


rst_df <- NULL

for (i in 1:length(filenames)){
	df <- read.table(filenames[i],h=T,sep = '\t')

	rst <- power_fdr(df)

	rst_df <- rbind.data.frame(rst_df,rst)
}


rst_df %>% 
	group_by(method) %>%
	summarise(fdr = mean(FDR),
		  power = mean(Power),
		  FDR_SE = sd(FDR)/sqrt(100),
		  Power_SE = sd(Power)/sqrt(100)) %>%
	dplyr::rename("Power" = "power","FDR" = 'fdr') %>%
	write.table(file = output_file,row.names = F, col.names = T, sep = '\t',quote = F)

#Table of GPPC and truth

rst_df <- NULL

for (i in 1:length(filenames)){
        df <- read.table(filenames[i],h=T,sep = '\t')

	df$cgene <- 0
	
	df$cgene[df$gamma != 0 | df$delta != 0] <- 1

	coloc_rst <- df[,c("gene","glcp_expr","glcp_protein")]
        multi_intact_rst <- multi_intact_scan(chisq_vec_dof = 2,
                                      chisq_vec = df$chisq,
                                      coloc_rst = coloc_rst)

	multi_intact_rst <- multi_intact_rst[,c(1,9)]

	df <- merge(df,multi_intact_rst,by="gene")

	tmp_rst <- df[,c("gene","cgene","posterior","expr_r2","protein_r2","ct_r2","ct_r2_med_expr","ct_r2_med_protein")]

	colnames(tmp_rst)[3] <- "GPPC"

	tmp_rst$sim_num <- i

        rst_df <- rbind.data.frame(rst_df,tmp_rst)
}

#calculate power stratified by PVE


write.table(rst_df,file="sim_rst/gppc_summary_table.txt",quote=F,sep='\t',row.names=F,col.names=T)


d = read.table("sim_rst/gppc_summary_table.txt", head=T)
attach(d)


lfdr = sort(1-GPPC)
FDR = cumsum(lfdr)/1:length(lfdr)
thresh = 1-lfdr[max(which(FDR<=0.05))]

powers <- c()

rej = cgene[GPPC>=thresh&expr_r2<=.05]
powers <- c(powers,sum(rej)/sum(cgene[expr_r2<=0.05]))
rej = cgene[GPPC>=thresh&protein_r2<=.05]
powers <- c(powers,sum(rej)/sum(cgene[protein_r2<=0.05]))
rej = cgene[GPPC>=thresh&ct_r2<=.05]
powers <- c(powers,sum(rej)/sum(cgene[ct_r2<=0.05]))

rej = cgene[GPPC>=thresh&expr_r2>0.05&expr_r2<=0.15]
powers <- c(powers,sum(rej)/sum(cgene[expr_r2>0.05&expr_r2<=0.15]))
rej = cgene[GPPC>=thresh&protein_r2>0.05&protein_r2<=0.15]
powers <- c(powers,sum(rej)/sum(cgene[protein_r2>0.05&protein_r2<=0.15]))
rej = cgene[GPPC>=thresh&ct_r2>0.05&ct_r2<=0.15]
powers <- c(powers,sum(rej)/sum(cgene[ct_r2>0.05&ct_r2<=0.15]))

rej = cgene[GPPC>=thresh&expr_r2>0.15&expr_r2<=0.25]
powers <- c(powers,sum(rej)/sum(cgene[expr_r2>0.15&expr_r2<=0.25]))
rej = cgene[GPPC>=thresh&protein_r2>0.15&protein_r2<=0.25]
powers <- c(powers,sum(rej)/sum(cgene[protein_r2>0.15&protein_r2<=0.25]))
rej = cgene[GPPC>=thresh&ct_r2>0.15&ct_r2<=0.25]
powers <- c(powers,sum(rej)/sum(cgene[ct_r2>0.15&ct_r2<=0.25]))

rej = cgene[GPPC>=thresh&expr_r2>0.25]
powers <- c(powers,sum(rej)/sum(cgene[expr_r2>0.25]))
rej = cgene[GPPC>=thresh&protein_r2>0.25]
powers <- c(powers,sum(rej)/sum(cgene[protein_r2>0.25]))
rej = cgene[GPPC>=thresh&ct_r2>0.25]
powers <- c(powers,sum(rej)/sum(cgene[ct_r2>0.25]))

plotdf <- data.frame("pve" = rep(c("Low\n(0, 0.05]",
				   "Low-Medium\n(0.05, 0.15]",
				   "Medium-High\n(0.15, 0.25]",
				   "High\n(0.25, 1]"),each=3),
		     "type" = rep(c("Expression","Protein","Complex Trait"),4),
		     "power" = powers)

pdf("sim_rst/power_stratified_by_pve.pdf")
plotdf %>%
	mutate(pve = factor(pve,levels = c("Low\n(0, 0.05]",
					    "Low-Medium\n(0.05, 0.15]",
					    "Medium-High\n(0.15, 0.25]",
					    "High\n(0.25, 1]"))) %>%
	mutate(type = factor(type,levels = c("Expression","Protein","Complex Trait"))) %>%
	ggplot(aes(x = type,y=power,fill=type)) + 
	geom_bar(position="dodge",stat="identity") +
	facet_wrap(~pve,nrow=1) +
	scale_fill_manual(values = c("blue","red","dark green")) +
	xlab("") +
	ylab("Power") + 
	theme_bw() + 
	theme(text = element_text(size = 10,face="bold"),legend.position="none",axis.text.x = element_text(angle = 45, hjust=1),aspect.ratio = 1)
dev.off()










