#source("scripts/int_enrich/.R")

library(INTACT)
library(qvalue)
source("scripts/multi_intact.R")



file_prefixes = c("e_to_y_and_p_to_y_effects/p_to_e_effects/",
		  "e_to_y_and_p_to_y_effects/no_effects_between_e_and_p/",
		  "e_to_y_and_p_to_y_effects/e_to_p_effects/",
		  "p_to_y_effects_only/p_to_e_effects/",
		  "p_to_y_effects_only/no_effects_between_e_and_p/",
		  "p_to_y_effects_only/e_to_p_effects/",
		  "e_to_y_effects_only/p_to_e_effects/",
		  "e_to_y_effects_only/no_effects_between_e_and_p/",
		  "e_to_y_effects_only/e_to_p_effects/")

filename = "sim.summary.v4.rst"

dag_seq = rep(seq(1,3),3)

for(i in 1:length(file_prefixes)){

d = read.table(paste0("sim_data/",file_prefixes[i],filename),head=T)
attach(d)

#List of causal genes (via expression, protein, or both)
cgene = gene[causal_protein==1 | causal_expr == 1]

#SMR pvals
p_smr_expr = 2*pnorm(-abs(z_smr_expr))
p_smr_protein = 2*pnorm(-abs(z_smr_protein))

#SMR FDR correction
smr_qrst_expr = qvalue(p_smr_expr,fdr.level=0.05)
smr_qrst_protein = qvalue(p_smr_protein,fdr.level=0.05)

rej_smr_expr = gene[which(smr_qrst_expr$sig)]
rej_smr_protein = gene[which(smr_qrst_protein$sig)]


#PTWAS pvals and FDR correction
p_ptwas_expr = 2*pnorm(-abs(z_ptwas_expr))
p_ptwas_protein = 2*pnorm(-abs(z_ptwas_protein))
ptwas_qrst_expr = qvalue(p_ptwas_expr,fdr.level=0.05)
ptwas_qrst_protein = qvalue(p_ptwas_protein,fdr.level=0.05)
rej_ptwas_expr=gene[which(ptwas_qrst_expr$sig)]
rej_ptwas_protein=gene[which(ptwas_qrst_protein$sig)]


#FOCUS model pvals
p_focus_expr = 2*pnorm(-abs(z_focus_expr))
p_focus_protein = 2*pnorm(-abs(z_focus_protein))
focus_qrst_expr = qvalue(p_focus_expr, fdr.level=0.05)
focus_qrst_protein = qvalue(p_focus_protein, fdr.level=0.05)
rej_focus_expr = gene[which(focus_qrst_expr$sig)]
rej_focus_protein = gene[which(focus_qrst_protein$sig)]


#GLCP FDR correction (GWAS trait - protein)

glcp_lfdr_protein = sort(1-glcp_protein)
FDR_LCP_protein = cumsum(glcp_lfdr_protein)/(1:length(glcp_lfdr_protein))
glcp_thresh_protein = 1-glcp_lfdr_protein[max(which(FDR_LCP_protein<=0.05))]
rej_glcp_protein = gene[glcp_protein>=  glcp_thresh_protein]

#GLCP FDR correction (GWAS trait - expr)

glcp_lfdr_expr = sort(1-glcp_expr)
FDR_LCP_expr = cumsum(glcp_lfdr_expr)/(1:length(glcp_lfdr_expr))
glcp_thresh_expr = 1-glcp_lfdr_expr[max(which(FDR_LCP_expr<=0.05))]
rej_glcp_expr = gene[glcp_expr>=  glcp_thresh_expr]


#INTACT FDR correction (protein, default settings)

intact_pprobs_protein = intact(z_vec = z_ptwas_protein,prior = linear, GLCP_vec = glcp_protein)
intact_lfdr_protein = sort(1-intact_pprobs_protein)
FDR_intact_protein = cumsum(intact_lfdr_protein)/(1:length(intact_lfdr_protein))
intact_thresh_protein = 1-intact_lfdr_protein[max(which(FDR_intact_protein<=0.05))]
rej_intact_protein = gene[intact_pprobs_protein >= intact_thresh_protein]


#INTACT FDR correction (expr)

intact_pprobs_expr = intact(z_vec = z_ptwas_expr,prior = linear, GLCP_vec = glcp_expr)
intact_lfdr_expr = sort(1-intact_pprobs_expr)
FDR_intact_expr = cumsum(intact_lfdr_expr)/(1:length(intact_lfdr_expr))
intact_thresh_expr = 1-intact_lfdr_expr[max(which(FDR_intact_expr<=0.05))]
rej_intact_expr = gene[intact_pprobs_expr >= intact_thresh_expr]



#F test pval and FDR correction

p_f = pf(f_stat,df1 = 2, df2 = 497,lower.tail = F)
f_qrst = qvalue(p_f, fdr.level = 0.05)
rej_f = gene[which(f_qrst$sig)]



#chisq test pval and FDR correction
p_chisq = pchisq(chisq,df = 2, lower.tail = F)
chisq_qrst = qvalue(p_chisq,fdr.level=0.05)
rej_chisq = gene[which(chisq_qrst$sig)]

#post hoc approach using F test and both GLCPs
#F test FDR5% and max(glcp_expr,glcp_protein) > 0.5)

rej_post_hoc = gene[which(chisq_qrst$sig & (glcp_protein > 0.5 | glcp_expr > 0.5))]


#new Multi-INTACT approach
#m_intact_pprobs = multi_intact(f_vec = f_stat,df1 = 2, df2 = 497,
#			    GLCP_mat = cbind(glcp_expr,glcp_protein),t = 0.05)
coloc_rst = d[,c("gene","glcp_expr","glcp_protein")]
multi_intact_rst = multi_intact_scan(chisq_vec_dof = 2,
				      chisq_vec = chisq,
				      coloc_rst = coloc_rst)
multi_rst = multi_intact_rst[order(match(multi_intact_rst$gene,d$gene)),]
m_intact_pprobs = multi_rst$posterior
m_intact_lfdr = sort(1-m_intact_pprobs)
FDR_m_intact = cumsum(m_intact_lfdr)/(1:length(m_intact_lfdr))
m_intact_thresh = 1-m_intact_lfdr[max(which(FDR_m_intact<=0.05))]
rej_m_intact = gene[m_intact_pprobs >= m_intact_thresh]




assess<-function(rej){

    fp_set = setdiff(rej, cgene)
    tp_set = intersect(rej, cgene)

    FDR = length(fp_set)/length(rej)

    power = length(tp_set)/length(cgene)

    return(data.frame(total_rej = length(rej), 
		      FDR = round(FDR,digits = 3),
		      Power = round(power,digits = 3)))

#    cat("total rej: ", length(rej), "  FDR: ", FDR, "  Power: ", power, "\n")
}


out = rbind.data.frame(assess(rej_ptwas_protein),
		 assess(rej_ptwas_expr),
		 assess(rej_smr_protein),
		 assess(rej_smr_expr),
		 assess(rej_focus_protein),
		 assess(rej_focus_expr),
		 assess(rej_glcp_protein),
		 assess(rej_glcp_expr),
		 assess(rej_intact_protein),
		 assess(rej_intact_expr),
		 assess(rej_chisq),
		 assess(rej_post_hoc),
		 assess(rej_m_intact))

out$method = c("PTWAS protein",
	       "PTWAS expression",
	       "SMR protein",
	       "SMR expression",
	       "FOCUS protein",
	       "FOCUS expression",
	       "GLCP protein",
	       "GLCP expression",
	       "INTACT protein",
	       "INTACT expression",
	       "Chi-square test",
	       "Post hoc",
	       "Multi-INTACT")

print(out[,c(4,1,2,3)])


out <- out[,c(4,1,2,3)]

out$dag = paste0("dag",dag_seq[i])

write.table(out,file=paste0("sim_rst/",file_prefixes[i],"power_fdr.rst"),row.names = F, col.names = T, sep = '\t')

}
