library(INTACT)
library(qvalue)

source("scripts/multi_intact.R")

filename = "sim_data/vary_prediction_models/sim.summary.v7.rst"

d = read.table(filename,head=T)

cgene = d$gene[d$delta != 0 | d$gamma != 0]

coloc_rst = d[,c("gene","glcp_expr","glcp_protein")]

multi_intact_posterior_ptwas = multi_intact_scan(chisq_vec_dof = 2,
						 chisq_vec = d$chisq,
						 coloc_rst = coloc_rst)

multi_intact_posterior_smr = multi_intact_scan(chisq_vec_dof = d$X2_df_smr,
                                                 chisq_vec = d$X2_smr,
                                                 coloc_rst = coloc_rst)

multi_intact_posterior_elastic_net = multi_intact_scan(chisq_vec_dof = d$X2_df_elastic_net,
                                                 chisq_vec = d$X2_elastic_net,
                                                 coloc_rst = coloc_rst)



rej_gene_ptwas = multi_intact_posterior_ptwas$gene[fdr_rst2(multi_intact_posterior_ptwas$posterior)$sig]

rej_gene_smr = multi_intact_posterior_smr$gene[fdr_rst2(multi_intact_posterior_smr$posterior)$sig]

rej_gene_elastic_net = multi_intact_posterior_elastic_net$gene[fdr_rst2(multi_intact_posterior_elastic_net$posterior)$sig]


assess<-function(rej){

    fp_set = setdiff(rej, cgene)
    tp_set = intersect(rej, cgene)

    FDR = length(fp_set)/length(rej)

    power = length(tp_set)/length(cgene)

    return(data.frame(total_rej = length(rej),
                      FDR = round(FDR,digits = 3),
                      Power = round(power,digits = 3)))

}

out = rbind.data.frame(assess(rej_gene_ptwas),
                 assess(rej_gene_smr),
                 assess(rej_gene_elastic_net))


out$Multi_INTACT_Prediction_Method = c("PTWAS",
                                	"SMR",
                                	"Elastic net")

print(out[,c(4,1,2,3)])

out <- out[c(4,3,2)]
write.table(out,"sim_rst/multi_intact_methods_power_fdr.rst",row.names = F, col.names = T, sep = '\t')

