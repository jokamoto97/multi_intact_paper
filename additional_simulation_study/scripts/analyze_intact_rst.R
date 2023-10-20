library(INTACT)
library(qvalue)

source("scripts/multi_intact.R")

filename = "sim_data/vary_prediction_models/sim.summary.v6.rst"

d = read.table(filename,head=T)

cgene = d$gene[d$delta != 0 | d$gamma != 0]

d$posterior_ptwas_expr = intact(z_vec = d$z_ptwas_expr,GLCP_vec = d$glcp_expr)
d$posterior_ptwas_protein = intact(z_vec = d$z_ptwas_protein,GLCP_vec = d$glcp_protein)

d$posterior_smr_expr = intact(z_vec = d$z_smr_expr,GLCP_vec = d$glcp_expr)
d$posterior_smr_protein = intact(z_vec = d$z_smr_protein,GLCP_vec = d$glcp_protein)

d$posterior_elastic_net_expr = intact(z_vec = d$elastic_net_expr_stat,GLCP_vec = d$glcp_expr)
d$posterior_elastic_net_protein = intact(z_vec = d$elastic_net_protein_stat,GLCP_vec = d$glcp_protein)

rej_gene_ptwas_expr = d$gene[fdr_rst2(d$posterior_ptwas_expr)$sig]
rej_gene_ptwas_protein = d$gene[fdr_rst2(d$posterior_ptwas_protein)$sig]

rej_gene_smr_expr = d$gene[fdr_rst2(d$posterior_smr_expr)$sig]
rej_gene_smr_protein = d$gene[fdr_rst2(d$posterior_smr_protein)$sig]

rej_gene_elastic_net_expr = d$gene[fdr_rst2(d$posterior_elastic_net_expr)$sig]
rej_gene_elastic_net_protein = d$gene[fdr_rst2(d$posterior_elastic_net_protein)$sig]



assess<-function(rej){

    fp_set = setdiff(rej, cgene)
    tp_set = intersect(rej, cgene)

    FDR = length(fp_set)/length(rej)

    power = length(tp_set)/length(cgene)

    return(data.frame(total_rej = length(rej),
                      FDR = round(FDR,digits = 3),
                      Power = round(power,digits = 3)))

}

out = rbind.data.frame(assess(rej_gene_ptwas_protein),
                 assess(rej_gene_ptwas_expr),
                 assess(rej_gene_smr_protein),
                 assess(rej_gene_smr_expr),
                 assess(rej_gene_elastic_net_protein),
                 assess(rej_gene_elastic_net_expr))


out$INTACT_Prediction_Method = c("PTWAS protein",
               			"PTWAS expression",
               			"SMR protein",
               			"SMR expression",
		                "Elastic net protein",
		                "Elastic net expression")

print(out[,c(4,1,2,3)])

out <- out[c(4,3,2)]
write.table(out,"sim_rst/intact_methods_power_fdr.rst",row.names = F, col.names = T, sep = '\t')


