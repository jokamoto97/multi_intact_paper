source('scripts/multi_intact.R')


sumstat_rst = read.table("sim_data/sum_stats_approaches/compare_x2_estimators.tab",sep = '\t',header=T)

filename = "sim_data/sum_stats_approaches/sim.summary.v3.rst"

d = read.table(filename,head=T)
d = merge(sumstat_rst,d,by = 'gene')


#List of causal genes (via expression, protein, or both)
cgene = d$gene[d$causal_protein==1 | d$causal_expr == 1]



coloc_rst = d[,c("gene","glcp_expr","glcp_protein")]
multi_intact_rst_indiv = multi_intact_scan(chisq_vec_dof = 2,
                                      chisq_vec = d$X2_indiv,
                                      coloc_rst = coloc_rst)

multi_intact_rst_dapg = multi_intact_scan(chisq_vec_dof = 2,
                                      chisq_vec = d$X2_dapg,
                                      coloc_rst = coloc_rst)

multi_intact_rst_multiX = multi_intact_scan(chisq_vec_dof = 2,
                                      chisq_vec = d$X2_multi,
                                      coloc_rst = coloc_rst)

rej_gene_indiv = multi_intact_rst_indiv$gene[fdr_rst2(posterior = multi_intact_rst_indiv$posterior)$sig]
rej_gene_dapg = multi_intact_rst_dapg$gene[fdr_rst2(posterior = multi_intact_rst_dapg$posterior)$sig]
rej_gene_multiX = multi_intact_rst_multiX$gene[fdr_rst2(posterior = multi_intact_rst_multiX$posterior)$sig]



assess<-function(rej){

    fp_set = setdiff(rej, cgene)
    tp_set = intersect(rej, cgene)

    FDR = length(fp_set)/length(rej)

    power = length(tp_set)/length(cgene)

    return(data.frame(total_rej = length(rej),
                      Multi_INTACT_FDR = round(FDR,digits = 3),
                      Multi_INTACT_Power = round(power,digits = 3)))
}


out = rbind.data.frame(assess(rej_gene_indiv),
                 assess(rej_gene_dapg),
                 assess(rej_gene_multiX))

out$Joint_Test_Method = c("Individual-level",
               "Summary-level (DAP-G formula)",
               "Summary-level (MultiXcan formula)")


out = out[,c(4,3,2)]


write.table(out,file = 'sim_rst/indiv_vs_summary_data_power_fdr.txt')
