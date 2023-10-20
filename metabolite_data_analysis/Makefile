main: setup run
setup:

	if [ ! -d "results/" ]; then mkdir results; fi
	if [ ! -d "results/venn_diagrams" ]; then mkdir results/venn_diagrams; fi
	if [ ! -d "results/ukbb_gtex_multi_intact" ]; then mkdir results/ukbb_gtex_multi_intact; fi
	if [ ! -d "results/ukbb_gtex_multi_intact/gppc_rst" ]; then mkdir results/ukbb_gtex_multi_intact/gppc_rst; fi
	if [ ! -d "results/ukbb_gtex_multi_intact/gppc_rst/gppc_summary" ]; then mkdir results/ukbb_gtex_multi_intact/gppc_rst/gppc_summary; fi
	if [ ! -d "results/ukbb_gtex_multi_intact/em_algorithm_rst" ]; then mkdir results/ukbb_gtex_multi_intact/em_algorithm_rst; fi

run:
	Rscript scripts/intact_venn_diagrams.R 
	#sbatch scripts/multi_intact_prior_estimation/multi_intact_prior_estimation.sh  
	openmp_wrapper -d scripts/multi_intact_gppc/multi_intact_gppc_ow.cmd -t 12 2>/dev/null
	openmp_wrapper -d scripts/multi_intact_model_posteriors/multi_intact_compute_posteriors_ow.cmd -t 12 2>/dev/null
	Rscript scripts/concat_gppc_rst.R
	Rscript scripts/plot_gppc_summary.R
	openmp_wrapper -d scripts/sig_gppc_triplets/sig_gppc_triplets_ow.cmd -t 12 2>/dev/null
	Rscript scripts/kba_gppc_venn_diagram.R
	Rscript scripts/plot_gppc_proportions.R
	Rscript scripts/n_genes_tested_per_tissue.R
	Rscript scripts/n_genes_tested_per_tissue_pairwise.R
	Rscript scripts/tissue_concordance_heatmap.R
	Rscript scripts/tissue_concordance_histogram.R
	openmp_wrapper -d scripts/get_marginal_e_and_p_posterior_0_5_triplets/get_marginal_e_and_p_posterior_0_5_triplets_ow.cmd -t 12 2>/dev/null
	Rscript scripts/analyze_high_marginal_e_and_p_posteriors.R
	Rscript scripts/analyze_high_marginal_e_and_p_posteriors_kba.R
	Rscript scripts/analyze_high_marginal_e_and_p_posteriors_facet.R
	Rscript scripts/analyze_high_marginal_e_and_p_posteriors_sign_pval.R
	openmp_wrapper -d scripts/twas_pwas_spearman_correlation/twas_pwas_spearman_correlation_ow.cmd -t 12 2>/dev/null
	Rscript scripts/twas_pwas_spearman_correlation/assess_twas_pwas_spearman_correlation.R
	Rscript scripts/lipid_gsea_rst.R

