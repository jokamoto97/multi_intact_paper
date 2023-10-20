main: setup run
setup:

	if [ ! -d "sim_rst/" ]; then mkdir sim_rst; fi

run:
	Rscript scripts/power_fdr.R 
	Rscript scripts/plot_power_fdr.R 
	Rscript scripts/bayes_factor_em_posterior_distributions.R
	Rscript scripts/bayes_factor_em_prior_estimates.R 
	Rscript	scripts/bayes_factor_em_roc.R
	Rscript scripts/bayes_factor_em_tp_fp_curve.R
	Rscript scripts/analyze_r2.R
	Rscript scripts/bayes_factor_em_posterior_distributions_separate_by_prior.R
