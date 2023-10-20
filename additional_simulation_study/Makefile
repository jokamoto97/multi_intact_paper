main: setup run
setup:

	if [ ! -d "sim_rst/" ]; then mkdir sim_rst; fi
	if [ ! -d "sim_rst/e_to_y_and_p_to_y_effects" ]; then mkdir sim_rst/e_to_y_and_p_to_y_effects; fi
	if [ ! -d "sim_rst/p_to_y_effects_only" ]; then mkdir sim_rst/p_to_y_effects_only; fi
	if [ ! -d "sim_rst/e_to_y_effects_only" ]; then mkdir sim_rst/e_to_y_effects_only; fi
	if [ ! -d "sim_rst/e_to_y_and_p_to_y_effects/p_to_e_effects" ]; then mkdir sim_rst/e_to_y_and_p_to_y_effects/p_to_e_effects; fi
	if [ ! -d "sim_rst/e_to_y_and_p_to_y_effects/no_effects_between_e_and_p" ]; then mkdir sim_rst/e_to_y_and_p_to_y_effects/no_effects_between_e_and_p; fi
	if [ ! -d "sim_rst/e_to_y_and_p_to_y_effects/e_to_p_effects" ]; then mkdir sim_rst/e_to_y_and_p_to_y_effects/e_to_p_effects; fi
	if [ ! -d "sim_rst/p_to_y_effects_only/p_to_e_effects" ]; then mkdir sim_rst/p_to_y_effects_only/p_to_e_effects; fi
	if [ ! -d "sim_rst/p_to_y_effects_only/no_effects_between_e_and_p" ]; then mkdir sim_rst/p_to_y_effects_only/no_effects_between_e_and_p; fi
	if [ ! -d "sim_rst/p_to_y_effects_only/e_to_p_effects" ]; then mkdir sim_rst/p_to_y_effects_only/e_to_p_effects; fi
	if [ ! -d "sim_rst/e_to_y_effects_only/p_to_e_effects" ]; then mkdir sim_rst/e_to_y_effects_only/p_to_e_effects; fi
	if [ ! -d "sim_rst/e_to_y_effects_only/no_effects_between_e_and_p" ]; then mkdir sim_rst/e_to_y_effects_only/no_effects_between_e_and_p; fi
	if [ ! -d "sim_rst/e_to_y_effects_only/e_to_p_effects" ]; then mkdir sim_rst/e_to_y_effects_only/e_to_p_effects; fi

run:
	Rscript scripts/analyze2.R
	Rscript scripts/analyze_across_classes.R 
	Rscript scripts/analyze_sum_stats_approaches.R
	Rscript scripts/check_equivalence_of_estimators.R
	Rscript scripts/analyze_intact_rst.R
	Rscript scripts/analyze_multi_intact_rst.R

