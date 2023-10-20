#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --job-name=pm_pmn
#SBATCH --time=2-0:0:0
#SBATCH --mail-user=[YOUR EMAIL]
#SBATCH --mail-type=end,fail
#SBATCH --partition="main"
#SBATCH --output=~/slurm-%A_%a.out
#SBATCH --array=1-49

Rscript --vanilla scripts/multi_intact_prior_estimation/multi_intact_prior_estimation.R

