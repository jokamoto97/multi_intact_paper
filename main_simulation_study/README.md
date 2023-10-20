# Main simulation study


## Simulation summary

In this set of simulations, we consider 1198 genes, each with 1500 common cis variants. The genotype data are obtained from the GTEx project. Due to the NIH policy regarding privacy of genetic data, we are not allowed to distribute real genotype data from the GTEx project. For this set of simulations, we do not generate effects between gene expression and protein levels. See the manuscript Methods for details on how phenotypes were generated.  

## Required R packages

data.table (1.14.8), dplyr (1.1.3), ggplot2 (3.4.2), tidyr (1.3.0), stringr (1.5.0), INTACT (1.0.2), SQUARE
M (2021.1), aod (1.3.2), pROC (1.18.4), ggpattern (1.1.0-0), ggh4x (0.2.5), qvalue (1.5.0)


## Running analysis

Before running the analysis, copy the sim_data/ directory from the Main simulation study Google drive folder at URL to this directory.

Run ```make``` to start the analysis process.
