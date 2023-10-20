# METSIM plasma metabolite data analysis


## Analysis details

In this analysis, we integrate the UK Biobank plasma pQTL data (Sun et al. (2023)) and the multi-tissue GTEx eQTL data (v.8) to implicate putative causal genes and mechanisms underlying plasma metabolite levels. We cannot distribute the individual-level genotype and phenotype data from the METSIM study due to privacy policy. See the manuscript Methods for details on the analysis.

We use ```openmp_wrapper``` to run some of the analysis steps in parallel. We recommend using this tool to speed up computation.

We estimate tissue-specific model priors using and EM algorithm. This step is computationally intensive and may take hours to run. If you wish to replicate this step, we have included a script to run on a high-performance cluster (commented out in the Makefile). If not, we have included the prior estimates in the data/ directory.

## Required R packages

data.table (1.14.8), dplyr (1.1.3), ggplot2 (3.4.2), tidyr (1.3.0), tidyverse (2.0.0), ggpattern (1.1.0-0), stringr (1.5.0), qvalue (1.5.0), INTACT (1.0.2), eulerr (7.0.0), org.Hs.eg.db (3.17.0), biomaRt (2.56.1), SQUAREM (2021.1), aod (1.3.2), pheatmap (1.0.12)

## Running analysis

Before running the analysis, copy the data/ directory from the Metabolite data analysis Google drive folder at https://tinyurl.com/ye2aysvf to this directory.

Run ```make``` to start the analysis process.

