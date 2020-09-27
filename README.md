# GWAS_composition

The GWAS_compositon directory contains scripts used to perform GAPIT and FarmCPU for NIR compositional data. 

1_pheno_pipeline.R: Plots raw data and gives summary statistics as well as run a random effects model. Make many plots and also outputs a BLUP data file.

2_GAPIT_FarmCPU.R: Runs GAPIT to generate PCAsa and kinship matrix. FarmCPU runs the GWAS and make QQ plots, Manhattan plots, and significant SNPs results table.

3_plots.R: Generates the plots that are used in the manuscript. 
