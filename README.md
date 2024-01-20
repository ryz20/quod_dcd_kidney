# quod_dcd_kidney
Repository for bulk RNA-sequencing analysis of n=271 human kidney biopsies. This is not a supported package.

# Analysis pipeline
This repository contains the following files which illustrate the methodology in our paper and accompanies supplementary material and sequencing data available on GEO and referenced in our paper.

The following files include:
explore.R (key data exploration including clustering samples and PCA)
de_analysis.R (differential expression analysis)
de_vis.R (data visualization)
pathway.R (pathway analysis)
WGCNA.R (weighted gene co-expression network analysis of cortical samples, and subsequent module analysis)

Additionally, files for supporting analysis include:
cortex_medulla_JRF.R (analysis of n=5 paired cortex and medulla biopsies undertaken by JR Ferdinand)
cortex_medulla.R (downstream analysis of n=5 paired cortex and medulla biopsies)
paired_kidneys.R (analysis of kidney pairs from n=5 individuals)

GSEA output files which are used in our paper are also uploaded under Additional Files.

# Genesets
This folder contains genesets used from external sources or other analysis for use in the analysis. 

# Additional files
This folder contains additional files relevant to the analyses above.
