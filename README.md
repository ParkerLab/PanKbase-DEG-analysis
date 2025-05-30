# PanKbase DEG analysis
This directory contains Code and docs for PanKbase islet DEG analysis. Data underlying this analysis is from the single cell map which was generated using code at https://github.com/ParkerLab/PanKbase-scRNA-seq

## Content
This directory contains Jupyter notebooks with details about how to execute code, along with scripts (R/Snakemake). <br>
- `scripts`: directory contains scripts to run analyses. <br>
- `1_pseudobulk_counts.R.ipynb`: Jupyter notebook with instructions on how pseudo-bulk counts were obtained. <br>
- `2_DEG-analysis.R.ipynb`: Jupyter notebook with instructions on how to run DEG analysis using latent variables with RUVSeq and DESeq2.

## Container
A Singularity container with R packages can be found at

## Contact
Ha Vu (vthihong at umich.edu) and Stephen Parker (scjp at umich.edu)

