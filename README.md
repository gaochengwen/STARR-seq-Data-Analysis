STARR-seq Data Analysis Pipeline
This repository contains a set of scripts designed to process, analyze, and visualize data related to variant analysis. The pipeline includes steps for filtering data, extracting alleles, performing quality control (QC), and conducting downstream analyses like MPRA (Massively Parallel Reporter Assay).

Prerequisites
R (version 4.2.0 or higher)
Python 3.x with Bio and numpy libraries
Required R libraries:
tidyverse
dplyr
stringr
readr
purrr
tibble
RColorBrewer
ape
reshape2
mpra (for MPRA analysis)
Ensure you have the necessary dependencies installed before running the pipeline.

Python Scripts
extract_alleles.py
Extracts alleles from sequence data using pairwise sequence alignment, and outputs allele information along with sequence comparisons.
R Scripts
analysis.R
Analyzes variant and experimental data, calculating ratios of gene expression between different experimental conditions. It generates QQ plots for data visualization and assesses allele-specific effects.

count_qc.R
Performs quality control on count data by filtering out low-expressed genes and outlier samples. It generates various QC plots, including hierarchical clustering and D-statistics.

sum.R
Aggregates QC-filtered data by allele and SNP, producing a summary file for downstream analysis.

mpra_analysis.R
Conducts MPRA analysis on DNA and RNA expression data, focusing on allele-specific effects. It outputs significant results, including allele-specific expression and statistical analysis.

Usage
Prepare Input Files

Variant data (e.g., SNP information)
Sequence data (e.g., DNA and RNA expression data)
Run the Scripts

Start with 01_add_seqinfo.sh to process the sequence and allele data.
Continue through each script in sequence, following the steps to filter, summarize, and analyze the data.
Final MPRA analysis can be conducted using 08_mpra_analysis.sh.

