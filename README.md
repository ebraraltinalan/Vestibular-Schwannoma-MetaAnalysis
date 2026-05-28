Transcriptomic Meta-Analysis of Vestibular Schwannoma

This repository contains the R scripts used in the study "Transcriptomic Meta-Analysis Identifies Dysregulated Pathways and Potential Therapeutic Targets in Vestibular Schwannoma."

Workflow
The repository includes scripts for:

Gene-level random-effects meta-analysis using the Sidik–Jonkman estimator.
Gene Ontology (GO) enrichment analysis.
Volcano plot visualization.
DrugBank-based drug–gene interaction analyses and network visualization.

Both covariate-adjusted and non-adjusted meta-analysis workflows are provided. In the covariate-adjusted analyses, dataset-specific differential expression estimates were generated using available clinical covariates prior to meta-analysis.

Repository Structure
The analysis pipeline is organized as follows:

01_metaanalysis.R / 01b_metaanalysis_withcov.R: Main meta-analysis scripts.

02_volcanoplots.R: Visualization of significant differentially expressed genes.

03_enrichment_analysis.R: GO Biological Process enrichment analysis.

04_drugbank_.R / 04b_drugbank_withcov.R: Drug–gene interaction filtering and network generation.

Prerequisites
To run these scripts, you will need R and the following key packages:
metafor, dplyr, clusterProfiler, org.Hs.eg.db, dbparser, ggplot2, ggraph, igraph, ggrepel, readxl, writexl.

Input Data
The analyses are based on the following publicly available GEO datasets:

GSE108237

GSE141801

GSE108524

GSE39645

The repository does not include raw microarray preprocessing scripts. Users should obtain the original datasets from GEO and generate the corresponding dataset-specific differential expression results prior to running the downstream analyses provided here.

Notes
The scripts contain local file paths used during the original analyses and may require modification (updating to your local machine's directory structure) before execution on another system.
