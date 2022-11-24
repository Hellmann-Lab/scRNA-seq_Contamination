# The effect of background noise and its removal on the analysis of single-cell expression data

This repository contains the necessary code to reproduce the analysis for the manuscript **"The effect of background noise and its removal on the analysis of single-cell expression data"**.

## Snakemake

We process multiple single-cell and single-nucleus RNA-seq datasets of mouse kidney with a common pipeline to estimate levels of background RNA contamination and compare different methods for correction of background RNA. This includes:  
- **cell calling**  
- **cell type assignment:** Reference based classification of single-cells using a publicly available annotated [dataset](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02048-6).
- **genotyping and strain assignment:** Using a list of genetic variants between mouse strains downloaded from the [Mouse Genomes Project](https://www.sanger.ac.uk/data/mouse-genomes-project//) we determine the strain identity of each cell by genotyping single cells with [cellsnp-lite](https://cellsnp-lite.readthedocs.io/) and demultiplexing with [vireo](https://vireosnp.readthedocs.io/). Importantly, this also provides matrices of allele counts per cell barcode. 
- **genotype estimation:** We identify cross-strain contamination based on genetic variants and estimate background RNA levels per single cell from this.
- **filtering and preprocessing:** Basic filtering, processing and clustering steps to prepare the count matrix for further analysis.
- **applying different correction methods:** We compare three methods that are designed to remove noise originating from ambient/ background RNA: [CellBender](https://www.biorxiv.org/content/10.1101/791699v1), [SoupX](https://doi.org/10.1093/gigascience/giaa151) and [DecontX](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1950-6) in a range of different parameter settings. 
- **evaluation of correction performance:** We evaluated the output of each method by comparing to our genotype based estimations and calculating metrics to assess denoising performance.

## Analysis

Beyond the standardized pipeline, we perform further analysis to compare empty droplet, contamination and endogenous profiles (**Deconvolution**) and summarize evaluation metrics of the method benchmark (**Benchmark**).

## Figure scripts

- **01_dataset_description**: Strain and cell type composition of input datasets, as well as some additional statistics about the genotype estimation strategy (related to Figure 1, Suppl. Figure S2). 
- **02_backgroundRNA_estimates**: Visualization of backround noise fractions per cell (related to Figure 2, Suppl. Figures S1,S3).
- **03_origin_of_background RNA**: Comparison of endogenous expression, background noise contamination and empty droplet profiles (related to Figure 3, Suppl. Figures S4,S5).
- **04_effect_on_downstream_analysis**: Impact of background noise on specificity and detectability of marker genes (related to Figure 4, Suppl. Figures S6,S7).
- **05_benchmark_estimation**: Comparison of background noise estimation accuracy of different computional methods (related to Figure 5, Suppl. Figure S8)
- **06_benchmark_downstream_analysis**: Summarizing the method comparison results on downstream analysis effect based on the Snakemake pipeline across datasets and parameter settings (related to Figure 6, Suppl. Figures S9-13). 

## Benchmark Data availability
We analysed 5 mouse 10X experiments. Each is a mix of kidney cells from 3 mouse strains (BL6, SVLMJ, CAST). The data can be downloaded at [zenodo](https://zenodo.org/record/7328632#.Y3YMtOzML0s).

We provide files with cell type, strain and contamination information for each replicate in a zip-folder, where each contains 5 files:

- filtered_feature_bc_matrix.h5 - CellRanger output, filtered count matrix
- raw_feature_bc_matrix.h5 - CellRanger output, raw count matrix
- seurat_CAST.RDS - Processed Seurat object with cell type annotations, *M.m. castaneus* cells only
- seurat.RDS - Processed Seurat object with cell type annotations, all cells
- perCell_noMito_CAST_binom.RDS - Estimated background noise levels per cell in *M.m. castaneus* cells



