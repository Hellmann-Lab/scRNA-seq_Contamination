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

## Figure scripts

Beyond the standardized pipeline, we summarize results and perform further analysis and create figures with regards to  
- **01_dataset_description**: Strain and cell type composition of input datasets, as well as some additional statistics about the genotype estimation strategy (related to Figure 1). 
- **02_per_cell_estimates**: Visualization of backround noise fractions per cell (related to Figure 2).
- **03_origin_of_background RNA**: Comparison of endogenous expression, background noise contamination and empty droplet profiles (related to Figure 3).
- **04_effect_on_downstream_analysis**: Impact of background noise on specificity and detectability of marker genes (related to Figure 4).
- **05_06_benchmark_decontamination_methods**: Summarizing the method comparison results based on the Snakemake pipeline across datasets and parameter settings (related to Figure 5 & 6). 

## Benchmark Data availability
We analysed 5 mouse 10X experiments. Each is a mix of kidney cells from 3 mouse strains (BL6, SVLMJ, CAST).

- single cell RNA-seq [rep1](https://zenodo.org/api/files/85576185-ca53-4d69-8205-0271a1f52150/rep1.zip?versionId=ea3978b8-87b9-4c07-83f4-60003d0ebd22)
- single cell RNA-seq [rep2](https://zenodo.org/api/files/85576185-ca53-4d69-8205-0271a1f52150/rep2.zip?versionId=be135ce4-d5ef-49bb-b43c-e6e1397fc78d)
- single cell RNA-seq [rep3](https://zenodo.org/api/files/85576185-ca53-4d69-8205-0271a1f52150/rep3.zip?versionId=12e2205c-62e0-4b24-bac8-f5761475a677)
- single nucleus RNA-seq [nuc2](https://zenodo.org/api/files/85576185-ca53-4d69-8205-0271a1f52150/nuc2.zip?versionId=36b4b26e-1084-4c8e-94c3-3b4c9f216652)
- single nucleus RNA-seq [nuc3](https://zenodo.org/api/files/85576185-ca53-4d69-8205-0271a1f52150/nuc3.zip?versionId=8fe8e9d4-b8c1-4937-8d27-b3ac4692ebc9)

We provide files with cell type, strain and contaminaiton information in zip-folder, where each contains 5 files:

- **filtered_feature_bc_matrix.h5**
- **perCell_noMito_CAST_binom.RDS**
- **raw_feature_bc_matrix.h5**
- **seurat_CAST.RDS**
- **seurat.RDS**

