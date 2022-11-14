# scRNA-seq background RNA

This repository contains the necessary code to reproduce the analysis for the manuscript **"The effect of background noise and its removal on
the analysis of single-cell expression data"**.

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
- **dataset composition**: Strain and cell type composition of input datasets, as well as some additional statistics about the genotype estimation strategy. 
- **origin of background RNA**: Comparison of endogenous expression, background RNA contamination and empty droplet profiles.
- **method benchmark**: Summarizing the method comparison results across datasets and parameter settings. 
