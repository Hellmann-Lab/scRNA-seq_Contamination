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

## Snakemake_benchmark

Since you might not be interested in running the whole pipeline from start to finish, we provide a reduced version of the workflow here that only covers **benchmarking: application of different methods & performance evaluation**.   
To complete the input for this pipeline, some bigger files have to downloaded from [zenodo](https://zenodo.org/record/7328632#.Y3YMtOzML0s) (see next chapter): Please copy for each dataset the files `seurat.RDS` and `seurat_CAST.RDS` into the folder `input/{dataset}` and the files `filtered_feature_bc_matrix.h5` and `raw_feature_bc_matrix.h5` into a subfolder `input/{dataset}/cellranger`.  
The `config.yml` file can be modified to select benchmarking datasets, methods and parameter settings. Each method is applied to the selected benchmark datasets for background noise corrections and the outputs are evaluated with several evaluation metrics described in the manuscript: 

![workflow](https://github.com/Hellmann-Lab/scRNA-seq_Contamination/blob/main/Snakemake_benchmark/rulegraph.svg?raw=true)

If you want to add and evaluate a new method, this can be achieved by adding a new script and rule to the Snakefile that produces as output a denoised count matrix (`benchmark/corrected/{method}/{dataset}/{parameter_setting}_cormat.RDS`) and a table with estimated background noise levels per cell (`benchmark/corrected/{method}/{dataset}/{parameter_setting}_contPerCell.RDS`), which are required for all evaluation steps. 


## Benchmark Data availability
We analysed 5 mouse 10X experiments. Each is a mix of kidney cells from 3 mouse strains (BL6, SVLMJ, CAST). The data can be downloaded at [zenodo](https://zenodo.org/record/7328632#.Y3YMtOzML0s).

We provide files with cell type, strain and contamination information for each replicate in a zip-folder, where each contains 5 files:

- **filtered_feature_bc_matrix.h5** - CellRanger output, filtered count matrix
- **raw_feature_bc_matrix.h5** - CellRanger output, raw count matrix
- **seurat_CAST.RDS** - Processed Seurat object with cell type annotations, *M.m. castaneus* cells only
- **seurat.RDS** - Processed Seurat object with cell type annotations, all cells
- **perCell_noMito_CAST_binom.RDS** - Estimated background noise levels per cell in *M.m. castaneus* cells

## Analysis

Beyond the standardized pipeline, we perform further analysis to compare empty droplet, contamination and endogenous profiles (**Deconvolution**) and summarize evaluation metrics of the method benchmark (**Benchmark**).  
This folder also contains some files that are necessary to reproduce the analysis and figures:   
- **cell_metadata.RDS**: cell-wise metadata information including *replicate*, *celltype*, *Strain*, background noise level (*contPerCell_binom*) and some cell QC metrics.  
- **benchmark_metrics.RDS**: For each combination of method (CellBender, SoupX, DecontX, raw), parameter setting and replicate this table contains a collection of metrics to evaluate the performance in estimating background noise levels and improving downstream analysis after correction. 
- Proximal tubule cell markers: 1) Downloaded from [PanglaoDB](https://panglaodb.se/markers.html?cell_type=%27Proximal%20tubule%20cells%27) that were detected (**panglao_markers_Mm.RDS**) and the top10 markers with the highest average expression in PT cells (**top10_PT_markers.RDS**). 2) Genes that were detected as DE between PT and other cells after correction with CellBender/SoupX/DecontX (**DE_seurat_sigUP.RDS**)
- Statistics related to informative variants and their coverage per cell (**per_cell_stats_CAST_variants.RDS, position_stats_summary.RDS**)

## Figure scripts

- **01_dataset_description**: Strain and cell type composition of input datasets, as well as some additional statistics about the genotype estimation strategy (related to Figure 1, Suppl. Figure S3). 
- **02_backgroundRNA_estimates**: Visualization of background noise fractions per cell (related to Figure 2, Suppl. Figures S1,S2,S4).
- **03_origin_of_background RNA**: Comparison of endogenous expression, background noise contamination and empty droplet profiles (related to Figure 3, Suppl. Figures S5,S6).
- **03_02_barcode_swapping**: Identification and quantification of barcode swapping events originating from PCR chimera (related to Suppl. Figure S7)
- **04_effect_on_downstream_analysis**: Impact of background noise on specificity and detectability of marker genes (related to Figure 4, Suppl. Figures S8,S9).
- **05_benchmark_estimation**: Comparison of background noise estimation accuracy of different computional methods (related to Figure 5, Suppl. Figure S10)
- **06_benchmark_downstream_analysis**: Summarizing the method comparison results on downstream analysis effect based on the Snakemake pipeline across datasets and parameter settings (related to Figure 6, Suppl. Figures S11-15). 
