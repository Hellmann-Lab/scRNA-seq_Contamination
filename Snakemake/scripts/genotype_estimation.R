args <- commandArgs(TRUE)

library(tidyverse)
library(vcfR)
library(cowplot)

source("/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/genotype_estimation/genotype_estimation_functions.R")
ENSEMBL_to_SYMBOL <- readRDS("/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/genotype_estimation/ENSEMBL_to_SYMBOL.RDS")

genotype_estimation <- function(cellSNP_out, vireo_out, vireo_vcf, outdir){
  system(paste0("mkdir ",outdir))
  
  # load vcf files:
  strain_vcf <- read.vcfR(vireo_vcf)
  cellSNP_vcf <- read.vcfR(paste0(cellSNP_out,"/cellSNP.base.vcf.gz"))
  
  ## Prepare reference vcf file
  #   Input: Reference vcf file (subsetted for variants with coverage to reduce size)  
  #   Output: Dataframe: genotype of each strain per position + position and gene annotations
  strain_gt <- prepare_reference_vcf(strain_vcf,ENSEMBL_to_SYMBOL)
  
  ## Check if variants in cellSNP and vireo vcf are matching:
  check_SNP_matching(strain_vcf, cellSNP_vcf, strain_gt)
  
  
  ## Get allele counts per barcode
  csc_intermediate <- get_allele_counts_per_cell(strain_gt, cellSNP_out, vireo_out,cellSNP_vcf, coverage_filter = 100)
  saveRDS(csc_intermediate, paste0(outdir,"/csc_intermediate.RDS"))
  gc()
  
  ## Calculate cross-genotype contamination
  csc_gt <- calculate_crossGT_contamination(csc_intermediate)
  
  ## Remove strain specific variants
  csc_filt <- remove_strain_specific_variants(csc_gt)
  # Detectable contamination from musculus strains in CAST cells:
  p.CAST_cont <- plot_contAllele_fraction(csc_gt, contaminating_strains = "SvImJ_BL6", quantile_thresh = 0.01)
  
  ## Estimate contamination fraction per cell: 
  # CAST cells and CAST specific SNPs only, no mitochondrial genes
  perCell_noMito_CAST_binom <- csc_filt %>% 
    filter(!grepl("mt-", SYMBOL), Strain == "CAST", SNPisALT == "CAST") %>% 
    calculate_cont_perCell_binom()
  saveRDS(perCell_noMito_CAST_binom,paste0(outdir, "/perCell_noMito_CAST_binom.RDS"))
  
  # No mitochondrial genes, all Cells
  perCell_noMito_binom <- csc_filt %>% 
    filter(!grepl("mt-", SYMBOL)) %>% 
    calculate_cont_perCell_binom()
  saveRDS(perCell_noMito_binom, paste0(outdir,"/perCell_noMito_allCells_binom.RDS"))
  
  ## some diagnostic plots:
  p.checkup <- plot_genotype_estimation(perCell_noMito_binom, perCell_noMito_CAST_binom)
  p.summary <- plot_grid(p.CAST_cont, p.checkup, nrow = 1, rel_widths = c(1,3))
  ggsave(paste0(outdir, "/plots.pdf"), p.summary, width=15, height=3.5)
}

genotype_estimation(args[1],args[2],args[3],args[4])
