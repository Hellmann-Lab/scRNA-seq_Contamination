args <- commandArgs(TRUE)

library(tidyverse)
library(Matrix)

ENSEMBL_to_SYMBOL <- readRDS("input/ENSEMBL_to_SYMBOL.RDS")

calculate_kNN <- function(mat, ndim = 20, k = 20){
  mat <- as.matrix(mat)
  mat <- mat[rowSums(mat) > 0,]
  #mat <- scater::logNormCounts(mat)
  red_dat <- BiocSingular::runPCA(t(mat), rank = ndim, get.rotation = FALSE)$x
  kNN <- BiocNeighbors::findAnnoy(red_dat, k = k, warn.ties = FALSE)$index
  return(kNN)
}

match_matrices <- function(mat,reference_mat){mat[rownames(mat) %in% rownames(reference_mat), colnames(reference_mat)]}

calculate_kNN_overlap <- function(ref_kNN, test_kNN){
  mean(sapply(seq_len(nrow(test_kNN)), function(cell_idx){
    length(intersect(ref_kNN[cell_idx,], test_kNN[cell_idx,]))}))
}

evaluate_kNN_overlap <- function(corrected_dir, reference_kNN, kNN_overlap){
  corrected_dir <- gsub("/dummy.txt","",corrected_dir)
  ref <- readRDS(reference_kNN)
  kNN_list <- list()
  for(mat in list.files(corrected_dir, pattern = "_cormat.RDS")){
    cor_mat <- readRDS(paste0(corrected_dir,"/",mat))
    
    kNN <- calculate_kNN(match_matrices(cor_mat, ref$mat),ndim = 30, k = 50)
    
    ovlp <- calculate_kNN_overlap(ref$ref_kNN, kNN)
    
    kNN_list[[gsub("_cormat.RDS","",mat)]] <- ovlp
  }
  
  kNN_overlap_df <- bind_rows(kNN_list, .id = "param")
  
  saveRDS(kNN_overlap_df, kNN_overlap)
}
  
evaluate_kNN_overlap(args[1],args[2],args[3])
      
