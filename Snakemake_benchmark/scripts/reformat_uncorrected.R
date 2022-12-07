args <- commandArgs(TRUE)

reformat_uncorrected <- function(seurat, cormat){
  # read seurat object
  seu <- readRDS(seurat)
  # save unocorrected count matrix 
  saveRDS(seu@assays$RNA@counts, cormat)
}

reformat_uncorrected(args[1],args[2])