# run DecontX 
args <- commandArgs(TRUE)

suppressPackageStartupMessages(
  library(celda, lib.loc = "/data/home/janssen/R/x86_64-pc-linux-gnu-library/4.1/") 
)
  
run_DecontX <- function(cellranger, seurat, cluster_res, use_empty, cormat, perCell){
  # read seurat object as input
  seu <- readRDS(seurat)
  
  # Clustering
  if(cluster_res == "Default"){
    z = NULL
  } else {
    z = seu[[paste0("RNA_snn_res.", cluster_res)]][,1]
  }
  
  # Background profile
  if(use_empty == "True"){
    print("Using empty droplet profile")
    raw <- Seurat::Read10X_h5(paste0(cellranger, "/raw_feature_bc_matrix.h5"))
    raw <- raw[rownames(seu),]
  } else {
    print(paste0("use_empty = ",use_empty, ": empty droplets not used"))
    raw = NULL
  }
  
  # Remove background
  res <- decontX(x = as.matrix(seu@assays$RNA@counts),
                 z = z,
                 background = raw)
  
  perCell_cont <- data.frame(cell = colnames(seu), cont = res$contamination)
  
  saveRDS(res$decontXcounts, cormat)
  saveRDS(perCell_cont, perCell)
}

run_DecontX(args[1],args[2],args[3],args[4],args[5],args[6])