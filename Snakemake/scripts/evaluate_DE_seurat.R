args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(tidyseurat)
})

evaluate_DE_seurat <- function(corrected_dir, seurat, differential_expression_seurat){
  seu_uncorrected <- readRDS(seurat)
  
  DE_list <- list()
  for(mat in list.files(corrected_dir, pattern = "_cormat.RDS")){
    cor_mat <- readRDS(paste0(corrected_dir,"/",mat))
    
    print("performing DE analysis with Seurat")
    seu <- CreateSeuratObject(counts = cor_mat, meta.data = seu_uncorrected@meta.data)
    seu <- filter(seu, Strain == "CAST")
    seu <- NormalizeData(seu)
    seu$inf <- ifelse(seu$celltype == "PT", "PT", "other")
    PT.markers <- FindMarkers(seu, group.by = "inf", ident.1 = "PT", logfc.threshold = 0)
    
    DE_list[[gsub("_cormat.RDS","",mat)]] <- PT.markers
  }
  
  DE_df <- bind_rows(DE_list, .id = "param")
  saveRDS(DE_df, differential_expression_seurat)
}

evaluate_DE_seurat(args[1],args[2],args[3])