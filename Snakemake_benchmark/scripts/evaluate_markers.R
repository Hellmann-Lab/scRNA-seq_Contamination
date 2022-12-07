args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(tidyseurat)
})

evaluate_markers <- function(corrected_dir, seurat, differential_expression_seurat, expression_fraction_PT_markers){
  seu_uncorrected <- readRDS(seurat)
  corrected_dir <- gsub("/dummy.txt","",corrected_dir)
  panglao_markers <- readRDS("input/panglao_markers_Mm.RDS")
  
  DE_list <- list()
  expr_fraction_list <- list()
  for(mat in list.files(corrected_dir, pattern = "_cormat.RDS")){
    cor_mat <- readRDS(paste0(corrected_dir,"/",mat))
    
    ### DIFFERENTIAL EXPRESSION WITH SEURAT #####
    print("performing DE analysis with Seurat")
    seu <- CreateSeuratObject(counts = cor_mat, meta.data = seu_uncorrected@meta.data)
    seu <- filter(seu, Strain == "CAST")
    seu <- NormalizeData(seu)
    seu$inf <- ifelse(seu$celltype == "PT", "PT", "other")
    PT.markers <- FindMarkers(seu, group.by = "inf", ident.1 = "PT", logfc.threshold = 0)
    
    DE_list[[gsub("_cormat.RDS","",mat)]] <- PT.markers
    
    ### EXPRESSION FRACTION ####
    cnts_round <- cor_mat[panglao_markers, ] %>% round
    if(sum(colnames(cnts_round) != colnames(seu_uncorrected)) == 0){
      df_percent <- data.frame(gene = rownames(cnts_round),
                               expr_within_celltype = Matrix::rowSums(cnts_round[,seu_uncorrected$celltype == "PT"] > 0) / ncol(cnts_round[,seu_uncorrected$celltype == "PT"]),
                               expr_within_others = Matrix::rowSums(cnts_round[,seu_uncorrected$celltype != "PT"] > 0) / ncol(cnts_round[,seu_uncorrected$celltype != "PT"]))
    } else {
      stop("corrected matrix doesn't match original input matrix")
    }
    expr_fraction_list[[gsub("_cormat.RDS","",mat)]] <- df_percent
  }
  
  DE_df <- bind_rows(DE_list, .id = "param")
  saveRDS(DE_df, differential_expression_seurat)
  
  expr_fraction_df <- bind_rows(expr_fraction_list, .id = "param")
  saveRDS(expr_fraction_df, expression_fraction_PT_markers)
}

evaluate_markers(args[1],args[2],args[3], args[4])