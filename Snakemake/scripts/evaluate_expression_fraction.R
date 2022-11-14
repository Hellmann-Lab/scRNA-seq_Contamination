args <- commandArgs(TRUE)

library(tidyverse)

evaluate_expression_fraction <- function(corrected_dir, seurat, expression_fraction_PT_markers){
  
  seu_uncorrected <- readRDS(seurat)
  panglao_markers <- readRDS("/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/files/panglao_markers_Mm.RDS")
  
  expr_fraction_list <- list()
  for(mat in list.files(corrected_dir, pattern = "_cormat.RDS")){
    cor_mat <- readRDS(paste0(corrected_dir,"/",mat))
    
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
  
  expr_fraction_df <- bind_rows(expr_fraction_list, .id = "param")
  
  saveRDS(expr_fraction_df, expression_fraction_PT_markers)
} 


evaluate_expression_fraction(args[1],args[2],args[3])
