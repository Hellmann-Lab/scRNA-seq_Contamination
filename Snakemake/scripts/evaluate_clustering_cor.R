# external and internal cluster evaluation
args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(cluster)
  library(tidyseurat)
  library(ClusterR)
  library(SingleR)
})

cluster_evaluation <- function(corrected_dir, seurat,classification_evaluation, cluster_evaluation_external_cor, cluster_evaluation_internal_cor, silhouette_per_cell_cor){
  # read in classification result for cell type labels
  print("reading in classification result")
  class <- readRDS(classification_evaluation)
  
  print("performing clustering")
  seu_list_corrected <- list()
  for(param in names(class)){
    celltype_annotation <- column_to_rownames(class[[param]]$celltype_labels, "cell")
    
    cor_mat <- readRDS(paste0(corrected_dir,"/",param,"_cormat.RDS"))
    filt_mat <- as.matrix(cor_mat)[,rownames(celltype_annotation)]
    seu <- CreateSeuratObject(counts = filt_mat, meta.data = celltype_annotation)
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    #seu <- ScaleData(seu) #only performs scaling on top 2000 variable features
    seu <- ScaleData(seu,features = rownames(seu)) #for all features
    seu <- RunPCA(seu, features = VariableFeatures(object = seu), verbose = F)
    # more PCs
    seu <- FindNeighbors(seu, dims = 1:30)
    seu <- FindClusters(seu, resolution = 1)
    seu_list_corrected[[param]] <- seu
  }
  
  # read in seurat object for uncorrected cell type labels in order to filter cell types to be comparable across methods
  seu_uncor <- readRDS(seurat)
  
  # only consider cell types with sufficient number of cells in CAST (>10) in uncorrected data
  celltypes_keep <- seu_uncor %>% 
    filter(Strain == "CAST") %>% 
    group_by(celltype) %>% 
    summarise(n_cells = n()) %>% 
    filter(n_cells > 10) %>% 
    pull(celltype)
  
  # filter corrected seurat objects
  seu_list_corrected <- lapply(seu_list_corrected, function(x){
    filter(x, celltype %in% celltypes_keep)
  })
  
  
  #### EXTERNAL CLUSTER EVALUATION ####
  external_eval <- lapply(seu_list_corrected, function(seu){
    ext_df <- capture.output(ClusterR::external_validation(true_labels = as.numeric(as.factor(seu$celltype)),
                                                           clusters = as.numeric(seu$seurat_clusters), 
                                                           summary_stats = TRUE)) %>% 
      strsplit(split = ":") %>% 
      .[lapply(., length) == 2] %>% 
      Reduce(rbind, .) %>% 
      apply(2, trimws) %>% 
      data.frame(row.names = NULL) %>% 
      dplyr::rename("metric" = "X1", "value" = "X2") %>% 
      mutate(value = as.numeric(value))
    
    return(ext_df)
  }) %>% bind_rows(.id = "param")
  
  saveRDS(external_eval, cluster_evaluation_external_cor)
  
  
  #### INTERNAL CLUSTER EVALUATION ####
  internal_eval <- lapply(seu_list_corrected, function(seu){
    # silhouette on 30 PCs and euclidean distance
    dist <- dist(seu@reductions$pca@cell.embeddings[,1:30], method = "euclidean")
    sil <- summary(silhouette(as.numeric(as.factor(seu$celltype)), dist))
    sil_df <- data.frame(metric = paste0("silhouette_",levels(as.factor(seu$celltype))),
                         value = sil$clus.avg.widths) %>% 
      add_row(metric = "silhouette", value = sil$avg.width)
    return(sil_df)
  }) %>% bind_rows(.id = "param")

  
  sil_per_cell <- lapply(seu_list_corrected, function(seu){
    # silhouette on 30 PCs and euclidean distance
    dist <- dist(seu@reductions$pca@cell.embeddings[,1:30], method = "euclidean")
    sil_allCells <- silhouette(as.numeric(as.factor(seu$celltype)), dist)
    sil_allCells_df <- data.frame(cell = colnames(seu),
                                  celltype = seu$celltype,
                                  sil_width = as.data.frame.matrix(sil_allCells)$sil_width)
    return(sil_allCells_df)
  })
  
  saveRDS(internal_eval, cluster_evaluation_internal_cor)
  saveRDS(sil_per_cell, silhouette_per_cell_cor)
  
}
  
cluster_evaluation(args[1],args[2],args[3], args[4], args[5], args[6])
  