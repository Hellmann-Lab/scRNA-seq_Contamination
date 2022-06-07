# silhouette evaluation
args <- commandArgs(TRUE)

library(tidyverse)
library(cluster)
library(tidyseurat)
library(ClusterR)

cluster_evaluation <- function(corrected_dir, seurat, cluster_evaluation_external, cluster_evaluation_internal){
  # read in annotated seurat object for celltype labels
  seu <- readRDS(seurat)
  # only consider cell types with sufficient number of cells
  seu_filt <- filter(seu, 
                     celltype %in% c("PT", "DCT", "aLOH"),
                     Strain == "CAST")
  
  
  # read output of different parameter settings
  seu_list_corrected <- list()
  for(mat in list.files(corrected_dir, pattern = "_cormat.RDS")){
    cor_mat <- readRDS(paste0(corrected_dir,"/",mat))
    filt_mat <- as.matrix(cor_mat)[,colnames(seu_filt)]
    seu <- CreateSeuratObject(counts = filt_mat, meta.data = seu_filt@meta.data)
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    #seu <- ScaleData(seu) #only performs scaling on top 2000 variable features
    seu <- ScaleData(seu,features = rownames(seu)) #for all features
    seu <- RunPCA(seu, features = VariableFeatures(object = seu), verbose = F)
    # more PCs
    seu <- FindNeighbors(seu, dims = 1:30)
    seu <- FindClusters(seu, resolution = 1)
    seu_list_corrected[[gsub("_cormat.RDS","",mat)]] <- seu
  }
  
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
  
  saveRDS(external_eval, cluster_evaluation_external)
  
  
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

  
  saveRDS(internal_eval, cluster_evaluation_internal)
}
  
cluster_evaluation(args[1],args[2],args[3], args[4])
  