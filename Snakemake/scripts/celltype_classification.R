args <- commandArgs(TRUE)

### cell type classification with pre-trained reference
suppressPackageStartupMessages({
  library(SingleR)
  library(Seurat)
  library(tidyverse)
})


celltype_classification <- function(cellranger_dir, cell_BC_path, vireo_out, prediction_results, celltype_assignment, UMAP_celltype){
  # load data: cellranger count matrix, cell calling results, strain assignments:
  mat <- Read10X_h5(paste0(cellranger_dir,"/filtered_feature_bc_matrix.h5"))
  cell_BC <- read.table(cell_BC_path)$V1
  donor_id <- read.table(paste0(vireo_out, "/donor_ids.tsv"), header = T)
  
  # filter for cells with strain assignment, that met cell calling criteria
  donor_id_filt <- donor_id %>% 
    filter(cell %in% cell_BC,
           !(donor_id %in% c("doublet", "unassigned"))) %>% 
    column_to_rownames("cell")
  
  # cell filtering
  mat_filt <- mat[,colnames(mat) %in% rownames(donor_id_filt)]
  
  # Normalization and basic processing:
  print("############ Pre-processing with Seurat ###############")
  seu <- CreateSeuratObject(mat_filt, meta.data = donor_id_filt)
  seu <- NormalizeData(seu)
  seu <- ScaleData(seu, verbose = F)
  seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = F)
  seu <- RunPCA(seu, verbose = F)
  seu <- RunUMAP(seu, dims = 1:30, verbose = F)
  
  # classify cells with pre-trained reference
  print("############ Reference based classification ###############")
  train_denisenko <- readRDS("/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/files/reference_data/trainSingleR_denisenko.RDS")
  
  pred <- classifySingleR(test = GetAssayData(seu, "data"), 
                          trained = train_denisenko, 
                          fine.tune = F, BPPARAM=BiocParallel::MulticoreParam(8))
  
  ct <- data.frame(cell = colnames(seu), celltype = pred$pruned.labels)
  
  saveRDS(pred, prediction_results)
  saveRDS(ct, celltype_assignment)
  
  # Add results to Seurat object
  seu$celltype <- pred$pruned.labels
  p.UMAP <- DimPlot(seu, group.by = c("donor_id", "celltype"))
  ggsave(UMAP_celltype, p.UMAP, width = 10, height = 3.5)
}

celltype_classification(args[1],args[2],args[3],args[4],args[5],args[6])