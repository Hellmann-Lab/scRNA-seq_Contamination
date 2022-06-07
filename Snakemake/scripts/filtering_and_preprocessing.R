args <- commandArgs(TRUE)

library(tidyseurat)
library(tidyverse)

filtering_and_preprocessing <- function(cellranger_dir, vireo_out, celltype_assignment, genotype_estimation_dir, cluster_res, seurat_out, seurat_CAST_out, raw_matrix){
  mat <- Read10X_h5(paste0(cellranger_dir,"/filtered_feature_bc_matrix.h5"))
  donor_id <- read.table(paste0(vireo_out, "/donor_ids.tsv"), header = T)
  celltype_assignment <- readRDS(celltype_assignment)
  genotype_estimation <- readRDS(paste0(genotype_estimation_dir,"/perCell_noMito_CAST_binom.RDS"))
  cluster_res <- as.numeric(strsplit(cluster_res, split = ",")[[1]])
  
  ### Prepare anntation
  # reformat strain assignment:
  strain_assignment <- donor_id %>% 
    mutate(Strain = case_when(donor_id == "C57BL_6NJ" ~ "BL6",
                              donor_id == "129S1_SvImJ" ~ "SvImJ",
                              donor_id == "CAST_EiJ" ~ "CAST",
                              T ~ donor_id)) %>% 
    filter(!Strain %in% c("unassigned", "doublet")) %>% 
    select(cell, Strain)
  
  # reformat celltype assignment: 
  ct <- celltype_assignment %>% 
    mutate(celltype = case_when(celltype %in% c("CD_IC","CD_IC_A","CD_IC_B") ~ "CD_IC",
                                #celltype %in% c("aLOH","dLOH") ~ "LOH",
                                celltype %in% c("T","NK","B","MPH") ~ "Immune",
                                T ~ celltype)) %>% 
    filter(!is.na(celltype))
  
  # reformat genotype estimates:
  background <- genotype_estimation %>% 
    dplyr::rename("bRNA" = "contPerCell_binom")
  
  # combine into one table
  cell_metadata <- ct %>% inner_join(strain_assignment) %>% left_join(background) %>% 
    column_to_rownames("cell")
  
  ### Filter and process count matrix
  # filtering cells & genes
  mat_filt <- mat[!grepl("mt-", rownames(mat)), rownames(cell_metadata)]
  
  # basic processing
  seu <- CreateSeuratObject(mat_filt, meta.data = cell_metadata)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, nfeatures = 5000)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
  seu <- RunUMAP(seu, dims = 1:30)
  seu <- FindNeighbors(seu, dims = 1:30)
  
  # Clustering at different resolutions
  for(res in cluster_res[!is.na(cluster_res)]){
    seu <- FindClusters(seu, resolution = res)
  }
  
  # CAST cells only
  seu_CAST <- filter(seu, Strain == "CAST")
  seu_CAST <- NormalizeData(seu_CAST)
  seu_CAST <- FindVariableFeatures(seu_CAST, nfeatures = 5000)
  seu_CAST <- ScaleData(seu_CAST)
  seu_CAST <- RunPCA(seu_CAST)
  seu_CAST <- RunUMAP(seu_CAST, dims = 1:30)
  
  saveRDS(seu, seurat_out)
  saveRDS(seu_CAST, seurat_CAST_out)
  saveRDS(mat_filt, file = raw_matrix)
}

filtering_and_preprocessing(args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8])
