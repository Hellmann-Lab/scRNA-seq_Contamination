args <- commandArgs(TRUE)

library(tidyverse)
library(scran)
library(SingleCellExperiment)
library(limma)

evaluate_DE <- function(corrected_dir, seurat, differential_expression){
  seu_uncorrected <- readRDS(seurat)
  
  DE_list <- list()
  for(mat in list.files(corrected_dir, pattern = "_cormat.RDS")){
    cor_mat <- readRDS(paste0(corrected_dir,"/",mat))
    
    print("performing normalization")
    sce <- SingleCellExperiment(assays = list(counts = cor_mat), colData = seu_uncorrected@meta.data)
    sce <- sce[,sce$Strain == "CAST"]
    sce <- sce[rowSums(counts(sce) > 0) > 5, colSums(counts(sce)) > 0]
    clusters <- quickCluster(sce, min.size=50)
    sce <- computeSumFactors(sce, cluster=clusters)
    sce <- logNormCounts(sce)
    cnts_norm <- logcounts(sce)
    sf <- calculateSumFactors(sce, cluster=clusters)
    
    print("performing DE analysis with limma")
    inf <- as.factor(ifelse(sce$celltype == "PT", "PT", "other"))
    # calculate normalisation factors
    nsf <- log(sf/Matrix::colSums(counts(sce)))
    nsf <- exp(nsf - mean(nsf, na.rm=T))
    # construct input object
    dge <- edgeR::DGEList(counts = counts(sce),
                          lib.size = Matrix::colSums(counts(sce)),
                          norm.factors = nsf,
                          group = factor(inf),
                          remove.zeros = FALSE)
    
    design.mat <- stats::model.matrix(~inf)
    y <- new("EList")
    y$E <- edgeR::cpm(dge, log = T)
    fit <- limma::lmFit(object = y, design = design.mat)
    fit <- eBayes(fit, trend = T, robust = T)
    tt <- topTable(fit, n = Inf, adjust.method = "BH", confint = T) %>% 
      rownames_to_column("gene")
    
    
    DE_list[[gsub("_cormat.RDS","",mat)]] <- tt
  }
  
  DE_df <- bind_rows(DE_list, .id = "param")
  saveRDS(DE_df, differential_expression)
}



evaluate_DE(args[1],args[2],args[3])
