# run SoupX
args <- commandArgs(TRUE)

library(SoupX)

run_SoupX <- function(cellranger, seurat, cluster_res, set_cont, cormat, perCell){
  # read input
  seu <- readRDS(seurat)
  cells <- colnames(seu)

  if(!(paste0("RNA_snn_res.", cluster_res) %in% colnames(seu@meta.data))){
    stop("SoupX requires pre-computed clustering!")
  }
  
  sc <- load10X(cellranger, cellIDs = cells) 
  
  # define clusters
  cl = seu[[paste0("RNA_snn_res.", cluster_res)]][,1]
  sc <- setClusters(sc, cl)
  
  # filter genes to match seurat object
  sc$toc <- sc$toc[rownames(sc$toc) %in% rownames(seu),]
  sc$soupProfile <- sc$soupProfile[rownames(sc$soupProfile) %in% rownames(seu),]
  
  # estimate rho
  if(set_cont == "Auto"){
    sc <- autoEstCont(sc)
  } else {
    sc <- setContaminationFraction(sc, contFrac = as.numeric(set_cont))
  }
  
  # correct counts
  corrected_counts <- adjustCounts(sc)
  
  perCell_cont <- data.frame(cell = colnames(corrected_counts),
                             cont = 1-(Matrix::colSums(corrected_counts) / Matrix::colSums(sc$toc)))
  
  saveRDS(corrected_counts, cormat)
  saveRDS(perCell_cont, perCell)
}

run_SoupX(args[1],args[2],args[3],args[4],args[5],args[6])

