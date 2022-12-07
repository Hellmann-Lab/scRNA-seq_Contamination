# run CellBender
args <- commandArgs(TRUE)

# required inputs: 
#   raw cellranger matrix
#   estimated number of cells
#   total droplets included (default 25000)

# parameters: 
#   number of epochs
#   false positive rate - higher value = more background removed

run_CellBender <- function(cellranger, seurat,total_droplets_included, fpr, CB_out,  cormat, perCell){
  
  outdir <- gsub("cormat.h5","",CB_out)
  print(outdir)
  seu <- readRDS(seurat)
  cells <- colnames(seu)
  genes <- rownames(seu)
  n_cells <- length(cells)
  
  system(paste("mkdir -p ", outdir, " ; sbatch --wait scripts/run_CellBender.sh",cellranger, CB_out, n_cells, total_droplets_included, fpr), wait = T)
  
  ### import old result, delete later #####
  #cormat_old <- readRDS(CB_out)
  #cormat_out <- cormat_old[genes,cells]
  
  ############
  
  # save output
  cormat_out <- Seurat::Read10X_h5(CB_out)[genes,cells]
  
  # calculate contamination per cell
  uncor <- Seurat::Read10X(paste0(cellranger, "/filtered_feature_bc_matrix"))[genes,cells] 
  perCell_cont <- data.frame(cell = colnames(cormat_out),
                             cont = 1-(Matrix::colSums(cormat_out) / Matrix::colSums(uncor)))
  
  saveRDS(cormat_out, cormat)
  saveRDS(perCell_cont, perCell)
}
  
run_CellBender(args[1],args[2],args[3],args[4],args[5],args[6],args[7])

# run_CellBender(cellranger = "/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/method_comparison/input/rep1/cellranger",
#                CB_out = "/home/janssen/scribble/cellbender/fpr_cormat.h5",
#                cell_BCs = "/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/method_comparison/input/rep1/cell_BC.txt",
#                total_droplets_included = 25000,
#                fpr = 0.01)

