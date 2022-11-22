
# GOALS -------------------------------------------------------------------

#' empty droplet analysis focusing on deconvolution
#' extended to deconvolution benchmarking of endogenous good cells, contamination profile of cells AND empty droplets


# SETTINGS ----------------------------------------------------------------

library(tidyverse)
library(cowplot)
# library(Seurat)
library(SingleCellExperiment)


# INPUT -------------------------------------------------------------------

# path to input files
workdir<-"/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/deconvolution_new/"
setwd(workdir)
source("scripts/Rscript/helper_functions.R")

# files
# count matrices and cell assignment

exp1_files<-load_files(workdir, BCxGene.path = "raw_inputs_mat/rep1_BCxGene_matrices_CAST.RDS",
                       assignment.path = "raw_inputs/cell_assignment_rep1.RDS", 
                       numb.exp = 1,freq_cutoff = 0.0075, type = "CAST")

exp2_files<-load_files(workdir, BCxGene.path = "raw_inputs_mat/rep2_BCxGene_matrices_CAST.RDS",
                       assignment.path = "raw_inputs/cell_assignment_rep2.RDS", 
                       numb.exp = 2,freq_cutoff = 0.0075, type = "CAST")

exp3_files<-load_files(workdir, BCxGene.path = "raw_inputs_mat/rep3_BCxGene_matrices_CAST.RDS",
                       assignment.path = "raw_inputs/cell_assignment_rep3.RDS", 
                       numb.exp = 3,freq_cutoff = 0.0075, type = "CAST")

nuc2_files<-load_files(workdir, BCxGene.path = "raw_inputs_mat/nuc2_BCxGene_matrices_CAST.RDS",
                       assignment.path = "raw_inputs/cell_assignment_nuc2.RDS", 
                       numb.exp = 4,freq_cutoff = 0.0075, type = "CAST")

nuc3_files<-load_files(workdir, BCxGene.path = "raw_inputs_mat/nuc3_BCxGene_matrices_CAST.RDS",
                       assignment.path = "raw_inputs/cell_assignment_nuc3.RDS", 
                       numb.exp = 5,freq_cutoff = 0.0075, type = "CAST")


# STEPS -------------------------------------------------------------------

# 1. generate input files
# 2. generate reference in SCDC using good endo cells
# 3. run deconvolution of pseudobulk per experiment based on ...
# a. endogenous expression profile: evaluate quality of reference dataset by generating pseudobulk consisting of cells with known proportions, 
# - all cells pseudobulk, random subsampling, random sampling of 50% per cell type
# b. contaminated expression profile: quantify the contribution of ambient/barcode swapping to the droplets with cells captured
# - all contaminated in one pseudobulk, random subsampling, random sampling of 50% per cell type
# c. empty droplets: quantify the contribution of cell types to the ambient profile in droplets w/o cells
# - all empty droplets in one pseudobulk, random samplings of 25%, 50%, 75% from all empty droplets 

# i noticed while combining the deconvolution results that the names of cell types are not consistent (see make.names function), and adapted them in the input generation!
# i added the number of umi/cell and detected genes/cell to the singlecellexperiment, given that they are now computed per assay type (endogenous, contaminated) they differ from nUMI

# 1 Generate Input -----------------------------------------------------------------------

sharedgenes <- Reduce(union, list(rownames(exp1_files$BCxGene$endo),
                                  rownames(exp2_files$BCxGene$endo),
                                  rownames(exp3_files$BCxGene$endo),
                                  rownames(nuc2_files$BCxGene$endo),
                                  rownames(nuc3_files$BCxGene$endo),
                                  rownames(exp1_files$BCxGene$cont),
                                  rownames(exp2_files$BCxGene$cont),
                                  rownames(exp3_files$BCxGene$cont),
                                  rownames(nuc2_files$BCxGene$cont),
                                  rownames(nuc3_files$BCxGene$cont),
                                  rownames(exp1_files$BCxGene$empty),
                                  rownames(exp2_files$BCxGene$empty),
                                  rownames(exp3_files$BCxGene$empty),
                                  rownames(nuc2_files$BCxGene$empty),
                                  rownames(nuc3_files$BCxGene$empty)))

cell_assignment <- data.table::rbindlist(list('1'=exp1_files$cell_assignment,
                                              '2'=exp2_files$cell_assignment,
                                              '3'=exp3_files$cell_assignment,
                                              '4'=nuc2_files$cell_assignment,
                                              '5'=nuc3_files$cell_assignment), 
                                         idcol = 'experiment') %>%
  mutate(ID = paste(BC,experiment,sep="_"))



#let's make the matrix per type like this to speed it up a bit
sparse_endo<-merge.sparse(listMatrixes = list(exp1_files$BCxGene$endo,exp2_files$BCxGene$endo,exp3_files$BCxGene$endo,nuc2_files$BCxGene$endo,nuc3_files$BCxGene$endo), 
                          allRownames=sharedgenes, 
                          allColnames=cell_assignment$ID)
sparse_cont<-merge.sparse(listMatrixes = list(exp1_files$BCxGene$cont,exp2_files$BCxGene$cont,exp3_files$BCxGene$cont,nuc2_files$BCxGene$cont,nuc3_files$BCxGene$cont), 
                          allRownames=sharedgenes, 
                          allColnames=cell_assignment$ID)
sparse_empty<-merge.sparse(listMatrixes = list(exp1_files$BCxGene$empty,exp2_files$BCxGene$empty,exp3_files$BCxGene$empty,nuc2_files$BCxGene$empty,nuc3_files$BCxGene$empty), 
                           allRownames=sharedgenes, 
                           allColnames=cell_assignment$ID)

#need to add empty "padding" to each of the matrices because they have to be of the same length.. not sure though why even generate 3 separate matrices? maybe there is some convenience aspect..
empty_endo<-empty_matrix(sparse_endo)
empty_cont<-empty_matrix(sparse_cont)
empty_empty<-empty_matrix(sparse_empty)

endo = cbind2(sparse_endo, empty_cont) %>% cbind2(empty_empty)
cont = cbind2(empty_endo, sparse_cont) %>% cbind2(empty_empty)
empty = cbind2(empty_endo, empty_cont) %>% cbind2(sparse_empty)

CAnnot <- DataFrame(cell_assignment, row.names = cell_assignment$ID)
CAnnot<-CAnnot[match(colnames(endo),CAnnot$ID),]


#this same empty matrix is set for endo, cont and empty and then filled in afterwards by the real numbers
sce_all <- SingleCellExperiment::SingleCellExperiment(list(endo = endo,cont = cont,empty = empty),
                                                      colData=CAnnot,
                                                      metadata=list(study="CZI_kidney", 
                                                                    organism = "mouse",
                                                                    strain = "all", 
                                                                    protocol = "singlecell"))


# add cell qc metrics
# for endogenous counts
endo.stats <- scuttle::perCellQCMetrics(sce_all,
                                        assay.type = 'endo',
                                        percent_top = NULL,
                                        flatten = TRUE)
colnames(endo.stats) <- paste('endo', colnames(endo.stats), sep = '_')
# for contaminated fraction
cont.stats <- scuttle::perCellQCMetrics(sce_all,
                                        assay.type = 'cont',
                                        percent_top = NULL,
                                        flatten = TRUE)
colnames(cont.stats) <- paste('cont', colnames(cont.stats), sep = '_')

# fraction of cont/endo 
cell.stats <- cbind(endo.stats, cont.stats) %>% 
  data.frame() %>% 
  dplyr::mutate(endofrac_sum = endo_sum/(cont_sum+endo_sum),
                endofrac_detected = endo_detected/(cont_detected+endo_detected),
                contfrac_sum = cont_sum/(cont_sum+endo_sum),
                contfrac_detected = cont_detected/(cont_detected+endo_detected)) %>% 
  dplyr::select(-endo_total, -cont_total)

colData(sce_all) <- cbind(colData(sce_all), cell.stats)
saveRDS(sce_all, file = paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_sce.rds"))


# 2 Generate Single Cell Reference -----------------------------------------------------------------------

library(SCDC)
library(RColorBrewer)

sce_all <- readRDS(file = paste0(workdir, "input/CZI_kidney_mouse_10XChromium_sc_CAST_sce.rds"))

# all experiments together
keep <- (colData(sce_all)[,"cell_assignment_new"] == "endo")
sce.red <- subset(sce_all, , keep)
eset<-generate_ESET(sce.red, type="endo")
saveRDS(object = eset, file = paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_endo_eset.rds"))

ct <- unique(phenoData(eset)$cluster)
n <- length(ct)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ct2 <- sample(col_vector, n)


#assignment overview
demo.p <- DemoPlot(eset, 
                   cluster = "cluster", 
                   sample = "experiment", 
                   select.ct = ct, 
                   Palette = ct2)

saveRDS(object = demo.p , file = paste0(workdir, "input/CZI_kidney_mouse_10XChromium_sc_CAST_endo_demoplot.rds"))


#To make SCDC robust to single-cell clustering, a quality control procedure is performed as a first step to remove cells with questionable cell-type assignments, as well as cells with low library preparation and sequencing quality. Specifically, each single cell is treated as a "bulk" sample and its cell-type composition can be derived by a first-pass run of SCDC. The threshold of keeping a single cell from an assigned cluster is chosen as 0.7 in this illustration example. This user-defined threshold should be selected carefully, to achieve the effect of quality control of single cells without obstructing the performance on the downstream analysis. An unproper threshold could potentially lead to filtering out a large proportion of cells from rare cell types. To be conservative, this can be chosen from (0.5,0.7).

# We have now less cells to choose from (only non-CAST as endo).. I might change qc treshold to 0.5

eset.qc <- SCDC_qc(eset, 
                   ct.varname = "cluster", 
                   sample = "experiment", 
                   scsetname = 'endo_reference',
                   ct.sub = ct,
                   qcthreshold = 0.6) 

#i guess this is the number of cells that made it past the quality threshold
length(colnames(eset.qc$sc.eset.qc)) #23461

saveRDS(object = eset.qc, file = paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_endo_eset_qc.rds"))

# extract ids of retained cells for reference and add them to the singlecellexperiment object
qc.id <- DataFrame('SCDC_Reference' = colnames(sce_all) %in% colnames(eset.qc$sc.eset.qc) ,
                   row.names = colnames(sce_all))
# append estimates to sce
sum(qc.id$SCDC_Reference)
colData(sce_all)[, 'SCDC_Reference'] <- qc.id$SCDC_Reference

saveRDS(sce_all, file = paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_sce.rds"))


# per experiment
# Similarly, we can perform clustering QC for single cells from only one subject using function SCDC_qc_ONE().
# hmmm, but why would we?
#forgot to change saved rds file names - next time instead of "allstrains" put CAST

lapply(unique(sce_all$experiment), function(i){
  print(i)
  keep <- (colData(sce_all)[,"cell_assignment_new"] == "endo" & colData(sce_all)[,"experiment"] == i)
  print(sum(keep))
  sce.red <- subset(sce_all, , keep)
  eset<-generate_ESET(sce.red, "endo")
  
  saveRDS(object = eset, file = paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_", i, "_endo_eset.rds"))
  
  ct <- unique(phenoData(eset)$cluster)
  n <- length(ct)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  ct2 <- sample(col_vector, n)
  
  demo.p <- DemoPlot(eset, 
                     cluster = "cluster", 
                     sample = "experiment", 
                     select.ct = ct, 
                     Palette = ct2)
  
  saveRDS(object = demo.p , file = paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_", i, "_endo_demoplot.rds"))
  
  eset.qc <- SCDC_qc_ONE(eset, 
                         ct.varname = "cluster", 
                         sample = "experiment", 
                         scsetname = 'endo_reference',
                         ct.sub = ct,
                         qcthreshold = 0.6) 
  
  saveRDS(object = eset.qc, file = paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_", i,"_endo_eset_qc.rds"))
})



# 3 Deconvolute pseudobulk -----------------------------------------------------------------------

sce_all <- readRDS(file = paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_sce.rds"))

# construct pseudobulk of 25 random subsamplings of cells at 50% of its original cell numbers per cell type
# for endogenous counts of good cells and contaminated counts of good cells

sce_good_endo <- subset(sce_all, , (cell_assignment_new == "endo"))# 
ID_sampling_good<-subsample_cells(sce_subset=sce_good_endo, rand=0.5, ndraw=25)

saveRDS(ID_sampling_good, file = paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_ID_good.rds"))


#now we need the cont profiles
sce_good_cont <- subset(sce_all, , (cell_assignment_new == "cont"))
ID_sampling_cont<-subsample_cells(sce_subset=sce_good_cont, rand=0.5, ndraw=25)

saveRDS(ID_sampling_cont, file = paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_ID_cont.rds"))


# construct pseudobulk of 25 random subsamplings of empty droplets at 25%, 50%, 75% and once 100%
sce_empty <- subset(sce_all, , cell_assignment_new == "empty")
ID_sampling_empty<-subsample_empty(sce_subset=sce_empty, rand=c(0.25,0.5,0.75), ndraw=25)
ID_sampling_empty <- do.call('c', ID_sampling_empty)
ID_sampling_empty[["1_1"]] <- colnames(sce_empty)

saveRDS(ID_sampling_empty, file = paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_ID_empty.rds"))




# PREPARE FOR SLURM -------------------------------------------------------------------

# input sce
input.sce <- paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_sce.rds")

# input qc-filtered ESET SCDC reference
# over all experiments
input.eset <- paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_endo_eset_qc.rds")
# per experiment 
input.esets <- c(paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_1_endo_eset_qc.rds"), 
                 paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_2_endo_eset_qc.rds"), 
                 paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_3_endo_eset_qc.rds"),
                 paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_4_endo_eset_qc.rds"),
                 paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_5_endo_eset_qc.rds"))


# input ID sampling objects
input.good.ID <- paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_ID_good.rds")
input.cont.ID <- paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_ID_cont.rds")
input.empty.ID <- paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_ID_empty.rds")

single <- TRUE

# experiment
experiment <- c("1", "2", "3", "4", "5")
# gene expression 
assay <- c('endo', 'cont', 'empty')
# pseudobulk
pseudobulk <- c('all', 'random', 'subsample')

# output path
outpath <- paste0(workdir,"output/scdc/single_ref/")
# output name
outname <- "CZI_kidney_mouse_10XChromium_sc_CAST"


# parameter combination
param.slurm <- base::expand.grid(experiment = experiment,
                                 assay = assay, 
                                 pseudobulk = pseudobulk) %>% 
  dplyr::filter(!c(assay == 'empty' & pseudobulk %in% c('all', 'random')))

slurm.input <- cbind(sce = input.sce,
                     eset = input.esets,
                     single = single,
                     goodID = input.good.ID, 
                     contID = input.cont.ID, 
                     emptyID = input.empty.ID,
                     param.slurm,
                     outpath = outpath,
                     outname = outname)

write.table(format(slurm.input, scientific=FALSE), 
            file =paste0(workdir,"scripts/bash/slurm_deconv_CAST_singleref.txt"),
            sep = "\t", row.names = T, quote = F, col.names = F, na='NULL')


# input sce
input.sce <- paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_sce.rds")

# input qc-filtered ESET SCDC reference
# over all experiments
input.eset <- paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_endo_eset_qc.rds")

# input ID sampling objects
input.good.ID <- paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_ID_good.rds")
input.cont.ID <- paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_ID_cont.rds")
input.empty.ID <- paste0(workdir,"input/CZI_kidney_mouse_10XChromium_sc_CAST_ID_empty.rds")

single <- FALSE

# experiment
experiment <- c("1", "2", "3", "4", "5")
# gene expression 
assay <- c('endo', 'cont', 'empty')
# pseudobulk
pseudobulk <- c('all', 'random', 'subsample')

# output path
outpath <- paste0(workdir,"output/scdc/multi_ref/")
# output name
outname <- "CZI_kidney_mouse_10XChromium_sc_CAST"


# parameter combination
param.slurm <- base::expand.grid(experiment = experiment,
                                 assay = assay, 
                                 pseudobulk = pseudobulk) %>% 
  dplyr::filter(!c(assay == 'empty' & pseudobulk %in% c('all', 'random')))

slurm.input <- cbind(sce = input.sce,
                     eset = input.eset,
                     single = single,
                     goodID = input.good.ID,
                     contID = input.cont.ID,
                     emptyID = input.empty.ID,
                     param.slurm,
                     outpath = outpath,
                     outname = outname)

write.table(format(slurm.input, scientific=FALSE), 
            file =paste0(workdir,"scripts/bash/slurm_deconv_CAST_multiref.txt"),
            sep = "\t", row.names = T, quote = F, col.names = F, na='NULL')
