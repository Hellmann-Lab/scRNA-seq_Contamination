# OPTPARSE ----------------------------------------------------------------

#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("--sce"), type="character",
              help="Input Single Cell Experiment with full path.", 
              metavar="character"),
  make_option(c("--eset"), type="character",
              help="Input SCDC QC filtering with ESet and figure with full path.", 
              metavar="character"),
  make_option(c("--single"), type="logical",
              help="Has the reference ESet one sample or more?", 
              metavar="character",
              default = TRUE),
  make_option(c("--goodID"), type="character",
              help="Input Object of subsampled good cell barcode ID.", 
              metavar="character",
              default = NULL),
  make_option(c("--contID"), type="character",
              help="Input Object of subsampled contamination barcode ID.", 
              metavar="character",
              default = NULL),
  make_option(c("--emptyID"), type="character",
              help="Input Object of subsampled empty droplet barcode ID.", 
              metavar="character",
              default = NULL),
  make_option(c("--experiment"), type="character",
              help="Experiment Run.", 
              metavar="character"),
  make_option(c("--assay"), type="character",
              help="Assay Type.", 
              metavar="character"),
  make_option(c("--pseudobulk"), type="character",
              help="Pseudobulk Type.", 
              metavar="character"),
  make_option(c("--outpath"), type="character", 
              help="full output path", 
              metavar="character"),
  make_option(c("--outname"), type="character", 
              help="project name used for output rds file name", 
              metavar="character",
              default = "results")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

# opt <- list("sce" = "/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/deconvolution/input/CZI_kidney_mouse_10XChromium_sc_allstrains_sce.rds",
# "eset" = "/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/deconvolution/input/CZI_kidney_mouse_10XChromium_sc_allstrains_1_endo_eset_qc.rds",
#             "goodID" = "/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/deconvolution/input/CZI_kidney_mouse_10XChromium_sc_allstrains_ID_good.rds",
#             "emptyID" = "/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/deconvolution/input/CZI_kidney_mouse_10XChromium_sc_allstrains_ID_empty.rds",
#             "experiment" = "1",
#             "assay" = "endo",
#             "pseudobulk" = "all",
#             "outpath" = "/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/deconvolution/output/scdc/",
#             "outname" = "CZI_kidney_mouse_10XChromium_sc_allstrains")


# SETTINGS ----------------------------------------------------------------

library(methods)
library(tidyverse)#, lib.loc = "/data/home/vieth/R/x86_64-pc-linux-gnu-library/4.0")
library(SingleCellExperiment)#, lib.loc = "/data/home/vieth/R/x86_64-pc-linux-gnu-library/4.0")
library(SCDC)
library(RColorBrewer)


# INPUT -------------------------------------------------------------------

sce <- readRDS(file = opt$sce)
eset <- readRDS(file = opt$eset)
ref_eset <- eset$sc.eset.qc

goodID <- readRDS(file = opt$goodID)
contID <- readRDS(file = opt$contID)
emptyID <- readRDS(file = opt$emptyID)

# PSEUDOBULK --------------------------------------------------------------

# select the assay and experiment
if(opt$assay %in% c('endo')){
  if(opt$pseudobulk=='all'){
    keep <- colData(sce)[,"cell_assignment_new"] == "endo" 
  } else {
    keep <- colData(sce)[,"cell_assignment_new"] == "endo" & colData(sce)[,"experiment"] == opt$experiment
  }
  sce.red <- subset(sce, , keep)
}


if(opt$assay %in% c('cont')){
  if(opt$pseudobulk=='all'){
    keep <- colData(sce)[,"cell_assignment_new"] == "cont" 
  } else {
    keep <- colData(sce)[,"cell_assignment_new"] == "cont" & colData(sce)[,"experiment"] == opt$experiment
  }
  sce.red <- subset(sce, , keep)
}

if(opt$assay %in% c('empty')){
  keep <- colData(sce)[,"cell_assignment_new"] == "empty" & colData(sce)[,"experiment"] == opt$experiment
  sce.red <- subset(sce, , keep)
  colData(sce.red)[,"celltype"] <- ifelse(is.na(colData(sce.red)[,"celltype"]), 'empty', 'unknown')
}

# combine counts and annotation to ExpressionSet
cnts <- as.matrix(assay(sce.red, opt$assay))
colnames(cnts) <- colnames(sce.red)
rownames(cnts) <- rownames(sce.red)
pdata <- cbind(cellname = colnames(cnts), 
               cluster = colData(sce.red)[,"celltype"], 
               experiment = colData(sce.red)[,"experiment"])
fdata <- rownames(cnts)
deconv.eset <- getESET(cnts, fdata = fdata, pdata = pdata)
ct <- unique(phenoData(deconv.eset)$cluster) 

# generate pseudobulk samples ...

# ... using all cells per experiment
if(opt$pseudobulk == 'all'){
  pseudo.eset <- generateBulk_allcells(deconv.eset, 
                                       ct.varname = "cluster", 
                                       sample = "experiment", 
                                       ct.sub = ct)
  pseudo_eset <- pseudo.eset$pseudo_eset
  true_p <- pseudo.eset$truep
  ct_ref <- unique(phenoData(ref_eset)$cluster) 
}

# ... using random sampling of cells per experiment
if(opt$pseudobulk == 'random'){
  pseudo.eset <- generateBulk_norep(deconv.eset, 
                                    ct.varname = "cluster", 
                                    sample = "experiment", 
                                    ct.sub = ct,
                                    nbulk = 25)
  pseudo_eset <- pseudo.eset$pseudo_eset
  true_p <- pseudo.eset$true_p
  ct_ref <- unique(phenoData(ref_eset)$cluster) 
}

# ... using subsampling of cells per experiment
if(opt$pseudobulk == 'subsample'){
  if(opt$assay %in% c('endo')){
    
    pseudo.eset.L <- lapply(1:length(goodID), function(i){
      print(i)
      deconv.eset.tmp <- deconv.eset[, colnames(deconv.eset) %in% goodID[[i]]]
      pseudo.eset.tmp <- generateBulk_allcells(deconv.eset.tmp, 
                                               ct.varname = "cluster", 
                                               sample = "experiment", 
                                               ct.sub = ct)
      colnames(pseudo.eset.tmp$pseudo_eset) <- paste(opt$experiment, as.character(i), sep = "_")
      pseudo.eset.tmp$true_p <- t(as.matrix(pseudo.eset.tmp$truep))
      rownames(pseudo.eset.tmp$true_p) <- paste(opt$experiment, as.character(i), sep = "_")
      
      res <- list('assayData' = assayData(pseudo.eset.tmp$pseudo_eset)[['exprs']],
                  'true_p' = pseudo.eset.tmp$true_p)
      res
    })
    pseudobulk.L <-  lapply(1:length(goodID), function(i) {pseudo.eset.L[[i]]$assayData})
    pseudobulk <- do.call('cbind', pseudobulk.L)
    true_p.L <-  lapply(1:length(goodID), function(i) {pseudo.eset.L[[i]]$true_p})
    true_p <- do.call('rbind', true_p.L)
    
    pdata <- data.frame(sample = colnames(pseudobulk), 
                        disease = NA, 
                        stringsAsFactors = F)
    fdata <- rownames(pseudobulk)
    pseudo_eset <- getESET(pseudobulk, 
                           fdata = fdata, 
                           pdata = pdata)
    ct_ref <- unique(phenoData(ref_eset)$cluster) 
    
  }
  if(opt$assay %in% c('cont')){
    
    pseudo.eset.L <- lapply(1:length(contID), function(i){
      print(i)
      deconv.eset.tmp <- deconv.eset[, colnames(deconv.eset) %in% contID[[i]]]
      pseudo.eset.tmp <- generateBulk_allcells(deconv.eset.tmp, 
                                               ct.varname = "cluster", 
                                               sample = "experiment", 
                                               ct.sub = ct)
      colnames(pseudo.eset.tmp$pseudo_eset) <- paste(opt$experiment, as.character(i), sep = "_")
      pseudo.eset.tmp$true_p <- t(as.matrix(pseudo.eset.tmp$truep))
      rownames(pseudo.eset.tmp$true_p) <- paste(opt$experiment, as.character(i), sep = "_")
      
      res <- list('assayData' = assayData(pseudo.eset.tmp$pseudo_eset)[['exprs']],
                  'true_p' = pseudo.eset.tmp$true_p)
      res
    })
    pseudobulk.L <-  lapply(1:length(contID), function(i) {pseudo.eset.L[[i]]$assayData})
    pseudobulk <- do.call('cbind', pseudobulk.L)
    true_p.L <-  lapply(1:length(contID), function(i) {pseudo.eset.L[[i]]$true_p})
    true_p <- do.call('rbind', true_p.L)
    
    pdata <- data.frame(sample = colnames(pseudobulk), 
                        disease = NA, 
                        stringsAsFactors = F)
    fdata <- rownames(pseudobulk)
    pseudo_eset <- getESET(pseudobulk, 
                           fdata = fdata, 
                           pdata = pdata)
    ct_ref <- unique(phenoData(ref_eset)$cluster) 
    
  }
  if(opt$assay %in% c('empty')){
    
    pseudo.eset.L <- lapply(names(emptyID), function(i){
      print(i)
      deconv.eset.tmp <- deconv.eset[, colnames(deconv.eset) %in% emptyID[[i]]]
      pseudo.eset.tmp <- generateBulk_allcells(deconv.eset.tmp, 
                                               ct.varname = "cluster", 
                                               sample = "experiment", 
                                               ct.sub = ct)
      colnames(pseudo.eset.tmp$pseudo_eset) <- paste(opt$experiment, 
                                                     as.character(i), sep = "_")
      res <- list('assayData' = assayData(pseudo.eset.tmp$pseudo_eset)[['exprs']])
      res
    })
    pseudobulk.L <-  lapply(1:length(emptyID), function(i) {
      pseudo.eset.L[[i]]$assayData
    })
    pseudobulk <- do.call('cbind', pseudobulk.L)

    pdata <- data.frame(sample = colnames(pseudobulk), 
                        disease = NA, 
                        stringsAsFactors = F)
    fdata <- rownames(pseudobulk)
    pseudo_eset <- getESET(pseudobulk,
                           fdata = fdata,
                           pdata = pdata)
    true_p <- NULL
    ct_ref <- unique(phenoData(ref_eset)$cluster) 
  }
}


# setting the true porportions for empty and contaminated reads to NULL
# this might be more appropriate if we think that we cannot know the true identity of contamination just like with droplets
if(any(opt$assay %in% c("empty", "cont"))){
  true_p <- NULL
}

# save the input for SCDC
saveRDS(object = pseudo_eset, file = paste0(opt$outpath, 
                                            opt$outname, "_",
                                            opt$experiment, "_",
                                            opt$assay, "_",
                                            opt$pseudobulk, "_",
                                            "pseudobulk_eset.rds"))
if(!is.null(true_p)){
  saveRDS(object = true_p, file = paste0(opt$outpath, 
                                         opt$outname, "_",
                                         opt$experiment, "_",
                                         opt$assay, "_",
                                         opt$pseudobulk, "_",
                                         "pseudobulk_truep.rds"))
}



# SCDC DECONVOLUTION ------------------------------------------------------

if(isTRUE(opt$single)){
  deconv_res <-  SCDC_prop_ONE(bulk.eset = pseudo_eset, 
                               sc.eset = ref_eset, 
                               ct.varname = "cluster", 
                               sample = "experiment", 
                               ct.sub = ct_ref,
                               truep = true_p)
  
}
if(isFALSE(opt$single)){
  deconv_res <-  SCDC_prop(bulk.eset = pseudo_eset, 
                           sc.eset = ref_eset, 
                           ct.varname = "cluster", 
                           sample = "experiment", 
                           ct.sub = ct_ref,
                           truep = true_p)
  
}

saveRDS(object = deconv_res, file = paste0(opt$outpath, 
                                            opt$outname, "_",
                                            opt$experiment, "_",
                                            opt$assay, "_",
                                            opt$pseudobulk, "_",
                                            "pseudobulk_deconv.rds"))


# FIN ---------------------------------------------------------------------

sessionInfo()

gc()

q()


