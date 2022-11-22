
# collection of helper functions for deconvolution_scdc.R
library(tidyverse)
library(cowplot)
library(SingleCellExperiment)
library(SCDC)
library(RColorBrewer)


#calculate cell type frequency
ct_freq<-function(cell_assignment_df, type ="allstrains"){
  
  if (type == "CAST"){
    cell_assignment_df<-cell_assignment_df %>%
      filter(cell_assignment_new =="endo") 
  }
  #need to filter for good_cell 
  goodcell_distinct<-cell_assignment_df %>%
    filter(cell_assignment == "good_cell") %>%
    #dplyr::group_by(experiment) %>%
    dplyr::mutate(total=length(BC)) %>%
    dplyr::group_by(celltype) %>%
    dplyr::mutate(n=length(celltype)) %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(prop=n/total) %>%
    distinct(celltype, prop) #%>%
    # mutate(above1perc=ifelse(prop > 0.007, T, F),
    #        combi=paste0(experiment,"_",celltype))
  
  return(goodcell_distinct)
}


# new cell assignment: change "good cell" to "endo" or "cont", filter for the 3 relevant categories
load_files<-function(path.files=path.files,BCxGene.path, assignment.path, numb.exp, freq_cutoff = 0.007,
                     type = "allstrains"){
  
  BCxGene <- readRDS(paste0(path.files,BCxGene.path))
  BCs<-unlist(lapply(BCxGene, colnames), use.names = F)
  
  cell_assignment<- readRDS(paste0(path.files,assignment.path))
  
  #cell_freqs<-ct_freq(cell_assignment, type = type)
  
  cell_assignment<-cell_assignment %>%
    #this one is being used only in the case of CAST subset
    mutate(cell_assignment_new=case_when(BC %in% colnames(BCxGene$endo) ~ "endo",
                                         BC %in% colnames(BCxGene$cont) ~ "cont",
                                         T ~ cell_assignment)) %>%
    left_join(ct_freq(., type = type)) %>%
    filter(cell_assignment =="empty" | prop > freq_cutoff) %>%
    filter(BC %in% BCs)
  
  BCxGene<-lapply(BCxGene, function(x){
    #x<-x[colnames(x) %in% cell_assignment$BC]
    colnames(x)<-paste(colnames(x),numb.exp,sep="_")
    return(x)
  })

  return(list(BCxGene=BCxGene,
              cell_assignment=cell_assignment,
              BCs=BCs))
  
}
#cell_assignment_ex1 <- readRDS(paste0(path.files,"for_deconvolution/cell_assignment_ex1.RDS")) %>%
#   mutate(endo = BC %in% colnames(experiment1_BCxGene$endo),
#         cont = BC %in% colnames(experiment1_BCxGene$cont))


#https://stackoverflow.com/questions/43117608/r-binding-sparse-matrices-of-different-sizes-on-rows
merge.sparse = function(listMatrixes, allRownames=sharedgenes, allColnames=cell_assignment$ID) {
  # takes a list of sparse matrixes with different columns and adds them row wise
  suppressWarnings(rm(matrixToReturn))
  for (currentMatrix in listMatrixes) {
    
    #filter relevant cells
    currentMatrix=currentMatrix[,colnames(currentMatrix) %in% allColnames]

    #add the missing rows
    misR=allRownames[!allRownames %in% rownames(currentMatrix)]
    misRmat=Matrix::Matrix(matrix(0,
                                nrow=length(misR),
                                ncol=ncol(currentMatrix),
                                dimnames = list(misR,colnames(currentMatrix))))
    newMatrix = do.call(rbind2, c(currentMatrix, misRmat))
    newMatrix = newMatrix[match(allRownames, rownames(newMatrix)),]
    
    if (!exists("matrixToReturn")) {
      matrixToReturn <- newMatrix
    }
    else {
      matrixToReturn <- cbind2(matrixToReturn,newMatrix)
    }}
  matrixToReturn  
}



empty_matrix<-function(example_mat){
  Matrix::Matrix(matrix(0, nrow=nrow(example_mat), ncol=ncol(example_mat), dimnames = dimnames(example_mat)))
}


#generate ESET
generate_ESET<-function(sce.red, type){
  cnts <- as.matrix(assay(sce.red, type))
  colnames(cnts) <- colnames(sce.red)
  rownames(cnts) <- rownames(sce.red)
  cnts[is.na(cnts)]<-0
  fdata <- rownames(cnts)
  pdata <- cbind(cellname = colnames(cnts), 
                 cluster = colData(sce.red)[,"celltype"], 
                 experiment = colData(sce.red)[,"experiment"])
  eset <- getESET(cnts, fdata = fdata, pdata = pdata)
  return(eset)
}



#subsampling of cells for pseudobulk
subsample_cells<-function(sce_subset, rand, ndraw){
  
  samp_no <- colData(sce_subset) %>% 
    data.frame() %>% 
    group_by(experiment, celltype) %>% 
    dplyr::summarise(n = n()) %>% 
    mutate(rno = floor(n * rand)) %>% 
    dplyr::select(-n)
  
  ID_sampling<- lapply(1:ndraw, function(i){
    BC_samp <- colData(sce_subset) %>% 
      data.frame() %>% 
      group_by(experiment, celltype) %>% 
      nest() %>%  
      ungroup() %>% 
      dplyr::left_join(samp_no, by = c('celltype' = "celltype", "experiment" = "experiment")) %>% 
      mutate(samp = map2(data, rno, sample_n)) %>% 
      dplyr::select(-data) %>%
      unnest(samp) %>% 
      pull(ID)
    BC_samp
  })
  return(ID_sampling)
}


subsample_empty<-function(sce_subset, rand, ndraw){
  # samp_no <- colData(sce_subset) %>% 
  #   data.frame() %>% 
  #   group_by(experiment) %>% 
  #   dplyr::summarise(n = n()) %>% 
  #   mutate(rno = floor(n * rand)) %>% 
  #   dplyr::select(-n)
  
  ID_sampling_empty <- lapply(rand, function(i){
    samp_no <- colData(sce_subset) %>% 
      data.frame() %>% 
      group_by(experiment) %>% 
      dplyr::summarise(n = n()) %>% 
      mutate(rno = floor(n * i)) %>% 
      dplyr::select(-n)
    IDs <- lapply(1:ndraw, function(j){
      BC_samp <- colData(sce_subset) %>% 
        data.frame() %>% 
        group_by(experiment) %>% 
        nest() %>%  
        ungroup() %>% 
        dplyr::left_join(samp_no, by = c("experiment" = "experiment")) %>% 
        mutate(samp = map2(data, rno, sample_n)) %>% 
        dplyr::select(-data) %>%
        unnest(samp) %>% 
        pull(ID)
      BC_samp
    })
    names(IDs) <- paste(i, 1:ndraw, sep = "_")
    IDs
  })
  return(ID_sampling_empty)
}



slide_theme <- function (base_size = 11, base_family = "") {
  theme_light(base_size = base_size, base_family = base_family) +
    theme(plot.background = element_blank(), #element_rect(colour = "black", fill=NA, size=1.25)
          panel.grid.major = element_line(colour = "grey75"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey25", size = 1),
          strip.background = element_rect(fill = NA, colour = NA),
          strip.text = element_text(colour = "black", size = rel(1.25), face = 'bold'),
          plot.title = element_text(size = rel(1.25), face = 'bold', colour = 'black',
                                    margin = margin(t = 5, r = 0,  b = 15, l = 0)),
          axis.text = element_text(colour = "black", size = rel(1.25)),
          axis.title = element_text(colour = "black", face="bold", size = rel(1.25)),
          axis.title.y = element_text(margin = margin(t = 0, r = 15,  b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15,  r = 0, b = 0, l = 0)),
          legend.title = element_text(colour = "black", size = rel(1), face = "bold"),
          legend.key.size = unit(0.75, "lines"),
          legend.text = element_text(size = rel(0.9), colour = "black"),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.background = element_rect(colour = NA, fill = NA))
}

stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data = fun, geom = geom, width = 0.1, ...)
}



