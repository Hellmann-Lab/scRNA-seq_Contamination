
library(tidyverse)
library(cowplot)
# library(Seurat)
library(SingleCellExperiment)



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

setwd("/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/deconvolution_new/")
workdir<-"/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/deconvolution_new/"

#COLORS

# cell type colors
celltype_colors <- setNames(
  c(RColorBrewer::brewer.pal(11, "Paired"),"darkgrey","grey"),
  c("PT","CD_IC" ,   "CD_PC", "CD_Trans" ,     "CNT" ,     "DCT",     "Endo",      "Fib","aLOH",   "dLOH",    "MC" ,    "Podo","Immune"))

# RESULTS  ---------------------------------------------------

res.deconv.files <- list.files(path = paste0(workdir,"output/scdc/single_ref/"), pattern = "pseudobulk_deconv.rds")

res.truep.files <- list.files(path = paste0(workdir,"output/scdc/single_ref/"), pattern = "pseudobulk_truep.rds")

# enodegenuous counts
res.deconv.endo.files <- res.deconv.files[grepl(pattern = "endo", res.deconv.files)]
names(res.deconv.endo.files) <- gsub(pattern = "_pseudobulk_deconv.rds", replacement = "", x = res.deconv.endo.files)
res.truep.endo.files <- res.truep.files[grepl(pattern = "endo", res.truep.files)]
names(res.truep.endo.files) <- gsub(pattern = "_pseudobulk_truep.rds", replacement = "", x = res.truep.endo.files)

# contaminated fraction
res.deconv.cont.files <- res.deconv.files[grepl(pattern = "cont", res.deconv.files)]
names(res.deconv.cont.files) <- gsub(pattern = "_pseudobulk_deconv.rds", replacement = "", x = res.deconv.cont.files)

# empty droplets 
res.deconv.empty.files <- res.deconv.files[grepl(pattern = "empty", res.deconv.files)]
names(res.deconv.empty.files) <- gsub(pattern = "_pseudobulk_deconv.rds", replacement = "", x = res.deconv.empty.files)

# reference cells per experiment
qc.eset.files <- list.files(path = paste0(workdir,"input/"), pattern = "_endo_eset_qc.rds")
names(qc.eset.files) <- c('1', '2', '3', '4', '5', 'all')

# number and proportion of cells per experiment in reference
qc.bc.L <- sapply(names(qc.eset.files), function(i){
  tmp <- readRDS(file = paste0(workdir,"input/", qc.eset.files[i]))
  data.frame('BC' = tmp$sc.eset.qc$cellname,
             'celltype' = tmp$sc.eset.qc$cluster,
             stringsAsFactors = F) %>% 
    dplyr::group_by(celltype) %>% 
    dplyr::count(name = 'cellnumber') %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(refproportion = cellnumber/sum(cellnumber))
}, USE.NAMES = TRUE, simplify = FALSE)

res.qc <- data.table::rbindlist(qc.bc.L, idcol = "Experiment") %>% 
  dplyr::filter(Experiment %in% c('1', '2', '3', '4', '5'))

saveRDS(res.qc, file=paste0("output/scdc/summaries/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_deconvref.rds"))

# number of cells in full set and subsampling of good cells
sce_all <- readRDS(file = paste0("input/CZI_kidney_mouse_10XChromium_sc_CAST_sce.rds"))
sce_good <- subset(sce_all, , cell_assignment_new == "endo")
rand <- 0.5

samp_no_endo_good <- colData(sce_good) %>% 
  data.frame() %>% 
  group_by(experiment, celltype) %>% 
  dplyr::summarise(no = n()) %>% 
  dplyr::ungroup() %>% 
  mutate(rno.50 = floor(no * rand)) %>% 
  dplyr::group_by(experiment) %>% 
  mutate(prop = no/sum(no), prop.50 = rno.50/sum(rno.50)) 

saveRDS(samp_no_endo_good, file=paste0(workdir,"output/scdc/summaries/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_sampling.rds"))

# cell qc metrics per experiment and cell type
good_cell_qc_stats <-  colData(sce_all) %>% 
  data.frame() %>% 
  dplyr::filter(cell_assignment_new == "endo") %>% 
  dplyr::select(experiment, celltype, 
                endo_sum, endo_detected, 
                cont_sum, cont_detected,  
                endofrac_sum, endofrac_detected,
                contfrac_sum, contfrac_detected) %>% 
  dplyr::group_by(experiment, celltype) %>% 
  dplyr::summarise(across(contains("_"), 
                          list(mean = mean, se=~sd(.)/sqrt(n())), .names = "{.col}.{.fn}"))

saveRDS(good_cell_qc_stats, file=paste0(workdir,"output/scdc/summaries/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_goodcell_qc.rds"))

ref_cell_qc_stats <-  colData(sce_all) %>% 
  data.frame() %>% 
  dplyr::filter(cell_assignment_new == "endo" & SCDC_Reference == TRUE) %>% 
  dplyr::select(experiment, celltype, 
                endo_sum, endo_detected, 
                cont_sum, cont_detected,  
                endofrac_sum, endofrac_detected,
                contfrac_sum, contfrac_detected) %>% 
  dplyr::group_by(experiment, celltype) %>% 
  dplyr::summarise(across(contains("_"), 
                          list(mean = mean, se=~sd(.)/sqrt(n())), .names = "{.col}.{.fn}"))

saveRDS(ref_cell_qc_stats, file=paste0(workdir,"output/scdc/summaries/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_refcell_qc.rds"))

# RESULTS ENDO ------------------------------------------------------------

res.deconv.endo.L <- sapply(names(res.deconv.endo.files), function(i){
  tmp <- readRDS(file = paste0(workdir,"output/scdc/single_ref/", res.deconv.endo.files[i]))
  dat <- tmp$prop.est.mvw %>% 
    data.frame() %>% 
    tibble::rownames_to_column(var='experiment') %>% 
    tidyr::pivot_longer(-experiment, names_to = "celltype", values_to = "EstimatedProportion")
  dat
}, USE.NAMES = TRUE, simplify = FALSE)

res.deconv.endo.dat <- data.table::rbindlist(res.deconv.endo.L, idcol = "setup") %>% 
  tidyr::separate(col="experiment", into = c("experiment", "run"), fill = "right", sep = "_")

res.eval.endo.L <- sapply(names(res.deconv.endo.files), function(i){
  tmp <- readRDS(file = paste0(workdir,"output/scdc/single_ref/", res.deconv.endo.files[i]))
  dat <- data.frame(RMSD_bysample = tmp$peval$evals[[1]]$RMSD_bysample,
                    mAD_bysample = tmp$peval$evals[[1]]$mAD_bysample,
                    Pearson_bysample = tmp$peval$evals[[1]]$Pearson_bysample) %>% 
    tibble::rownames_to_column(var='experiment')
  dat
}, USE.NAMES = TRUE, simplify = FALSE)

res.eval.endo.dat <- data.table::rbindlist(res.eval.endo.L, idcol = "setup") %>% 
  tidyr::separate(col="experiment", into = c("experiment", "run"), fill = "right", sep = "_")

res.truep.endo.L <- sapply(names(res.truep.endo.files), function(i){
  print(i)
  tmp <- readRDS(file = paste0(workdir,"output/scdc/single_ref/", res.truep.endo.files[i]))
  dat <- data.frame(tmp, stringsAsFactors = F)
  
  if(any(colnames(dat) %in% "sample.id")){
    dat <- dat   %>% 
      dplyr::rename(experiment = sample.id, 
                    celltype = cluster.id, 
                    TrueProportion = Freq)
  }else{
    dat<- dat %>% 
      tibble::rownames_to_column(var = "experiment") %>% 
      tidyr::pivot_longer(-experiment, names_to = "celltype", values_to = "TrueProportion")
  }
  dat$experiment <- as.character(dat$experiment)
  dat$celltype <- as.character(dat$celltype)
  
  dat
}, USE.NAMES = TRUE, simplify = FALSE)

res.truep.endo.dat <- data.table::rbindlist(res.truep.endo.L, idcol = "setup") %>% 
  tidyr::separate(col="experiment", into = c("experiment", "run"), fill = "right", sep = "_")

# combine true and estimated proportions
res.endo.dat <- dplyr::left_join(res.deconv.endo.dat,
                                 res.truep.endo.dat,
                                 by = c("setup" = "setup",
                                        "experiment" = "experiment",
                                        "run" = "run",
                                        "celltype"="celltype"))
# kick out wrong combinations of 'all'
# Explanation, SDCDC does not work with only one bulk sample to deconvolute so I just made either enough random pseudobulk samples or considered all experiments for pseudobulk samples.
res.endo.dat <-
  res.endo.dat %>% 
  tidyr::separate(setup, into = c("Project", "Tissue", "Species", "RNAseq", "Unit", "Strain", "Experiment", "Read", "Pseudobulk"), sep = "_", remove = FALSE) %>% 
  dplyr::mutate(Match = case_when(Experiment == experiment ~ 'correct',
                                  TRUE ~ 'kickout')) %>% 
  dplyr::filter(Match == 'correct') %>% 
  dplyr::select(-Match, -experiment)

res.eval.endo.dat <-
  res.eval.endo.dat %>% 
  tidyr::separate(setup, into = c("Project", "Tissue", "Species", "RNAseq", "Unit", "Strain", "Experiment", "Read", "Pseudobulk"), sep = "_", remove = FALSE) %>% 
  dplyr::mutate(Match = case_when(Experiment == experiment ~ 'correct',
                                  TRUE ~ 'kickout')) %>% 
  dplyr::filter(Match == 'correct') %>% 
  dplyr::select(-Match, -experiment)


# INES ####
saveRDS(res.endo.dat, file=paste0(workdir,"output/scdc/summaries/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_endo_prop.rds"))
saveRDS(res.eval.endo.dat, file=paste0(workdir,"output/scdc/summaries/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_endo_eval.rds"))

# table of evaluation metrics
res.eval.endo.table <- res.eval.endo.dat %>% 
  dplyr::group_by(setup) %>% 
  dplyr::summarise(RMSD = mean(RMSD_bysample),
                   SE.RMSD = plotrix::std.error(RMSD_bysample),
                   mAD = mean(mAD_bysample),
                   SE.mAD = plotrix::std.error(mAD_bysample),
                   Pearson = mean(Pearson_bysample),
                   SE.Pearson = plotrix::std.error(Pearson_bysample)) %>% 
  ungroup() %>% 
  tidyr::separate(setup, into = c("Project", "Tissue", "Species", "RNAseq", "Unit", "Strain", "Experiment", "Read", "Pseudobulk"), sep = "_", remove = TRUE) %>% 
  dplyr::select(Experiment, Pseudobulk, RMSD,SE.RMSD,mAD,SE.mAD,Pearson,SE.Pearson) %>% 
  data.frame()
res.eval.endo.ptable <- dplyr::select(res.eval.endo.table, Experiment, Pseudobulk, RMSD,mAD,Pearson)

# plotting
# New facet label names for pseudobulk variable
pseudobulk.labs <- c("all cells", "25 times random sampling", "25 times random sampling\n of downsampled cells") 
names(pseudobulk.labs) <- c("all", "random", "subsample")
experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3", "Experiment 4", "Experiment 5") 
names(experiment.labs) <- c("1", "2", "3", "4", "5")

# number of cells, true proportion, reference proportion per experiment
nocells.barplot <- res.endo.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion, TrueProportion) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  dplyr::left_join(samp_no_endo_good, 
                   by = c('Experiment' = 'experiment', 'celltype' = 'celltype')) %>% 
  dplyr::filter(Pseudobulk == 'all') %>% 
  ggplot(data = ., aes(x = reorder(celltype, no), y = no, fill = Experiment)) +
  geom_bar(stat = 'identity', color = 'black', position = 'dodge') +
  coord_flip() +
  scale_y_log10( labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides='b') +
  labs(x = NULL,
       y = NULL, 
       title  = 'Number of Cells per Experiment') +
  slide_theme(base_size = 10) +
  theme(legend.position = 'bottom')
nocells.barplot

allpropcells.plot <- res.endo.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion, TrueProportion) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  dplyr::left_join(samp_no_endo_good, 
                   by = c('Experiment' = 'experiment', 'celltype' = 'celltype')) %>% 
  dplyr::filter(Pseudobulk == 'all') %>% 
  ggplot(data = ., aes(x = reorder(celltype, no), y = prop, fill = Experiment)) +
  geom_bar(stat = 'identity', color = 'black', position = 'dodge') +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL,
       y = NULL, 
       title  = 'Proportion of Cells per Experiment') +
  slide_theme(base_size = 10) +
  theme(legend.position = 'bottom')

refpropcells.plot <- res.endo.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion, TrueProportion) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  dplyr::left_join(samp_no_endo_good, 
                   by = c('Experiment' = 'experiment', 'celltype' = 'celltype')) %>% 
  dplyr::filter(Pseudobulk == 'all') %>% 
  ggplot(data = ., aes(x = reorder(celltype, no), y = refproportion, fill = Experiment)) +
  geom_bar(stat = 'identity', color = 'black', position = 'dodge') +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL,
       y = NULL, 
       title  = 'Proportion of Cells in Deconvolution Reference per Experiment') +
  slide_theme(base_size = 10) +
  theme(legend.position = 'bottom')

qc.plot <- cowplot::plot_grid(nocells.barplot + theme(legend.position = 'none'),
                              allpropcells.plot + theme(legend.position = 'none'),
                              refpropcells.plot + theme(legend.position = 'none'),
                              cowplot::get_legend(nocells.barplot),
                              rel_heights = c(1,1,1,0.1),
                              labels = c("A", "B", "C", NULL),
                              ncol = 1)

cowplot::save_plot(paste0(workdir,'figures/SCDC/benchmark/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_cell_prop.pdf'),
                   qc.plot,
                   base_height = 12,
                   base_width = 7)
cowplot::save_plot(paste0(workdir,'figures/SCDC/benchmark/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_cell_prop.png'),
                   qc.plot,
                   base_height = 12,
                   base_width = 7)

# reference vs true proportion: scatterplot
endo.prop.plot <- res.endo.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion, TrueProportion) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  dplyr::filter(Pseudobulk == 'all') %>% 
  ggplot(data = ., aes(x = refproportion, y = TrueProportion)) +
  geom_point(aes(fill=celltype, shape = Experiment), color = 'black', size = 2, alpha = 0.75) +
  facet_wrap(~ Experiment, labeller = labeller(Experiment = experiment.labs)) +
  geom_abline(slope = 1, intercept = 0) +
  #scale_shape_manual(values=c(21,22,24)) +
  scale_fill_manual(values=celltype_colors)+
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  labs(x='Reference Proportion', y = 'Cellular Proportion', title  = 'Comparison of Proportions per Experiment',
       subtitle = 'Reference = all SCDC QC-Filtered cells; Cellular = all good cells') +
  slide_theme(base_size = 8) +
  theme(legend.position = 'bottom')
endo.prop.plot
cowplot::save_plot(paste0(workdir,'figures/SCDC/benchmark/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_endo_benchmark_eval_prop.pdf'),
                   endo.prop.plot,
                   base_height = 8,
                   base_width = 12)
cowplot::save_plot(paste0(workdir,'figures/SCDC/benchmark/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_endo_benchmark_eval_prop.png'),
                   endo.prop.plot,
                   base_height = 8,
                   base_width = 12)

# estimated vs true proportion: scatterplot
endo.eval.plot <- ggplot(data = res.endo.dat, aes(x = TrueProportion, y = EstimatedProportion)) +
  geom_point(aes(fill=celltype, shape = Experiment), color = 'black', size = 2, alpha = 0.75) +
  facet_wrap(~Pseudobulk, labeller = labeller(Pseudobulk = pseudobulk.labs)) +
  geom_abline(slope = 1, intercept = 0) +
  #scale_shape_manual(values=c(21,22,24)) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  scale_fill_manual(values=celltype_colors)+
  labs(x='True Proportion', y = 'Estimated Proportion', title  = 'SCDC Deconvolution of pseudobulk samples derived from good cells using endogenous expression counts of all strains',
       subtitle = 'Experiment-Specific Reference using SCDC QC-Filtered good cells') +
  slide_theme(base_size = 8) +
  theme(legend.position = 'bottom')

endo.eval.table <- ggpubr::ggtexttable(res.eval.endo.ptable,
                                       rows = NULL, 
                                       theme = ggpubr::ttheme("classic", base_size = 8))

endo.overview.plot <- cowplot::ggdraw(endo.eval.plot + 
                                        theme(plot.margin = unit(c(t=0, r=10, b=0, l=0), "cm"))) +
  cowplot::draw_plot(endo.eval.table, 0.39, 0.05)
endo.overview.plot
cowplot::save_plot(paste0(workdir,'figures/SCDC/benchmark/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_endo_benchmark_eval_overview.pdf'),
                   endo.overview.plot,
                   base_height = 5,
                   base_width = 18)
cowplot::save_plot(paste0(workdir,'figures/SCDC/benchmark/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_endo_benchmark_eval_overview.png'),
                   endo.overview.plot,
                   base_height = 5,
                   base_width = 18)

# estimated vs true proportion: bar plot of all cells
all.barplot <- res.endo.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion, TrueProportion) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  tidyr::pivot_longer(cols = -c(celltype, Experiment, Pseudobulk, run), names_to = "Type", values_to = 'Proportion') %>% 
  dplyr::mutate(run=ifelse(is.na(run), "1", run)) %>% 
  dplyr::mutate(Type=case_when(Type == 'EstimatedProportion' ~ 'Estimated',
                               Type == 'TrueProportion' ~ 'True',
                               Type == 'refproportion' ~ 'Reference')) %>% 
  dplyr::mutate(Type = factor(Type, levels = c("Reference", "True", "Estimated"))) %>% 
  dplyr::filter(Pseudobulk == 'all') %>% 
  ggplot(., aes(x = Type, y = Proportion, fill = celltype)) +
  geom_bar(stat = 'identity', color = 'black') +
  scale_fill_manual(values=celltype_colors)+
  coord_flip() +
  facet_wrap(~ Experiment, labeller = labeller(Experiment = experiment.labs)) +
  labs(x = NULL) +
  labs(x=NULL,
       y = 'Proportion', 
       title  = 'SCDC Deconvolution of pseudobulk samples derived from good cells using endogenous expression counts of all strains',
       subtitle = 'Using all good cells per experiment') +
  slide_theme(base_size = 10) +
  theme(legend.position = 'bottom')
all.barplot

cowplot::save_plot(paste0(workdir,'figures/SCDC/benchmark/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_endo_benchmark_eval_allcells_barplot.pdf'),
                   all.barplot,
                   base_height = 8,
                   base_width = 15)
cowplot::save_plot(paste0(workdir,'figures/SCDC/benchmark/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_endo_benchmark_eval_allcells_barplot.png'),
                   all.barplot,
                   base_height = 8,
                   base_width = 15)

#INES ####
saveRDS(res.endo.dat, "/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/Extra_analyses/SCDC_res_endo_dat.rds")
# difference between true and estimated vs number of cells: 'scatterplot'
random.error.plot <- res.endo.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion, TrueProportion) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  dplyr::left_join(samp_no_endo_good, 
                   by = c('Experiment' = 'experiment', 'celltype' = 'celltype')) %>% 
  dplyr::filter(Pseudobulk == 'random') %>% 
  dplyr::mutate(Deviance = TrueProportion-EstimatedProportion) %>% 
  ggplot(., aes(x = no, y = Deviance, color = celltype)) +
  geom_hline(yintercept = 0,linetype='dashed') +
  scale_color_manual(values=celltype_colors)+
  stat_sum_df("mean_se") +
  scale_x_log10( labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides='b') +
  labs(x='Number of cells', 
       y = 'Difference between \n True and Estimated Proportion', 
       title = 'SCDC Deconvolution of pseudobulk samples derived from good cells using endogenous expression counts of all strains',
       subtitle = 'Using 25 times random sampling of good cells per experiment') +
  facet_wrap(~Experiment, labeller = labeller(Experiment = experiment.labs)) +
  slide_theme(base_size = 10) +
  theme(legend.position = 'bottom', 
        plot.margin = unit(c(t=0, r=1, b=0, l=0), "cm"))
random.error.plot
cowplot::save_plot(paste0(workdir,'figures/SCDC/benchmark/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_endo_benchmark_eval_randomcells_meandiff.pdf'),
                   random.error.plot,
                   base_height = 8,
                   base_width = 12)
cowplot::save_plot(paste0(workdir,'figures/SCDC/benchmark/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_endo_benchmark_eval_randomcells_meandiff.png'),
                   random.error.plot,
                   base_height = 8,
                   base_width = 12)

# difference between true and estimated of all cells sorted by cell numbers: barplot
all.error.barplot <- res.endo.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion, TrueProportion) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  dplyr::filter(Pseudobulk == 'all') %>% 
  dplyr::mutate(Deviance = TrueProportion-EstimatedProportion) %>% 
  ggplot(., aes(x = fct_reorder(celltype, refproportion), y = Deviance, fill = Experiment)) +
  geom_hline(yintercept = 0,linetype='dashed') +
  geom_bar(stat="identity",position='dodge', color = 'black') +
  coord_flip() + 
  labs(x='Cell type sorted by proportion in reference', 
       y = 'Difference between True and Estimated Proportion', 
       title = 'SCDC Deconvolution of pseudobulk samples derived from good cells \nusing endogenous expression counts of all strains',
       subtitle = 'Using all good cells per experiment') +
  slide_theme(base_size = 10) +
  theme(legend.position = 'right', 
        plot.margin = unit(c(t=0, r=1, b=0, l=0), "cm"))
all.error.barplot
cowplot::save_plot(paste0(workdir,'figures/SCDC/benchmark/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_endo_benchmark_eval_allcells_celltype_meandiff.pdf'),
                   all.error.barplot,
                   base_height = 10,
                   base_width = 10)
cowplot::save_plot(paste0(workdir,'figures/SCDC/benchmark/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_endo_benchmark_eval_allcells_celltype_meandiff.png'),
                   all.error.barplot,
                   base_height = 10,
                   base_width = 10)

# estimated vs true proportion of random samplings: scatterplot facet per cellltype
endo.random.plot <- res.endo.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion, TrueProportion) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  dplyr::left_join(samp_no_endo_good, 
                   by = c('Experiment' = 'experiment', 'celltype' = 'celltype')) %>% 
  dplyr::filter(Pseudobulk == 'random') %>% 
  dplyr::mutate(Deviance = TrueProportion-EstimatedProportion) %>% 
  ggplot(data = ., aes(x = TrueProportion, y = EstimatedProportion)) +
  geom_point(aes(fill=celltype),size = 2, color = 'black', shape = 21) +
  scale_fill_manual(values=celltype_colors)+
  facet_grid(fct_reorder(celltype, no)~Experiment,scales='free', labeller = labeller(Experiment = experiment.labs)) +
  geom_abline(slope = 1, intercept = 0) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  labs(x='True Proportion', 
       y = 'Estimated Proportion', 
       title  = 'SCDC Deconvolution of pseudobulk samples derived from\n good cells using endogenous expression counts of all strains',
       subtitle = '25 times random samplings of good cells per experiment') +
  slide_theme(base_size = 8) +
  theme(legend.position = 'bottom',
        strip.text.y = element_blank())
endo.random.plot
cowplot::save_plot(paste0(workdir,'figures/SCDC/benchmark/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_endo_benchmark_eval_randomcells_scatter.pdf'),
                   endo.random.plot,
                   base_height = 12,
                   base_width = 10)
cowplot::save_plot(paste0(workdir,'figures/SCDC/benchmark/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_endo_benchmark_eval_randomcells_scatter.png'),
                   endo.random.plot,
                   base_height = 12,
                   base_width = 10)

# difference between true and estimated vs number of cells: 'scatterplot'
subsample.error.plot <- res.endo.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion, TrueProportion) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  dplyr::left_join(samp_no_endo_good, 
                   by = c('Experiment' = 'experiment', 'celltype' = 'celltype')) %>% 
  dplyr::filter(Pseudobulk == 'subsample') %>% 
  dplyr::mutate(Deviance = TrueProportion-EstimatedProportion) %>% 
  ggplot(., aes(x = rno.50, y = Deviance, color = celltype)) +
  scale_color_manual(values=celltype_colors)+
  geom_hline(yintercept = 0,linetype='dashed') +
  stat_sum_df("mean_se") +
  scale_x_log10( labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides='b') +
  labs(x='Number of cells', 
       y = 'Difference between \n True and Estimated Proportion', 
       title = 'SCDC Deconvolution of pseudobulk samples derived from good cells using endogenous expression counts of all strains',
       subtitle = '25 times random sampling of downsampled cells') +
  facet_wrap(~Experiment, labeller = labeller(Experiment = experiment.labs)) +
  slide_theme(base_size = 10) +
  theme(legend.position = 'bottom', 
        plot.margin = unit(c(t=0, r=1, b=0, l=0), "cm"))
subsample.error.plot
cowplot::save_plot(paste0(workdir,'figures/SCDC/benchmark/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_endo_benchmark_eval_subsamplecells_meandiff.pdf'),
                   subsample.error.plot,
                   base_height = 8,
                   base_width = 12)
cowplot::save_plot(paste0(workdir,'figures/SCDC/benchmark/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_endo_benchmark_eval_subsamplecells_meandiff.png'),
                   subsample.error.plot,
                   base_height = 8,
                   base_width = 12)

# RESULTS CONTAMINATION ---------------------------------------------------

res.deconv.cont.L <- sapply(names(res.deconv.cont.files), function(i){
  tmp <- readRDS(file = paste0(workdir,"output/scdc/single_ref/", res.deconv.cont.files[i]))
  dat <- tmp$prop.est.mvw %>% 
    data.frame() %>% 
    tibble::rownames_to_column(var='experiment') %>% 
    tidyr::pivot_longer(-experiment, names_to = "celltype", values_to = "EstimatedProportion")
  dat
}, USE.NAMES = TRUE, simplify = FALSE)

res.deconv.cont.dat <- data.table::rbindlist(res.deconv.cont.L, idcol = "setup") %>% 
  tidyr::separate(col="experiment", into = c("experiment", "run"), fill = "right", sep = "_")

res.eval.cont.L <- sapply(names(res.deconv.cont.files), function(i){
  tmp <- readRDS(file = paste0(workdir,"output/scdc/single_ref/", res.deconv.cont.files[i]))
  dat <- data.frame(RMSD_bysample = tmp$yeval$evals[[1]]$RMSDy_bysample,
                    mAD_bysample = tmp$yeval$evals[[1]]$mADy_bysample,
                    Spearman_bysample = tmp$yeval$evals[[1]]$spearmany_bysample) %>% 
    tibble::rownames_to_column(var='experiment')
  dat
}, USE.NAMES = TRUE, simplify = FALSE)

res.eval.cont.dat <- data.table::rbindlist(res.eval.cont.L, idcol = "setup") %>% 
  tidyr::separate(col="experiment", into = c("experiment", "run"), fill = "right", sep = "_")

# kick out wrong combinations of 'all'
# Explanation, SDCDC does not work with only one bulk sample to deconvolute so I just made either enough random pseudobulk samples or considered all experiments for pseudobulk samples.
res.cont.dat <-
  res.deconv.cont.dat %>% 
  tidyr::separate(setup, into = c("Project", "Tissue", "Species", "RNAseq", "Unit", "Strain", "Experiment", "Read", "Pseudobulk"), sep = "_", remove = FALSE) %>% 
  dplyr::mutate(Match = case_when(Experiment == experiment ~ 'correct',
                                  TRUE ~ 'kickout')) %>% 
  dplyr::filter(Match == 'correct') %>% 
  dplyr::select(-Match, -experiment)

res.eval.cont.dat <-
  res.eval.cont.dat %>% 
  tidyr::separate(setup, into = c("Project", "Tissue", "Species", "RNAseq", "Unit", "Strain", "Experiment", "Read", "Pseudobulk"), sep = "_", remove = FALSE) %>% 
  dplyr::mutate(Match = case_when(Experiment == experiment ~ 'correct',
                                  TRUE ~ 'kickout')) %>% 
  dplyr::filter(Match == 'correct') %>% 
  dplyr::select(-Match, -experiment)

saveRDS(res.cont.dat, file=paste0(workdir,"output/scdc/summaries/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_cont_prop.rds"))
saveRDS(res.eval.cont.dat, file=paste0(workdir,"output/scdc/summaries/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_cont_eval.rds"))

# table of evaluation metrics
res.eval.cont.table <- res.eval.cont.dat %>% 
  dplyr::group_by(setup) %>% 
  dplyr::summarise(RMSD = mean(RMSD_bysample),
                   SE.RMSD = plotrix::std.error(RMSD_bysample),
                   mAD = mean(mAD_bysample),
                   SE.mAD = plotrix::std.error(mAD_bysample),
                   Spearman = mean(Spearman_bysample),
                   SE.Spearman = plotrix::std.error(Spearman_bysample)) %>% 
  ungroup() %>% 
  tidyr::separate(setup, into = c("Project", "Tissue", "Species", "RNAseq", "Unit", "Strain", "Experiment", "Read", "Pseudobulk"), sep = "_", remove = TRUE) %>% 
  dplyr::select(Experiment, Pseudobulk, RMSD,SE.RMSD,mAD,SE.mAD,Spearman,SE.Spearman) %>% 
  data.frame()
res.eval.cont.ptable <- dplyr::select(res.eval.cont.table, Experiment, Pseudobulk, RMSD,mAD,Spearman)

# plotting
# New facet label names for pseudobulk variable
pseudobulk.labs <- c("all cells", "25 times random sampling", "25 times random sampling\n of downsampled cells") 
names(pseudobulk.labs) <- c("all", "random", "subsample")
experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3", "Experiment 4") 
names(experiment.labs) <- c("1", "2", "3", "4")

# contaminated vs endogenous counts qc scatterplot
qc.contno.plot <- good_cell_qc_stats %>% 
  ggplot(., aes(y = cont_sum.mean, x = endo_sum.mean)) +
  geom_errorbar(aes(color = celltype, ymin=cont_sum.mean-cont_sum.se, ymax=cont_sum.mean+cont_sum.se)) +
  geom_errorbarh(aes(color = celltype, xmin=endo_sum.mean-endo_sum.se, xmax=endo_sum.mean+endo_sum.se)) +
  geom_point(aes(fill=celltype, shape = experiment), color = 'black', size = 2) +
  scale_fill_manual(values = celltype_colors)+
  scale_color_manual(values = celltype_colors)+
  #scale_shape_manual(values=c(21,22,24)) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  labs(x='Endogenous Counts',
       y = 'Contaminated Counts', 
       title  = 'Relationship between contaminated counts and endogenous counts of all strains per cell type',
       subtitle = 'Using all good cells per experiment') +
  slide_theme(base_size = 10) +
  theme(legend.position = 'bottom')
qc.contno.exp.plot <- good_cell_qc_stats %>% 
  ggplot(., aes(y = cont_sum.mean, x = endo_sum.mean)) +
  geom_errorbar(aes(color = celltype, ymin=cont_sum.mean-cont_sum.se, ymax=cont_sum.mean+cont_sum.se)) +
  geom_errorbarh(aes(color = celltype, xmin=endo_sum.mean-endo_sum.se, xmax=endo_sum.mean+endo_sum.se)) +
  geom_point(aes(fill=celltype, shape = experiment), color = 'black', size = 2) +
  scale_color_manual(values = celltype_colors)+
  scale_fill_manual(values = celltype_colors)+
  facet_wrap(~ experiment, labeller = labeller(experiment = experiment.labs)) +
  #scale_shape_manual(values=c(21,22,24)) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  labs(x='Endogenous Counts',
       y = 'Contaminated Counts', 
       title  = 'Relationship between contaminated counts and endogenous counts of all strains per cell type',
       subtitle = 'Using all good cells per experiment') +
  slide_theme(base_size = 10) +
  theme(legend.position = 'bottom')

qc.contno.plots <- cowplot::plot_grid(qc.contno.plot +theme(legend.position = 'none'),
                                      qc.contno.exp.plot+theme(legend.position = 'none'),
                                      get_legend(qc.contno.exp.plot),
                                      ncol=1,
                                      rel_heights = c(1,1.25,0.1),
                                      labels=c('A', 'B', NULL))

cowplot::save_plot(paste0(workdir,'figures/SCDC/contamination/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_cont_endo_counts.pdf'),
                   qc.contno.plots,
                   base_height = 10,
                   base_width = 10)
cowplot::save_plot(paste0(workdir,'figures/SCDC/contamination/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_cont_endo_counts.png'),
                   qc.contno.plots,
                   base_height = 10,
                   base_width = 10)

# contaminated fraction vs endogenous counts qc scatterplot
qc.contfrac.plot <- good_cell_qc_stats %>% 
  ggplot(., aes(y = contfrac_sum.mean, x = endofrac_sum.mean)) +
  geom_errorbar(aes(color = celltype, ymin=contfrac_sum.mean-contfrac_sum.se, ymax=contfrac_sum.mean+contfrac_sum.se)) +
  geom_errorbarh(aes(color = celltype, xmin=endofrac_sum.mean-endofrac_sum.se, xmax=endofrac_sum.mean+endofrac_sum.se)) +
  geom_point(aes(fill=celltype, shape = experiment), color = 'black', size = 2) +
  scale_fill_manual(values = celltype_colors)+
  scale_color_manual(values = celltype_colors)+
  #scale_shape_manual(values=c(21,22,24)) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  labs(x='Endogenous Fraction',
       y = 'Contaminated Fraction', 
       title  = 'Comparing contaminated fraction to endogenous fraction of all strains per cell type',
       subtitle = 'Using all good cells per experiment') +
  slide_theme(base_size = 10) +
  theme(legend.position = 'bottom')

qc.contfrac.exp.plot <- good_cell_qc_stats %>% 
  ggplot(., aes(y = contfrac_sum.mean, x = endofrac_sum.mean)) +
  geom_errorbar(aes(color = celltype, ymin=contfrac_sum.mean-contfrac_sum.se, ymax=contfrac_sum.mean+contfrac_sum.se)) +
  geom_errorbarh(aes(color = celltype, xmin=endofrac_sum.mean-endofrac_sum.se, xmax=endofrac_sum.mean+endofrac_sum.se)) +
  geom_point(aes(fill=celltype, shape = experiment), color = 'black', size = 2) +
  scale_fill_manual(values = celltype_colors)+
  scale_color_manual(values = celltype_colors)+
  facet_wrap(~ experiment, labeller = labeller(experiment = experiment.labs)) +
  #scale_shape_manual(values=c(21,22,24)) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  labs(x='Endogenous Counts',
       y = 'Contaminated Fraction', 
       title  = 'Comparing contaminated fraction to endogenous counts of all strains per cell type',
       subtitle = 'Using all good cells per experiment') +
  slide_theme(base_size = 10) +
  theme(legend.position = 'bottom')

qc.contfrac.plots <- cowplot::plot_grid(qc.contfrac.plot +theme(legend.position = 'none'),
                                        qc.contfrac.exp.plot+theme(legend.position = 'none'),
                                        get_legend(qc.contfrac.exp.plot),
                                        ncol=1,
                                        rel_heights = c(1,1.25,0.1),
                                        labels=c('A', 'B', NULL))

cowplot::save_plot(paste0(workdir,'figures/SCDC/contamination/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_cont_endo_frac.pdf'),
                   qc.contfrac.plots,
                   base_height = 10,
                   base_width = 10)
cowplot::save_plot(paste0(workdir,'figures/SCDC/contamination/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_cont_endo_frac.png'),
                   qc.contfrac.plots,
                   base_height = 10,
                   base_width = 10)

# contaminated fraction per experiment and celltype barplot
qc.contfrac.barplot <- good_cell_qc_stats %>% 
  dplyr::left_join(samp_no_endo_good %>%  dplyr::select(experiment, celltype, no), 
                   by = c('experiment' = 'experiment', 'celltype' = 'celltype')) %>% 
  ggplot(., aes(x = reorder(celltype, no), y = contfrac_sum.mean, fill=experiment)) +
  geom_bar(stat='identity', position=position_dodge(),color='black') +
  geom_errorbar(width=.2,
                position=position_dodge(.9),
                color = 'black', aes(ymin=contfrac_sum.mean-contfrac_sum.se, 
                                     ymax=contfrac_sum.mean+contfrac_sum.se)) +
  coord_flip() +
  labs(x=NULL,
       y = 'Contaminated Fraction', 
       title  = 'Contaminated fraction per cell type and experiment',
       subtitle = 'Using all good cells per experiment') +
  slide_theme(base_size = 10) +
  theme(legend.position = 'bottom')

qc.contfrac.barplot

cowplot::save_plot(paste0(workdir,'figures/SCDC/contamination/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_cont_randomcells_barplot.pdf'),
                   qc.contfrac.barplot,
                   base_height = 10,
                   base_width = 10)
cowplot::save_plot(paste0(workdir,'figures/SCDC/contamination/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_cont_randomcells_barplot.png'),
                   qc.contfrac.barplot,
                   base_height = 10,
                   base_width = 10)

# estimated proportion: bar plot of all cells
all.cont.barplot <- res.cont.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  dplyr::left_join(samp_no_endo_good %>%  dplyr::select(experiment, celltype, prop), 
                   by = c('Experiment' = 'experiment', 'celltype' = 'celltype')) %>% 
  tidyr::pivot_longer(cols = -c(celltype, Experiment, Pseudobulk, run), names_to = "Type", values_to = 'Proportion') %>% 
  dplyr::mutate(run=ifelse(is.na(run), "1", run)) %>% 
  dplyr::mutate(Type=case_when(Type == 'EstimatedProportion' ~ 'Contamination',
                               Type == 'refproportion' ~ 'Reference',
                               Type == 'prop' ~ 'Good cells')) %>% 
  dplyr::mutate(Type = factor(Type, levels = c("Reference", 'Good cells', "Contamination"))) %>% 
  dplyr::filter(Pseudobulk == 'all') %>% 
  ggplot(., aes(x = Type, y = Proportion, fill = celltype)) +
  scale_fill_manual(values = celltype_colors)+
  geom_bar(stat = 'identity', color = 'black') +
  coord_flip() +
  facet_wrap(~ Experiment, labeller = labeller(Experiment = experiment.labs)) +
  labs(x=NULL,
       y = 'Proportion', 
       title  = 'SCDC Deconvolution of pseudobulk samples derived from good cells using contamination expression counts of all strains',
       subtitle = 'Using all good cells per experiment') +
  slide_theme(base_size = 10) +
  theme(legend.position = 'bottom')
all.cont.barplot

cont.eval.table <- ggpubr::ggtexttable(res.eval.cont.ptable %>% 
                                         dplyr::filter(Pseudobulk == 'all') %>% 
                                         dplyr::select(-Pseudobulk),
                                       rows = NULL, 
                                       theme = ggpubr::ttheme("classic", base_size = 8))

cont.overview.barplot <- cowplot::ggdraw(all.cont.barplot + 
                                           theme(plot.margin = unit(c(t=0, r=8, b=0, l=0), "cm"))) +
  cowplot::draw_plot(cont.eval.table, 0.39, 0.05)

cowplot::save_plot(paste0(workdir,'figures/SCDC/contamination/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_cont_eval_allcells_barplot.pdf'),
                   cont.overview.barplot,
                   base_height = 8,
                   base_width = 15)
cowplot::save_plot(paste0(workdir,'figures/SCDC/contamination/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_cont_eval_allcells_barplot.png'),
                   cont.overview.barplot,
                   base_height = 8,
                   base_width = 15)

#  estimated proportion of all cells sorted by cell numbers: barplot
all.error.barplot <- res.cont.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  dplyr::left_join(samp_no_endo_good %>%  dplyr::select(experiment, celltype, prop), 
                   by = c('Experiment' = 'experiment', 'celltype' = 'celltype')) %>% 
  dplyr::filter(Pseudobulk == 'all') %>% 
  ggplot(., aes(x = fct_reorder(celltype, prop), y = EstimatedProportion, fill = Experiment)) +
  geom_hline(yintercept = 0,linetype='dashed') +
  geom_bar(stat="identity",position='dodge', color = 'black') +
  coord_flip() + 
  labs(x=NULL, 
       y = 'Estimated Proportion in Contaminated Fraction', 
       title = 'SCDC Deconvolution of pseudobulk samples derived from good cells \nusing contaminated expression counts of all strains',
       subtitle = 'Using all good cells per experiment sorted by numbers per cell type') +
  slide_theme(base_size = 10) +
  theme(legend.position = 'right', 
        plot.margin = unit(c(t=0, r=1, b=0, l=0), "cm"))
all.error.barplot

cowplot::save_plot(paste0(workdir,'figures/SCDC/contamination/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_cont_eval_allcells_celltype_barplot.pdf'),
                   all.error.barplot,
                   base_height = 10,
                   base_width = 10)
cowplot::save_plot(paste0(workdir,'figures/SCDC/contamination/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_cont_eval_allcells_celltype_barplot.png'),
                   all.error.barplot,
                   base_height = 10,
                   base_width = 10)

# estimated proportion vs number of cells for random sampling, mean +- se
random.prop.plot <- res.cont.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion,) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  dplyr::left_join(samp_no_endo_good, 
                   by = c('Experiment' = 'experiment', 'celltype' = 'celltype')) %>% 
  dplyr::filter(Pseudobulk == 'random') %>% 
  ggplot(., aes(x = no, y = EstimatedProportion, color = celltype)) +
  scale_color_manual(values = celltype_colors)+
  geom_hline(yintercept = 0,linetype='dashed') +
  stat_sum_df("mean_se") +
  scale_x_log10( labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides='b') +
  labs(x='Number of cells', 
       y = 'Estimated Proportion', 
       title = 'SCDC Deconvolution of pseudobulk samples derived from good cells \nusing contaminated expression counts of all strains',
       subtitle = 'Using 25 times random sampling of good cells per experiment') +
  facet_wrap(~Experiment, labeller = labeller(Experiment = experiment.labs)) +
  slide_theme(base_size = 10) +
  theme(legend.position = 'bottom', 
        plot.margin = unit(c(t=0, r=1, b=0, l=0), "cm"))

random.error.plot <- res.cont.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion,) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  dplyr::left_join(samp_no_endo_good, 
                   by = c('Experiment' = 'experiment', 'celltype' = 'celltype')) %>% 
  dplyr::filter(Pseudobulk == 'random') %>% 
  ggplot(., aes(x = no, y = prop-EstimatedProportion, color = celltype)) +
  scale_color_manual(values = celltype_colors)+
  geom_hline(yintercept = 0,linetype='dashed') +
  stat_sum_df("mean_se") +
  scale_x_log10( labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides='b') +
  labs(x='Number of cells', 
       y = 'Difference between \n Cellular and Estimated Proportion', 
       title = 'SCDC Deconvolution of pseudobulk samples derived from good cells using contaminated expression counts of all strains',
       subtitle = 'Using 25 times random sampling of good cells per experiment') +
  facet_wrap(~Experiment, labeller = labeller(Experiment = experiment.labs)) +
  slide_theme(base_size = 10) +
  theme(legend.position = 'bottom', 
        plot.margin = unit(c(t=0, r=1, b=0, l=0), "cm"))
random.error.plot


# RESULTS EMPTY DROPLETS --------------------------------------------------

res.deconv.empty.L <- sapply(names(res.deconv.empty.files), function(i){
  tmp <- readRDS(file = paste0(workdir,"output/scdc/single_ref/", res.deconv.empty.files[i]))
  dat <- tmp$prop.est.mvw %>% 
    data.frame() %>% 
    tibble::rownames_to_column(var='experiment') %>% 
    tidyr::pivot_longer(-experiment, names_to = "celltype", values_to = "EstimatedProportion")
  dat
}, USE.NAMES = TRUE, simplify = FALSE)

res.deconv.empty.dat <- data.table::rbindlist(res.deconv.empty.L, idcol = "setup") %>% 
  tidyr::separate(col="experiment", into = c("experiment", "rprop", "run"), fill = "right", sep = "_")

res.eval.empty.L <- sapply(names(res.deconv.empty.files), function(i){
  tmp <- readRDS(file = paste0(workdir,"output/scdc/single_ref/", res.deconv.empty.files[i]))
  dat <- data.frame(RMSD_bysample = tmp$yeval$evals[[1]]$RMSDy_bysample,
                    mAD_bysample = tmp$yeval$evals[[1]]$mADy_bysample,
                    Spearman_bysample = tmp$yeval$evals[[1]]$spearmany_bysample) %>% 
    tibble::rownames_to_column(var='experiment')
  dat
}, USE.NAMES = TRUE, simplify = FALSE)

res.eval.empty.dat <- data.table::rbindlist(res.eval.empty.L, idcol = "setup") %>% 
  tidyr::separate(col="experiment", into = c("experiment", 'rprop', "run"), fill = "right", sep = "_")

# kick out wrong combinations of 'all'
# Explanation, SDCDC does not work with only one bulk sample to deconvolute so I just made either enough random pseudobulk samples or considered all experiments for pseudobulk samples.
res.empty.dat <-
  res.deconv.empty.dat %>% 
  tidyr::separate(setup, into = c("Project", "Tissue", "Species", "RNAseq", "Unit", "Strain", "Experiment", "Read", "Pseudobulk"), sep = "_", remove = FALSE) %>% 
  dplyr::mutate(Match = case_when(Experiment == experiment ~ 'correct',
                                  TRUE ~ 'kickout')) %>% 
  dplyr::filter(Match == 'correct') %>% 
  dplyr::select(-Match, -experiment)

res.eval.empty.dat <-res.eval.empty.dat %>% 
  tidyr::separate(setup, into = c("Project", "Tissue", "Species", "RNAseq", "Unit", "Strain", "Experiment", "Read", "Pseudobulk"), sep = "_", remove = FALSE) %>% 
  dplyr::mutate(Match = case_when(Experiment == experiment ~ 'correct',
                                  TRUE ~ 'kickout')) %>% 
  dplyr::filter(Match == 'correct') %>% 
  dplyr::select(-Match, -experiment)

saveRDS(res.empty.dat, file=paste0(workdir,"output/scdc/summaries/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_empty_prop.rds"))
saveRDS(res.eval.empty.dat, file=paste0(workdir,"output/scdc/summaries/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_empty_eval.rds"))

# table of evaluation metrics
res.eval.empty.table <- res.eval.empty.dat %>% 
  dplyr::group_by(setup, rprop) %>% 
  dplyr::summarise(RMSD = mean(RMSD_bysample),
                   SE.RMSD = plotrix::std.error(RMSD_bysample),
                   mAD = mean(mAD_bysample),
                   SE.mAD = plotrix::std.error(mAD_bysample),
                   Spearman = mean(Spearman_bysample),
                   SE.Spearman = plotrix::std.error(Spearman_bysample)) %>% 
  ungroup() %>% 
  tidyr::separate(setup, into = c("Project", "Tissue", "Species", "RNAseq", "Unit", "Strain", "Experiment", "Read", "Pseudobulk"), sep = "_", remove = TRUE) %>% 
  dplyr::select(Experiment, Pseudobulk, rprop, RMSD,SE.RMSD,mAD,SE.mAD,Spearman,SE.Spearman) %>% 
  data.frame()
res.eval.empty.ptable <- dplyr::select(res.eval.empty.table, Experiment, rprop, RMSD,mAD,Spearman)

# plotting
# New facet label names for pseudobulk variable
pseudobulk.labs <- c("all cells", "25 times random sampling", "25 times random sampling\n of downsampled cells") 
names(pseudobulk.labs) <- c("all", "random", "subsample")
rprop.labs <- c("25%", "75%", "100%") 
names(rprop.labs) <- c("0.25", "0.75", "1")
experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3","Experiment 4") 
names(experiment.labs) <- c("1", "2", "3", "4", "5")

# estimated proportion: bar plot of all cells
all.empty.barplot <- res.empty.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, rprop, run, EstimatedProportion) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  dplyr::left_join(samp_no_endo_good %>%  dplyr::select(experiment, celltype, prop), 
                   by = c('Experiment' = 'experiment', 'celltype' = 'celltype')) %>% 
  #the step where all cell setting is selected
  dplyr::filter(Pseudobulk == 'subsample' & rprop == "1") %>% 
  dplyr::select(-Pseudobulk, -rprop, -run) %>% 
  tidyr::pivot_longer(cols = -c(celltype, Experiment), names_to = "Type", values_to = 'Proportion') %>% 
  dplyr::mutate(Type=case_when(Type == 'EstimatedProportion' ~ 'Empty',
                               Type == 'refproportion' ~ 'Reference',
                               Type == 'prop' ~ 'Good cells')) %>% 
  dplyr::mutate(Type = factor(Type, levels = c("Reference", 'Good cells', "Empty"))) %>% 
  ggplot(., aes(x = Type, y = Proportion, fill = celltype)) +
  geom_bar(stat = 'identity', color = 'black') +
  scale_fill_manual(values = celltype_colors)+
  coord_flip() +
  facet_wrap(~ Experiment, labeller = labeller(Experiment = experiment.labs)) +
  labs(x=NULL,
       y = 'Proportion', 
       title  = 'SCDC Deconvolution of pseudobulk samples derived from empty droplets using expression counts of all strains',
       subtitle = 'Using all empty droplets per experiment') +
  slide_theme(base_size = 10) +
  theme(legend.position = 'bottom')
all.empty.barplot

empty.eval.table <- ggpubr::ggtexttable(res.eval.empty.ptable %>% 
                                          dplyr::filter( rprop == "1") %>% 
                                          dplyr::select( -rprop),
                                        rows = NULL, 
                                        theme = ggpubr::ttheme("classic", base_size = 8))

empty.overview.barplot <- cowplot::ggdraw(all.empty.barplot + 
                                            theme(plot.margin = unit(c(t=0, r=8, b=0, l=0), "cm"))) +
  cowplot::draw_plot(empty.eval.table, 0.39, 0.05)

cowplot::save_plot('figures/SCDC/emptydroplets/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_empty_eval_alldroplets_barplot.pdf',
                   empty.overview.barplot,
                   base_height = 8,
                   base_width = 15)
cowplot::save_plot('figures/SCDC/emptydroplets/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_empty_eval_alldroplets_barplot.png',
                   empty.overview.barplot,
                   base_height = 8,
                   base_width = 15)

#  estimated proportion of celltypes in droplets sorted by cell numbers: barplot
all.error.barplot <- res.empty.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, rprop, run, EstimatedProportion) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  dplyr::left_join(samp_no_endo_good %>%  dplyr::select(experiment, celltype, prop), 
                   by = c('Experiment' = 'experiment', 'celltype' = 'celltype')) %>% 
  dplyr::filter(Pseudobulk == 'subsample' & rprop == "1") %>% 
  ggplot(., aes(x = fct_reorder(celltype, prop), y = EstimatedProportion, fill = Experiment)) +
  geom_hline(yintercept = 0,linetype='dashed') +
  geom_bar(stat="identity",position='dodge', color = 'black') +
  coord_flip() + 
  labs(x=NULL, 
       y = 'Estimated Proportion in Empty Droplets', 
       title = 'SCDC Deconvolution of pseudobulk samples derived from empty droplets \nusing expression counts of all strains',
       subtitle = 'Using all empty droplets per experiment sorted by numbers per cell type') +
  slide_theme(base_size = 10) +
  theme(legend.position = 'right', 
        plot.margin = unit(c(t=0, r=1, b=0, l=0), "cm"))
all.error.barplot

cowplot::save_plot(paste0(workdir,'figures/SCDC/emptydroplets/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_empty_eval_alldroplets_celltype_barplot.pdf'),
                   all.error.barplot,
                   base_height = 10,
                   base_width = 10)
cowplot::save_plot(paste0(workdir,'figures/SCDC/emptydroplets/CZI_kidney_mouse_10XChromium_sc_CAST_singleref_empty_eval_alldroplets_celltype_barplot.png'),
                   all.error.barplot,
                   base_height = 10,
                   base_width = 10)

# estimated proportion vs number of cells for random sampling, mean +- se
random.prop.plot <- res.cont.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion,) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  dplyr::left_join(samp_no_endo_good, 
                   by = c('Experiment' = 'experiment', 'celltype' = 'celltype')) %>% 
  dplyr::filter(Pseudobulk == 'random') %>% 
  ggplot(., aes(x = no, y = EstimatedProportion, color = celltype)) +
  geom_hline(yintercept = 0,linetype='dashed') +
  scale_color_manual(values = celltype_colors)+
  stat_sum_df("mean_se") +
  scale_x_log10( labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides='b') +
  labs(x='Number of cells', 
       y = 'Estimated Proportion', 
       title = 'SCDC Deconvolution of pseudobulk samples derived from good cells \nusing contaminated expression counts of all strains',
       subtitle = 'Using 25 times random sampling of good cells per experiment') +
  facet_wrap(~Experiment, labeller = labeller(Experiment = experiment.labs)) +
  slide_theme(base_size = 10) +
  theme(legend.position = 'bottom', 
        plot.margin = unit(c(t=0, r=1, b=0, l=0), "cm"))

random.error.plot <- res.cont.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion,) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  dplyr::left_join(samp_no_endo_good, 
                   by = c('Experiment' = 'experiment', 'celltype' = 'celltype')) %>% 
  dplyr::filter(Pseudobulk == 'random') %>% 
  ggplot(., aes(x = no, y = prop-EstimatedProportion, color = celltype)) +
  geom_hline(yintercept = 0,linetype='dashed') +
  scale_color_manual(values = celltype_colors)+
  stat_sum_df("mean_se") +
  scale_x_log10( labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides='b') +
  labs(x='Number of cells', 
       y = 'Difference between \n Cellular and Estimated Proportion', 
       title = 'SCDC Deconvolution of pseudobulk samples derived from good cells using contaminated expression counts of all strains',
       subtitle = 'Using 25 times random sampling of good cells per experiment') +
  facet_wrap(~Experiment, labeller = labeller(Experiment = experiment.labs)) +
  slide_theme(base_size = 10) +
  theme(legend.position = 'bottom', 
        plot.margin = unit(c(t=0, r=1, b=0, l=0), "cm"))
random.error.plot









# now plot all together (main plot) ####

# estimated vs true proportion: bar plot of all cells
all.endo<- res.endo.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion, TrueProportion) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  tidyr::pivot_longer(cols = -c(celltype, Experiment, Pseudobulk, run), names_to = "Type", values_to = 'Proportion') %>% 
  dplyr::mutate(run=ifelse(is.na(run), "1", run)) %>% 
  dplyr::mutate(Type=case_when(Type == 'EstimatedProportion' ~ 'Estimated',
                               Type == 'TrueProportion' ~ 'True',
                               Type == 'refproportion' ~ 'Reference')) %>% 
  dplyr::filter(Pseudobulk == 'all') 



all.cont<- res.cont.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion) %>%
  tidyr::pivot_longer(cols = -c(celltype, Experiment, Pseudobulk, run), names_to = "Type", values_to = 'Proportion') %>% 
  dplyr::mutate(run=ifelse(is.na(run), "1", run)) %>% 
  dplyr::mutate(Type=case_when(Type == 'EstimatedProportion' ~ 'Contamination')) %>% 
  dplyr::filter(Pseudobulk == 'all') 



# estimated proportion: bar plot of all cells
all.empty <- res.empty.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, rprop, run, EstimatedProportion) %>%
  dplyr::filter(Pseudobulk == 'subsample' & rprop == "1") %>% 
  dplyr::select(-Pseudobulk, -rprop, -run) %>% 
  tidyr::pivot_longer(cols = -c(celltype, Experiment), names_to = "Type", values_to = 'Proportion') %>% 
  dplyr::mutate(Type=case_when(Type == 'EstimatedProportion' ~ 'Empty'))

all.props<-bind_rows(all.endo,all.cont, all.empty) %>%
  dplyr::mutate(Type = factor(Type, levels = c("Reference","True", "Estimated", "Contamination", "Empty"))) 
  
  
  
prop.plot <- ggplot(data = all.props, 
                    aes(x = Type, y = Proportion, fill = forcats::fct_rev(celltype))) + 
  geom_col(color = 'black', lwd=0.1) + 
  scale_fill_manual(values = celltype_colors) +
  guides(fill = guide_legend(nrow = 2)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  coord_flip() +
  facet_wrap(~Experiment, scales='free') +
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face='bold'),
        strip.text = element_text(size = 12, face = 'bold'),
        strip.background = element_blank()) 

ggsave("figures/SCDC/paper/Deconvolution_Proportion.pdf",
       prop.plot, 
       width = 12, height = 8)  

# deconv_eval_table %>%
#   kable() %>%
#   kable_styling() %>%
#   save_kable("figures/SCDC/paper/Deconvolution_Evaluation.html")





# cell type colors
celltype_colors <- setNames(
  c(RColorBrewer::brewer.pal(11, "Paired"),"darkgrey","grey"),
  c("PT","CD_IC" ,   "CD_PC", "CD_Trans" ,     "CNT" ,     "DCT",     "Endo",      "Fib","aLOH",   "dLOH",    "MC" ,    "Podo","Immune"))


# estimated vs true proportion: bar plot of all cells
all.endo<- res.endo.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion, TrueProportion) %>%
  dplyr::left_join(res.qc %>% dplyr::select(-cellnumber), 
                   by = c('Experiment' = 'Experiment', 'celltype' = 'celltype')) %>% 
  tidyr::pivot_longer(cols = -c(celltype, Experiment, Pseudobulk, run), names_to = "Type", values_to = 'Proportion') %>% 
  dplyr::mutate(run=ifelse(is.na(run), "1", run)) %>% 
  dplyr::mutate(Type=case_when(Type == 'EstimatedProportion' ~ 'Estimated',
                               Type == 'TrueProportion' ~ 'True',
                               Type == 'refproportion' ~ 'Reference')) %>% 
  dplyr::filter(Pseudobulk == 'all') 



all.cont<- res.cont.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, run, EstimatedProportion) %>%
  tidyr::pivot_longer(cols = -c(celltype, Experiment, Pseudobulk, run), names_to = "Type", values_to = 'Proportion') %>% 
  dplyr::mutate(run=ifelse(is.na(run), "1", run)) %>% 
  dplyr::mutate(Type=case_when(Type == 'EstimatedProportion' ~ 'Contamination')) %>% 
  dplyr::filter(Pseudobulk == 'all') 



# estimated proportion: bar plot of all cells
all.empty <- res.empty.dat %>% 
  dplyr::select(Experiment, celltype, Pseudobulk, rprop, run, EstimatedProportion) %>%
  dplyr::filter(Pseudobulk == 'subsample' & rprop == "1") %>% 
  dplyr::select(-Pseudobulk, -rprop, -run) %>% 
  tidyr::pivot_longer(cols = -c(celltype, Experiment), names_to = "Type", values_to = 'Proportion') %>% 
  dplyr::mutate(Type=case_when(Type == 'EstimatedProportion' ~ 'Empty'))

all.props<-bind_rows(all.endo,all.cont, all.empty) %>%
  dplyr::mutate(Type = factor(Type, levels = c("Reference","True", "Estimated", "Contamination", "Empty"))) 



prop.plot <- ggplot(data = all.props, 
                    aes(x = Type, y = Proportion, fill = forcats::fct_rev(celltype))) + 
  geom_col(color = 'black', lwd=0.1) + 
  scale_fill_manual(values = celltype_colors) +
  guides(fill = guide_legend(nrow = 2)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  coord_flip() +
  facet_wrap(~Experiment, scales='free') +
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face='bold'),
        strip.text = element_text(size = 12, face = 'bold'),
        strip.background = element_blank()) 

ggsave("figures/SCDC/paper/Deconvolution_Proportion_CAST.pdf",
       prop.plot, 
       width = 10, height = 8)  




saveRDS(all.props, "/data/share/htp/CZI_kidney/mouse_experiments/backgroundRNA/Extra_analyses/SCDC_CT_proportions.rds")







