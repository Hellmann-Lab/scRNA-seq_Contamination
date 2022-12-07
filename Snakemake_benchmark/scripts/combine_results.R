args <- commandArgs(TRUE)

library(tidyverse)

combine_results <- function(estimation_accuracy, cluster_internal, cluster_external, kNN_overlap, expr_fraction, marker_lfc, out){
  # Estimation accuracy
  l_est <- list()
  for(res in strsplit(estimation_accuracy, ",")[[1]]){
    l_est[[res]] <- readRDS(res) %>% data.frame() %>% 
      pivot_longer(cols = c("tau","rmsle"), names_to = "metric") %>% 
      mutate(method = word(res,3,3,"/"),
             replicate = word(res,4,4,"/"))
  }
  df_est <-  bind_rows(l_est) %>% mutate(evaluation_category = "estimation_accuracy")
  
  # Cluster evaluation internal
  l_cl_int <- list()
  for(res in strsplit(cluster_internal, ",")[[1]]){
    l_cl_int[[res]] <- readRDS(res) %>% data.frame() %>% 
      mutate(method = word(res,3,3,"/"),
             replicate = word(res,4,4,"/"))
  }
  df_cl_int <-  bind_rows(l_cl_int) %>% mutate(evaluation_category = "cluster_evaluation")
  
  # Silhouette per cell type
  df_silhouette_ct <- df_cl_int %>%
    dplyr::filter(grepl("silhouette_", metric)) %>%
    group_by(method, replicate, param) %>%
    dplyr::summarize(metric = "avg_silhouette",
                     value = mean(value),
                     evaluation_category = "cluster_evaluation")
  
  # Cluster evaluation external
  l_cl_ext <- list()
  for(res in strsplit(cluster_external, ",")[[1]]){
    l_cl_ext[[res]] <- readRDS(res) %>% data.frame() %>% 
      mutate(method = word(res,3,3,"/"),
             replicate = word(res,4,4,"/"))
  }
  df_cl_ext <-  bind_rows(l_cl_ext) %>% mutate(evaluation_category = "cluster_evaluation")
  
  # kNN overlap 
  l_kNN <- list()
  for(res in strsplit(kNN_overlap, ",")[[1]]){
    l_kNN[[res]] <- readRDS(res) %>% data.frame() %>% 
      pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>% 
      mutate(method = word(res,3,3,"/"),
             replicate = word(res,4,4,"/"),
             metric = "kNN_overlap") %>% 
      relocate(metric, .after = param)
  }
  df_kNN <-  bind_rows(l_kNN) %>% mutate(evaluation_category = "cluster_evaluation")
  
  #### Marker genes ####
  topGenes <- readRDS("input/top10_PT_markers.RDS")
  
  # expression fraction
  l_expr <- list()
  for(res in strsplit(expr_fraction, ",")[[1]]){
    l_expr[[res]] <- readRDS(res) %>% 
      dplyr::filter(expr_within_celltype > 0, gene %in% topGenes) %>% 
      group_by(param) %>% 
      summarise(expression_fraction = mean(expr_within_others),
                log_ratio_expression = mean(-log((expr_within_others+1)/(expr_within_celltype + 1)))) %>% 
      pivot_longer(cols = c(expression_fraction, log_ratio_expression), names_to = "metric", values_to = "value") %>% 
      mutate(method = word(res,3,3,"/"), replicate = word(res,4,4,"/"))
    }
  df_expr <-  bind_rows(l_expr) %>% mutate(evaluation_category = "marker_evaluation")
  
  # log2 fold changes
  l_lfc <- list()
  for(res in strsplit(marker_lfc, ",")[[1]]){
    l_lfc[[res]] <- readRDS(res) %>% data.frame %>% 
        rownames_to_column("gene") %>% 
        mutate(gene = gsub("\\..*", "", gene)) %>% 
        dplyr::filter(gene %in% topGenes) %>% 
        group_by(param) %>% 
        summarise(metric = "lfc", value = mean(avg_log2FC)) %>% 
        mutate(method = word(res,3,3,"/"), replicate = word(res,4,4,"/"))
    }
  
  df_lfc <-  bind_rows(l_lfc) %>% mutate(evaluation_category = "marker_evaluation")
  
  all_metrics <- Reduce(bind_rows, list(df_est, df_cl_int, df_silhouette_ct, df_cl_ext, df_kNN, df_expr, df_lfc))
  
  saveRDS(all_metrics, out)
}

combine_results(args[1],args[2],args[3],args[4],args[5],args[6],args[7])
