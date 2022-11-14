args <- commandArgs(TRUE)

library(tidyverse)

evaluate_estimation_accuracy <- function(corrected_dir, genotype_estimate, estimation_accuracy){
  CpC_gt <- readRDS(paste0(genotype_estimate, "/perCell_noMito_CAST_binom.RDS"))
  
  eval_list <- list()
  for(param in list.files(corrected_dir, pattern = "_contPerCell.RDS")){
    CpC <- readRDS(paste0(corrected_dir,"/",param))
    df <- inner_join(CpC, CpC_gt)
    eval_list[[gsub("_contPerCell.RDS","",param)]] <- data.frame(
      tau = cor(df$contPerCell_binom,df$cont, method = "kendall"),
      rmsle = sqrt(mean((log1p(df$contPerCell_binom) - log1p(df$cont))^2))
    )
  }
  
  eval_df <- bind_rows(eval_list, .id="param")
  
  saveRDS(eval_df, estimation_accuracy)
}
  
evaluate_estimation_accuracy(args[1],args[2],args[3])
