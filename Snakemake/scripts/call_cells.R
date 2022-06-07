args <- commandArgs(TRUE)

library(DropletUtils)
library(tidyverse)

call_cells <- function(in.cellranger, out.cell_BC, out.top75k_BC, out.UMI_curve){
  #order cells with barcodeRanks function:
  raw <- Seurat::Read10X(paste0(in.cellranger, "/raw_feature_bc_matrix"))
  cellranger_calls <- read.table(paste0(in.cellranger, "/filtered_feature_bc_matrix/barcodes.tsv.gz"))$V1
  
  # order cells by UMI count and identify inflection point for more stringent filtering
  br.out <- barcodeRanks(raw, lower = 700)
  
  # Call cells: Above inflection point and called as cell by cellranger v3
  cell_assignment <- br.out %>% data.frame %>% 
    rownames_to_column("BC") %>% 
    dplyr::filter(rank <= 100000) %>%
    mutate(cellranger_call = ifelse(BC %in% cellranger_calls, T,F),
           above_inflection = ifelse(total > br.out@metadata$inflection,T,F),
           isCell = ifelse(cellranger_call == T & above_inflection == T, T, F))
  
  # Define cell containing droplets:
  cell_BC <- cell_assignment$BC[cell_assignment$isCell == T]
  # Define top 75k barcodes with highest total UMI:
  top75k_BC <- arrange(cell_assignment, -total)$BC[1:75000]
  
  # Plot UMI curve:
  p.UMI_curve <- ggplot(cell_assignment, aes(x = rank, y = total, color = isCell))+
      scale_color_manual(values = c("grey","red"))+
      geom_point(size = 0.2)+
      scale_y_log10(limits = c(1, 100000))+
      scale_x_log10()+
      geom_hline(yintercept = br.out@metadata$inflection, linetype = "dashed", color = "green")+
      #geom_vline(xintercept = c(nUMI_low,nUMI_high), color = "#68C3D4") +
      labs(title = paste0("nCells =  ", sum(cell_assignment$isCell)), x = "cell rank", y = "total UMI count")+
      theme_bw()+
      theme(legend.position = c(0.2,0.2))+ 
      guides(colour = guide_legend(override.aes = list(size=4)))
  
  write.table(cell_BC, out.cell_BC, quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(top75k_BC, out.top75k_BC,  quote = F, row.names = F, col.names = F)
  ggsave(out.UMI_curve, p.UMI_curve,  width = 4, height = 3.5, device = "pdf")
}

call_cells(args[1],args[2],args[3], args[4])
