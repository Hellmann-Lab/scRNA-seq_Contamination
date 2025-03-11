#### Genotype based background RNA estimation functions ####
library(tidyverse)
library(data.table)
library(stats)

prepare_reference_vcf <- function(strain_vcf, ENSEMBL_to_SYMBOL){
  # Filtering for the relevant strains
  strain_gt <- strain_vcf@gt %>% data.frame %>% 
    dplyr::select(c("X129S1_SvImJ", "C57BL_6NJ", "CAST_EiJ"))
  
  # Getting the genotype per strain
  strain_gt <- sapply(strain_gt, substring, 1,3) %>% data.frame() %>% 
    
    # Add chromosome,position and gene information
    mutate(position = paste0(data.frame(strain_vcf@fix)$CHROM,"_",data.frame(strain_vcf@fix)$POS),
           ENSEMBL = str_extract(strain_vcf@fix[,"INFO"], pattern = "ENS.+?(?=\\|)")) %>% 
    # Convert to symbols
    left_join(ENSEMBL_to_SYMBOL) %>% 
  
    # Categorize the SNPs based on which allele is considered as ALT in the different strains. 
    # This summarizes the 3 genotype columns and makes sure, that only homozygous SNPs are considered.
    # SNPs that do not vary between strains or are heterozygous in either of them are removed
      mutate(SNPisALT = case_when(
        CAST_EiJ == "1/1" & X129S1_SvImJ == "0/0" & C57BL_6NJ == "0/0" ~ "CAST",
        CAST_EiJ == "1/1" & X129S1_SvImJ == "1/1" & C57BL_6NJ == "0/0" ~ "CAST_SvImJ",
        CAST_EiJ == "1/1" & X129S1_SvImJ == "0/0" & C57BL_6NJ == "1/1" ~ "CAST_BL6",
        CAST_EiJ == "0/0" & X129S1_SvImJ == "1/1" & C57BL_6NJ == "0/0" ~ "SvImJ",
        CAST_EiJ == "0/0" & X129S1_SvImJ == "1/1" & C57BL_6NJ == "1/1" ~ "SvImJ_BL6",
        CAST_EiJ == "0/0" & X129S1_SvImJ == "0/0" & C57BL_6NJ == "1/1" ~ "BL6"
      )) %>%
      filter(!is.na(SNPisALT))
  
  return(strain_gt)
}

check_SNP_matching <- function(strain_vcf, cellSNP_vcf, strain_gt){
  cellSNP <- data.frame(position =  paste(data.frame(cellSNP_vcf@fix)[,1], data.frame(cellSNP_vcf@fix)[,2], sep = "_"),
                       REFk = data.frame(cellSNP_vcf@fix)[,"REF"],
                       ALTk = data.frame(cellSNP_vcf@fix)[,"ALT"])
  
  strain <- data.frame(position =  paste(data.frame(strain_vcf@fix)[,1], data.frame(strain_vcf@fix)[,2], sep = "_"),
                       REFs = data.frame(strain_vcf@fix)[,"REF"],
                       ALTs = data.frame(strain_vcf@fix)[,"ALT"])
  
  combined <- inner_join(strain, cellSNP, "position") %>% 
    filter(position %in% strain_gt$position)
  
  if(identical(combined$REFk, combined$REFs) & identical(combined$ALTk, combined$ALTs)){
    print("All good! SNP annotations are matching between reference and cellSNP vcf.")
  } else {
    print("STOP! SNP ANNOTATIONS ARE NOT MATCHING PERFECTLY! Calculations will be off, check input vcfs.")
  }
}

get_allele_counts_per_cell <- function(strain_gt, cellSNP_out, vireo_out, cellSNP_vcf, coverage_filter = 100){
  # Cell assignment data from Vireo
  donor_id <- read.table(paste0(vireo_out,"/donor_ids.tsv"), header = T)
  # ALT (alternative allele counts) and DP (ALT+REF coverage per position) matrices from cellSNP
  print("loading AD and DP matrices")
  AD <- Matrix::readMM(paste0(cellSNP_out,"/cellSNP.tag.AD.mtx"))
  DP <- Matrix::readMM(paste0(cellSNP_out,"/cellSNP.tag.DP.mtx"))
  
  print(paste0("filtering positions with coverage > ", coverage_filter))
  pos <- paste0(data.frame(cellSNP_vcf@fix)$CHROM,"_",data.frame(cellSNP_vcf@fix)$POS)
  keepPos<- Matrix::rowSums(DP)>coverage_filter & pos %in% strain_gt$position
  DP_filt<-DP[keepPos,]
  AD_filt<-AD[keepPos,]
  colnames(AD_filt) = colnames(DP_filt) = gsub(".1", "-1",donor_id$cell)
  rownames(AD_filt) = rownames(DP_filt) = pos[keepPos]
  
  print("Combine AD and DP matrices")
  # code copied from https://stackoverflow.com/questions/52662748/from-sparsematrix-to-dataframe
  AD_long <- data.frame(position=rownames(AD_filt)[AD_filt@i + 1], cell=colnames(AD_filt)[AD_filt@j + 1], AD=AD_filt@x)
  DP_long <- data.frame(position=rownames(DP_filt)[DP_filt@i + 1], cell=colnames(DP_filt)[DP_filt@j + 1], DP=DP_filt@x)
  csc <- left_join(DP_long, AD_long)
  csc$AD[is.na(csc$AD)] <- 0
  
  
  # AD_long<-melt(as.matrix(AD_filt), id=rownames(AD_filt))
  # DP_long<-melt(as.matrix(DP_filt), id=rownames(AD_filt))
  # csc<-cbind(AD_long,DP_long[,"value"])
  # colnames(csc)<-c("position","cell","AD","DP")
  # csc <- csc[csc$DP != 0,]
  
  print("Combine with strain assignments (vireo) and reference genotype table (strain_gt)")
  csc <- csc %>% 
    left_join(dplyr::select(donor_id, c(cell, donor_id)), "cell") %>% 
    left_join(dplyr::select(strain_gt, c(position, ENSEMBL, SYMBOL, SNPisALT)), "position") %>%
    dplyr::rename("Strain" = "donor_id") %>% 
    mutate(Strain = case_when(Strain == "C57BL_6NJ" ~ "BL6",
                              Strain == "129S1_SvImJ" ~ "SvImJ",
                              Strain == "CAST_EiJ" ~ "CAST",
                              T ~ Strain))
  
  print("Calculate ALT allele fraction per SNP & remove BL6 ALT SNPs")
  csc <- csc %>% 
    group_by(position) %>% 
    mutate(ALT_fraction = sum(AD) / sum(DP)) %>% 
    ungroup() %>% 
    filter(!(grepl("BL6", SNPisALT)))

  return(csc)
}

plot_ALT_allele_frac <- function(csc_intermediate){
  # Within assigned cells only 
  csc_filt <- csc_intermediate %>% 
    filter(!(Strain %in% c("unassigned", "doublet"))) %>% 
    # recalculate ALT fractions within good cells only
    group_by(position) %>% 
    mutate(ALT_fraction_cells = sum(AD) / sum(DP)) %>% 
    ungroup() 
  
  # Plot fraction of ALT alleles per SNP
  p <- plot_grid(plotlist = 
              lapply(unique(csc_filt$SNPisALT), function(i){
                filt <- csc_filt %>% filter(SNPisALT == i)
                nSNP <- length(unique(filt$position))
                nGeneCoverage <- length(unique(filt$ENSEMBL))
                
                ALT_strains <- as.vector(strsplit(i, split = "_")[[1]])
                cell_ass <- csc_filt %>% dplyr::select(cell, Strain) %>% distinct %>% .$Strain# / nrow(.)
                fraction <- sum(cell_ass %in% ALT_strains) / length(cell_ass)
                
                filt %>% dplyr::select(position,SNPisALT, ALT_fraction, ALT_fraction_cells) %>% 
                  distinct() %>% 
                  ggplot()+
                  geom_density(aes(x = ALT_fraction))+
                  geom_density(aes(x = ALT_fraction_cells), linetype = "dashed", color = "grey")+
                  geom_vline(xintercept = fraction, color = "darkgrey")+
                  ggtitle(paste0(i,", variants:",nSNP, " (in ", nGeneCoverage," genes)"))+
                  scale_x_continuous(limits = c(0,1))+
                  labs(x = "% ALT allele/ SNP")+
                  theme_bw()+
                  theme(axis.title.y = element_blank())
              }), nrow = 1)
  
  return(p)
}

calculate_crossGT_contamination <- function(csc){
  strains <- c("CAST", "SvImJ", "BL6")
  csc <- csc %>% 
    filter(!(Strain %in% c("unassigned", "doublet"))) %>% 
    group_by(SNPisALT) %>% #faster than rowwise, necessary because strsplit doesn't go rowwise in mutate
    mutate(genotype = ifelse(Strain %in% unlist(strsplit(SNPisALT,"_")), "1/1", "0/0"),
           genoNormAD = case_when(genotype == "1/1" ~ DP - AD,
                                  genotype == "0/0" ~ AD),
           contaminator = case_when(genotype == "1/1" ~ paste(strains[!(strains %in% unlist(str_split(SNPisALT  , "_")))], collapse = "_"),
                                    genotype == "0/0" ~ SNPisALT),
           contAllele_fraction =  ifelse(genotype == "1/1", 1-ALT_fraction, ALT_fraction)) %>% 
    ungroup()
  return(csc)
}

remove_strain_specific_variants <- function(csc, quantile_thresh = 0.01){
  filt_thresh <- csc %>% 
    dplyr::select(position, genotype, SNPisALT, contAllele_fraction, contaminator) %>% 
    distinct() %>% 
    group_by(genotype, SNPisALT, contaminator) %>% 
    dplyr::summarise(contAllele_thresh = quantile(contAllele_fraction, quantile_thresh))
  
  print(filt_thresh)
  
  csc_filt <- left_join(csc, filt_thresh) %>% 
    mutate(contAllele_filter = ifelse(contAllele_fraction < contAllele_thresh, T, F)) %>% 
    filter(contAllele_filter == F)
  
  # csc <- csc %>% 
  #   group_by(genotype, SNPisALT) %>% 
  #   mutate(contAllele_thresh = quantile(contAllele_fraction, quantile_thresh),
  #          contAllele_filter = ifelse(contAllele_fraction < contAllele_thresh, T, F)) %>% 
  #   filter(contAllele_filter == F)
  
  return(csc_filt)
}


plot_contAllele_fraction <- function(csc, contaminating_strains, quantile_thresh = 0.01){
  csc_filt <- csc %>% 
    filter(contaminator == contaminating_strains) %>% 
    select(position, genotype, contAllele_fraction) %>% 
    distinct()
  
  thresh <-  quantile(csc_filt$contAllele_fraction, quantile_thresh)
  median_frac <- median(csc_filt$contAllele_fraction)
  nSNP <- nrow(csc_filt)
  
  p <- csc_filt %>% 
    #filter(Strain == strain, SNPisALT == SNP) %>% 
    dplyr::select(position, contAllele_fraction) %>% 
    distinct() %>% 
    ggplot()+
    annotate("rect", xmin=0, xmax=thresh, ymin=0, ymax=Inf, alpha=0.2, fill="red")+
      #geom_rect(aes(xmin=0, xmax=thresh, ymin=0, ymax=Inf), fill.alpha = 0.4, fill = "red")+
      geom_density(aes(x = contAllele_fraction))+
      geom_vline(xintercept = thresh, color = "darkred")+
      ggtitle(label = paste0("Allele frequencies of SNPs indicating contamination coming from ", contaminating_strains, " cells"),
              subtitle = paste0("Total number of SNPs: ",nSNP," , Detectable contamination: ",round(median_frac,2)*100, "% (median), Filter: <", round(thresh,2)*100, "% ", contaminating_strains, " alleles"))+
      scale_x_continuous(limits = c(0,1))+
      labs(x = "% contaminating allele/ SNP")+
      theme_bw()+
      theme(axis.title.y = element_blank())
  
  return(p)
}

### Correction factors ####

calculate_correction2 <- function(csc){
  cor_df <- csc %>% 
    group_by(position, contaminator) %>% 
    summarise(cont_fraction = sum(genoNormAD) / sum(DP)) %>% 
    group_by(position) %>% 
    mutate(cont_ratio = case_when(contaminator == contaminator[1] ~ cont_fraction[contaminator == contaminator[2]] /  cont_fraction[contaminator == contaminator[1]],
                                  contaminator == contaminator[2] ~ cont_fraction[contaminator == contaminator[1]] /  cont_fraction[contaminator == contaminator[2]]),
           corrfac2 = 1+cont_ratio,
           corrfac2 = ifelse(is.na(corrfac2) | corrfac2 == Inf,1,corrfac2))
  
  csc_cor2 <- left_join(csc,
                        select(cor_df, position, contaminator, corrfac2)) %>% 
    mutate(AD_corrected2 = genoNormAD * corrfac2)
  
  return(csc_cor2)
}

### Summarise per cell ####
calculate_cont_perCell <- function(csc){
  csc %>% 
    group_by(Strain, cell) %>% 
    summarise(contPerCell_uncorrected = sum(genoNormAD) / sum(DP),
              contPerCell_corrected1 = sum(AD_corrected1) / sum(DP),
              contPerCell_corrected2 = sum(AD_corrected2) / sum(DP),
              contPerCell_corrected1 = ifelse(contPerCell_corrected1>1,1,contPerCell_corrected1),
              contPerCell_corrected2 = ifelse(contPerCell_corrected2>1,1,contPerCell_corrected2))
}

calculate_perCell_bootstrapping <- function(csc, nboot = 100){
  calc_cont = function(data,index){
    d = data[index,]  #create bootstrap sample of all columns of original data?
    return(sum(d$AD_corrected1) / sum(d$DP))  #calculate weighted mean using 'counts' and 'weights' columns
    
  }
  
  library(boot)
  csc %>% 
    group_split(cell) %>% 
    purrr::map_dfr(
      function(x){
        cell = unique(x$cell)
        Strain = unique(x$Strain)
        contPerCell_uncorrected = sum(x$genoNormAD) / sum(x$DP)
        contPerCell_corrected1 = sum(x$AD_corrected1) / sum(x$DP)
        contPerCell_corrected1 = ifelse(contPerCell_corrected1>1,1,contPerCell_corrected1)
        #cont_uncor = sum(x$genoNormAD) /sum(x$DP)
        CI.LL.cont_cor1 = boot.ci(boot(x, calc_cont, R = nboot), type = "basic")$basic[4]
        CI.UL.cont_cor1 = boot.ci(boot(x, calc_cont, R = nboot), type = "basic")$basic[5]
        data.frame(Strain, cell, contPerCell_uncorrected, contPerCell_corrected1, CI.LL.cont_cor1, CI.UL.cont_cor1)
      }
    )
}

# binomial way
likelihood_function <- function(rho,genoNormAD,DP,contAllele_fraction){
  density <- c()
  for(l in 1:length(DP)){
    density[l] <- dbinom(genoNormAD[l],DP[l], rho*contAllele_fraction[l])
  }
  -sum(log(density))
}

# binomial way
likelihood_function_weighted <- function(rho,genoNormAD,DP,contAllele_fraction){
  -sum(sapply( 1:length(DP), function(l){ 
    d <-dbinom(genoNormAD[l],DP[l], rho*contAllele_fraction[l])  
    return(log(DP[l]/sum(DP)) * log(d))
  }))
}

CI_function <- function(alpha = 0.05, genoNormAD, DP, contAllele_fraction, objective, cont_binom, direction){
  cut <- qchisq(p = 1-alpha, df = 1)/2
  if(objective < cut){
    res <- NA
  } 
  if(direction == "lower" & objective > cut){
    if(likelihood_function(0,genoNormAD, DP, contAllele_fraction)-objective-cut < 0){
      res <- 0
    } else {
      res <- uniroot(f = function(x)likelihood_function(x,genoNormAD, DP, contAllele_fraction)-objective-cut,
                     interval = c(0,cont_binom))$root
    }}
  if(direction == "upper" & objective > cut){
    if(likelihood_function(1,genoNormAD, DP, contAllele_fraction)-objective-cut < 0){
      res <- 1
    } else {
      res <- uniroot(f = function(x)likelihood_function(x,genoNormAD, DP, contAllele_fraction)-objective-cut,
                     interval = c(cont_binom,1))$root
    }}
  return(res)}

calculate_cont_perCell_binom <- function(csc){
  csc %>%
    group_by(cell,Strain) %>% 
    summarise(optimize_res = optimize(likelihood_function, interval = c(0,1), genoNormAD, DP, contAllele_fraction),
              contPerCell_binom = optimize_res$minimum,
              CI_low = CI_function(alpha=0.05, genoNormAD, DP, contAllele_fraction, objective = optimize_res$objective, cont_binom = contPerCell_binom, direction = "lower"),
              CI_high = CI_function(alpha=0.05,genoNormAD, DP, contAllele_fraction, objective = optimize_res$objective, cont_binom = contPerCell_binom, direction = "upper"),
              nSNP = length(DP)) %>% 
    select(-optimize_res) %>% distinct
}

calculate_cont_perCell_binom_weighted <- function(csc){
  csc %>%
    group_by(cell,Strain) %>% 
    summarise(optimize_res = optimize(likelihood_function_weighted, interval = c(0,1), genoNormAD, DP, contAllele_fraction),
              contPerCell_binom_weighted = optimize_res$minimum,
              #CI_low = CI_function(alpha=0.05, genoNormAD, DP, contAllele_fraction, objective = optimize_res$objective, cont_binom = contPerCell_binom, direction = "lower"),
              #CI_high = CI_function(alpha=0.05,genoNormAD, DP, contAllele_fraction, objective = optimize_res$objective, cont_binom = contPerCell_binom, direction = "upper"),
              nSNP = length(DP)) %>% 
    select(-optimize_res) %>% distinct
}


plot_genotype_estimation <- function(perCell_noMito_binom, perCell_noMito_CAST_binom){
  p.perCell_allStrains <- perCell_noMito_binom %>% 
    data.frame() %>% 
    mutate(rank = dense_rank(contPerCell_binom)) %>% 
    ggplot(aes(x = rank, y = contPerCell_binom))+
    geom_errorbar(aes(ymin=CI_low, ymax=CI_high), colour="grey", width=.1)+
    geom_point(size = 0.4)+
    facet_grid(~Strain)+
    scale_y_continuous(limits = c(0,1))+
    labs(x = "cell rank", y = "% background RNA (corrected)")+
    theme_bw()
  p.perCell_allStrains
  
  
  p.CAST_conf <- perCell_noMito_CAST_binom %>% 
    data.frame() %>% 
    mutate(rank = dense_rank(contPerCell_binom)) %>% 
    filter(rank < 0.75*max(rank)) %>% 
    ggplot(aes(x = rank, y = contPerCell_binom))+
    geom_errorbar(aes(ymin=CI_low, ymax=CI_high), colour="grey", width=.1)+
    geom_point(size = 0.4)+
    #scale_y_log10()+
    scale_y_continuous(limits = c(0,0.3))+
    ggtitle("First 3 quantiles of CAST cells")+
    labs(x = "cell rank", y = "% background RNA (corrected)")+
    theme_bw()
  
  p.CI_width <- perCell_noMito_CAST_binom %>% 
    data.frame() %>% 
    mutate(CI_width = CI_high - CI_low) %>% 
    ggplot(aes(x = contPerCell_binom, y = CI_width))+
    geom_point(size = 0.3)+
    stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.3) +
    scale_fill_viridis_b()+
    labs(x = "% background RNA (corrected)", y = "CI width")+
    ggtitle("Confidence vs contamination level")+
    theme_bw()+
    theme(legend.position = "none")
  
  plot_grid(p.CAST_conf, p.CI_width, nrow = 1)
}
