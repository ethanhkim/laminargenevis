library(ggplot2)
library(tidyverse)
library(magrittr)
library(plotly)
library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)
library(here)
library(data.table)

MTG_matrix <- fread(here("data", "processed", "Allen_scRNA", "MTG_matrix.csv")) %>%
  select(-V1)
test_genelist <- c("RELN", "CXCL14", "DISC1", "CNR1", "CHRNA7", "NDNF")

na_MTG_matrix <- MTG_matrix %>%
  filter(is.na(median_exp_value))
MTG_NA_genes <- unique(na_MTG_matrix$gene)

test_scRNA_auroc <- function(source_dataset, multiple_genelist) {
  
  #Function for ranking datasets
  rank_dataset <- function(source_dataset) {
    dataset_ranked <- source_dataset %>%
      select(Layer_1:WM) %>%
      map_df(rank)
    dataset_ranked$gene_symbol <- source_dataset$gene_symbol
    dataset_ranked$class_label <- source_dataset$class_label
    dataset_ranked$region_label <- source_dataset$region_label
    return(dataset_ranked)
  }
  
  #Function to return indices 
  return_indices <- function(ranked_dataset, multiple_genelist) {
    Indices <- as_tibble(ranked_dataset$gene_symbol)
    names(Indices) <- "gene_symbol"
    Indices %<>% mutate(isTargetGene = gene_symbol %in% multiple_genelist)
    targetIndices <- Indices$isTargetGene
    return(targetIndices)
  }
  
  ## Function for AUROC
  # from https://github.com/sarbal/EGAD/blob/master/EGAD/R/auroc_analytic.R
  # by Sara Ballouz
  auroc_analytic <- function(scores, labels) {
    negatives <- which(labels == 0, arr.ind = TRUE)
    scores[negatives] <- 0
    #Calculate coefficients
    p <- sum(scores, na.rm = TRUE)
    nL <- length(labels)
    np <- sum(labels, na.rm = TRUE)
    nn <- nL - np
    #Calculate AUROC score
    auroc <- (p/np - (np + 1)/2)/nn
    return(auroc)
  } 
  
  #Function to determine AUROC values for each region, per class
  auroc_region <- function(dataset, multiple_genelist) {
    GABA <- dataset %>%
      filter(class_label == "GABAergic")
    rank_GABA <- rank_dataset(GABA)
    GABA_indices <- return_indices(GABA, multiple_genelist)
    rank_GABA %<>% select(-gene_symbol, -region_label, -class_label, -WM)
    AUROC_GABA <- map_df(rank_GABA, auroc_analytic, GABA_indices) %>%
      add_column(class_label = "GABAergic")
    
    GLUT <- dataset %>%
      filter(class_label == "Glutamatergic")
    rank_GLUT <- rank_dataset(GLUT)
    GLUT_indices <- return_indices(GLUT, multiple_genelist)
    rank_GLUT %<>% select(-gene_symbol, -region_label, -class_label, -WM)
    AUROC_GLUT <- map_df(rank_GLUT, auroc_analytic, GLUT_indices) %>%
      add_column(class_label = "Glutamatergic")
    
    NON_N <- dataset %>%
      filter(class_label == "Non-neuronal")
    rank_NON_N <- rank_dataset(NON_N)
    NON_N_indices <- return_indices(NON_N, multiple_genelist)
    rank_NON_N %<>% select(-gene_symbol, -region_label, -class_label, -WM)
    AUROC_NON_N <- map_df(rank_NON_N, auroc_analytic, NON_N_indices) %>%
      add_column(class_label = "Non-neuronal")
    
    AUROC_table <- rbind(AUROC_GABA, AUROC_GLUT, AUROC_NON_N)
    return(AUROC_table)
  }
  
  #Create list of tibbles, separating the scRNA matrix into region-specific tibbles
  region_list <- list(A1C = tibble(), V1C = tibble(), CgG = tibble(), 
                      MTG = tibble(), S1 = tibble(), M1 = tibble())
  for (i in unique(scRNA_matrix$region_label)) {
    region_list[[i]] <- scRNA_matrix %>%
      filter(region_label == i) %>%
      as_tibble()
  }
  
  #Create AUROC tables per region
  A1C <- auroc_region(region_list$A1C, multiple_genelist) %>%
    add_column(region_label = "A1C")
  V1C <- auroc_region(region_list$V1C, multiple_genelist) %>%
    add_column(region_label = "V1C")
  CgG <- auroc_region(region_list$CgG, multiple_genelist) %>%
    add_column(region_label = "CgG")
  MTG <- auroc_region(region_list$MTG, multiple_genelist) %>%
    add_column(region_label = "MTG")
  S1 <- auroc_region(region_list$S1, multiple_genelist) %>%
    add_column(region_label = "S1")
  M1 <- auroc_region(region_list$M1, multiple_genelist) %>%
    add_column(region_label = "M1")
  
  #Bind region-specific AUROC tables into one
  scRNA_auroc <- rbind(A1C, V1C, CgG, MTG, S1, M1)
  return(scRNA_auroc)
}
  


test <- MTG_matrix %>%
  filter(gene %in% test_genelist)


  ggplot(data = He_heatmap_data, mapping = aes(x = layer, y = gene_symbol, fill = Z_score)) +
    geom_tile() +
    scale_fill_distiller(palette = "RdYlBu", limits = c(-1,1)*max(abs(He_heatmap_data$Z_score))) +
    scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
    labs(y = "", x = "", title = "He et al Heatmap") +
    #Puts stars on layer marker annotations
    geom_text(aes(label = layer_label), size = 7, vjust = 1)

test %>%
  filter(class_label == "GABAergic") %>%
  ggplot(mapping = aes(x = cortical_layer_label, y = gene, fill = median_exp_value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu", limits = c(-1,1)*max(abs(test$median_exp_value))) +
  scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0))

test %>%
  filter(class_label == "Non-neuronal") %>%
  ggplot(mapping = aes(x = cortical_layer_label, y = gene, fill = median_exp_value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu", limits = c(-1,1)*max(abs(test$median_exp_value))) +
  scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0))

test %>%
  ggplot(mapping = aes(x = region_label, y = AUC, colour = class_label)) +
  geom_point() +
  facet_wrap( ~ Layer, ncol = 2)

test_process_heatmap_function <- function(source_dataset, input_genelist){
  processed_heatmap_data <- source_dataset %>%
    dplyr::filter(gene_symbol %in% input_genelist) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    as.matrix()
  processed_heatmap_data[is.nan(processed_heatmap_data)] <- 0
  processed_heatmap_data %<>% 
    as_tibble() %>%
    filter_at(vars(Layer_1, Layer_2, Layer_3, Layer_4, Layer_5, Layer_6, WM),all_vars(!is.na(.))) %>%
    mutate_if(is.numeric,as.character, is.factor, as.character) %>%
    pivot_longer(
      Layer_1:Layer_6,
      names_to = "layer"
    ) %>%
    select(-"marker_label") %>%
    # Create significance markers depending on if gene is marked in a specific layer
    mutate(layer_label = 
             case_when(
               Layer_1_marker == layer ~ "*",
               Layer_2_marker == layer ~ "*",
               Layer_3_marker == layer ~ "*",
               Layer_4_marker == layer ~ "*",
               Layer_5_marker == layer ~ "*",
               Layer_6_marker == layer ~ "*"
             )) %>%
    select(gene_symbol, layer, value, layer_label) %>%
    rename(Z_score = "value") %>%
    mutate_at("Z_score", as.numeric) %>%
    mutate_at("layer_label", ~replace(., is.na(.), ""))
  #Order data according to similarity in expression profile
  ordered_data <- source_dataset %>%
    dplyr::filter(gene_symbol %in% input_genelist) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    replace_na(list(
      Layer_1 = 0,
      Layer_2 = 0,
      Layer_3 = 0,
      Layer_4 = 0,
      Layer_5 = 0,
      Layer_6 = 0,
      WM = 0
    )) %>%
    mutate_if(is.numeric,as.character, is.factor, as.character) %>%
    select("gene_symbol":"WM")
  ordered_data_matrix <- as.matrix(ordered_data)
  rownames(ordered_data) <- ordered_data$gene_symbol
  data_dendro <- as.dendrogram(hclust(d = dist(x = ordered_data_matrix)))
  data_order <- order.dendrogram(data_dendro)
  processed_heatmap_data$gene_symbol <- factor(x = processed_heatmap_data$gene_symbol,
                                               levels = ordered_data$gene_symbol[data_order], 
                                               ordered = TRUE)
  processed_heatmap_data %<>% ungroup()
  return(processed_heatmap_data)
}



test_He_heatmap_data <- test_process_heatmap_function(He_DS1_Human_averaged, test_genelist)


test_He_heatmap_data %>%
  ggplot(mapping = aes(x = layer, y = gene_symbol, fill = Z_score)) +
    geom_tile() +
    scale_fill_distiller(palette = "RdYlBu", limits = c(-1,1)*max(abs(test_He_heatmap_data$Z_score))) +
    scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
    labs(y = "", x = "", title = "He et al Heatmap") +
    #Puts stars on layer marker annotations
    geom_text(aes(label = layer_label), size = 7, vjust = 1) +
    theme(ifelse((length(unique(test_He_heatmap_data$gene_symbol)) >= 50),
                 axis.text.y = element_blank(),
                 axis.text.y = test_He_heatmap_data$gene_symbol
    ))
 


ggplot(data = test_He_heatmap_data, mapping = aes(x = layer, y = Z_score, names = gene_symbol)) +
  geom_jitter(width = 0.1)






