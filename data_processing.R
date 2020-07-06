library(reshape2)
library(dplyr)

## Function to process heatmap data ##
process_heatmap_function <- function(source_dataset, input_genelist){
  
  source_dataset %>%
    dplyr::filter(gene_symbol %in% input_genelist) %>%
    distinct(gene_symbol, .keep_all = TRUE)  %>%
    mutate_if(is.numeric,as.character, is.factor, as.character) %>%
    pivot_longer(
      Layer_1:WM,
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
               Layer_6_marker == layer ~ "*",
               Layer_WM_marker == layer ~ "*"
             )) %>%
    select(gene_symbol, layer, value, layer_label) %>%
    rename(Z_score = "value") %>%
    mutate_at("Z_score", as.numeric) %>%
    mutate_at("layer_label", ~replace(., is.na(.), ""))
}

## Function to process barplot data ##

process_barplot_data <- function(input_genelist, dataset1, dataset2) {
  
  He_barplot_data <- dataset1 %>%
    dplyr::filter(gene_symbol %in% input_genelist) %>%
    t() %>%
    as_tibble(rownames = NA)
  Maynard_barplot_data <- dataset2 %>%
    dplyr::filter(gene_symbol %in% input_genelist) %>%
    t() %>%
    as_tibble(rownames = NA)
  
  Barplot_data <- He_barplot_data %>%
    rownames_to_column(var = "variables") %>%
    add_column("Maynard_data" = Maynard_barplot_data$V1) %>%
    column_to_rownames(var = "variables") %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    pivot_longer(cols = c("Layer_1", "Layer_2", "Layer_3", "Layer_4", "Layer_5", "Layer_6", "WM"),
                 names_to = "layer") %>%
    dplyr::rename(Z_score = "value") %>%
    add_column("layers" = c("He_Layer_1", "He_Layer_2", "He_Layer_3", "He_Layer_4", 
                            "He_Layer_5", "He_Layer_6", "He_WM", "Maynard_Layer_1",
                            "Maynard_Layer_2", "Maynard_Layer_3", "Maynard_Layer_4", 
                            "Maynard_Layer_5", "Maynard_Layer_6","Maynard_WM")) %>%
    mutate(layer_label = 
             case_when(
               Layer_1_marker == layer ~ "*",
               Layer_2_marker == layer ~ "*",
               Layer_3_marker == layer ~ "*",
               Layer_4_marker == layer ~ "*",
               Layer_5_marker == layer ~ "*",
               Layer_6_marker == layer ~ "*",
               Layer_WM_marker == layer ~ "*"
             )) %>%
    select(gene_symbol, layers, Z_score, layer_label) %>%
    mutate_at("Z_score", as.numeric) %>%
    mutate_at("layer_label", ~replace(., is.na(.), "")) %>%
    mutate(Dataset = ifelse(test = str_detect(layers, "He"),
                            yes = "He",
                            no = "Maynard")) %>%
    mutate(Layer = case_when(
      str_detect(layers, "1") ~ "1", str_detect(layers, "2") ~ "2", str_detect(layers, "3") ~ "3",
      str_detect(layers, "4") ~ "4", str_detect(layers, "5") ~ "5", str_detect(layers, "6") ~ "6",
      str_detect(layers, "WM") ~ "WM"
    ))
}

separate_layers <- function(input_table, input_genelist) {
  
  list_of_layers <- list()
  
  for(i in c(1:7)) {
    layer <- input_table %>%
      dplyr::filter(str_detect(source_dataset, 'Zeng')) %>%
      dplyr::filter(gene_symbol %in% input_genelist, 
             layer_marker == i) %>% 
      select(-"layer_marker", -"source_dataset")
    layer_gene_symbol <- as.vector(layer$gene_symbol)
    list_of_layers[[i]] <- layer_gene_symbol
  }
  
 return(list_of_layers)
}


## Function for generating single-gene correlation across He & Maynard datasets

single_gene_correlation <- function(input_gene, He_dataset, Maynard_dataset) {
  
  He_gene <- He_dataset %>%
    select(gene_symbol:WM) %>%
    filter(gene_symbol %in% input_gene) %>%
    column_to_rownames(var = "gene_symbol") %>%
    t()
  
  Maynard_gene <- Maynard_dataset %>%
    select(gene_symbol:WM) %>%
    filter(gene_symbol %in% input_gene) %>%
    column_to_rownames(var = "gene_symbol") %>%
    t()
  
  He_Maynard_cormatrix <- cor(He_gene, Maynard_gene, method = "pearson") %>%
    as_tibble(rownames = NA)
  
  He_Maynard_vector <- He_Maynard_cormatrix %>%
    pull(var = 1) %>%
    as.numeric()
  
  return(format(round(He_Maynard_vector, 2), nsmall = 2))
}


## Functions for generating multi-gene correlation across He & Maynard datasets

dataset_correlation <- function(He_dataset, Maynard_dataset) {
  
  He_Maynard_common_gene_list <- intersect(He_dataset$gene_symbol, Maynard_dataset$gene_symbol)
  
  He_matrix <- He_dataset %>%
    select(gene_symbol:WM) %>%
    filter(gene_symbol %in% He_Maynard_common_gene_list) %>%
    column_to_rownames(var = "gene_symbol") %>%
    t()
  
  Maynard_matrix <- Maynard_dataset %>%
    select(gene_symbol:WM) %>%
    filter(gene_symbol %in% He_Maynard_common_gene_list) %>%
    column_to_rownames(var = "gene_symbol") %>%
    t()
  
  He_Maynard_cormatrix <- cor(He_matrix, Maynard_matrix, method = "pearson")
  He_Maynard_diag <- diag(He_Maynard_cormatrix, names = TRUE)
  
  return(He_Maynard_diag)
}

multi_gene_correlation <- function(input_genes, He_dataset, Maynard_dataset) {
  
  He_genes <- He_dataset %>%
    select(gene_symbol:WM) %>%
    filter(gene_symbol %in% input_genes) %>%
    column_to_rownames(var = "gene_symbol") %>%
    t()
  
  Maynard_genes <- Maynard_dataset %>%
    select(gene_symbol:WM) %>%
    filter(gene_symbol %in% input_genes) %>%
    column_to_rownames(var = "gene_symbol") %>%
    t() 
 
  He_Maynard_genes_cormatrix <- cor(He_genes, Maynard_genes, method = "pearson")
  He_Maynard_diag_genes <- diag(He_Maynard_genes_cormatrix, names = TRUE)
  He_Maynard_diag_genes <- format(round(He_Maynard_diag_genes, 2), nsmall = 2) %>%
    as_tibble() %>%
    pull(var = 1) %>%
    as.numeric() %>%
    mean() 
  
  return(format(round(He_Maynard_diag_genes, 2), nsmall = 2))
}


## Function for generating 

quantile_distribution <- function(dataset_diagonal, genes_diagonal) {
  
  quantile_distribution <- ecdf(dataset_diagonal)
  multiple_genes_quantile <- quantile_distribution(genes_diagonal)
  
  multiple_genes_quantile <- multiple_genes_quantile * 100
  
  return(format(floor(multiple_genes_quantile)))
}
