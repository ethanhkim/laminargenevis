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

process_barplot_data <- function(input_genelist) {
  
  He_barplot_data <- He_DS1_Human_averaged %>%
    dplyr::filter(gene_symbol %in% input_genelist) %>%
    t() %>%
    as_tibble(rownames = NA)
  Maynard_barplot_data <- Maynard_dataset_average %>%
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


