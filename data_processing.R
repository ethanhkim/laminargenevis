library(reshape2)
library(dplyr)

## Function to process heatmap data ##
process_heatmap_function <- function(source_dataset, input_genelist){
  
  source_dataset %>%
    dplyr::filter(gene_symbol %in% input_genelist) %>%
    distinct(gene_symbol, .keep_all = TRUE)  %>%
    arrange(gene_symbol) %>%
    column_to_rownames(var = "gene_symbol") %>%
    t() %>%
    melt() %>%
    as_tibble() %>%
    dplyr::rename(Cuts = "Var1", Gene_symbol = "Var2", Z_score = "value")
}

## Function to process barplot data ##

process_barplot_data <- function(input_genelist) {
  
  He_barplot_data <- He_DS1_Human_averaged %>%
    dplyr::filter(gene_symbol %in% input_genelist)
  Maynard_barplot_data <- Maynard_dataset_average %>%
    dplyr::filter(gene_symbol %in% input_genelist)
  
  Barplot_data <- tibble(
    "gene_symbol" = He_barplot_data$gene_symbol,
    "He_Layer_1" = He_barplot_data$Layer_1, "He_Layer_2" = He_barplot_data$Layer_2,
    "He_Layer_3" = He_barplot_data$Layer_3, "He_Layer_4" = He_barplot_data$Layer_4,
    "He_Layer_5" = He_barplot_data$Layer_5, "He_Layer_6" = He_barplot_data$Layer_6,
    "He_WM" = He_barplot_data$WM, 
    "Maynard_Layer_1" = Maynard_barplot_data$Layer_1, "Maynard_Layer_2" = Maynard_barplot_data$Layer_2,
    "Maynard_Layer_3" = Maynard_barplot_data$Layer_3, "Maynard_Layer_4" = Maynard_barplot_data$Layer_4,
    "Maynard_Layer_5" = Maynard_barplot_data$Layer_5, "Maynard_Layer_6" = Maynard_barplot_data$Layer_6,
    "Maynard_WM" = Maynard_barplot_data$WM) %>%
    pivot_longer(cols = -gene_symbol, names_to = "Layers") %>%
    dplyr::rename(Z_score = "value") %>%
    mutate(Dataset = ifelse(test = str_detect(Layers, "He"),
                            yes = "He",
                            no = "Maynard")) %>%
    mutate(Layer = case_when(
      str_detect(Layers, "1") ~ "1", str_detect(Layers, "2") ~ "2", str_detect(Layers, "3") ~ "3",
      str_detect(Layers, "4") ~ "4", str_detect(Layers, "5") ~ "5", str_detect(Layers, "6") ~ "6",
      str_detect(Layers, "WM") ~ "WM"
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


