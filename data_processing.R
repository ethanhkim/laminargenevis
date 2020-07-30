## Function to format Zeng dataset
process_Zeng_dataset <- function(Zeng_dataset) {
  Zeng_dataset %>%
    filter(region == "V1") %>%
    select("gene_symbol", "marker_annotation", "expression_level")
}

## Function to generate height for heatmap ##
heatmap_height <- function(genelist) {
  if (length(genelist) <= 10) {
    300
  } else {
    400
  }
}

## Function to process heatmap data ##
process_heatmap_function <- function(source_dataset, input_genelist){
  
  processed_heatmap_data <- source_dataset %>%
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
  
  #Order data according to similarity in expression profile
  ordered_data <- source_dataset %>%
    dplyr::filter(gene_symbol %in% input_genelist) %>%
    distinct(gene_symbol, .keep_all = TRUE)  %>%
    mutate_if(is.numeric,as.character, is.factor, as.character) %>%
    select("gene_symbol":"WM")
  
  ordered_data_matrix <- as.matrix(ordered_data)
  rownames(ordered_data) <- ordered_data$gene_symbol
  data_dendro <- as.dendrogram(hclust(d = dist(x = ordered_data_matrix)))
  
  data_order <- order.dendrogram(data_dendro)
  
  processed_heatmap_data$gene_symbol <- factor(x = processed_heatmap_data$gene_symbol,
                                         levels = ordered_data$gene_symbol[data_order], 
                                         ordered = TRUE)
  return(processed_heatmap_data)
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


## Function for generating gene correlations across He & Maynard datasets

single_gene_correlation <- function(input_genes, He_dataset, Maynard_dataset) {
  
  He_gene <- He_dataset %>%
    select(gene_symbol:WM) %>%
    filter(gene_symbol %in% input_genes) %>%
    column_to_rownames(var = "gene_symbol") %>%
    t()
  
  Maynard_gene <- Maynard_dataset %>%
    select(gene_symbol:WM) %>%
    filter(gene_symbol %in% input_genes) %>%
    column_to_rownames(var = "gene_symbol") %>%
    t() 
  
  He_Maynard_diag_gene <- cor(He_gene, Maynard_gene, method = "pearson") %>%
    as_tibble() %>%
    pull(var = 1) %>%
    as.numeric()
  
  return(format(signif(He_Maynard_diag_gene, digits = 3)))
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


## Function for generating quantile

quantile_distribution <- function(dataset_diagonal, genes_diagonal) {
  
  quantile_distribution <- ecdf(dataset_diagonal)
  multiple_genes_quantile <- quantile_distribution(genes_diagonal)
  multiple_genes_quantile <- multiple_genes_quantile * 100
  
  return(format(floor(multiple_genes_quantile)))
}

## Function for Wilcoxon test
wilcoxtest <- function(input_genes, He_dataset, Maynard_dataset, He_Maynard_diagonal) {

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
  
  wilcoxtest <- wilcox.test(He_Maynard_diag_genes, He_Maynard_diagonal)$p.value
  
  return(format(signif(wilcoxtest, digits = 4)))
}

## Overall function for AUROC and AUROC table:
AUROC_function <- function(He_dataset, Maynard_dataset, multiple_genelist) {
  
  #Function for ranking datasets
  rank_dataset <- function(source_dataset) {
    dataset_ranked <- source_dataset %>%
      select(Layer_1:WM) %>%
      map_df(rank)
    dataset_ranked$gene_symbol <- source_dataset$gene_symbol
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
    
    p <- sum(scores, na.rm = TRUE)
    nL <- length(labels)
    np <- sum(labels, na.rm = TRUE)
    nn <- nL - np
    
    auroc <- (p/np - (np + 1)/2)/nn
    
    return(auroc)
  } 
  
  ## Function for MWU for AUROC
  apply_MWU <- function(column, targetIndices) {
    wilcox.test(column[targetIndices], column[!targetIndices], conf.int = F)$p.value
  }
  
  He_df <- rank_dataset(He_dataset)
  Maynard_df <- rank_dataset(Maynard_dataset)
  
  He_indices <- return_indices(He_df, multiple_genelist)
  Maynard_indices <- return_indices(Maynard_df, multiple_genelist)
  
  He_df %<>% select(-gene_symbol)
  Maynard_df %<>% select(-gene_symbol)
  
  AUROC_He <- map_df(He_df, auroc_analytic, He_indices)
  wilcox_AUROC_He <- map_df(He_df, apply_MWU, He_indices)
  AUROC_Maynard <- map_df(Maynard_df, auroc_analytic, Maynard_indices)
  wilcox_AUROC_Maynard <- map_df(Maynard_df, apply_MWU, Maynard_indices)
  
  # Create AUROC table for display
  AUROC_table <- bind_cols(gather(AUROC_He, key = Layers, value = AUROC_He),
                           gather(wilcox_AUROC_He, value = pValue),
                           gather(AUROC_Maynard, key = Layers, value = AUROC_Maynard),
                           gather(wilcox_AUROC_Maynard, value = pValue)) %>%
    rename(Layers = Layers...1) %>%
    rename(key = key...3) %>%
    rename(pValue_He = pValue...4) %>%
    rename(pValue_Maynard = pValue...8) %>%
    select(Layers, AUROC_He, pValue_He, AUROC_Maynard, pValue_Maynard) %>% 
    mutate(pValue_He = signif(pValue_He, digits = 3),
           pValue_Maynard = signif(pValue_Maynard, digits = 3),
           AUROC_He = signif(AUROC_He, digits = 3),
           AUROC_Maynard = signif(AUROC_Maynard, digits = 3),
           adjusted_P_He = signif(p.adjust(pValue_He), digits = 3),
           adjusted_P_Maynard = signif(p.adjust(pValue_Maynard), digits = 3)) %>%
    select(Layers, AUROC_He, pValue_He, adjusted_P_He, AUROC_Maynard, pValue_Maynard, adjusted_P_Maynard)
}


