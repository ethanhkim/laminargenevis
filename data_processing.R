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
    distinct(gene_symbol, .keep_all = TRUE) %>%
    as.matrix()
  
  processed_heatmap_data[is.nan(processed_heatmap_data)] <- 0
  
  processed_heatmap_data %<>% 
    as_tibble() %>%
    #Filter out genes that have NA for all layers
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
  
  return(processed_heatmap_data)
}

## Function to process barplot data ##

process_barplot_data <- function(input_genelist, He_dataset, Maynard_dataset, scRNA_dataset) {
  
  He_barplot_data <- He_dataset %>%
    dplyr::filter(gene_symbol %in% input_genelist) %>%
    t() %>%
    as_tibble(rownames = NA)
  Maynard_barplot_data <- Maynard_dataset %>%
    dplyr::filter(gene_symbol %in% input_genelist) %>%
    t() %>%
    as_tibble(rownames = NA)
  scRNA_barplpot_data <- MTG_matrix_scaled %>%
    dplyr::filter(gene_symbol %in% input_genelist) %>%
    rename("Z_score" = mean_expression_scaled, layer = cortical_layer_label) %>%
    mutate(Layer = gsub("L", "Layer_", layer)) %>%
    unite(layers, c("class_label", "Layer"), sep = "_", remove = F) %>%
    mutate(Layer = gsub("Layer_", "", Layer)) %>%
    rename(Dataset = class_label) %>%
    add_column(layer_label = "") %>%
    select(gene_symbol, layers, Z_score, layer_label, Dataset, Layer)
  
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
  
  Barplot_data %<>%
    rbind(scRNA_barplpot_data) %>%
    mutate(Dataset = factor(Dataset, levels = c("He", "Maynard", "GABAergic", "Glutamatergic", "Non-neuronal")))
  
  return(Barplot_data)
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


## Functions for generating singe or multi-gene correlations across He & Maynard datasets

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
    filter(across(everything(), ~ !is.na(.))) %>%
    t()
  
  Maynard_genes <- Maynard_dataset %>%
    select(gene_symbol:WM) %>%
    filter(gene_symbol %in% input_genes) %>%
    column_to_rownames(var = "gene_symbol") %>%
    filter(across(everything(), ~ !is.na(.))) %>%
    t() 
  
  He_Maynard_genes_cormatrix <- cor(He_genes, Maynard_genes, method = "pearson")
  He_Maynard_diag_genes <- diag(He_Maynard_genes_cormatrix, names = TRUE) %>%
    #mean correlation value of all genes included
    mean() 
  
  return(format(round(He_Maynard_diag_genes, 3), nsmall = 3))
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

#Function for ranking bulk tissue datasets
rank_bulk_dataset <- function(source_dataset) {
  dataset_ranked <- source_dataset %>%
    select(Layer_1:WM) %>%
    map_df(rank)
  dataset_ranked$gene_symbol <- source_dataset$gene_symbol
  return(dataset_ranked)
}

#Function for ranking scRNA datasets
rank_dataset_scRNA <- function(source_dataset) {
  dataset_ranked <- source_dataset %>%
    ungroup() %>%
    pivot_wider(
      names_from = cortical_layer_label,
      values_from = mean_expression_scaled
    )
  gene_symbol <- dataset_ranked$gene_symbol
  dataset_ranked %<>% select(L1:L6) %>% map_df(rank, ties.method = "min") %>%
    add_column(gene_symbol = gene_symbol, class_label = "GABAergic")
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

## Function for performing AUROC on cell types in scRNA data
auroc_cell_type <- function(dataset, multiple_genelist, cell_type) {
  cell_type_df <- dataset %>%
    filter(class_label == cell_type)
  ranked_cell_type_df <- rank_dataset_scRNA(cell_type_df)
  indices <- return_indices(ranked_cell_type_df, multiple_genelist)
  ranked_cell_type_df %<>% select(-gene_symbol, -class_label)
  AUROC_cell_type <- map_df(ranked_cell_type_df, auroc_analytic, indices)
  AUROC_MWU <- map_df(ranked_cell_type_df, apply_MWU, indices)
  AUROC_cell_type %<>% rbind(AUROC_MWU) %>% as.data.frame()
  row.names(AUROC_cell_type) <- c(paste("AUROC", cell_type, sep = "_"), "p-Value")
  return(AUROC_cell_type)
}


## Function for generating bulk tissue AUROC and p-value data
AUROC_bulk <- function(He_dataset, Maynard_dataset, multiple_genelist) {
  He_df <- rank_bulk_dataset(He_dataset)
  Maynard_df <- rank_bulk_dataset(Maynard_dataset)
  
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
    select(Layers, AUROC_He, pValue_He, adjusted_P_He, AUROC_Maynard, pValue_Maynard, adjusted_P_Maynard) %>%
    rename("AUROC (He)" = AUROC_He,
           "p-value (He)" = pValue_He,
           "Adjusted p-value (He)" = adjusted_P_He,
           "AUROC (Maynard)" = AUROC_Maynard,
           "p-value (Maynard)" = pValue_Maynard,
           "Adjusted p-value (Maynard)" = adjusted_P_Maynard)
  
  AUROC <- AUROC_table %>% select(starts_with("AUROC"), Layers) %>%
    mutate(Layers = gsub("Layer_", "L", Layers)) %>%
    pivot_longer(cols = starts_with("AUROC"),
                 names_to = c("class_label"),
                 values_to = "AUROC_values") %>%
    mutate(class_label = gsub("AUROC_", "", class_label)) %>%
    rename(AUROC = AUROC_values, Layer = Layers)
  bulk_table_pValue <- AUROC_table %>% select(starts_with("Adjusted"), Layers) %>%
    mutate(Layers = gsub("Layer_", "L", Layers)) %>%
    pivot_longer(cols = starts_with("Adjusted"),
                 names_to = c("class_label"),
                 values_to = "adj.pvalues") %>%
    mutate(adj.pvalues = ifelse(adj.pvalues > 0.05, NA, adj.pvalues))
  AUROC %<>%
    add_column(pValue = bulk_table_pValue$adj.pvalues) %>%
    mutate(signif_marker = ifelse(!is.na(pValue), "*", ""))
}

## Function for generating scRNA AUROC and p-value data
AUROC_scRNA <- function(source_dataset, multiple_genelist) {
  scRNA_AUROC_list <- list()
  for (i in unique(source_dataset$class_label)) {
    scRNA_AUROC_list[[i]] <- auroc_cell_type(source_dataset, multiple_genelist, i)
  }
  scRNA_AUROC_table <- rbindlist(scRNA_AUROC_list, idcol = "rownames") %>%
    mutate(rownames = c("AUROC_GABA", "pValue_GABA", "AUROC_GLUT", "pValue_GLUT", "AUROC_NON", "pValue_NON")) %>%
    column_to_rownames(var = "rownames") %>% 
    add_column(WM = NA) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Layer")
  scRNA_pValue_table <- scRNA_AUROC_table %>% select(starts_with("pValue"), Layer) %>%
    pivot_longer(cols = starts_with("pValue"),
                 names_to = c("p_Value"),
                 values_to = "pValue") %>%
    mutate(p_Value = gsub("pValue_", "", p_Value)) %>%
    rename(class_label = p_Value) %>%
    mutate(pValue = p.adjust(pValue)) %>%
    mutate(pValue = ifelse(pValue > 0.05, NA, pValue))
  scRNA_AUROC_table %<>% select(starts_with("AUROC"), Layer) %>%
    pivot_longer(cols = starts_with("AUROC"),
                 names_to = c("class_label"),
                 values_to = "AUROC_values") %>%
    mutate(class_label = gsub("AUROC_", "", class_label)) %>%
    rename(AUROC = AUROC_values) %>%
    add_column(pValue = scRNA_pValue_table$pValue) %>%
    mutate(signif_marker = ifelse(!is.na(pValue), "*", ""))
  return(scRNA_AUROC_table)
}

## Function for binding scRNA and bulk tissue AUROC data to view as heatmap
AUROC_data <- function(scRNA_AUROC_table, bulk_AUROC_table) {
  AUROC_table <- rbind(scRNA_AUROC_table, bulk_AUROC_table) %>%
    rename(dataset = class_label) %>%
    mutate(dataset = gsub("GABA", "scRNA_GABA", dataset)) %>%
    mutate(dataset = gsub("GLUT", "scRNA_GLUT", dataset)) %>%
    mutate(dataset = gsub("NON", "scRNA_Non-neuronal", dataset)) %>%
    mutate(dataset = gsub("AUROC ", "", dataset)) %>%
    mutate(dataset = gsub("\\(", "", dataset)) %>%
    mutate(dataset = gsub("\\)", "", dataset))
  AUROC_table$Layer <- factor(AUROC_table$Layer, levels = c("WM", "L6", "L5", "L4", "L3", "L2", "L1"))
  return(AUROC_table)
}


#scRNA Heatmaps ---- 

scRNA_Heatmap_data <- function(data, genelist, cellType) {
  heatmap_data <- data %>%
    filter(gene_symbol %in% genelist) %>%
    filter(class_label == cellType) %>%
    rename(Layer = cortical_layer_label, Mean_Expression = mean_expression_scaled, Cell_Type = class_label)
}


