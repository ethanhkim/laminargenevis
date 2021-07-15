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

## Function to process barplot data ##
process_barplot_data <- function(input_genelist, He_dataset, Maynard_dataset, 
                                 Allen_dataset) {
  
  he_data <- He_dataset %>%
    filter(gene_symbol %in% input_genelist) %>%
    add_column(source_dataset = "He")
  maynard_data <- Maynard_dataset %>%
    filter(gene_symbol %in% input_genelist) %>%
    add_column(source_dataset = "Maynard")
  allen_MTG_data <- Allen_dataset %>%
    filter(gene_symbol %in% input_genelist) %>% 
    mutate(source_dataset = paste("ABI", class_label, sep = "_")) %>%
    select(-class_label) %>% 
    select(gene_symbol, L1, L2, L3, L4, L5, L6, WM, source_dataset)
  
  barplot_data <- rbind(he_data, maynard_data, 
                        allen_MTG_data) %>%
    pivot_longer(cols = L1:WM,
                 names_to = "layer",
                 values_to = "expression") %>%
    mutate(source_dataset = factor(source_dataset, 
                                   levels = c("He", "Maynard", "ABI_GABAergic", 
                                              "ABI_Glutamatergic", 
                                              "ABI_Non-neuronal")))
  return(barplot_data)
}

separate_layers <- function(input_table, input_genelist, source) {
  
  layer_markers <- list()
  for(i in c(1:7)) {
    layer <- input_table %>%
      filter(stri_detect_fixed(source_dataset, source)) %>%
      filter(gene_symbol %in% input_genelist, layer_marker == i) %>%
      pull(gene_symbol)
    layer_markers[[i]] <- layer
  }
  
  names(layer_markers) <- c("L1", "L2", "L3", "L4", "L5", "L6", "WM")

  return(layer_markers)
}

## Function to process heatmap data for bulk-tissue and sn-RNAseq ##
process_heatmap_data <- function(source, source_dataset, input_genelist, 
                                 cell_type = NA) {
  
  # If processing Allen data:
  if (source == "Allen") {
    # Filter source data to only take specific cell_type
    source_dataset %<>% filter(class_label == cell_type) %>%
      select(-class_label)
    # Create heatmap data
    processed_heatmap_data <- source_dataset %>%
      # Filter out duplicate genes
      distinct(gene_symbol, .keep_all = TRUE) %>%
      column_to_rownames(var = "gene_symbol") %>%
      t() %>% scale() %>% t() %>% as.data.frame() %>%
      rownames_to_column(var = "gene_symbol") %>%
      # Filter for genes in inputted gene list
      filter(gene_symbol %in% input_genelist) %>%
      # Lengthen wide data
      pivot_longer(cols = L1:WM, names_to = "cortical_layer_label",
                   values_to = "expression") %>%
      # Rename cortical_layer_label column to layer
      rename(layer = cortical_layer_label)
    # If processing bulk-tissue (i.e., He or Maynard)
  } else if (source == "He" | source == "Maynard") {
    # Create heatmap data
    processed_heatmap_data <- source_dataset %>%
      column_to_rownames(var = "gene_symbol") %>%
      t() %>% scale() %>% t() %>% as.data.frame() %>%
      rownames_to_column(var = "gene_symbol") %>%
      # Filter for genes in inputted gene list
      filter(gene_symbol %in% input_genelist) %>%
      # Filter out duplicate genes
      distinct(gene_symbol, .keep_all = TRUE) %>%
      as.matrix()
    
    # Set NaN values to 0
    processed_heatmap_data[is.nan(processed_heatmap_data)] <- 0
    
    processed_heatmap_data %<>% 
      as_tibble() %>%
      #Filter out genes that have NA for all layers
      filter_at(vars(L1, L2, L3, L4, L5, L6, WM),
                all_vars(!is.na(.))) %>%
      mutate_if(is.numeric,as.character, is.factor, as.character) %>%
      # Lengthen wide data
      pivot_longer(cols = L1:WM, names_to = "layer", values_to = "expression"
      ) %>%
      # Select gene_symbol, layer, expression columns in order
      select(gene_symbol, layer, expression) %>%
      # Mutate expression column to be numeric (just in case it's chr)
      mutate_at("expression", as.numeric)
  }
  
  #Order data according to similarity in expression profile
  ordered_data <- source_dataset %>%
    # Filter out for genes in gene list
    filter(gene_symbol %in% input_genelist) %>%
    # Replace NA values with 0
    replace_na(list(
      L1 = 0, L2 = 0, L3 = 0, L4 = 0, L5 = 0,
      L6 = 0, WM = 0)) %>%
    mutate_if(is.numeric,as.character, is.factor, as.character) %>%
    select("gene_symbol":"WM")
  
  # Code to order data according to similarity in expression profile
  #ordered_data_matrix <- as.matrix(ordered_data)
  #rownames(ordered_data) <- ordered_data$gene_symbol
  #data_dendro <- as.dendrogram(hclust(d = dist(x = ordered_data_matrix)))
  
  #data_order <- order.dendrogram(data_dendro)
  
  # Set order in processed heatmap data to match the order
  #processed_heatmap_data$gene_symbol <- 
  #  factor(x = processed_heatmap_data$gene_symbol,
  #         levels = ordered_data$gene_symbol[data_order], 
  #         ordered = TRUE)
  
  processed_heatmap_data$layer <- with(processed_heatmap_data, 
                                       factor(layer, 
                                              levels = rev(sort(unique(layer)))))
  
  # Return ordered heatmap data in long format
  return(processed_heatmap_data)
}

## Functions for generating singe or multi-gene correlations across He & Maynard datasets

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
  
  He_Maynard_diag_gene <- cor(He_gene, Maynard_gene, method = "pearson") %>%
    as_tibble() %>%
    pull(var = 1) %>%
    as.numeric()
  
  return(format(signif(He_Maynard_diag_gene, digits = 3)))
}

multi_gene_correlation <- function(input_genes, He_dataset, Maynard_dataset) {
  
  He_genes <- He_dataset %>%
    filter(gene_symbol %in% input_genes) %>%
    arrange(gene_symbol) %>%
    column_to_rownames(var = "gene_symbol") %>%
    filter(across(everything(), ~ !is.na(.))) %>%
    t()
  
  Maynard_genes <- Maynard_dataset %>%
    filter(gene_symbol %in% input_genes) %>%
    arrange(gene_symbol) %>%
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
quantile_distribution <- function(dataset_correlation, selected_gene_correlation) {
  
  # Create quantiles using ecdf()
  quantile_distribution <- ecdf(dataset_correlation)
  # Query quantile of the selected gene's correlation
  query_quantile <- quantile_distribution(selected_gene_correlation)
  # Multiply by 100
  quantile <- query_quantile * 100
  
  return(format(floor(quantile)))
}

## Function for creating p-value for correlation of user-selected genes
wilcoxtest <- function(input_genelist, He_dataset, Maynard_dataset, He_Maynard_diagonal) {
  
  He_genes <- He_dataset %>%
    filter(gene_symbol %in% input_genelist) %>%
    arrange(gene_symbol) %>%
    column_to_rownames(var = "gene_symbol") %>%
    t()
  
  Maynard_genes <- Maynard_dataset %>%
    filter(gene_symbol %in% input_genelist) %>%
    arrange(gene_symbol) %>%
    column_to_rownames(var = "gene_symbol") %>%
    t() 
  if (length(input_genelist) == 1) {
    return(format(signif(cor.test(He_genes[,1], Maynard_genes[,1])$p.value), digits = 4))
  }
  # Pearson's correlation test for the user-selected genes
  He_Maynard_genes_cormatrix <- cor(He_genes, Maynard_genes, method = "pearson")
  # Get the diagonal (correlation of each gene against itself)
  He_Maynard_diag_genes <- diag(He_Maynard_genes_cormatrix, names = TRUE)
  
  # P-value on Pearson's
  wilcoxtest <- wilcox.test(He_Maynard_diag_genes, He_Maynard_diagonal)$p.value
  
  return(format(signif(wilcoxtest, digits = 4)))
}

#Function for ranking bulk tissue datasets
rank_bulk_dataset <- function(source_dataset) {
  dataset_ranked <- source_dataset %>%
    select(L1:WM) %>%
    map_df(rank, ties.method = "min")
  dataset_ranked$gene_symbol <- source_dataset$gene_symbol
  return(dataset_ranked)
}

#Function for ranking scRNA datasets
rank_AIBS_dataset <- function(source_dataset) {
  dataset_ranked <- source_dataset 
  gene_symbol <- dataset_ranked$gene_symbol
  dataset_ranked %<>% select(L1:L6) %>% map_df(rank, ties.method = "min") %>%
    add_column(gene_symbol = gene_symbol)
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
# from https://github.com/sarbal/EGAD/blob/master/R/auroc_analytic.R
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
  ranked_cell_type_df <- rank_AIBS_dataset(cell_type_df)
  indices <- return_indices(ranked_cell_type_df, multiple_genelist)
  ranked_cell_type_df %<>% select(-gene_symbol)
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
           # Multiple test correction using bonferroni
           adjusted_P_He = signif(p.adjust(pValue_He, method = "bonferroni"), digits = 3),
           adjusted_P_Maynard = signif(p.adjust(pValue_Maynard, method = "bonferroni"), digits = 3)) %>%
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
  return(AUROC)
}

## Function for generating scRNA AUROC and p-value data
AUROC_AIBS <- function(source_dataset, multiple_genelist) {
  scRNA_AUROC_list <- list()
  for (i in unique(source_dataset$class_label)) {
    scRNA_AUROC_list[[i]] <- auroc_cell_type(source_dataset, multiple_genelist, i)
  }
  scRNA_AUROC_table <- rbindlist(scRNA_AUROC_list, idcol = "rownames") %>%
    mutate(rownames = c("AUROC_GABA", "pValue_GABA", "AUROC_GLUT", 
                        "pValue_GLUT", "AUROC_NON", "pValue_NON")) %>%
    column_to_rownames(var = "rownames") %>% 
    add_column(WM = NA) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Layer")
  scRNA_pValue_table <- scRNA_AUROC_table %>% 
    select(starts_with("pValue"), Layer) %>%
    pivot_longer(cols = starts_with("pValue"),
                 names_to = c("p_Value"),
                 values_to = "pValue") %>%
    mutate(p_Value = gsub("pValue_", "", p_Value)) %>%
    rename(class_label = p_Value) %>%
    mutate(pValue = p.adjust(pValue, method = "bonferroni")) %>%
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
AUROC_data <- function(AIBS_AUROC_table, bulk_AUROC_table) {
  AUROC_table <- rbind(AIBS_AUROC_table, bulk_AUROC_table) %>%
    rename(dataset = class_label) %>%
    mutate(dataset = gsub("GABA", "scRNA_GABA", dataset)) %>%
    mutate(dataset = gsub("GLUT", "scRNA_GLUT", dataset)) %>%
    mutate(dataset = gsub("NON", "scRNA_Non-neuronal", dataset)) %>%
    mutate(dataset = gsub("AUROC ", "", dataset)) %>%
    mutate(dataset = gsub("\\(", "", dataset)) %>%
    mutate(dataset = gsub("\\)", "", dataset))
  AUROC_table$Layer <- factor(AUROC_table$Layer, 
                              levels = c("WM", "L6", "L5", "L4", "L3", "L2", "L1"))
  return(AUROC_table)
}

