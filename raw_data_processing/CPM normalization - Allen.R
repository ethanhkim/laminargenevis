## Allen SMART-seq Multiple Regions CPM normalization ##

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(purrr)
library(data.table)
library(moments)
library(here)
library(magrittr)
library(edgeR)
library(conflicted)

# Set conflicts
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("transpose", "data.table")
conflict_prefer("cpm", "edgeR")

# If working on SCC ----
if (getwd() == "/external/mgmt3/genome/scratch/Neuroinformatics/ekim/transcriptome_project") {
  # Check if Allen count matrix and metadata exist
  if (!file.exists(here("Data", "raw_data", "Allen",
                        "singlecellMatrix.csv")) == TRUE) {
    singlecellMatrix_dest <- here("Data", "raw_data", "Allen", "singlecellMatrix.csv")
    singlecellMatrix_url <- "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/matrix.csv"
    allen_singlecellMatrix <- download.file(singlecellMatrix_url, singlecellMatrix_dest)
    
    rm(singlecellMatrix_dest, singlecellMatrix_url, allen_singlecellMatrix)
  }
  
  if (!file.exists(here("Data", "raw_data", "Allen",
                        "singlecellMetadata.csv")) == TRUE) {
    # Download Allen sc-RNAseq matrix and metadata for CAMH SCC
    singlecellMetadata_dest <- here("Data", "raw_data", "Allen", "singlecellMetadata.csv")
    singlecellMetadata_url <- "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/metadata.csv"
    allen_singlecellMetadata <- download.file(singlecellMetadata_url, singlecellMetadata_dest)
    
    rm(singlecellMetadata_dest, singlecellMetadata_url, allen_singlecellMetadata)
  }
  
  # Functions to use #
  
  # Create AIBS data with merged metadata and remove outliers ----
  clean_AIBS_data <- function(AIBS_countMatrix, AIBS_metadata) {
    metadata <- AIBS_metadata
    matrix <- AIBS_countMatrix
    
    matrix$outlier_call <- metadata$outlier_call
    matrix$class_label <- metadata$class_label
    matrix$region_label <- metadata$region_label
    matrix$cortical_layer_label <- metadata$cortical_layer_label
    matrix$external_donor_name_label <- metadata$external_donor_name_label
    
    matrix %<>% filter(outlier_call == FALSE) %>%
      select(-outlier_call) 
    return(matrix)
  }
  
  # Downsample
  sample_nuclei <- function(df, layer_label, n_sample) {
    
    set.seed(0)
    df %<>% filter(cortical_layer_label == layer_label)
    sampled_layer <- df[sample(nrow(df), n_sample),]
    
    return(sampled_layer)
  }
  
  # Summarize gene counts by layer and cell-type
  sum_gene_count <- function(region_separated_df, aggregate_level) {
    
    df <- region_separated_df
    
    # Lengthen dataframe
    df %<>%
      select(sample_name, class_label, cortical_layer_label,
             external_donor_name_label, everything()) %>%
      # Lengthen dataframe to:
      # gene_symbol | count | everything()
      data.table::melt(id.vars = c("sample_name", "class_label", 
                                   "cortical_layer_label",
                                   "external_donor_name_label"),
                       variable.factor = FALSE,
                       variable.name = "gene_symbol",
                       value.name = "count")
    
    # If aggregating at the gene level:
    if (aggregate_level == "gene") {
      df %<>% group_by(cortical_layer_label, gene_symbol) %>%
        # Sum up all gene counts at a layer and gene level
        summarize(countSum = sum(count)) %>%
        # Widen dataframe to gene (row) by layer (column)
        spread(cortical_layer_label, countSum)
    } else if (aggregate_level == "cell_type") {
      df %<>%
        # Group by layer, cell type (class_label) and gene 
        group_by(cortical_layer_label, class_label,
                 gene_symbol) %>%
        # Sum up all gene counts at a layer and cell type level
        summarize(countSum = sum(count)) %>%
        # Widen dataframe to gene (row) by layer (column)
        spread(cortical_layer_label, countSum)
    }
    
    return(df)
  }
  
  # Read in data and metadata ----
  AIBS_metadata <- fread(here("Data", "raw_data", "Allen", 
                              "singlecellMetadata.csv"), header = T) %>%
    select(sample_name, external_donor_name_label, class_label, 
           subclass_label, region_label, cortical_layer_label, 
           outlier_call) %>%
    as_tibble()
  AIBS_matrix <- fread(here("Data", "raw_data", "Allen",
                            "singlecellMatrix.csv"), header = T)
  
  # Create cleaned dataframe and remove source data ----
  AIBS_cleaned_df <- clean_AIBS_data(AIBS_matrix, AIBS_metadata)
  rm(AIBS_matrix, AIBS_metadata)
  
  # Separate by region (MTG) ----
  MTG <- AIBS_cleaned_df %>%
    filter(region_label == "MTG") %>%
    select(-region_label)
  rm(AIBS_cleaned_df)
  
  # Downsample to lowest number of nuclei per cell type
  cell_type_specific_data <- list()
  cell_type_nuclei_count <- list()
  for (cell_type in c('GABAergic', 'Glutamatergic', 'Non-neuronal')) {
    cell_type_specific_data[[cell_type]] <- MTG %>%
      filter(class_label == cell_type)
    cell_type_nuclei_count[[cell_type]] <- cell_type_specific_data[[cell_type]] %>%
      group_by(cortical_layer_label) %>%
      count()
  }
  
  cell_type_nuc_count <- tibble(
    GABAergic = cell_type_nuclei_count[['GABAergic']] %>% pull(),
    Glutamatergic = cell_type_nuclei_count[['Glutamatergic']] %>% pull(),
    Non_neuronal = cell_type_nuclei_count[['Non-neuronal']] %>% pull()
  )
  
  GABA_sampled <- list()
  GLUT_sampled <- list()
  NONN_sampled <- list()
  for (layer in c("L1", "L2", "L3", "L4", "L5", "L6")) {
    GABA_sampled[[layer]] <- sample_nuclei(cell_type_specific_data[['GABAergic']], layer, 381)
    GLUT_sampled[[layer]] <- sample_nuclei(cell_type_specific_data[['Glutamatergic']], layer, 274)
    NONN_sampled[[layer]] <- sample_nuclei(cell_type_specific_data[['Non-neuronal']], layer, 125)
  }
  
  GABA_sampled_df <- rbindlist(GABA_sampled)
  GLUT_sampled_df <- rbindlist(GLUT_sampled)
  NONN_sampled_df <- rbindlist(NONN_sampled)
  
  Allen_downsampled_df <- rbind(GABA_sampled_df, GLUT_sampled_df, NONN_sampled_df)
  Allen_downsampled_cell_sum_count <- sum_gene_count(Allen_downsampled_df, "cell_type")
  write.csv(Allen_downsampled_cell_sum_count, 
            file = here("Data", "raw_data", "Allen", "Allen_downsampled_cell_sum_count.csv"))
  
  # Aggregate MTG data ----
  ## Aggregate at layer level across each cell types
  MTG_cell_type_sum_count <- sum_gene_count(MTG, "cell_type")
  ## Aggregate at layer level across all cell types
  MTG_gene_sum_count <- sum_gene_count(MTG, "gene")
  
  # Write aggregate data as .csv
  write.csv(MTG_cell_type_sum_count, 
            file = here("Data", "raw_data", "Allen", "MTG_cell_type_sum_count.csv"))
  write.csv(MTG_gene_sum_count, 
            file = here("Data", "raw_data", "Allen", "MTG_gene_sum_count.csv"))

  # If working on local ----
} else if (getwd() == "/Users/ethankim/Google Drive/Desk_Laptop/U of T/Grad School/French Lab/transcriptome_project") {
  
  # Load data aggregated at cell type AND gene level ----
  MTG_cell_type_sum_count <- fread(here("Data", "raw_data", "Allen", "MTG_cell_type_sum_count.csv")) %>%
    select(-V1)
  MTG_gene_sum_count <- fread(here("Data", "raw_data", "Allen", "MTG_gene_sum_count.csv")) %>%
    select(-V1)
  Allen_downsampled_df <- fread(here("Data", "raw_data", "Allen", "Allen_downsampled_cell_sum_count.csv"))
  
}


## CPM normalize gene-level aggregate data ----

# Non-filtered data:
Allen_gene_logCPM_dataset <- MTG_gene_sum_count %>%
  # Remove gene_symbol column for cpm()
  select(-gene_symbol) %>% 
  # Add +1 to remove 0's for log transformation
  mutate(across(where(is.numeric), ~. +1)) %>%
  # CPM normalize with log = T
  cpm(log = T) %>% as.data.frame() %>%
  # Add back in gene symbols
  add_column(gene_symbol = MTG_gene_sum_count$gene_symbol)

# Filtered data: CPM > 0.1
Allen_gene_logCPM_filtered_dataset <- MTG_gene_sum_count %>%
  # Add one to counts to avoid taking cpm of 0
  mutate_at(c("L1", "L2", "L3", "L4", "L5", "L6"), ~. +1) %>%
  column_to_rownames(var = "gene_symbol") %>%
  cpm() %>%
  as.data.frame() %>%
  # Retain gene symbols, as filter removes rownames
  rownames_to_column(var = "gene_symbol") %>%
  # Filter out samples of CPM < 0.1
  filter_at(vars(-gene_symbol), all_vars(. > .1)) %>%
  column_to_rownames(var = "gene_symbol")
names <- rownames(Allen_gene_logCPM_filtered_dataset)
Allen_gene_logCPM_filtered_dataset %<>%
  # Take log2 of CPM
  map_df(log2) %>%
  select(L1, L2, L3, L4, L5, L6) %>%
  as.data.frame() %>%
  add_column(gene_symbol = names) %>%
  select(gene_symbol, everything())

## CPM normalize cell-type-level aggregate data ----

# Non-filtered data:
Allen_logCPM_dataset <- MTG_cell_type_sum_count %>%
  # Add one to counts to avoid taking cpm of 0
  mutate_at(c("L1", "L2", "L3", "L4", "L5", "L6"), ~. +1) %>%
  unite(gene_class, c("gene_symbol", "class_label")) %>%
  column_to_rownames(var = "gene_class") %>%
  cpm(log = T) %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_class") %>%
  separate(gene_class, into = c("gene_symbol", "class_label"),
           sep = "_") %>%
  add_column(WM = NA) %>%
  select(gene_symbol, class_label, L1, L2, L3, L4, L5, L6, WM)

# Filtered data: CPM > 0.1
Allen_logCPM_filtered_dataset <- MTG_cell_type_sum_count %>%
  # Add one to counts to avoid taking cpm of 0
  mutate_at(c("L1", "L2", "L3", "L4", "L5", "L6"), ~. +1) %>%
  unite(gene_class, c("gene_symbol", "class_label")) %>%
  column_to_rownames(var = "gene_class") %>%
  cpm() %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_class") %>%
  separate(gene_class, into = c("gene_symbol", "class_label"),
           sep = "_") %>%
  # Filter out samples of CPM < 0.1
  filter_at(vars(-gene_symbol, -class_label), all_vars(. > .1)) %>%
  unite(gene_class, c("gene_symbol", "class_label")) %>%
  column_to_rownames(var = "gene_class")
names <- rownames(Allen_logCPM_filtered_dataset)
Allen_logCPM_filtered_dataset %<>%
  # Take log2 of CPM
  map_df(log2) %>%
  select(L1, L2, L3, L4, L5, L6) %>%
  # Take z-score (for app)
  t() %>% scale() %>% t() %>%
  as.data.frame() %>%
  add_column(gene_class = names, WM = NA) %>%
  separate(gene_class, into = c("gene_symbol", "class_label"),
           sep = "_") %>%
  select(gene_symbol, class_label, L1, L2, L3, L4, L5, L6, WM)

## CPM normalize downsampled data ----
Allen_downsampled_logCPM_dataset <- Allen_downsampled_df %>%
  # Add one to counts to avoid taking cpm of 0
  mutate_at(c("L1", "L2", "L3", "L4", "L5", "L6"), ~. +1) %>%
  unite(gene_class, c("gene_symbol", "class_label")) %>%
  column_to_rownames(var = "gene_class") %>%
  cpm(log = T) %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_class") %>%
  separate(gene_class, into = c("gene_symbol", "class_label"),
           sep = "_") %>%
  add_column(WM = NA) %>%
  select(gene_symbol, class_label, L1, L2, L3, L4, L5, L6, WM)


# Filtered data: CPM > 0.1
Allen_downsampled_logCPM_filtered_dataset <- Allen_downsampled_df %>%
  # Add one to counts to avoid taking cpm of 0
  mutate_at(c("L1", "L2", "L3", "L4", "L5", "L6"), ~. +1) %>%
  unite(gene_class, c("gene_symbol", "class_label")) %>%
  column_to_rownames(var = "gene_class") %>%
  cpm() %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_class") %>%
  separate(gene_class, into = c("gene_symbol", "class_label"),
           sep = "_") %>%
  # Filter out samples of CPM < 0.1
  filter_at(vars(-gene_symbol, -class_label), all_vars(. > .1)) %>%
  unite(gene_class, c("gene_symbol", "class_label")) %>%
  column_to_rownames(var = "gene_class")
names <- rownames(Allen_downsampled_logCPM_filtered_dataset)
Allen_downsampled_logCPM_filtered_dataset %<>%
  # Take log2 of CPM
  map_df(log2) %>%
  select(L1, L2, L3, L4, L5, L6) %>%
  # Take z-score (for app)
  t() %>% scale() %>% t() %>%
  as.data.frame() %>%
  add_column(gene_class = names, WM = NA) %>%
  separate(gene_class, into = c("gene_symbol", "class_label"),
           sep = "_") %>%
  select(gene_symbol, class_label, L1, L2, L3, L4, L5, L6, WM)


# Save logCPM data ----
save(Allen_gene_logCPM_dataset, file = here("Data", "processed_data", "Allen_gene_logCPM_dataset.Rdata"))
save(Allen_gene_logCPM_filtered_dataset, 
     file = here("Data", "processed_data", "Allen_gene_logCPM_filtered_dataset.Rdata"))
save(Allen_logCPM_dataset, file = here("Data", "processed_data", "Allen_logCPM_dataset.Rdata"))
save(Allen_logCPM_filtered_dataset, 
     file = here("Data", "processed_data", "Allen_logCPM_filtered_dataset.Rdata"))
save(Allen_downsampled_logCPM_dataset, file = here("Data", "processed_data", "Allen_downsampled_logCPM_dataset.Rdata"))
save(Allen_downsampled_logCPM_filtered_dataset, 
     file = here("Data", "processed_data", "Allen_downsampled_logCPM_filtered_dataset.Rdata"))

