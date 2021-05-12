## Maynard et al CPM normalization ##

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("LieberInstitute/spatialLIBD")
library(spatialLIBD)
library(data.table)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
library(purrr)
library(magrittr)
library(stringr)
library(conflicted)
library(biomaRt)
library(HGNChelper)
library(edgeR)
library(here)

# Set conflicts #
conflict_prefer('select', 'dplyr')
conflict_prefer('filter', 'dplyr')
conflict_prefer('cpm', 'edgeR')
conflict_prefer('rename', 'dplyr')

# Fetch data from spatialLIBD
sce_layer <- fetch_data(type = "sce_layer")
# Get raw UMI counts per gene
Maynard_dataset <- sce_layer@assays@data@listData$counts %>%
  as_tibble(rownames=NA)
# Get HGNC gene symbol from ensembl ID's
Maynard_ensembl_list <- sce_layer@rowRanges@ranges@NAMES
Maynard_dataset$gene_symbol <- mapIds(org.Hs.eg.db, 
                                      keys = Maynard_ensembl_list, keytype = "ENSEMBL", column="SYMBOL")
Maynard_dataset %<>%
  select(gene_symbol, everything()) %>% filter(!is.na(gene_symbol)) %>%
  group_by(gene_symbol) %>%
  summarise(across(.fns = mean))

# Write layer-level data for working on SCC
write.csv(Maynard_dataset, file = here("Data", "raw_data", "Maynard et al", 
                                       "layer_level_data.csv"))


# Select only individuals with all cortical layers (n = 2)
Maynard_dataset_subset <- Maynard_dataset %>%
  select("gene_symbol", contains('151507'), contains('151508'), contains('151509'),
         contains('151510'), contains('151673'), contains('151674'), 
         contains('151675'), contains('151676')) %>%
  column_to_rownames(var = "gene_symbol")

# Function to create columns of mean values across individuals for Maynard
create_sum_col <- function(Maynard_subset, layer_label) {
  df <- Maynard_subset
  layer <- layer_label
  
  df %<>% select(contains(layer)) %>%
    mutate(countSum = rowSums(.)) %>%
    pull(countSum)
  
  return(df)
}

# Create list of mean columns
Maynard_sum_col_list <- list()
for (label in c("Layer1", "Layer2", "Layer3", "Layer4", "Layer5", "Layer6", "WM")) {
  Maynard_sum_col_list[[label]] <- create_sum_col(Maynard_dataset_subset, label)
}
# Bind list and convert to dataframe
Maynard_dataset_averaged <- as.data.frame(do.call(rbind, Maynard_sum_col_list)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = 'gene_symbol') %>%
  rename(L1 = Layer1, L2 = Layer2, L3 = Layer3, L4 = Layer4, L5 = Layer5, 
         L6 = Layer6)

# Normalize Maynard UMI counts with CPM, log = T
Maynard_logCPM_dataset <- Maynard_dataset_averaged %>%
  # Remove gene_symbol column for cpm()
  select(-gene_symbol) %>% 
  # Add +1 to remove 0's for log transformation
  mutate(across(where(is.numeric), ~. +1)) %>%
  # CPM normalize with log = T
  cpm(log = T) %>% as.data.frame() %>%
  # Add back in gene symbols
  add_column(gene_symbol = Maynard_dataset_averaged$gene_symbol) %>%
  select(gene_symbol, everything())

Maynard_logCPM_filtered_dataset <- Maynard_dataset_averaged %>%
  # Add one to counts to avoid taking cpm of 0
  mutate_at(c("L1", "L2", "L3", "L4", "L5", "L6", "WM"), ~. +1) %>%
  column_to_rownames(var = "gene_symbol") %>%
  cpm() %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_symbol") %>%
  # Filter out samples of CPM < 0.1
  filter_at(vars(-gene_symbol), all_vars(. > .1)) %>%
  column_to_rownames(var = "gene_symbol")
names <- rownames(Maynard_logCPM_filtered_dataset)
Maynard_logCPM_filtered_dataset %<>%
  # Take log2 of CPM
  map_df(log2) %>%
  select(L1, L2, L3, L4, L5, L6, WM) %>%
  # Take z-score (for app)
  t() %>% scale() %>% t() %>%
  as.data.frame() %>%
  add_column(gene_symbol = names) %>%
  select(gene_symbol, everything())

# Clean up workspace
rm(Maynard_dataset, Maynard_dataset_averaged, Maynard_dataset_subset,
   Maynard_sum_col_list, label, Maynard_ensembl_list,
   create_sum_col, sce_layer)

# Write normalized data as .Rdata
save(Maynard_logCPM_dataset, file = here("Data", "processed_data", 
                                         "Maynard_logCPM_dataset.Rdata"))
save(Maynard_logCPM_filtered_dataset, file = here("Data", "processed_data", 
                                         "Maynard_logCPM_filtered_dataset.Rdata"))

# Write normalized data as .csv
write.csv(Maynard_logCPM_dataset, file = here("Data", "processed_Data", 
                                               "Maynard_logCPM_dataset.csv"))


