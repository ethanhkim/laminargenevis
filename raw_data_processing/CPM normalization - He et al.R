## He et al Normalization ##

# Load in libraries
library(edgeR)
library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(here)

# Read in raw count matrix
He_count_matrix <- fread(file = here("Data", "raw_data", "He et al", "he_human_102_counts_matrix.csv"))

# Clean up column names
names(He_count_matrix) <- gsub(x = names(He_count_matrix), pattern = "\\-", replacement = "\\_") 

# Duplicate gene symbols exist - take the average across each section
He_count_matrix %<>%
  group_by(gene_symbol) %>%
  summarise(across(.fns = mean))

# Separate by DS1
He_DS1_matrix <- He_count_matrix %>%
  select(gene_symbol, contains("DS1")) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "sample_metadata")
  
# Clean up sample metadata column
He_DS1_matrix %<>%
  mutate(cleaned_metadata = str_replace(sample_metadata, ".*1_", ""),
         cleaned_metadata = str_replace(cleaned_metadata, "S1", "S01"),
         cleaned_metadata = str_replace(cleaned_metadata, "S010", "S10"),
         cleaned_metadata = str_replace(cleaned_metadata, "S011", "S11"),
         cleaned_metadata = str_replace(cleaned_metadata, "S012", "S12"),
         cleaned_metadata = str_replace(cleaned_metadata, "S013", "S13"),
         cleaned_metadata = str_replace(cleaned_metadata, "S014", "S14"),
         cleaned_metadata = str_replace(cleaned_metadata, "S015", "S15"),
         cleaned_metadata = str_replace(cleaned_metadata, "S016", "S16"),
         cleaned_metadata = str_replace(cleaned_metadata, "S017", "S17"),
         cleaned_metadata = str_replace(cleaned_metadata, "S018", "S18")) %>%
  select(-sample_metadata) %>% select(cleaned_metadata, everything())
# Transpose back to original dataframe (sample [col] by gene [row])
He_DS1_matrix %<>% 
  data.table::transpose(
    keep.names = "gene_symbol",
    make.names = names(He_DS1_matrix[1])
    ) %>% as.data.frame()

# Average across slices for participants
He_DS1_matrix_slice_sum <- data.frame(matrix(NA, ncol = 59453, nrow = 1)[-1,])
for (slice in c("S01", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10",
                   "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18")) {
  He_DS1_matrix_slice_sum[slice,] <- rowSums(select(He_DS1_matrix, contains(slice)), na.rm = T)
}
rm(slice)

# Sanity check - did it average correctly? Check using S01-containing columns
test <- He_DS1_matrix %>% select(gene_symbol, contains("S01")) %>%
  column_to_rownames(var = "gene_symbol") 
test %<>%
  mutate(count = rowSums(select(test, contains("S01")), na.rm = T)) %>%
  select(count)
# Select first row - should be the product of rowMeans of S01
test2 <- He_DS1_matrix_slice_sum[1,] %>%
  t() %>%
  as_tibble()
# Check to see if all numbers in vector are equal
all(test == test2)
# Passes. Remove test items
rm(list = c("test", "test2"))

# Create layer representations from He et al figure
He_DS1_transposed <- He_DS1_matrix_slice_sum %>% 
  t() %>% as_tibble() %>% 
  # Add back in gene symbol column
  add_column(gene_symbol = He_DS1_matrix$gene_symbol) %>%
  # Rename S01 to S1, and reorder columns
  rename(S1 = S01) %>% select(gene_symbol, everything()) %>%
  # Create weighted columns, according to He et al figure
  mutate(S4_weighted = S4 * 0.5) %>%
  mutate(S9_weighted = S9 * 0.5)

# Function to create a weighted average column
wt_sum_col_fn <- function(transposed_df, cols_to_sum) {
  df <- transposed_df
  # Create wt.mean column - averages across columns
  sum_col <- df %>%
    select(all_of(cols_to_sum)) %>%
    data.table::transpose() %>%
    colSums()
  # Remove name attributes
  names(sum_col) <- NULL
  return(sum_col)
}

# Create averaged He et al. tibble from the scaled values. Don't include gene_symbol
# as it needs to go through CPM normalization (requires no rownames)
He_DS1_sum_layer <- tibble(
  gene_symbol = He_DS1_transposed$gene_symbol,
  L1 = He_DS1_transposed$S1,
  L2 = wt_sum_col_fn(He_DS1_transposed, c("S2", "S3", "S4_weighted")),
  L3 = wt_sum_col_fn(He_DS1_transposed, c("S4_weighted", "S5", "S6")),
  L4 = wt_sum_col_fn(He_DS1_transposed, c("S7", "S8", "S9_weighted")),
  L5 = wt_sum_col_fn(He_DS1_transposed, c("S9_weighted", "S10", "S11", "S12")),
  L6 = wt_sum_col_fn(He_DS1_transposed, c("S13","S14", "S15", "S16")),
  WM = He_DS1_transposed$S17
)

# CPM normalize with log = T
He_DS1_logCPM_dataset <- He_DS1_sum_layer %>%
  # Remove gene_symbol column for cpm()
  select(-gene_symbol) %>% 
  # Add +1 to remove 0's for log transformation
  mutate(across(where(is.numeric), ~. +1)) %>%
  # CPM normalize with log = T
  cpm(log = T) %>% as.data.frame() %>%
  # Add back in gene symbols
  add_column(gene_symbol = He_DS1_sum_layer$gene_symbol) %>%
  select(gene_symbol, everything())

# CPM normalize, filter out CPM < 0.1
He_DS1_logCPM_filtered_dataset <- He_DS1_sum_layer %>%
  # Add one to counts to avoid taking cpm of 0
  mutate_at(c("L1", "L2", "L3", "L4", "L5", "L6", "WM"), ~. +1) %>%
  column_to_rownames(var = "gene_symbol") %>%
  cpm() %>%
  as.data.frame() %>%
  # Retain gene symbols, as filter removes rownames
  rownames_to_column(var = "gene_symbol") %>%
  # Filter out samples of CPM < 0.1
  filter_at(vars(-gene_symbol), all_vars(. > .1)) %>%
  column_to_rownames(var = "gene_symbol")
names <- rownames(He_DS1_logCPM_filtered_dataset)
He_DS1_logCPM_filtered_dataset %<>%
  # Take log2 of CPM
  map_df(log2) %>%
  select(L1, L2, L3, L4, L5, L6, WM) %>%
  # Take z-score (for app)
  t() %>% scale() %>% t() %>%
  as.data.frame() %>%
  add_column(gene_symbol = names) %>%
  select(gene_symbol, everything())


# Clean up remaining DS1 data
rm(He_DS1_matrix, He_DS1_sum_layer, He_DS1_transposed,
   He_DS1_matrix_slice_sum, He_count_matrix, wt_sum_col_fn,
   slice, WM)

# Write normalized data as .Rdata
save(He_DS1_logCPM_dataset, file = here("Data", "processed_data", 
                                         "He_DS1_logCPM_dataset.Rdata"))
save(He_DS1_logCPM_filtered_dataset, file = here("Data", "processed_data",
                                                 "He_DS1_logCPM_filtered_dataset.Rdata"))

# Write normalized data as .csv
write.csv(He_DS1_logCPM_dataset, file = here("Data", "processed_Data", 
                                              "He_DS1_logCPM_dataset.csv"))

