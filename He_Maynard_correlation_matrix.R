library(here)
library(tidyverse)
library(magrittr)

load(here("data", "processed", "He_DS1_Human_averaged.Rdata"), verbose = TRUE)
load(here("data", "processed", "Maynard_dataset_average.Rdata"), verbose = TRUE)

He_Maynard_common_gene_list <- intersect(He_DS1_averaged_by_layer$gene_symbol, Maynard_dataset_average$gene_symbol)

He_cormatrix <- He_DS1_averaged_by_layer %>%
  filter(gene_symbol %in% He_Maynard_common_gene_list) %>%
  arrange(gene_symbol) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t()

Maynard_cormatrix <- Maynard_dataset_average %>%
  select(gene_symbol:WM) %>%
  filter(gene_symbol %in% He_Maynard_common_gene_list) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t()

He_Maynard_cormatrix <- cor(He_cormatrix, Maynard_cormatrix, method = "pearson") 
He_Maynard_diag <- diag(He_Maynard_cormatrix, names = TRUE)

quantile(He_Maynard_diag, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)