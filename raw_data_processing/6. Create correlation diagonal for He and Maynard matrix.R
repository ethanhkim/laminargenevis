## Create diagonal for multiple gene correlation

library(tidyverse)
library(magrittr)

load(here("data", "processed", "He_DS1_Human_averaged.Rdata"), verbose = TRUE)
load(here("data", "processed", "Maynard_dataset_average.Rdata"), verbose = TRUE)

He_Maynard_common_gene_list <- intersect(He_DS1_Human_averaged$gene_symbol, Maynard_dataset_average$gene_symbol)

He_genes <- He_DS1_Human_averaged %>%
  select(gene_symbol:WM) %>%
  filter(gene_symbol %in% He_Maynard_common_gene_list) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t()

Maynard_genes <- Maynard_dataset_average %>%
  select(gene_symbol:WM) %>%
  filter(gene_symbol %in% He_Maynard_common_gene_list) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t() 

He_Maynard_genes_cormatrix <- cor(He_genes, Maynard_genes, method = "pearson")
He_Maynard_diag_genes <- diag(He_Maynard_genes_cormatrix, names = TRUE) %>%
  as_tibble()

write_csv(He_Maynard_diag_genes, here("data", "processed", "He_Maynard_diag_genes.csv"))
save(He_Maynard_diag_genes, file = here("data", "processed", "He_Maynard_diag_genes.Rdata"))
