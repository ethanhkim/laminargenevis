# Create correlation diagonal values for He vs. Maynard

# Load libraries
library(tidyverse)
library(magrittr)
library(WGCNA)
library(conflicted)

# Set conflicts
conflict_prefer("cor", "WGCNA")

# Load data
load(here("Data", "processed_data", "He_DS1_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_dataset.Rdata"))

# Create vector of genes commom between studies
common_genelist <- intersect(He_DS1_logCPM_dataset$gene_symbol, 
                             Maynard_logCPM_dataset$gene_symbol)

He_cormatrix <- He_DS1_logCPM_dataset %>%
  filter(gene_symbol %in% common_genelist) %>%
  arrange(gene_symbol) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t()

Maynard_cormatrix <- Maynard_logCPM_dataset %>%
  filter(gene_symbol %in% common_genelist) %>%
  arrange(gene_symbol) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t()

He_Maynard_cormatrix <- cor(He_cormatrix, Maynard_cormatrix, 
                            method = "pearson",
                            nThreads = 12) 

He_Maynard_gene_correlation <- diag(He_Maynard_cormatrix)

# Save all gene-to-gene correlation vector
save(He_Maynard_gene_correlation, 
     file = here("Data", "processed_data", "He_Maynard_gene_correlation.Rdata"))

# Clean workspace
rm(He_Maynard_cormatrix, Maynard_cormatrix, He_cormatrix, 
   common_genelist)


