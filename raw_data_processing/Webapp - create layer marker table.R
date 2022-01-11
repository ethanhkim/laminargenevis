# Script to create laminar marker lookup table

library(tidyverse)
library(here)
library(magrittr)
library(spatialLIBD)
library(limma)
library(purrr)
library(conflicted)
library(readxl)

# Set conflicts
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")

## He Layer markers ##
He_LM_Path <- here("Data", "raw", "He et al", "Supplementary Table 2.xlsx")
He_layer_markers <- read_xlsx(path = He_LM_Path) %>%
  select("Gene symbol", "Layer marker in human") %>%
  mutate_at(vars("Layer marker in human"), na_if, "No") %>%
  mutate_at(vars("Layer marker in human"), na_if, "NA") %>%
  rename(gene_symbol = "Gene symbol",  layer_marker = "Layer marker in human") %>%
  mutate(layer_marker = str_replace(layer_marker, "L", "")) %>%
  distinct(gene_symbol, .keep_all = T)

## Maynard layer markers ##
modeling_results <- fetch_data(type = "modeling_results")
Maynard_layer_markers <- as_tibble(modeling_results$enrichment) %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::rename(tstat_WM = "t_stat_WM", tstat_Layer1 = "t_stat_Layer1", tstat_Layer2 = "t_stat_Layer2", tstat_Layer3 = "t_stat_Layer3",
                tstat_Layer4 = "t_stat_Layer4", tstat_Layer5 = "t_stat_Layer5", tstat_Layer6 = "t_stat_Layer6", gene_symbol = "gene") %>%
  dplyr::select(-"ensembl", -starts_with("p")) %>%
  filter(gene_symbol %in% modeling_results$enrichment$gene) %>%
  pivot_longer(
    cols = -gene_symbol,
    names_to = c(".value", "initial_layer_marker"),
    names_sep = "_"
  ) %>%
  group_by(gene_symbol) %>% 
  slice(which.max(tstat)) %>%
  mutate(fdr_true_false = ifelse(test = (fdr < 0.1),
                                 yes = "yes",
                                 no = "no")) %>%
  mutate(layer_marker = ifelse(test = (fdr_true_false == "yes"),
                               yes = initial_layer_marker,
                               no = NA)) %>%
  select(gene_symbol, layer_marker) %>%
  mutate(layer_marker = str_replace(layer_marker, "Layer", "")) %>%
  mutate(layer_marker = str_replace(layer_marker, "WM", "7")) %>%
  distinct(gene_symbol, .keep_all = T)

## Zeng (Allen) layer markers ##

Zeng_Path <- here("Data", "raw", "Zeng et al", "Table S2.xlsx")
Zeng_dataset <- read_xlsx(path = Zeng_Path, sheet = "Final1000New", skip=1) %>%
  select("Gene symbol", "Entrez Gene ID", "Cortical marker (human)", "Level...20", "Pattern...21",
         "Pattern...23", "Pattern...25") %>% 
  rename(gene_symbol = "Gene symbol", entrez_id = "Entrez Gene ID", marker_annotation = "Cortical marker (human)", 
         expression_level = "Level...20", V1_pattern = "Pattern...21", V2_pattern = "Pattern...23", Temporal_pattern = "Pattern...25")

Zeng_marker_annotation <- Zeng_dataset %>%
  dplyr::select("gene_symbol", "marker_annotation") %>%
  rename(layer_marker_Zeng = "marker_annotation") %>%
  mutate(layer_marker_Zeng = gsub("layer", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("interneuron", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("astrocyte", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("astrocyte?", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("laminar", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub(" ", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("VEC", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("or", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("others", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("\\+", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("oligodendrocyte", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("5a", "5", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("6b", "6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("4c", "4", layer_marker_Zeng)) %>%
  na_if("") %>%
  mutate(layer_marker_Zeng = gsub("\\?", NA, layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("\\/1", "1", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("1\\/2\\/\\/", "1\\/2", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("\\/2\\/3", "2\\/3", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("34\\/5\\/6", "3\\/4\\/5\\/6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("\\/2\\/3", "2\\/3", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("6\\/\\/", "6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("6\\/6", "6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("2\\/3\\/5\\/6", "2,3,5,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("3\\/4\\/5\\/6", "3,4,5,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("3\\/5\\/6", "3,5,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("2\\/3\\/6", "2,3,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("2\\/3\\/4", "2,3,4", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("2\\/3\\/5", "2,3,5", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("4\\/5\\/6", "3,5,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("1\\/2\\/6", "1,2,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("\\/5\\/6", "5,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("\\/4\\/5", "4,5", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("5\\/6", "5,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("3\\/5", "3,5", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("2\\/3", "2,3", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("1\\/2", "1,2", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("4\\/6", "4,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("3\\/4", "3,4", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("/2", "2", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("/4", "4", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("/6", "6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("\\/", NA, layer_marker_Zeng))

# Create lookup table for all layer markers
layer_marker_table <- left_join(Maynard_layer_markers, He_layer_markers, by = "gene_symbol") %>%
  rename(He = layer_marker.x) %>%
  rename(Maynard = layer_marker.y)
layer_marker_table %<>% left_join(Zeng_marker_annotation, by = "gene_symbol") %>%
  rename(Zeng = layer_marker_Zeng)

# Export lookup table
save(layer_marker_table, file = here("Data", "processed_data", "layer_marker_table.Rdata"))
