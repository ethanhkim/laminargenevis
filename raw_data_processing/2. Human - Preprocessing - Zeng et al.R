#Adapted from https://github.com/jritch/mri_transcriptomics/blob/master/R%20Code/RunSingleGO.AUROC.Analysis.R

#Load required libraries

library(tidyverse)
library(readxl)
library(dplyr)
library(magrittr)
library(HGNChelper)

Zeng_Path <- here("Data", "Zeng et al", "Table S2.xlsx")

Zeng_dataset <- read_xlsx(path = Zeng_Path, sheet = "Final1000New", skip=1)

Zeng_dataset %<>% 
  dplyr::select("Gene symbol", "Cortical marker (human)", "Level...20") %>% 
  dplyr::rename(Gene.symbol = "Gene symbol", Cortical.marker..human. = "Cortical marker (human)", Expression.level = "Level...20")

#Separate by layer marker
Zeng_dataset_expanded <- Zeng_dataset %>% 
  mutate(Cortical.marker..human. = strsplit(as.character(Cortical.marker..human.), "[/+]|( or )")) %>% 
  unnest(Cortical.marker..human.) %>%
  as_tibble()

#Assimilate layers
Zeng_dataset_expanded %<>% 
  mutate(Cortical.marker..human. = gsub("layer( )?","", Cortical.marker..human.)) %<>% 
  mutate(Cortical.marker..human. = gsub("[?]","", Cortical.marker..human.)) %<>% 
  mutate(Cortical.marker..human. = gsub("4c","4", Cortical.marker..human.)) %<>% 
  mutate(Cortical.marker..human. = gsub("5a","5", Cortical.marker..human.)) %<>% 
  mutate(Cortical.marker..human. = gsub("6b","6", Cortical.marker..human.)) %<>% 
  mutate(Cortical.marker..human. = gsub("([0-6])","layer \\1", Cortical.marker..human.)) %<>% 
  mutate(Cortical.marker..human. = gsub("VEC","vascular endothelial cell", Cortical.marker..human.))

Zeng_dataset_cleaned <- Zeng_dataset_expanded %>% 
  dplyr::filter(Cortical.marker..human.  != 'others' & Cortical.marker..human. != "laminar" | is.na(NA)) %>%
  dplyr::rename(gene_symbol = "Gene.symbol")

#Check if gene_symbol list is updated

check_ZengList <- checkGeneSymbols(Zeng_dataset_cleaned$gene_symbol, unmapped.as.na = TRUE, map = NULL, 
                                      species = "human")

#Update gene_symbol list
Zeng_dataset_updated <- Zeng_dataset_cleaned
Zeng_dataset_updated$gene_symbol <- check_ZengList$Suggested.Symbol

#Remove alternative names
Zeng_dataset_updated$gene_symbol<- gsub("///.*","", Zeng_dataset_updated$gene_symbol)

#Remove rows with null values for gene_symbol
Zeng_dataset_updated %<>%
  filter(!is.na(gene_symbol))

#Summarize results
Zeng_dataset_updated %>% 
  group_by(Cortical.marker..human.) %>% 
  summarise(n())

write_csv(Zeng_dataset_updated, './R Scripts/export_data/Cleaned_Zeng_dataset.csv')
save(Zeng_dataset_updated, file = './R Scripts/export_data/Zeng_dataset_updated.Rdata')
