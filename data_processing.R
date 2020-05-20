library(org.Hs.eg.db)
library(tidyverse)
library(spatialLIBD)
library(here)
library(readxl)
library(dplyr)


## He et al, DS1 ##
here("data", "He et al", "tab.average_RPKM_section_human.tsv")
He_DS1_Human <- read.table(here("Data", "He et al", "tab.average_RPKM_section_human.tsv"))

#Truncate the version numbers on ENSEMBL ID
He_DS1_Human$V1 <- gsub("\\.\\d+$", "", He_DS1_Human$V1)

#Change ENSEMBL ID to Gene symbol
He_DS1_Human$gene_symbol <- mapIds(org.Hs.eg.db, keys = He_DS1_Human$V1, keytype = "ENSEMBL", column="SYMBOL")

#Rename
He_DS1_Human <- He_DS1_Human %>%
  dplyr::rename(Ensembl_ID = "V1",
                V1 = "V2", V2 = "V3", V3 = "V4", V4 = "V5", V5 = "V6", V6 = "V7", V7 = "V8",
                V8 = "V9", V9 = "V10", V10 = "V11", V11 = "V12", V12 = "V13", V13 = "V14", 
                V14 = "V15", V15 = "V16", V16 = "V17", V17 = "V18") 

#Remove unnecessary columns, rows and reorder so gene_symbol is at the front
He_DS1_Human <- He_DS1_Human %>%
  dplyr::select("gene_symbol", everything()) %>%
  filter(!is.na(gene_symbol)) %>%
  dplyr::select(-"Ensembl_ID") %>%
  #Remove duplicates
  distinct(gene_symbol, .keep_all = TRUE)

#Average cuts into layers
He_Layer1_averaged <- He_DS1_Human %>% 
  dplyr::select(gene_symbol, V1)

He_Layer2_averaged <- He_DS1_Human %>% 
  dplyr::select(gene_symbol, V2, V3, V4) 
He_Layer2_averaged$wt.mean <- He_Layer2_averaged %>%
  rowwise() %>%
  do(data.frame(
    wt.mean = weighted.mean(
      x = c(.$V2, .$V3, .$V4),
      w = c(1, 1, 0.5)
    )
  )) %>%
  ungroup() %>%
  use_series("wt.mean") 

He_Layer3_averaged <- He_DS1_Human %>% 
  dplyr::select(gene_symbol, V4, V5, V6)
He_Layer3_averaged$wt.mean <- He_Layer3_averaged %>%
  rowwise() %>%
  do(data.frame(
    wt.mean = weighted.mean(
      x = c(.$V4, .$V5, .$V6),
      w = c(0.5, 1, 1)
    )
  )) %>%
  ungroup() %>%
  use_series("wt.mean") 

He_Layer4_averaged <- He_DS1_Human %>% 
  dplyr::select(gene_symbol, V7, V8, V9)
He_Layer4_averaged$wt.mean <- He_Layer4_averaged %>%
  rowwise() %>%
  do(data.frame(
    wt.mean = weighted.mean(
      x = c(.$V7, .$V8, .$V9),
      w = c(1, 1, 0.5)
    )
  )) %>%
  ungroup() %>%
  use_series("wt.mean") 

He_Layer5_averaged <- He_DS1_Human %>% 
  dplyr::select(gene_symbol, V9, V10, V11, V12) 
He_Layer5_averaged$wt.mean <- He_Layer5_averaged %>%
  rowwise() %>%
  do(data.frame(
    wt.mean = weighted.mean(
      x = c(.$V9, .$V10, .$V11, .$V12),
      w = c(0.5, 1, 1, 1)
    )
  )) %>%
  ungroup() %>%
  use_series("wt.mean")

He_Layer6_averaged <- He_DS1_Human %>% 
  dplyr::select(gene_symbol, V13, V14, V15, V16)
He_Layer6_averaged$mean <- He_Layer6_averaged %>%
  rowwise() %>%
  do(data.frame(
    mean = mean(
      x = c(.$V13, .$V14, .$V15, .$V16),
      w = c(1, 1, 1, 1)
    )
  )) %>%
  ungroup() %>%
  use_series("mean")

He_WM_averaged <- He_DS1_Human  %>%
  dplyr::select(gene_symbol, V17)

He_DS1_Human_averaged <- tibble(
  gene_symbol = He_DS1_Human$gene_symbol,
  Layer_1 = He_Layer1_averaged$V1,
  Layer_2 = He_Layer2_averaged$wt.mean,
  Layer_3 = He_Layer3_averaged$wt.mean,
  Layer_4 = He_Layer4_averaged$wt.mean,
  Layer_5 = He_Layer5_averaged$wt.mean,
  Layer_6 = He_Layer6_averaged$mean,
  WM = He_WM_averaged$V17
) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t() %>%
  scale(center = TRUE) %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  add_column("gene_symbol" = He_DS1_Human$gene_symbol)


## Maynard et al ##

sce_layer <- fetch_data(type = "sce_layer")
Maynard_dataset <- as_tibble(sce_layer@assays@data@listData$logcounts)
Maynard_ensembl_list <- sce_layer@rowRanges@ranges@NAMES
Maynard_dataset$Ensembl_ID <- Maynard_ensembl_list

#Change Ensembl ID to gene symbol
Maynard_dataset$gene_symbol <- mapIds(org.Hs.eg.db, keys = Maynard_dataset$Ensembl_ID, keytype = "ENSEMBL", column="SYMBOL")
Maynard_dataset <- Maynard_dataset %>%
  dplyr::select(-"Ensembl_ID") %>% 
  dplyr::select("gene_symbol", everything()) %>%
  filter(!is.na(gene_symbol)) %>%
  distinct(gene_symbol, .keep_all = TRUE)

## Separate by layers and add average column
Maynard_dataset_Layer1 <- Maynard_dataset %>%
  dplyr::select("gene_symbol", contains("Layer1")) %>%
  mutate(dataset_mean = rowMeans(dplyr::select(., `151507_Layer1`:`151676_Layer1`)))

Maynard_dataset_Layer2 <- Maynard_dataset %>%
  dplyr::select("gene_symbol", contains("Layer2")) %>%
  mutate(dataset_mean = rowMeans(dplyr::select(., `151507_Layer2`:`151676_Layer2`)))

Maynard_dataset_Layer3 <- Maynard_dataset %>%
  dplyr::select("gene_symbol", contains("Layer3")) %>%
  dplyr::select(-"151669_Layer3", -"151670_Layer3", -"151671_Layer3", -"151672_Layer3") %>%
  mutate(dataset_mean = rowMeans(dplyr::select(., `151507_Layer3`:`151676_Layer3`)))

Maynard_dataset_Layer4 <- Maynard_dataset %>%
  dplyr::select("gene_symbol", contains("Layer4")) %>%
  dplyr::select(-"151669_Layer4", -"151670_Layer4", -"151671_Layer4", -"151672_Layer4") %>%
  mutate(dataset_mean = rowMeans(dplyr::select(., `151507_Layer4`:`151676_Layer4`)))

Maynard_dataset_Layer5 <- Maynard_dataset %>%
  dplyr::select("gene_symbol", contains("Layer5")) %>%
  dplyr::select(-"151669_Layer5", -"151670_Layer5", -"151671_Layer5", -"151672_Layer5") %>%
  mutate(dataset_mean = rowMeans(dplyr::select(., `151507_Layer5`:`151676_Layer5`)))

Maynard_dataset_Layer6 <- Maynard_dataset %>%
  dplyr::select("gene_symbol", contains("Layer6")) %>%
  dplyr::select(-"151669_Layer6", -"151670_Layer6", -"151671_Layer6", -"151672_Layer6") %>%
  mutate(dataset_mean = rowMeans(dplyr::select(., `151507_Layer6`:`151676_Layer6`)))

Maynard_dataset_WM <- Maynard_dataset %>%
  dplyr::select("gene_symbol", contains("WM")) %>%
  dplyr::select(-"151669_WM", -"151670_WM", -"151671_WM", -"151672_WM") %>%
  mutate(dataset_mean = rowMeans(dplyr::select(., `151507_WM`:`151676_WM`)))


Maynard_dataset_average <- tibble(
  gene_symbol = Maynard_dataset$gene_symbol,
  Layer_1 = Maynard_dataset_Layer1$dataset_mean,
  Layer_2 = Maynard_dataset_Layer2$dataset_mean,
  Layer_3 = Maynard_dataset_Layer3$dataset_mean,
  Layer_4 = Maynard_dataset_Layer4$dataset_mean,
  Layer_5 = Maynard_dataset_Layer5$dataset_mean,
  Layer_6 = Maynard_dataset_Layer6$dataset_mean,
  WM = Maynard_dataset_WM$dataset_mean) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t() %>%
  scale(center = TRUE) %>%
  t() %>%
  as_tibble() %>%
  add_column("gene_symbol" = Maynard_dataset$gene_symbol)

## Create .RData files for loading into app

dir.create("./data/processed")
save(He_DS1_Human_averaged, file = './data/processed/He_DS1_Human_averaged.Rdata')
save(Maynard_dataset_average, file = './data/processed/Maynard_dataset_average.Rdata')



