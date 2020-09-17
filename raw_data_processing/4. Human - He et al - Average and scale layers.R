## Load required libraries ##

library(tidyverse)
library(magrittr)

## Average He et al ##

He_DS1_values_averaged <- He_DS1_Human_values

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
  select(gene_symbol, V4, V5, V6)
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
  select(gene_symbol, V7, V8, V9)
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
  select(gene_symbol, V9, V10, V11, V12) 
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
  select(gene_symbol, V13, V14, V15, V16) %>%
  rowwise() %>%
  mutate(mean = mean(V13, V14, V15, V16))

He_WM_averaged <- He_DS1_Human  %>%
  select(gene_symbol, V17)


## Create averaged He et al. data table from the scaled values
He_DS1_averaged_by_layer <- tibble(
  gene_symbol = He_DS1_Human$gene_symbol,
  Layer_1 = He_Layer1_averaged$V1,
  Layer_2 = He_Layer2_averaged$wt.mean,
  Layer_3 = He_Layer3_averaged$wt.mean,
  Layer_4 = He_Layer4_averaged$wt.mean,
  Layer_5 = He_Layer5_averaged$wt.mean,
  Layer_6 = He_Layer6_averaged$mean,
  Layer_7 = He_WM_averaged$V17
)
<<<<<<< HEAD
=======


## Scale Layers ##

He_values_scaled <- He_values %>%
  dplyr::select(-"gene_symbol") %>%
  t() %>%
  scale() %>%
  t() %>%
  as_tibble() %>%
  add_column(He_values$gene_symbol) %>%
  rename(gene_symbol= "He_values$gene_symbol") %>%
  dplyr::select("gene_symbol", everything())


## Transpose data ##

#Transpose table such that table represents the 17 separate cuts for each gene
rownames(He_values) <- common_gene_list 
He_values_transposed <- He_values %>%
  dplyr::select(-"gene_symbol") %>%
  t() %>%
  as_tibble()



>>>>>>> 4e70cfd78b048eb455440c3023a2e50d19d2e770
