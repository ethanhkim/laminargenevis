## Export data to Webapp ##

Maynard_dataset_average <- Maynard_dataset_average %>%
  arrange(gene_symbol) %>%
  filter(gene_symbol %in% Maynard_layer_enrichment$gene_symbol) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t() %>%
  scale() %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "gene_symbol") %>%
  add_column(
    marker_label = Maynard_layer_enrichment$layer_marker_label,
    Layer_1_marker = Maynard_layer_enrichment$layer_marker_label_1,
    Layer_2_marker = Maynard_layer_enrichment$layer_marker_label_2,
    Layer_3_marker = Maynard_layer_enrichment$layer_marker_label_3,
    Layer_4_marker = Maynard_layer_enrichment$layer_marker_label_4,
    Layer_5_marker = Maynard_layer_enrichment$layer_marker_label_5,
    Layer_6_marker = Maynard_layer_enrichment$layer_marker_label_6,
    Layer_WM_marker = Maynard_layer_enrichment$layer_marker_label_WM)

He_DS1_averaged_by_layer_export <- He_DS1_averaged_by_layer %>%
  column_to_rownames(var = "gene_symbol") %>%
  t() %>%
  scale() %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "gene_symbol")

He_DS1_Human_averaged <-  merge(x = He_DS1_averaged_by_layer_export, y = He_Layer_markers, 
                                       by = "gene_symbol", all.x = TRUE) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  rename(Layer_WM = "Layer_7") %>%
  mutate(layer_marker_label = ifelse(test = (is.na(layer_marker)),
                                     yes = NA,
                                     no = "*")) %>%
  mutate(Layer_1_marker = ifelse(test = (layer_marker == "L1"),
                                       yes = "Layer_1",
                                       no = NA)) %>%
  mutate(Layer_2_marker = ifelse(test = (layer_marker == "L2"),
                                       yes = "Layer_2",
                                       no = NA)) %>%
  mutate(Layer_3_marker = ifelse(test = (layer_marker == "L3"),
                                       yes = "Layer_3",
                                       no = NA)) %>%
  mutate(Layer_4_marker = ifelse(test = (layer_marker == "L4"),
                                       yes = "Layer_4",
                                       no = NA)) %>%
  mutate(Layer_5_marker = ifelse(test = (layer_marker == "L5"),
                                       yes = "Layer_5",
                                       no = NA)) %>%
  mutate(Layer_6_marker = ifelse(test = (layer_marker == "L6"),
                                       yes = "Layer_6",
                                       no = NA)) %>%
  mutate(Layer_WM_marker = ifelse(test = (layer_marker == "WM"),
                                   yes = "WM",
                                   no = NA)) %>%
  rename(marker_label = "layer_marker_label", WM = "Layer_WM") %>%
  select(-"layer_marker")



write_csv(He_DS1_Human_averaged, here("data", "raw", "He et al", "He_DS1_Human_averaged.csv"))
save(He_DS1_Human_averaged, file = here("data", "raw", "He et al", "He_DS1_Human_averaged.Rdata"))


write_csv(Maynard_dataset_average, here("data", "raw", "He et al", "Maynard_dataset_average.csv"))
save(Maynard_dataset_average, file = './R Scripts/export_data/Maynard_dataset_average.Rdata')
                                       
  