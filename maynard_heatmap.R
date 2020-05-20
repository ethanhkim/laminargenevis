
Maynard_dataset_Layer1 %<>%
  mutate(average = rowMeans(dplyr::select(., `151507_Layer1`:`151676_Layer1`))) %>%
  dplyr::select("gene_symbol", "average", everything())

Maynard_dataset_Layer2 %<>%
  mutate(average = rowMeans(dplyr::select(., `151507_Layer2`:`151676_Layer2`))) %>%
  dplyr::select("gene_symbol", "average", everything())

Maynard_dataset_Layer3 %<>%
  mutate(average = rowMeans(dplyr::select(., `151507_Layer3`:`151676_Layer3`))) %>%
  dplyr::select("gene_symbol", "average", everything())

Maynard_dataset_Layer4 %<>%
  mutate(average = rowMeans(dplyr::select(., `151507_Layer4`:`151676_Layer4`))) %>%
  dplyr::select("gene_symbol", "average", everything())

Maynard_dataset_Layer5 %<>%
  mutate(average = rowMeans(dplyr::select(., `151507_Layer5`:`151676_Layer5`))) %>%
  dplyr::select("gene_symbol", "average", everything())

Maynard_dataset_Layer6 %<>%
  mutate(average = rowMeans(dplyr::select(., `151507_Layer6`:`151676_Layer6`))) %>%
  dplyr::select("gene_symbol", "average", everything())

Maynard_dataset_WM %<>%
  mutate(average = rowMeans(dplyr::select(., `151507_WM`:`151676_WM`))) %>%
  dplyr::select("gene_symbol", "average", everything())

Maynard_dataset_average <- Maynard_dataset$gene_symbol %>%
  as_tibble() %>%
  rename(gene_symbol = "value")
Maynard_dataset_average$Layer1 <- Maynard_dataset_Layer1$average
Maynard_dataset_average$Layer2 <- Maynard_dataset_Layer2$average
Maynard_dataset_average$Layer3 <- Maynard_dataset_Layer3$average
Maynard_dataset_average$Layer4 <- Maynard_dataset_Layer4$average
Maynard_dataset_average$Layer5 <- Maynard_dataset_Layer5$average
Maynard_dataset_average$Layer6 <- Maynard_dataset_Layer6$average

Maynard_dataset_average %<>%
  dplyr::select(-"gene_symbol") %>%
  t() %>%
  scale() %>%
  t() %>%
  as_tibble() %>%
  add_column(Maynard_dataset$gene_symbol) %>%
  rename(gene_symbol = "Maynard_dataset$gene_symbol")

create_subset_Maynard_layers <- function(x) {
  data <- Maynard_dataset_average %>%
    filter(gene_symbol %in% x) %>%
    slice(match(x, gene_symbol)) %>%
    column_to_rownames(var = "gene_symbol") %>%
    as.data.frame()
}

Maynard_Layer1_subset_Zeng <- create_subset_Maynard_layers(Zeng_layer1_gene_list)
Maynard_Layer2_subset_Zeng <- create_subset_Maynard_layers(Zeng_layer2_gene_list)
Maynard_Layer3_subset_Zeng <- create_subset_Maynard_layers(Zeng_layer3_gene_list)
Maynard_Layer4_subset_Zeng <- create_subset_Maynard_layers(Zeng_layer4_gene_list)
Maynard_Layer5_subset_Zeng <- create_subset_Maynard_layers(Zeng_layer5_gene_list)
Maynard_Layer6_subset_Zeng <- create_subset_Maynard_layers(Zeng_layer6_gene_list)


transform_tibble <- function(x) {
  t(x) %>%
    melt() %>%
    as_tibble() %>%
    rename(Layer = "Var1", Gene_symbol = "Var2", Z_score = "value")
}

Maynard_Layer1_heatmap_data <- transform_tibble(Maynard_Layer1_subset_Zeng)
ggplot(data = Maynard_Layer1_heatmap_data, mapping = aes(x = Layer, y = Gene_symbol, fill = Z_score)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu")

Maynard_Layer2_heatmap_data <- transform_tibble(Maynard_Layer2_subset_Zeng)
ggplot(data = Maynard_Layer2_heatmap_data, mapping = aes(x = Layer, y = Gene_symbol, fill = Z_score)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu")

Maynard_Layer3_heatmap_data <- transform_tibble(Maynard_Layer3_subset_Zeng)
ggplot(data = Maynard_Layer3_heatmap_data, mapping = aes(x = Layer, y = Gene_symbol, fill = Z_score)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu")

Maynard_Layer4_heatmap_data <- transform_tibble(Maynard_Layer4_subset_Zeng)
ggplot(data = Maynard_Layer4_heatmap_data, mapping = aes(x = Layer, y = Gene_symbol, fill = Z_score)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu")

Maynard_Layer5_heatmap_data <- transform_tibble(Maynard_Layer5_subset_Zeng)
ggplot(data = Maynard_Layer5_heatmap_data, mapping = aes(x = Layer, y = Gene_symbol, fill = Z_score)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu")

Maynard_Layer6_heatmap_data <- transform_tibble(Maynard_Layer6_subset_Zeng)
ggplot(data = Maynard_Layer6_heatmap_data, mapping = aes(x = Layer, y = Gene_symbol, fill = Z_score)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu")

