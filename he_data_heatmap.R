#Load required libraries

library(reshape2)

#Create lists of layer-specific genes

Zeng_layer1_gene_list <- Zeng_Layer1$gene_symbol
Zeng_layer2_gene_list <- Zeng_Layer2$gene_symbol
Zeng_layer3_gene_list <- Zeng_Layer3$gene_symbol
Zeng_layer4_gene_list <- Zeng_Layer4$gene_symbol
Zeng_layer5_gene_list <- Zeng_Layer5$gene_symbol
Zeng_layer6_gene_list <- Zeng_Layer6$gene_symbol

#Subset scaled values of He dataset with above lists

create_subset_He_layers <- function(x) {
  data <- He_values_scaled %>%
    filter(gene_symbol %in% x) %>%
    slice(match(x, gene_symbol)) %>%
    column_to_rownames(var = "gene_symbol") %>%
    dplyr::select(-"V17") %>%
    as.data.frame() 
}

He_Layer1_subset_Zeng <- create_subset_He_layers(Zeng_layer1_gene_list)
He_Layer2_subset_Zeng <- create_subset_He_layers(Zeng_layer2_gene_list) 
He_Layer3_subset_Zeng <- create_subset_He_layers(Zeng_layer3_gene_list)
He_Layer4_subset_Zeng <- create_subset_He_layers(Zeng_layer4_gene_list)
He_Layer5_subset_Zeng <- create_subset_He_layers(Zeng_layer5_gene_list)
He_Layer6_subset_Zeng <- create_subset_He_layers(Zeng_layer6_gene_list)


#Heatmaps

transform_tibble <- function(x) {
  t(x) %>%
    melt() %>%
    as_tibble() %>%
    rename(Cuts = "Var1", Gene_symbol = "Var2", Z_score = "value")
}

#Heatmap of layer 1 specific genes
He_Layer1_heatmap_data <- transform_tibble(He_Layer1_subset_Zeng)
ggplot(data = He_Layer1_heatmap_data, mapping = aes(x = Cuts, y = Gene_symbol, fill = Z_score)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu")

#Heatmap of layer 2 specific genes
He_Layer2_heatmap_data <- transform_tibble(He_Layer2_subset_Zeng)
ggplot(data = He_Layer2_heatmap_data, mapping = aes(x = Cuts, y = Gene_symbol, fill = Z_score)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu")

#Heatmap of layer 3 specific genes
He_Layer3_heatmap_data <- transform_tibble(He_Layer3_subset_Zeng)
ggplot(data = He_Layer3_heatmap_data, mapping = aes(x = Cuts, y = Gene_symbol, fill = Z_score)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu")

#Heatmaps of layer 4 specific genes
He_Layer4_heatmap_data <- transform_tibble(He_Layer4_subset_Zeng)
ggplot(data = He_Layer4_heatmap_data, mapping = aes(x = Cuts, y = Gene_symbol, fill = Z_score)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu")

#Heatmap of layer 5 specific genes
He_Layer5_heatmap_data <- transform_tibble(He_Layer5_subset_Zeng)
ggplot(data = He_Layer5_heatmap_data, mapping = aes(x = Cuts, y = Gene_symbol, fill = Z_score)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu")

#Heatmap of layer 6 specific genes
He_Layer6_heatmap_data <- transform_tibble(He_Layer6_subset_Zeng)
ggplot(data = He_Layer6_heatmap_data, mapping = aes(x = Cuts, y = Gene_symbol, fill = Z_score)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu")

