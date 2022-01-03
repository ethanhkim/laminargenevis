# Load required libraries
library(tidyverse)
library(magrittr)
library(data.table)
library(gridExtra)

# Source data processing functions
source("src/data_processing.R")

# Load data
load(here("data", "processed", "He_DS1_logCPM_filtered_dataset.Rdata"), verbose = TRUE)
load(here("data", "processed", "Maynard_logCPM_filtered_dataset.Rdata"), verbose = TRUE)
load(here("data", "processed", "Allen_downsampled_logCPM_filtered_dataset.Rdata"), verbose = TRUE)
load(here("data", "processed", "He_Maynard_gene_correlation.Rdata"), verbose = TRUE)
load(here("data", "processed", "layer_marker_table.Rdata"))

## Create Zeng et al marker lists -----

# Read csv and format
Zeng_markers <- read_csv("data/raw/Zeng et al/Cleaned_Zeng_dataset.csv") %>%
  rename(cortical_marker_human = `Cortical.marker..human.`) %>%
  select(gene_symbol, cortical_marker_human) %>%
  filter(!is.na(cortical_marker_human))

Zeng_marker_lists <- list()
for (i in c("layer 1", "layer 2", "layer 3", "layer 4", "layer 5", "layer 6")) {
  Zeng_marker_lists[[i]] <- Zeng_markers %>%
    filter(cortical_marker_human == i) %>%
    select(gene_symbol) %>%
    pull()
} 

# Verify layer marker count (optional)
#Zeng_markers %>% 
#  group_by(cortical_marker_human) %>% 
#  summarise(n())

## Create AUC tables using Zeng et al marker lists ------

# Create AUC dataframe using bulk tissue and scRNA-seq data
create_AUC_df <- function(Zeng_marker_list, pivot = FALSE) {
  # Create AUC data for bulk tissue data
  AUC_bulk_df <- AUROC_bulk(He_DS1_logCPM_filtered_dataset,
                            Maynard_logCPM_filtered_dataset,
                            Zeng_marker_list)
  # Create AUC data for scRNA data
  AUC_scRNA_df <- AUROC_AIBS(Allen_downsampled_logCPM_filtered_dataset,
                             Zeng_marker_list)
  # Combine the two
  AUC_df <- AUROC_data(AUC_bulk_df, AUC_scRNA_df)
  
  # If no need to pivot (for creating combined dataframe)
  if (pivot == FALSE) {
    # Return overall AUC dataframe
    return(AUC_df)
  } else {
    AUC_df %<>%
      select(-signif_marker) %>%
      pivot_wider(names_from = dataset,
                  values_from = AUROC:pValue)
    return(AUC_df)
  }
}

## Create lists of dataframes:
# To use for AUC heatmap
AUC_df_list <- list()
# To use for creating combined table
AUC_combined_df <- list()
for (i in c("layer 1", "layer 2", "layer 3", "layer 4", "layer 5", "layer 6")) {
  Zeng_markers <- Zeng_marker_lists[[i]]
  AUC_df_list[[i]] <- create_AUC_df(Zeng_markers) %>% 
    arrange(desc(Layer), dataset)
  AUC_combined_df[[i]] <- create_AUC_df(Zeng_markers, pivot = TRUE) %>%
    add_column(Zeng_layer = i)
}
# Bind list to one dataframe
AUC_combined_df <- rbindlist(AUC_combined_df)
write_csv(AUC_combined_df, file = "data/processed/Zeng_list_AUC_values.csv")

## Create AUC heatmaps with Zeng et al marker lists -----
AUC_heatmap_list <- list()
for (i in c("layer 1", "layer 2", "layer 3", "layer 4", "layer 5", "layer 6")) {
  AUC_heatmap_list[[i]] <- ggplot(data = AUC_df_list[[i]], 
                 mapping = aes(x = dataset, y = Layer, fill = AUROC)) +
    geom_tile() +
    scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(0,1)) +
    scale_y_discrete(expand=c(0,0)) + 
    scale_x_discrete(expand=c(0,0),
                     breaks=c("He", "Maynard", "scRNA_GABA", 
                              "scRNA_GLUT", "scRNA_Non-neuronal"),
                     labels=c("He (DLPFC)", "Maynard (DLPFC)", 
                              "AIBS: GABA (MTG)", 
                              "AIBS: GLUT (MTG)",
                              "AIBS: Non-neuron (MTG)")) +
    geom_text(aes(label = signif_marker), size = 7, vjust = 1) +
    labs(x = "\nSource Dataset", y = "Cortical Layer", fill = "AUC",
         title = paste0("Zeng et. al ", str_to_title(i), " Markers")) +
    theme(axis.text.x = element_text(angle = 25, size = 13, hjust = 0.95),
          axis.text.y = element_text(size = 13),
          axis.title.x = element_text(size = 17),
          axis.title.y = element_text(size = 17),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(size = 19))
}

## Print AUC heatmaps -----
ggsave(
  filename = "figures/Zeng_AUC_heatmaps.pdf",
  plot = marrangeGrob(AUC_heatmap_list, nrow = 1, ncol = 1),
  width = 15, height = 9
)
