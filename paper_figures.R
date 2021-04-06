## Script to produce figures for paper

library(ggplot2)
library(magrittr)
library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)
library(here)
library(scales)
library(ggdendro)
library(data.table)
library(DT)

# Source functions
source(here("data_processing.R"))

# Load datasets
load(here("data", "processed", "He_DS1_logCPM_dataset.Rdata"), verbose = TRUE)
load(here("data", "processed", "Maynard_logCPM_dataset.Rdata"), verbose = TRUE)
load(here("data", "processed", "Allen_logCPM_dataset.Rdata"), verbose = TRUE)
load(here("data", "processed", "He_DS1_logCPM_filtered_dataset.Rdata"), verbose = TRUE)
load(here("data", "processed", "Maynard_logCPM_filtered_dataset.Rdata"), verbose = TRUE)
load(here("data", "processed", "Allen_logCPM_filtered_dataset.Rdata"), verbose = TRUE)
load(here("data", "processed", "He_Maynard_gene_correlation.Rdata"), verbose = TRUE)
load(here("data", "processed", "layer_marker_table.Rdata"))

# Single gene Barplot (B)

sample_gene <- "RELN"
sample_list <- c("RELN", "CUX2", "FOXP2", "RASGRF2")
sample_quantile <- top_and_bottom_quantile(
  Maynard_logCPM_dataset, 
  He_DS1_logCPM_dataset, 
  Allen_logCPM_dataset
)
sample_barplot_df <- process_barplot_data(
  sample_gene, 
  He_DS1_logCPM_dataset, 
  Maynard_logCPM_dataset, 
  Allen_logCPM_dataset)

sample_barplot_df %>%
  ggplot(aes(x = layer, y = expression, fill = source_dataset, 
        group = source_dataset)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.75) + 
  ggtitle(paste0('Expression of ', sample_gene, 
                 ' across the human neocortex')) +
  geom_hline(yintercept = sample_quantile$top_5,
             color="black", linetype="dashed") +
  geom_hline(yintercept = sample_quantile$bottom_5,
             color="black", linetype="dashed") +
  theme_bw() + 
  scale_fill_discrete(
    name="Source Dataset", 
    breaks=c("He", "Maynard", "ABI_GABAergic", "ABI_Glutamatergic", 
             "ABI_Non-neuronal"),
    labels=c("He (DLPFC)", "Maynard (DLPFC)", "AIBS: GABA (MTG)", 
             "AIBS: GLUT (MTG)", "AIBS: Non-neuron (MTG)")) +
  scale_x_discrete(name = "\nCortical Layer") +
  theme(axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        plot.title = element_text(size=21)) +
  xlab("\nCortical Layer") + ylab("Expression (log(CPM))")

# AUC heatmap

sample_AUC_bulk <- AUROC_bulk(He_DS1_logCPM_filtered_dataset, 
                              Maynard_logCPM_filtered_dataset, 
                              sample_list) 
sample_AUC_AIBS <- AUROC_AIBS(Allen_logCPM_filtered_dataset, sample_list)
sample_AUC_df <- AUROC_data(sample_AUC_bulk, sample_AUC_scRNA)

sample_AUC_df %>%
  ggplot(mapping = aes(x = dataset, y = Layer, fill = AUROC)) +
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
  labs(x = "\nSource Dataset", y = "Cortical Layer", fill = "AUC") +
  ggtitle('Layer-specific Gene Enrichment Heatmap') +
  theme(axis.text.x = element_text(angle = 25, size = 13, hjust = 0.95),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.key.size = unit(1, 'cm'),
        plot.title = element_text(size=21))

