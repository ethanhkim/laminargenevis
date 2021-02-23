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

# Single gene Barplot (B)

sample_gene <- "RELN"
sample_list <- c("RELN", "CUX2", "FOXP2", "RASGRF2")
sample_quantile <- top_and_bottom_quantile(
  Maynard_dataset_average, He_DS1_Human_averaged, MTG_matrix_scaled
)
sample_barplot_df <- process_barplot_data(
  sample_gene, He_DS1_Human_averaged, Maynard_dataset_average, 
  MTG_matrix_scaled)

sample_barplot_df %>%
  ggplot(aes(x = Layer, y = Z_score, fill = Source_Dataset, 
                                  group = Source_Dataset)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.75) + 
  ggtitle(paste0('Expression of ', sample_gene, 
                 '\nacross the human neocortex')) +
  geom_hline(yintercept = sample_quantile$top_5,
             color="black", linetype="dashed") +
  geom_hline(yintercept = sample_quantile$bottom_5,
             color="black", linetype="dashed") +
  theme_bw() + 
  scale_fill_discrete(name="Source Dataset",
                      breaks=c("He", "Maynard", "ABI_GABAergic", 
                               "ABI_Glutamatergic", "ABI_Non-neuronal"),
                      labels=c("He (DLPFC)", "Maynard (DLPFC)", 
                               "AIBS: GABA (MTG)", 
                               "AIBS: GLUT (MTG)",
                               "AIBS: Non-neuron (MTG)")) +
  scale_x_discrete(name = "\nCortical Layer",
                   breaks = c("1", "2", "3", "4", "5", "6", "WM"),
                   labels = c("L1", "L2", "L3", "L4", "L5", "L6", "WM")) +
  theme(axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 19),
        axis.title.y = element_text(size = 19),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=21)) +
  xlab("\nCortical Layer") + ylab("mRNA expression (normalized)")

# AUC heatmap

sample_AUC_bulk <- AUROC_bulk(He_DS1_Human_averaged, Maynard_dataset_average, 
                             sample_list) 
sample_AUC_scRNA <- AUROC_scRNA(MTG_matrix_scaled, sample_list)
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
  #Puts stars on layer marker annotations
  geom_text(aes(label = signif_marker), size = 7, vjust = 1) +
  labs(x = "\nSource Dataset", y = "Cortical Layer", fill = "AUC") +
  ggtitle('Layer-specific Gene Enrichment Heatmap') +
  theme(axis.text.x = element_text(angle = 25, size = 15, hjust = 0.95),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 19),
        axis.title.y = element_text(size = 19),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.key.size = unit(1, 'cm'),
        plot.title = element_text(size=21))
