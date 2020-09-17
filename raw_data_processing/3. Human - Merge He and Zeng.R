#Create merged data between human data

#Create gene lists containing unique values

He_gene_list <- He_DS1_Human$gene_symbol
Zeng_gene_list <- unique(Zeng_dataset_updated$gene_symbol)
Maynard_gene_list <- Maynard_dataset$gene_symbol

#Create list of common genes
common_gene_list <- dplyr::intersect(intersect(He_gene_list, Zeng_gene_list), Maynard_gene_list)

#Filter rows containing only gene symbols from above lists
He_filtered_dataset <- He_DS1_Human %>%
  filter(gene_symbol %in% common_gene_list)
Maynard_filtered_dataset <- Maynard_dataset %>%
  filter(gene_symbol %in% common_gene_list)
Zeng_filtered_dataset <- Zeng_dataset_updated %>%
  filter(gene_symbol %in% common_gene_list)
  
#Note - SERF1A, TUT1, CRHR1 contain duplicates

#Find duplicate rows
duplicated(He_filtered_dataset$gene_symbol) %>%
  which() %>%
  print()

duplicated(Maynard_filtered_dataset$gene_symbol) %>%
  which() %>%
  print()

#Remove duplicate rows
He_filtered_dataset <- He_filtered_dataset %>% slice(-913, -930, -934)
Maynard_filtered_dataset <- Maynard_filtered_dataset %>% slice(-265, -785)

#Separate by layer 

Zeng_Layer1 <- filter(Zeng_filtered_dataset, Cortical.marker..human. == "layer 1")
Zeng_Layer2 <- filter(Zeng_filtered_dataset, Cortical.marker..human. == "layer 2")
Zeng_Layer3 <- filter(Zeng_filtered_dataset, Cortical.marker..human. == "layer 3")
Zeng_Layer4 <- filter(Zeng_filtered_dataset, Cortical.marker..human. == "layer 4")
Zeng_Layer5 <- filter(Zeng_filtered_dataset, Cortical.marker..human. == "layer 5")
Zeng_Layer6 <- filter(Zeng_filtered_dataset, Cortical.marker..human. == "layer 6") %>%
  #remove duplicated rows
  slice(-22, -31)

