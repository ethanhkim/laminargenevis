## Process gene input from app ##

## From Derek Howard, Leon French's tool "Polygenic Layers" 

library(stringr)

process_gene_input <- function(input_genes) {
  processed_genes <- unlist(str_split(trimws(input_genes), "[, \n]+") )
  return(processed_genes)
}
