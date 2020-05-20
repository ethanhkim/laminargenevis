## Process gene input from app ##

# adapted from Derek Howard, Leon French

library(stringr)
library(homologene)

process_gene_input <- function(input_genes) {
  processed_genes <- unlist(str_split(trimws(input_genes), "[, \n]+") )
  return(processed_genes)
}

convert_to_human <- function(input_genes, in_species) {
  if (in_species == "Rhesus Macaque") {
    human_genes <- homologene(input_genes, inTax = 9544, outTax = 9606)
    return(unique(human_genes$humanGene))
  } else if (in_species == "Chimp") {
    human_genes <- homologene(input_genes, inTax = 9598, outTax = 9606)
  } else {
    return(unique(input_genes))
  }
  #return(unique(human_genes$humanGene))
}


#genes <- c('CADM2', 'ZNF704', 'NCAM1', 'RABEP2', 'ATP2A1')

#genes <- c('CYP3A5','CYP3A4','CYP3A7','CYP39A1','CYP3A43')

#convert_genes(genes)