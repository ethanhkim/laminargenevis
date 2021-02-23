## Scripts for processing strings and/or data ##

## Process gene input from app ##

## From Derek Howard, Leon French's tool "Polygenic Layers" 

library(stringr)

process_gene_input <- function(input_genes) {
  processed_genes <- unlist(str_split(trimws(input_genes), "[, \n]+") )
  return(processed_genes)
}

assayed_gene_string <- function(genelist, He_df, Maynard_df, single_or_multiple) {
  if (single_or_multiple == "single") {
    statement <- cat(paste0(
      "Your submitted gene (", genelist,
      ") was found to be ",
      if (length(genelist %in% He_df) == 0) {
        "not assayed by He et al., "
      } else "assayed by He et al.,",
      if (length(genelist %in% Maynard_df) == 0) {
        "not assayed by Maynard et al., "
      } else " and assayed by Maynard et al.\n"
      ))
  } else {
    statement <- cat(paste0(
      "You inputted ", length(genelist),
      " genes. Of those genes: ",
      sum(genelist %in% unique(He_df$gene_symbol)),
      " genes were assayed by He et al., and ",
      sum(genelist %in% unique(Maynard_df$gene_symbol)),
      " genes were assayed by Maynard et al.\n"
    ))
  }
  return(statement)
}

layer_marker_preempt_string <- function(source_df, single_or_multiple, genelist) {
  if (single_or_multiple == "single") {
    statement <- cat(paste0("\nIn the ", source_df,
                            " et al data, ", genelist, 
                            " was found to "))
  } else {
    statement <- cat(paste0("\nIn the ", source_df, 
                            " et al data, these genes marked these specific layers: \n"))
  }
  
  return(statement)
}

layer_marker_string <- function(layer_marker_list, single_or_multiple) {
  
  if (single_or_multiple == "single") {
    if (length(unlist(layer_marker_list$Layer_1)) == 0 &
        length(unlist(layer_marker_list$Layer_2)) == 0 &
        length(unlist(layer_marker_list$Layer_3)) == 0 &
        length(unlist(layer_marker_list$Layer_4)) == 0 & 
        length(unlist(layer_marker_list$Layer_5)) == 0 &
        length(unlist(layer_marker_list$Layer_6)) == 0 &
        length(unlist(layer_marker_list$WM)) == 0) {
      statement <- "not mark any layer."
    } else {
      statement <- cat(paste0(
        "mark layer ",
        if (length(unlist(layer_marker_list$Layer_1)) == 0) {
          ""
        } else "1.",
        if (length(unlist(layer_marker_list$Layer_2)) == 0) {
          ""
        } else "2.",
        if (length(unlist(layer_marker_list$Layer_3)) == 0) {
          ""
        } else "3.",
        if (length(unlist(layer_marker_list$Layer_4)) == 0) {
          ""
        } else "4.",
        if (length(unlist(layer_marker_list$Layer_5)) == 0) {
          ""
        } else "5.",
        if (length(unlist(layer_marker_list$Layer_6)) == 0) {
          ""
        } else "6."
      ))
    }
  } else {
    if (length(unlist(layer_marker_list$Layer_1)) == 0 &
        length(unlist(layer_marker_list$Layer_2)) == 0 &
        length(unlist(layer_marker_list$Layer_3)) == 0 &
        length(unlist(layer_marker_list$Layer_4)) == 0 & 
        length(unlist(layer_marker_list$Layer_5)) == 0 &
        length(unlist(layer_marker_list$Layer_6)) == 0 &
        length(unlist(layer_marker_list$WM)) == 0) {
      statement <- "There were no genes that marked any layer."
    } else {
      statement <- cat(paste0( 
        if (length(unlist(layer_marker_list$Layer_1)) == 0) {
          ""
        } else paste0("- ",length(unlist(layer_marker_list$Layer_1)),
                      " gene(s) marked layer 1 (", 
                      paste(layer_marker_list$Layer_1, collapse = ", "), ")\n"),
        if (length(unlist(layer_marker_list $Layer_2)) == 0) {
          ""
        } else paste0("- ",length(unlist(layer_marker_list$Layer_2)), 
                      " gene(s) marked layer 2 (", 
                      paste(layer_marker_list $Layer_2, collapse = ", "), ")\n"),
        if (length(unlist(layer_marker_list $Layer_3)) == 0) {
          ""
        } else paste0("- ",length(unlist(layer_marker_list $Layer_3)), 
                      " gene(s) marked layer 3 (",
                      paste(layer_marker_list $Layer_3, collapse = ", "), ")\n"),
        if (length(unlist(layer_marker_list $Layer_4)) == 0) {
          ""
        } else paste0("- ",length(unlist(layer_marker_list $Layer_4)), 
                      " gene(s) marked layer 4 (",
                      paste(layer_marker_list $Layer_4, collapse = ", "), ")\n"),
        if (length(unlist(layer_marker_list $Layer_5)) == 0) {
          ""
        } else paste0("- ",length(unlist(layer_marker_list $Layer_5)), 
                      " gene(s) marked layer 5 (",
                      paste(layer_marker_list $Layer_5, collapse = ", "), ")\n"),
        if (length(unlist(layer_marker_list $Layer_6)) == 0) {
          ""
        } else paste0("- ", length(unlist(layer_marker_list $Layer_6)), 
                      " gene(s) marked layer 6 (",
                      paste(layer_marker_list $Layer_6, collapse = ", "), ").\n")
      ))
    }
  }
  return(statement)
}

stats_string <- function(genelist, correlation, p_value, quantile_stat, single_or_multiple) {
  
  if (single_or_multiple == "single") {
    statement <- cat(paste0(
      "\nBetween the He and Maynard datasets, ",
      genelist,
      "'s expression across the layers had a Pearson correlation value of ", 
      correlation,
      " (p = ",
      p_value,
      "), which ranks in the ",
      quantile_stat,
      "th quantile.\n"
    ))
  } else {
    statement <- cat(paste0(
      "\nBetween the He and Maynard datasets, the input genes' expression had a ",
      "mean Pearson correlation value of ", 
      correlation,
      " (p = ",
      p_value,
      "), which ranks in the ",
      quantile_stat,
      "th quantile.\n"
    ))
  }
  return(statement)
}


