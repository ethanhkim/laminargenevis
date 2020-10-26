## App script ##

library(shiny)
library(ggplot2)
library(magrittr)
library(plotly)
library(tibble)
library(shinyjs)
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)
library(here)
library(scales)
library(ggdendro)
library(data.table)
library(shinycssloaders)
library(DT)
source("string_processing.R") 
source("data_processing.R")
load(here("data", "processed", "He_DS1_Human_averaged.Rdata"), verbose = TRUE)
load(here("data", "processed", "Maynard_dataset_average.Rdata"), verbose = TRUE)
load(here("data", "processed", "Zeng_dataset_long.Rdata"), verbose = TRUE)
load(here("data", "processed", "Compared_Layer_markers.Rdata"), verbose = TRUE)
load(here("data", "processed", "He_Maynard_diag_genes.Rdata"), verbose = TRUE)
load(here("data", "processed", "Allen_scRNA", "MTG_matrix_scaled.Rdata"))


# Define UI ----
ui <- fluidPage(
  shinyjs::useShinyjs(),
  tags$head(includeHTML("google-analytics.html")),
  navbarPage(title = "Gene Expression Profile Comparison", 
             
  # Visualize gene expression across layers through heatmap or barplot    
  tabPanel(title = "Gene Visualization",
    sidebarLayout(
      sidebarPanel(
        # Single or multiple gene selector
        radioButtons(
          inputId = "selector", label =  "Single or multiple genes?",
          choices = c("Single", "Multiple")
        ),
        # Single input side-panel  ----
        conditionalPanel(
          condition = "input.selector == 'Single'",
          selectizeInput(
            inputId = "genelist", label = "Input gene:", choices = NULL,
            selected = NULL, multiple = FALSE, options = NULL),
            actionButton(inputId = "submit_barplot", label = "Submit")
        ),
        # Multiple input side-panel ----     
        conditionalPanel(
          condition = "input.selector == 'Multiple'",
          textAreaInput(
            inputId = "multiple_genelist", 
            label = "Input your gene list:", 
            placeholder = "GAD1, CCK, GRIN1"),
            actionButton(inputId = "submit_heatmap", label = "Submit")
        ),
      ),
      # Main page                 
      mainPanel(
        tabsetPanel(type = "tabs", id = "tabset",
          #Page for displaying information about datasets ----
          tabPanel(title = "Dataset Overview", value = "overview",
            br(),
            h3("Welcome to the Gene Visualization project!"),
            br(),
            h4("This web application allows you to examine layer-specific gene expression across the cortex and determine layer annotations."),
            h4("The data for this application has been sourced from these following studies and institutions:"),
            br(),
            br(),
            h3("Bulk Tissue:"),
            h3(a("Zeng et al. (2012)", href = "https://pubmed.ncbi.nlm.nih.gov/22500809/", target = "_blank")),
            h4("This study assayed roughly 1000 genes through in-situ hybridization in samples from the midtemporal and primary visual cortices."),
            br(),
            h3(a("He et al. (2017)", href = "https://pubmed.ncbi.nlm.nih.gov/28414332/", target = "_blank")),
            h4("This study assayed the whole genome using high-throughput RNA-seq in samples from the DLPFC."),
            br(),
            h3(a("Maynard et al. (2020)*", href = "https://www.biorxiv.org/content/10.1101/2020.02.28.969931v1", target = "_blank")),
            h4("*This study is currently a pre-print; it assayed the whole genome through the 10X Genomics Visium Platform in samples from the DLPFC."),
            br(),
            h3("Single-cell:"),
            h3(a("Allen Brain Atlas", href = "https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq", target = "_blank")),
            h4("This dataset assayed multiple cortical regions and assayed the whole genome across roughly 49,000 single-cell nuclei.")
          ),
          #Page for displaying single or multiple gene visualizations ----
          tabPanel(title = "Gene Visualization", value = "visualization",
            # Single gene visualization
            br(),
            conditionalPanel(
              h3("Layer-specific gene expression"),
              #Only show when the input selector is Single 
              condition = "input.selector == 'Single'",
              #Output barplot visualization
              plotOutput("Barplot") %>% withSpinner(),
              br(),
              br(),
              h5(textOutput("barplot_caption")),
              br(),
              br(),
              h4(verbatimTextOutput("summary_single")),
            ),
            # Multiple gene visualization
            conditionalPanel(
              h4("Layer-specific Heatmaps"),
              #Only show when the input selector is Multiple
              condition = "input.selector == 'Multiple'",
              br(),
              #Output heatmap visualizations for bulk tissue RNA-seq
              plotOutput("He_figure", height = "auto") %>% withSpinner(),
              plotOutput("Maynard_figure", height = "auto") %>% withSpinner(),
              br(),
              h5(textOutput("bulk_figure_caption")),
              br(),
              h4(verbatimTextOutput("summary_multiple")),
              br(),
              br(),
              #Output heatmap visualizations for scRNA-seq data per cell class type
              plotOutput("scRNA_heatmap_GABA", height = "auto") %>% withSpinner(),
              plotOutput("scRNA_heatmap_GLUT", height = "auto") %>% withSpinner(),
              plotOutput("scRNA_heatmap_NON", height = "auto") %>% withSpinner(),
              h5(textOutput("scRNA_figure_caption")),
              br(),
              br(),
              plotlyOutput("AUROC_heatmap") %>% withSpinner(),
              h5(textOutput("AUROC_heatmap_caption")),
              br(),
            )
           )
          )
        )
      )
    )
  )
)


# Define server logic ----
server <- function(input, output, session) {
  
  # List of genes that are common through the He and Maynard datasets
  common_genelist <- intersect(He_DS1_Human_averaged$gene_symbol, Maynard_dataset_average$gene_symbol) %>%
    sort()
  
  # Table of layer-marker annotations provided by Zeng et al.
  layer_marker_table <- Compared_Layer_Markers %>%
    filter(!is.na(layer_marker)) %>%
    filter(str_detect(source_dataset, 'Zeng'))
  
  # The diagonal of the p-value correlation matrix of He and Maynard gene expression values 
  He_Maynard_cor_diagonal <- He_Maynard_diag_genes %>%
    pull(var = 1)
  
  updateSelectizeInput(session, inputId = "genelist", 
                       choices = common_genelist, server = TRUE)
  
  #Hide plots until when genes are submitted
  output$Barplot <- NULL
  output$He_figure <- NULL
  output$Maynard_figure <- NULL 
  output$heatmap_caption <- NULL
  output$AUROC_heatmap <- NULL
  output$scRNA_heatmap_GABA <- NULL
  output$scRNA_heatmap_GLUT <- NULL
  output$scRNA_heatmap_NON <- NULL
  
  # Single gene input ----
  observeEvent(input$submit_barplot, {
    #When the submit button is pressed, change to the Visualization page
    updateTabsetPanel(session, "tabset", selected = "visualization")
    
    # List of selected gene(s)
    selected_gene_list_single <- isolate(process_gene_input(input$genelist))
    # Process dataset to correct format for heatmap and barplot
    Barplot_data <- process_barplot_data(selected_gene_list_single, He_DS1_Human_averaged, Maynard_dataset_average, MTG_matrix_scaled)
    
    # Single gene visualization; barplot
    output$Barplot <- renderPlot({
      ggplot(data = Barplot_data, aes(x = Layer, y = Z_score, fill = Source_Dataset, group = Source_Dataset)) +
        geom_bar(stat = "identity", position = "dodge", width = 0.75) + theme(axis.text.x = element_text(angle = 45)) +
        geom_text(aes(label = layer_label, group = Source_Dataset, 
                      vjust = ifelse(Z_score >= -0.1, 0, 2.5)), 
                  position = position_dodge(width = 0.75)) +
        geom_hline(yintercept = 1.681020) +
        geom_hline(yintercept = -1.506113) +
        theme_bw() +
        labs(fill = "Source Dataset")
    }) 
    
    output$barplot_caption <- renderPrint({
      cat(paste("Fig 1. The barplots were created using the data from He et al and Maynard et al studies, as well as the 
                data from the Allen Brain Institute. Raw RNA-seq data was normalized using z-score normalization. 
                The horizontal lines represent the value of the top (95th) and bottom (5th) quantile of expression 
                levels across both datasets."))
    })
    
    #Filter for selected genes from table containing Zeng et al layer marker annotations
    layer_marker_table_single <- layer_marker_table %>%
      dplyr::filter(gene_symbol %in% selected_gene_list_single)
    
    layer_specific_gene_list_single <- separate_layers(layer_marker_table_single, selected_gene_list_single)
    names(layer_specific_gene_list_single) <- c("Layer 1", "Layer 2", "Layer 3", "Layer 4",
                                                "Layer 5", "Layer 6", "White_matter")
    
    #Generate correlation value for single gene, the quantile and associated p-value
    single_gene_cor <- single_gene_correlation(selected_gene_list_single, He_DS1_Human_averaged, Maynard_dataset_average)
    single_gene_quantile <- quantile_distribution(He_Maynard_cor_diagonal, single_gene_cor)
    p_value_single_gene <- wilcoxtest(selected_gene_list_single, He_DS1_Human_averaged, Maynard_dataset_average, He_Maynard_cor_diagonal)
    
    #Single gene summary table ----
    output$summary_single <- renderPrint({
      cat(paste0(
        "Your submitted gene (",
        selected_gene_list_single,
        ") was found to be ",
        if (sum(length(selected_gene_list_single %in% He_DS1_Human_averaged$gene_symbol)) == 0) {
          "not assayed by He et al., "
        } else "assayed by He et al., ",
        if (sum(length(selected_gene_list_single %in% Maynard_dataset_average$gene_symbol)) == 0) {
          "not assayed by Maynard et al., "
        } else "assayed by Maynard et al., ",
        if (sum(length(selected_gene_list_single %in% Zeng_dataset_long$gene_symbol)) == 0) {
          "and not assayed by Zeng et al."
        } else "assayed by Zeng et al.\n\n",
        "Between the He and Maynard datasets, ",
        selected_gene_list_single,
        " has a Pearson correlation value of ",
        single_gene_cor,
        ", ranking in the ",
        single_gene_quantile,
        "th quantile (p = ",
        p_value_single_gene,
        ").\n\n",
        if (sum(length(selected_gene_list_single %in% Zeng_dataset_long$gene_symbol)) == 0) {
          ""
        } else paste0(
          "Specifically in the Zeng dataset, ",
          selected_gene_list_single,
          " was found to ",
          if (length(unlist(layer_specific_gene_list_single $`Layer 1`)) == 0 &
              length(unlist(layer_specific_gene_list_single $`Layer 2`)) == 0 &
              length(unlist(layer_specific_gene_list_single $`Layer 3`)) == 0 &
              length(unlist(layer_specific_gene_list_single $`Layer 4`)) == 0 & 
              length(unlist(layer_specific_gene_list_single $`Layer 5`)) == 0 &
              length(unlist(layer_specific_gene_list_single $`Layer 6`)) == 0 &
              length(unlist(layer_specific_gene_list_single $White_matter)) == 0) {
            "not mark any layer."
          } else {
            paste0(
              "mark ",
              if (length(unlist(layer_specific_gene_list_single $`Layer 1`)) == 0) {
                ""
              } else "layer 1, ",
              if (length(unlist(layer_specific_gene_list_single $`Layer 2`)) == 0) {
                ""
              } else "layer 2, ",
              if (length(unlist(layer_specific_gene_list_single $`Layer 3`)) == 0) {
                ""
              } else "layer 3, ",
              if (length(unlist(layer_specific_gene_list_single $`Layer 4`)) == 0) {
                ""
              } else "layer 4, ",
              if (length(unlist(layer_specific_gene_list_single $`Layer 5`)) == 0) {
                ""
              } else "layer 5, ",
              if (length(unlist(layer_specific_gene_list_single $`Layer 6`)) == 0) {
                ""
              } else "layer 6."
            )
          }
        ),
        sep="\n"
      ))
    }) 
  })

  
  #Multiple gene input ----
  observeEvent(input$submit_heatmap, {
    updateTabsetPanel(session, "tabset", selected = "visualization")
    
    # Create a vector of selected genes
    selected_gene_list_multiple <- isolate(process_gene_input(input$multiple_genelist))
    # Process He & Maynard data according to input genes
    He_heatmap_data <- process_heatmap_function(He_DS1_Human_averaged, selected_gene_list_multiple)
    Maynard_heatmap_data <- process_heatmap_function(Maynard_dataset_average, selected_gene_list_multiple)
    # Process scRNA heatmap data and separate by class
    scRNA_GABA <- scRNA_Heatmap_data(MTG_matrix_scaled, selected_gene_list_multiple, "GABAergic")
    scRNA_GLUT <- scRNA_Heatmap_data(MTG_matrix_scaled, selected_gene_list_multiple, "Glutamatergic") 
    scRNA_NON <- scRNA_Heatmap_data(MTG_matrix_scaled, selected_gene_list_multiple, "Non-neuronal") 
    # Filter the layer marker table for the genes inputted
    layer_marker_table_multiple <- layer_marker_table %>%
      dplyr::filter(gene_symbol %in% selected_gene_list_multiple)
    layer_specific_gene_list_multiple <- separate_layers(layer_marker_table_multiple, selected_gene_list_multiple)
    names(layer_specific_gene_list_multiple) <- c("Layer 1", "Layer 2", "Layer 3", "Layer 4",
                                                  "Layer 5", "Layer 6")
    
    # Generate correlation value for multiple gene and the quantile that value belongs in
    multi_gene_cor <- multi_gene_correlation(selected_gene_list_multiple, He_DS1_Human_averaged, Maynard_dataset_average)
    multi_gene_quantile <- quantile_distribution(He_Maynard_cor_diagonal, multi_gene_cor)
    p_value_multiple_gene <- wilcoxtest(selected_gene_list_multiple, He_DS1_Human_averaged, Maynard_dataset_average, He_Maynard_cor_diagonal)
    
    ## Code for AUROC analysis adapted from Derek Howard & Leon French - refer to data_processing.R
    AUROC_bulk_data <- AUROC_bulk(He_DS1_Human_averaged, Maynard_dataset_average, selected_gene_list_multiple) 
    AUROC_scRNA_data <- AUROC_scRNA(MTG_matrix_scaled, selected_gene_list_multiple)
    AUROC_df <- AUROC_data(AUROC_bulk_data, AUROC_scRNA_data)
    
    #Bulk tissue RNA-seq heatmaps ---- 
    heatmapHeight <- heatmap_height(selected_gene_list_multiple)
    if (length(selected_gene_list_multiple) <= 30) {
      output$He_figure <- renderPlot({
          ggplot(data = He_heatmap_data, mapping = aes(x = layer, y = gene_symbol, fill = Z_score)) +
            geom_tile() +
            scale_fill_distiller(palette = "RdYlBu", limits = c(-3,3)) +
            scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
            labs(y = "", x = "", title = "He et al Heatmap") +
            #Puts stars on layer marker annotations
            geom_text(aes(label = layer_label), size = 7, vjust = 1)
      }, height = heatmapHeight)
      output$Maynard_figure <- renderPlot({
        ggplot(data = Maynard_heatmap_data, mapping = aes(x = layer, y = gene_symbol, fill = Z_score)) +
          geom_tile() +
          scale_fill_distiller(palette = "RdYlBu", limits = c(-3,3)) +
          scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
          labs(y = "", x = "", title = "Maynard et al Heatmap") +
          #Puts stars on layer marker annotations
          geom_text(aes(label = layer_label), size = 7, vjust = 1)
      }, height = heatmapHeight)
      
      output$bulk_figure_caption <- renderPrint({
        cat(paste("Fig 1. The heatmaps were created using the data from He et al and Maynard et al studies. The raw 
                RNA-seq data was normalized using z-score normalization. The stars (*) indicate if the
                source paper denoted the gene to be highly enriched within the specific cortical layer."))
      })
    } else {
      output$He_figure <- renderPlot({
        He_heatmap_data %<>% inner_join(He_heatmap_data %>% group_by(layer) %>% summarize(median_rank = median(Z_score)), by = "layer") %>%
          mutate(layer = factor(layer, levels = c("Layer_1", "Layer_2", "Layer_3", "Layer_4", "Layer_5", "Layer_6"))) %>%
          ggplot(aes(x = layer, y = Z_score, group = layer, names = gene_symbol, fill = layer)) +
          geom_errorbar(aes(ymax = median_rank, ymin = median_rank), colour = "black", linetype = 1) +
          geom_jitter(width = .05, alpha = 0.4) +
          ylim(-3, 3) +
          guides(fill = "none") +
          theme_bw() +
          labs(x = "Z score", y = "Layer", title = "He et al Scatter plot")
      }, height = heatmapHeight)
      output$Maynard_figure <- renderPlot({
        Maynard_heatmap_data %<>% inner_join(Maynard_heatmap_data %>% group_by(layer) %>% summarize(median_rank = median(Z_score)), by = "layer") %>%
          mutate(layer = factor(layer, levels = c("Layer_1", "Layer_2", "Layer_3", "Layer_4", "Layer_5", "Layer_6"))) %>%
          ggplot(aes(x = layer, y = Z_score, group = layer, names = gene_symbol, fill = layer)) +
          geom_errorbar(aes(ymax = median_rank, ymin = median_rank), colour = "black", linetype = 1) +
          geom_jitter(width = .05, alpha = 0.4) +
          ylim(-3, 3) +
          guides(fill = "none") +
          theme_bw() +
          labs(x = "Z score", y = "Layer", title = "Maynard et al Scatter plot")
      }, height = heatmapHeight)

      output$bulk_figure_caption <- renderPrint({
        cat(paste("Fig 1. The dot plots were created using the data from He et al and Maynard et al studies. The raw 
                  RNA-seq data was normalized using z-score normalization. The horizontal bars indicate the median of the 
                  gene expression values for that layer across all genes."))
        })
    }
    
    # AUROC table ----
    output$AUROC_heatmap <- renderPlotly({
      ggplot(data = AUROC_df, mapping = aes(x = dataset, y = Layer, fill = AUROC)) +
        geom_tile() +
        scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(0,1)) +
        scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
        #Puts stars on layer marker annotations
        geom_text(aes(label = signif_marker), size = 7, vjust = 1)
    })
    
    output$AUROC_heatmap_caption <- renderPrint({
      cat(paste("Table 1. The layers were ranked in each dataset with respect to the gene expression of the chosen genes, and
                the AUC score was calculated per layer. P-values were calcuated using the Mann-Whitney U test. Stars indicate that
                the p-values passed correction, and were > 0.05."))
    })
    
    # Summary textbox ----
    output$summary_multiple <- renderPrint({
      #count of intersection of submitted genes with total gene list
      cat(paste0(
        "You inputted ",
        length(selected_gene_list_multiple),
        " genes. Of those genes:\n\n",
        sum(selected_gene_list_multiple %in% unique(He_DS1_Human_averaged$gene_symbol)),
        " genes were assayed by He et al., ",
        sum(selected_gene_list_multiple %in% unique(Maynard_dataset_average$gene_symbol)),
        " were assayed by Maynard et al., and ",
        sum(selected_gene_list_multiple %in% unique(Zeng_dataset_long$gene_symbol)),
        " were assayed by Zeng et al.\n\n",
        #Compared genome-wide, the AUC value for the input genes is [?] (p = <as currently setup>).",
        "Between the He and Maynard datasets, across the layers, the input genes had a mean Pearson correlation value of", 
        multi_gene_cor,
        " (p = ",
        p_value_multiple_gene,
        "), which\nranks in the ",
        multi_gene_quantile,
        "th quantile.\n\n",
        "Specifically in the Zeng dataset: \n\n",
        if (length(unlist(layer_specific_gene_list_multiple $`Layer 1`)) == 0 &
            length(unlist(layer_specific_gene_list_multiple $`Layer 2`)) == 0 &
            length(unlist(layer_specific_gene_list_multiple $`Layer 3`)) == 0 &
            length(unlist(layer_specific_gene_list_multiple $`Layer 4`)) == 0 & 
            length(unlist(layer_specific_gene_list_multiple $`Layer 5`)) == 0 &
            length(unlist(layer_specific_gene_list_multiple $`Layer 6`)) == 0 &
            length(unlist(layer_specific_gene_list_multiple $White_matter)) == 0) {
          "There were no genes that marked any layer."
        } else {
          paste0( if (length(unlist(layer_specific_gene_list_multiple$`Layer 1`)) == 0) {
            ""
          } else paste0(length(unlist(layer_specific_gene_list_multiple $`Layer 1`))," gene(s) marked layer 1 (", 
                        paste(layer_specific_gene_list_multiple $`Layer 1`, collapse = ", "), ").\n"),
          if (length(unlist(layer_specific_gene_list_multiple $`Layer 2`)) == 0) {
            ""
          } else paste0(length(unlist(layer_specific_gene_list_multiple $`Layer 2`)), " gene(s) marked layer 2 (", 
                        paste(layer_specific_gene_list_multiple $`Layer 2`, collapse = ", "), ").\n"),
          if (length(unlist(layer_specific_gene_list_multiple $`Layer 3`)) == 0) {
            ""
          } else paste0(length(unlist(layer_specific_gene_list_multiple $`Layer 3`)), " gene(s) marked layer 3 (",
                        paste(layer_specific_gene_list_multiple $`Layer 3`, collapse = ", "), ").\n"),
          if (length(unlist(layer_specific_gene_list_multiple $`Layer 4`)) == 0) {
            ""
          } else paste0(length(unlist(layer_specific_gene_list_multiple $`Layer 4`)), " gene(s) marked layer 4 (",
                        paste(layer_specific_gene_list_multiple $`Layer 4`, collapse = ", "), ").\n"),
          if (length(unlist(layer_specific_gene_list_multiple $`Layer 5`)) == 0) {
            ""
          } else paste0(length(unlist(layer_specific_gene_list_multiple $`Layer 5`)), " gene(s) marked layer 5 (",
                        paste(layer_specific_gene_list_multiple $`Layer 5`, collapse = ", "), ").\n"),
          if (length(unlist(layer_specific_gene_list_multiple $`Layer 6`)) == 0) {
            ""
          } else paste0(length(unlist(layer_specific_gene_list_multiple $`Layer 6`)), " gene(s) marked layer 6 (",
                        paste(layer_specific_gene_list_multiple $`Layer 6`, collapse = ", "), ").\n"),
          sep = "\n")
        }
      ))
    })
    
    # scRNA Heatmaps ----
    
    #Heatmap for scRNA data for GABAergic 
    if (length(selected_gene_list_multiple) <= 30) {
      output$scRNA_heatmap_GABA <- renderPlot({
        ggplot(data = scRNA_GABA, mapping = aes(x = Layer, y = gene_symbol, fill = Mean_Expression)) +
          geom_tile() +
          scale_fill_distiller(palette = "RdYlBu", limits = c(-3,3)) +
          scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
          labs(y = "", x = "", title = "GABAergic expression") 
      }, height = heatmapHeight)
      
      output$scRNA_heatmap_GLUT <- renderPlot({
        ggplot(data = scRNA_GLUT, mapping = aes(x = Layer, y = gene_symbol, fill = Mean_Expression)) +
          geom_tile() +
          scale_fill_distiller(palette = "RdYlBu", limits = c(-3,3)) +
          scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
          labs(y = "", x = "", title = "Glutamatergic expression") 
      }, height = heatmapHeight)
      
      output$scRNA_heatmap_NON <- renderPlot({
        ggplot(data = scRNA_NON, mapping = aes(x = Layer, y = gene_symbol, fill = Mean_Expression)) +
          geom_tile() +
          scale_fill_distiller(palette = "RdYlBu", limits = c(-3,3)) +
          scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
          labs(y = "", x = "", title = "Non-neuronal expression") 
      }, height = heatmapHeight)
      
      #scRNA heatmap caption
      output$scRNA_figure_caption <- renderPrint({
        cat(paste("Fig 2. The heatmaps were created using the scRNA-seq data from the Allen Brain Institute which was
                taken from live tissue, specifically from the middle temporal gyrus (MTG). The data was first
                log(x+1) normalized on a per-gene basis and the mean of all samples were taken on a per-gene and layer basis. 
                The mean was then transformed using z-score normalization. As there were three cell types identified, the heatmaps 
                represent the expression of the selected genes on a per-cell type basis."))
      })
    } else {
      output$scRNA_heatmap_GABA <- renderPlot({
        scRNA_GABA %<>% inner_join(scRNA_GABA %>% group_by(Layer) %>% summarize(median_rank = median(Mean_Expression)), by = "Layer") %>%
          mutate(Layer = factor(Layer, levels = c("L1", "L2", "L3", "L4", "L5", "L6"))) %>%
          ggplot(aes(x = Layer, y = Mean_Expression, group = Layer, names = gene_symbol, fill = Layer)) +
          geom_errorbar(aes(ymax = median_rank, ymin = median_rank), colour = "black", linetype = 1) +
          geom_jitter(width = .05, alpha = 0.4) +
          ylim(-3, 3) +
          guides(fill = "none") +
          theme_bw() +
          labs(x = "Z score", y = "Layer", title = "GABAergic Expression")
      }, height = heatmapHeight)
      
      output$scRNA_heatmap_GLUT <- renderPlot({
        scRNA_GLUT %<>% inner_join(scRNA_GABA %>% group_by(Layer) %>% summarize(median_rank = median(Mean_Expression)), by = "Layer") %>%
          mutate(Layer = factor(Layer, levels = c("L1", "L2", "L3", "L4", "L5", "L6"))) %>%
          ggplot(aes(x = Layer, y = Mean_Expression, group = Layer, names = gene_symbol, fill = Layer)) +
          geom_errorbar(aes(ymax = median_rank, ymin = median_rank), colour = "black", linetype = 1) +
          geom_jitter(width = .05, alpha = 0.4) +
          ylim(-3, 3) +
          guides(fill = "none") +
          theme_bw() +
          labs(x = "Z score", y = "Layer", title = "Glutamatergic Expression")
      }, height = heatmapHeight)
      
      output$scRNA_heatmap_NON <- renderPlot({
        scRNA_NON %<>% inner_join(scRNA_GABA %>% group_by(Layer) %>% summarize(median_rank = median(Mean_Expression)), by = "Layer") %>%
          mutate(Layer = factor(Layer, levels = c("L1", "L2", "L3", "L4", "L5", "L6"))) %>%
          ggplot(aes(x = Layer, y = Mean_Expression, group = Layer, names = gene_symbol, fill = Layer)) +
          geom_errorbar(aes(ymax = median_rank, ymin = median_rank), colour = "black", linetype = 1) +
          geom_jitter(width = .05, alpha = 0.4) +
          ylim(-3, 3) +
          guides(fill = "none") +
          theme_bw() +
          labs(x = "Z score", y = "Layer", title = "Non-neuronal Expression")
      }, height = heatmapHeight)
      
      output$scRNA_figure_caption <- renderPrint({
        cat(paste("Fig 2. The dot plots were created using the scRNA-seq data from the Allen Brain Institute which was
                taken from live tissue, specifically from the middle temporal gyrus (MTG). The data was first
                log(x+1) normalized on a per-gene basis and the mean of all samples were taken on a per-gene and layer basis. 
                The mean was then transformed using z-score normalization. As there were three cell types identified, the dot plots 
                represent the expression of the selected genes on a per-cell type basis. The line represents the median of the gene 
                expression values for all genes in that layer."))
      })
    }
  })
}




# Run the app ----
shinyApp(ui = ui, server = server)
