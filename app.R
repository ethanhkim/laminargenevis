## App script ##


# Load required libraries
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
library(stringi)
library(here)
library(scales)
library(ggdendro)
library(data.table)
library(shinycssloaders)
library(DT)

# Source processing scripts #
source("string_processing.R") 
source("data_processing.R")

# Load in data #
load(here("data", "processed", "He_DS1_logCPM_dataset.Rdata"), verbose = TRUE)
load(here("data", "processed", "Maynard_logCPM_dataset.Rdata"), verbose = TRUE)
load(here("data", "processed", "Allen_logCPM_dataset.Rdata"), verbose = TRUE)
load(here("data", "processed", "He_DS1_logCPM_filtered_dataset.Rdata"), verbose = TRUE)
load(here("data", "processed", "Maynard_logCPM_filtered_dataset.Rdata"), verbose = TRUE)
load(here("data", "processed", "Allen_logCPM_filtered_dataset.Rdata"), verbose = TRUE)
load(here("data", "processed", "He_Maynard_gene_correlation.Rdata"), verbose = TRUE)
load(here("data", "processed", "layer_marker_table.Rdata"))


# Define UI ----
ui <- fluidPage(
  shinyjs::useShinyjs(),
  tags$head(
    includeHTML("google-analytics.html"),
    tags$style(
    HTML("
          #summary_statistics_single {
              font-family: 
                'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif;
              font-size: 
                14px; }
          #summary_markers_single {
              font-family: 
                'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif;
              font-size: 
                14px; }
          #summary_statistics_multiple {
              font-family: 
                'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif;
              font-size: 
                14px; }
          #summary_markers_multiple {
              font-family: 
                'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif;
              font-size: 
                14px; }
          ")),
    # Wrap text in long text boxes
    tags$style(type='text/css', '#summary_statistics_single {white-space: pre-wrap;}'),
    tags$style(type='text/css', '#summary_statistics_multiple {white-space: pre-wrap;}'),
  ),
  navbarPage(
    title = "LaminaRGeneVis", 
    # Visualize gene expression across layers through heatmap or barplot
    tabPanel(
      title = "Gene Visualization",
      # Side panel for gene selection
      sidebarLayout(
        sidebarPanel(
          # Single or multiple gene selector
          radioButtons(
            inputId = "selector",
            label =  "Choose to examine either a single gene, 
            or multiple genes: \n", 
            choices = c("Single Gene", "Multiple Genes")),
          # Single input side-panel  ----
          conditionalPanel(
            condition = "input.selector == 'Single Gene'",
            # Drop-down bar
            selectizeInput(
              inputId = "genelist", label = "Input gene:", 
              choices = NULL, selected = 'RELN', multiple = FALSE, 
              options = NULL),
            # Submit button
            actionButton(inputId = "submit_barplot", label = "Submit")),
          # Multiple input side-panel ----
          conditionalPanel(
            condition = "input.selector == 'Multiple Genes'",
            # Text box
            textAreaInput(
              inputId = "multiple_genelist", 
              label = "Input your gene list:", 
              value = 'RELN, CUX2, FOXP2, RASGRF2'),
            # Submit button
            actionButton(inputId = "submit_heatmap", label = "Submit")),
          ),
        # Main page
        mainPanel(
          tabsetPanel(
            type = "tabs", id = "tabset",
            #Page for displaying information about datasets ----
            tabPanel(
              title = "App Overview", value = "overview",
              br(),
              h2("Welcome to LaminaRGeneVis!"),
              br(),
              p("This web application allows you to examine layer-specific 
                 gene expression across the cortex and determine layer
                 annotations.", style='font-size:19px'),
               br(),
               h3("App Workflow:"),
               p("Here's a quick overview of the app! Multiple datasets of 
                  RNA-seq expression (described below) have been standardized 
                  for ease of comparison as shown in (A). Users have the option 
                  of choosing to examine either", strong("Single Gene"), "or", 
                  strong("Multiple Genes")," in the HGNC gene symbol 
                  format. Choosing Single Gene will give you a plot similar 
                  to (B), and choosing Multiple Genes will give you multiple 
                  plots, including the plot shown in (C).", 
                  style = 'font-size:17px'),
              br(),
              # Figure
              img(src = 'pageFigure.png', 
                  style = "display: block; margin-left: auto; 
                           margin-right: auto;"),
              br(),
              h3("Source Data:"),
              p("The data for this application has been sourced from 
                 these following studies and institutions:", 
                style='font-size:17px'),
              # Bullet points for describing datasets
              tags$ul(
                # He et al description
                tags$li(p(a("He et al. (2017):", 
                             href = "https://pubmed.ncbi.nlm.nih.gov/28414332/", 
                             target = "_blank"), 
                             p("This study assayed the whole genome using
                                high-throughput RNA-seq in samples from the 
                               DLPFC.")), style = 'font-size:17px'),
                # Maynard et al description
                tags$li(p(a("Maynard et al. (2021):", 
                             href = "https://www.nature.com/articles/s41593-020-00787-0", 
                             target = "_blank"), 
                             p("This study assayed the whole genome through 
                                the Visium Platform (10X Genomics) in samples 
                                from the DLPFC.")), style = 'font-size:17px'),
                # AIBS description
                tags$li(p(a("Allen Institute for Brain Science (AIBS): 
                             Cell-Type Database", 
                             href = "https://portal.brain-map.org/atlases-and-data/
                                     rnaseq/human-multiple-cortical-areas-smart-seq", 
                             target = "_blank"), 
                             p("This dataset contains multiple cortical regions 
                                and assays the whole genome across roughly 
                                49,000 single-cell nuclei.")), 
                                style = 'font-size:17px'),
                ),
              br(),
              # Github link
              tags$div(
                'The code for the application and the analyses used are available on ',
                tags$a(href = "href = 'https://github.com/ethanhkim/transcriptome_app",
                       "Github."), style = 'font-size:17px')),
            # Single or multiple gene visualizations ----
            tabPanel(
              title = "Gene Visualization", value = "visualization",
              ### Single gene visualization ----
              br(),
              conditionalPanel(
                h3("Layer-specific gene expression"),
                #Only show when the input selector is Single
                condition = "input.selector == 'Single Gene'",
                #Output barplot visualization
                br(),
                plotOutput("Barplot") %>% withSpinner(),
                br(),
                br(),
                p(textOutput("barplot_caption"),style='font-size:15px'),
                br(),
                br(),
                h4(textOutput('summary_statistics_title')),
                p(verbatimTextOutput("summary_statistics_single"),
                  style='font-size:15px'),
                br(),
                h4(textOutput('summary_markers_title')),
                p(verbatimTextOutput("summary_markers_single"),
                  style='font-size:15px'),
              ),
              ### Multiple gene visualization ----
              conditionalPanel(
                h3("Layer-specific Figures"),
                #Only show when the input selector is Multiple
                condition = "input.selector == 'Multiple Genes'",
                br(),
                # Show AUC heatmap
                h4(textOutput('AUC_heatmap_title')),
                plotOutput("AUC_heatmap") %>% withSpinner(),
                p(textOutput("AUC_heatmap_caption"),style='font-size:15px'),
                br(),
                #Output heatmap visualizations for bulk tissue RNA-seq
                h4(textOutput('bulk_figure_title')),
                br(),
                plotOutput("He_figure", height = "auto") %>% withSpinner(),
                plotOutput("Maynard_figure", height = "auto") %>% withSpinner(),
                br(),
                p(textOutput("bulk_figure_caption"),style='font-size:15px'),
                br(),
                #Output heatmap visualizations for Allen data per cell class type
                h4(textOutput('AIBS_figure_title')),
                br(),
                plotOutput("AIBS_figure_GABA", height = "auto") %>% withSpinner(),
                plotOutput("AIBS_figure_GLUT", height = "auto") %>% withSpinner(),
                plotOutput("AIBS_figure_NONN", height = "auto") %>% withSpinner(),
                p(textOutput("AIBS_figure_caption"),style='font-size:15px'),
                br(),
                h4(textOutput('summary_statistics_multiple_title')),
                p(verbatimTextOutput("summary_statistics_multiple"),style='font-size:15px'),
                br(),
                h4(textOutput('summary_markers_multiple_title')),
                p(verbatimTextOutput("summary_markers_multiple")),
              )
            )
          )
        )
      )
    )
  )
)

# Server ----
server <- function(input, output, session) {
  
  # List of genes that are common through the He and Maynard datasets
  common_genelist <- intersect(He_DS1_logCPM_dataset$gene_symbol, 
                               Maynard_logCPM_dataset$gene_symbol) %>%
    sort()
  
  # Table of layer-marker annotations
  layer_marker_lookup_tbl <- layer_marker_table %>%
    pivot_longer(
      cols=c("He", "Maynard"), 
      names_to = "source_dataset",
      values_to = "layer_marker")
  
  # Calculate top and bottom 5 percentile for values
  top_and_bottom_5th_perc <- top_and_bottom_quantile(
    Maynard_logCPM_dataset, 
    He_DS1_logCPM_dataset, 
    Allen_logCPM_dataset
  )
  
  # Give available values for drop-down bar for single gene
  updateSelectizeInput(session, inputId = "genelist", selected = 'RELN',
                       choices = common_genelist, server = TRUE)
  
  #Hide plot and titles until when genes are submitted
  output$Barplot <- NULL
  output$bulk_figure_title <- NULL
  output$He_figure <- NULL
  output$Maynard_figure <- NULL 
  output$AUC_heatmap_title <- NULL
  output$heatmap_caption <- NULL
  output$AUC_heatmap <- NULL
  output$scRNA_figure_title <- NULL
  output$AIBS_figure_GABA <- NULL
  output$AIBS_figure_GLUT <- NULL
  output$AIBS_figure_NONN <- NULL
  output$summary_statistics_title <- NULL
  output$summary_markers_title <- NULL
  
  # Single gene input ----
  observeEvent(input$submit_barplot, {
    #When the submit button is pressed, change to the Visualization page
    updateTabsetPanel(session, "tabset", selected = "visualization")
    
    # List of selected gene(s)
    selected_gene_list_single <- isolate(process_gene_input(input$genelist))
    # Process dataset to correct format for heatmap and barplot
    Barplot_data <- process_barplot_data(
      input_genelist = selected_gene_list_single,
      He_dataset = He_DS1_logCPM_dataset, 
      Maynard_dataset = Maynard_logCPM_dataset,
      Allen_dataset = Allen_logCPM_dataset)
    
    # Layer marker for input gene
    He_layer_marker <- separate_layers(
      input_table = layer_marker_lookup_tbl, 
      input_genelist = selected_gene_list_single,
      source = "He")
    Maynard_layer_marker <- separate_layers(
      input_table = layer_marker_lookup_tbl,
      input_genelist = selected_gene_list_single,
      source = "Maynard")
    
    # Single gene barplot ----
    output$Barplot <- renderPlot({
      ggplot(
        data = Barplot_data, 
        aes(x = layer, y = expression, fill = source_dataset, 
            group = source_dataset)) +
        geom_bar(stat = "identity", position = "dodge", width = 0.75) + 
        ggtitle(paste0('Expression of ', selected_gene_list_single, 
                       ' across the human neocortex')) +
        geom_hline(yintercept = top_and_bottom_5th_perc$top_5,
                   color="black", linetype="dashed") +
        geom_hline(yintercept = top_and_bottom_5th_perc$bottom_5,
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
    }) 
    
    # Barplot caption
    output$barplot_caption <- renderPrint({
      cat(paste("Fig 1. The barplots were created using data from He et al, 
                Maynard et al and the Allen Institute for Brain Science (AIBS). 
                Raw RNA-seq data was processed and normalized through counts 
                per million (CPM) and log2-transformed. The horizontal dashed 
                lines represent the value of the top (95th) and bottom (5th) 
                quantile of normalized expression values across all 
                normalized data."))
    })
    
    #Filter for selected genes from table containing Zeng et al layer marker annotations
    layer_marker_table_single <- layer_marker_table %>%
      filter(selected_gene_list_single %in% gene_symbol)
    
    # Generate correlation value for single gene
    single_gene_cor <- single_gene_correlation(
      input_gene = selected_gene_list_single, 
      He_dataset = He_DS1_logCPM_dataset, 
      Maynard_dataset = Maynard_logCPM_dataset)
    # Generate what quantile the gene-to-gene correlation would be in
    single_gene_quantile <- quantile_distribution(
      dataset_correlation = He_Maynard_gene_correlation, 
      selected_gene_correlation = single_gene_cor)
    # Generate p-value for gene-to-gene correlation
    p_value_single_gene <- wilcoxtest(
      input_genelist = selected_gene_list_single, 
      He_dataset = He_DS1_logCPM_dataset, 
      Maynard_dataset = Maynard_logCPM_dataset, 
      He_Maynard_gene_correlation)
    
    # Agreement Statistics textbox title
    output$summary_statistics_title <- renderPrint({
      cat(paste('Agreement Statistics'))
    })
    
    # Summary textbox
    output$summary_statistics_single <- renderPrint({
      cat(paste0(
        assayed_gene_string(
          genelist = selected_gene_list_single, 
          He_df = He_DS1_logCPM_dataset,
          Maynard_df = Maynard_logCPM_dataset, 
          single_or_multiple = "single"),
        stats_string(
          genelist = selected_gene_list_single, 
          correlation = single_gene_cor, 
          p_value = p_value_single_gene, 
          quantile_stat = single_gene_quantile, 
          single_or_multiple = "single")
      ))
    })
    
    # Layer marker textbox title
    output$summary_markers_title <- renderPrint({
      cat(paste('Layer Markers'))
    })
    
    # Summary textbox
    output$summary_markers_single <- renderPrint({
      cat(paste0(
        layer_marker_preempt_string("He", "single", selected_gene_list_single),
        paste0(layer_marker_string(He_layer_marker, "single")),
        layer_marker_preempt_string("Maynard", "single", selected_gene_list_single),
        paste0(layer_marker_string(Maynard_layer_marker,  "single"))
      ))
    })
  })
  
  
  ## Multiple genes ----
  observeEvent(input$submit_heatmap, {
    updateTabsetPanel(session, "tabset", selected = "visualization")
    
    # Create a vector of selected genes
    selected_gene_list_multiple <- isolate(process_gene_input(input$multiple_genelist))
    # Process data with input genes
    He_heatmap_data <- process_heatmap_data(
      source = "He", source_dataset = He_DS1_logCPM_dataset,
      input_genelist = selected_gene_list_multiple)
    Maynard_heatmap_data <- process_heatmap_data(
      source = "Maynard", source_dataset = Maynard_logCPM_dataset, 
      input_genelist = selected_gene_list_multiple)
    AIBS_GABA_heatmap_data <- process_heatmap_data(
      source = "Allen", source_dataset = Allen_logCPM_dataset,
      input_genelist = selected_gene_list_multiple, cell_type = "GABAergic")
    AIBS_GLUT_heatmap_data <- process_heatmap_data(
      source = "Allen", source_dataset = Allen_logCPM_dataset,
      input_genelist = selected_gene_list_multiple, cell_type = "Glutamatergic")
    AIBS_NONN_heatmap_data <- process_heatmap_data(
      source = "Allen", source_dataset = Allen_logCPM_dataset,
      input_genelist = selected_gene_list_multiple, cell_type = "Non-neuronal")

    He_layer_marker <- separate_layers(layer_marker_lookup_tbl, 
                                       selected_gene_list_multiple,
                                       "He")
    Maynard_layer_marker <- separate_layers(layer_marker_lookup_tbl, 
                                            selected_gene_list_multiple,
                                            "Maynard")
    
    # Generate correlation value for multiple gene and the quantile that value belongs in
    multi_gene_cor <- multi_gene_correlation(selected_gene_list_multiple, 
                                             He_DS1_logCPM_dataset, Maynard_logCPM_dataset)
    multi_gene_quantile <- quantile_distribution(He_Maynard_gene_correlation, multi_gene_cor)
    p_value_multiple_gene <- wilcoxtest(selected_gene_list_multiple, 
                                        He_DS1_logCPM_dataset, Maynard_logCPM_dataset, 
                                        He_Maynard_gene_correlation)
    
    # Code for AUROC analysis adapted from Derek Howard & Leon French - refer to data_processing.R
    AUROC_bulk_data <- AUROC_bulk(He_DS1_logCPM_filtered_dataset, 
                                  Maynard_logCPM_filtered_dataset, 
                                  selected_gene_list_multiple) 
    AUROC_AIBS_data <- AUROC_AIBS(Allen_logCPM_filtered_dataset, 
                                  selected_gene_list_multiple)
    AUROC_df <- AUROC_data(AUROC_bulk_data, AUROC_AIBS_data)
    
    # Set dynamic heatmap height
    heatmapHeight <- heatmap_height(selected_gene_list_multiple)
    
    ### AUC heatmap ----
    # Heatmap figure title
    output$AUC_heatmap_title <- renderPrint({
      cat(paste('Layer-specific Gene Enrichment Heatmap'))
    })
    
    # Heatmap figure
    output$AUC_heatmap <- renderPlot({
      ggplot(data = AUROC_df, mapping = aes(x = dataset, y = Layer, fill = AUROC)) +
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
        theme(axis.text.x = element_text(angle = 25, size = 13, hjust = 0.95),
              axis.text.y = element_text(size = 13),
              axis.title.x = element_text(size = 17),
              axis.title.y = element_text(size = 17),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 11),
              legend.key.size = unit(1, 'cm'))
    })
    
    # AUC heatmap figure caption
    output$AUC_heatmap_caption <- renderPrint({
      cat(paste("Fig 1. Each dataset was filtered for genes with CPM > 0.1 
                 across all layers prior to log2 and z-score normalization. 
                 Layers were ranked in each dataset with respect to the gene 
                 expression of the chosen genes. The AUC score was calculated per 
                 layer using the rankings in each respective dataset. P-values 
                 were calcuated using the Mann-Whitney U test and adjusted for
                 multiple test correction through Bonferroni correction; stars (*) 
                 indicate p < 0.05."))
    })
    
    # If less than 30 genes input, create heatmaps for bulk tissue
    if (length(selected_gene_list_multiple) <= 30) {
      
      # Bulk-tissue expression figure title
      output$bulk_figure_title <- renderPrint({
        cat(paste('Bulk-tissue Gene Expression Heatmaps'))
      })
      
      # Bulk tissue heatmaps ----
      output$He_figure <- renderPlot({
        He_heatmap_data %>%
          ggplot(mapping = aes(x = gene_symbol, y = layer, 
                                                     fill = expression)) +
          geom_tile() +
          scale_fill_distiller(palette = "RdYlBu", limits = c(-3, 3)) +
          scale_y_discrete(expand=c(0,0)) + 
          scale_x_discrete(expand=c(0,0)) +
          labs(y = "Cortical Layer", x = "", title = "He et al data", 
               fill = "Z-Score\nNormalized\nExpression") +
          theme(axis.text.x = element_text(size = 13, angle = 25, hjust = 0.95), 
                axis.text.y = element_text(size = 13),
                axis.title.x = element_text(size = 17),
                axis.title.y = element_text(size = 17),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 11),
                plot.title = element_text(size=17))
      }, height = heatmapHeight)
      output$Maynard_figure <- renderPlot({
        Maynard_heatmap_data %>%
        ggplot(data = Maynard_heatmap_data, 
               mapping = aes(x = gene_symbol, y = layer, fill = expression)) +
          geom_tile() +
          scale_fill_distiller(palette = "RdYlBu", limits = c(-3, 3)) +
          scale_y_discrete(expand=c(0,0)) + 
          scale_x_discrete(expand=c(0,0)) +
          labs(y = "Cortical Layer", x = "", title = "Maynard et al data",
               fill = "Z-Score\nNormalized\nExpression") +
          theme(axis.text.x = element_text(angle = 25, size = 13, hjust = 0.95), 
                axis.text.y = element_text(size = 13),
                axis.title.x = element_text(size = 17),
                axis.title.y = element_text(size = 17),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 11),
                plot.title = element_text(size=17))
      }, height = heatmapHeight)
      
      # Bulk-tissue expression figure caption
      output$bulk_figure_caption <- renderPrint({
        cat(paste("Fig 2. The heatmaps were created using data from He et al 
                   Maynard et al. Raw RNA-seq data was processed and normalized 
                   through counts per million (CPM), log2-transformed and
                   z-score normalized."))
      })
    } else {
      
      # Bulk-tissue expression figure title
      output$bulk_figure_title <- renderPrint({
        cat(paste('Bulk-tissue Gene Expression Scatterplots'))
      })
      
      # Bulk-tissue Scatterplots ----
      output$He_figure <- renderPlot({
        He_heatmap_data %<>% inner_join(He_heatmap_data %>% group_by(layer) %>% 
                                          summarize(median_rank = median(expression)), 
                                        by = "layer") %>%
          mutate(layer = factor(layer, levels = c("L1", "L2", "L3", "L4", "L5", 
                                                  "L6", "WM"))) %>%
          ggplot(aes(x = layer, y = expression, group = layer, names = gene_symbol, 
                     fill = layer)) +
          # Median line
          geom_errorbar(aes(ymax = median_rank, ymin = median_rank), 
                        colour = "black", linetype = 1) +
          geom_jitter(width = .05, alpha = 0.4) +
          # Limit set to min and max expression value across all datasets
          ylim(-3,3) +
          guides(fill = "none") +
          theme_bw() +
          labs(x = "\nCortical Layer", y = "Z-Score Normalized Expression", 
               title = "He et al Scatter plot") +
          theme(axis.text.x = element_text(size = 13), 
                axis.text.y = element_text(size = 13),
                axis.title.x = element_text(size = 17),
                axis.title.y = element_text(size = 17),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 11),
                plot.title = element_text(size=17))
      }, height = heatmapHeight)
      
      output$Maynard_figure <- renderPlot({
        Maynard_heatmap_data %<>% inner_join(Maynard_heatmap_data %>% 
                                               group_by(layer) %>% 
                                               summarize(median_rank = median(expression)), 
                                             by = "layer") %>%
          mutate(layer = factor(layer, levels = c("L1", "L2", "L3", "L4", "L5", 
                                                  "L6", "WM"))) %>%
          ggplot(aes(x = layer, y = expression, group = layer, names = gene_symbol, 
                     fill = layer)) +
          # Median line
          geom_errorbar(aes(ymax = median_rank, ymin = median_rank), 
                        colour = "black", linetype = 1) +
          geom_jitter(width = .05, alpha = 0.4) +
          # Limit set to min and max expression value across all datasets
          ylim(-3,3) +
          guides(fill = "none") +
          theme_bw() +
          labs(x = "\nCortical Layer", y = "Z-Score Normalized Expression", 
               title = "Maynard et al Scatter plot") +
          theme(axis.text.x = element_text(size = 13), 
                axis.text.y = element_text(size = 13),
                axis.title.x = element_text(size = 17),
                axis.title.y = element_text(size = 17),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 11),
                plot.title = element_text(size=17))
      }, height = heatmapHeight) 
      
      output$bulk_figure_caption <- renderPrint({
        cat(paste("Fig 2. The scatter plots were created using data from He et al 
                   and Maynard et al. Raw RNA-seq data was processed and 
                   normalized through counts per million (CPM), log2-transformed
                   and z-score normalized. The horizontal bars indicate the median 
                   of the gene expression values for that layer across all genes."))
      })
    }
    
    ## AIBS Figures ##
    
    # AIBS heatmaps ----
    if (length(selected_gene_list_multiple) <= 30) {
      
      # AIBS figure title
      output$AIBS_figure_title <- renderPrint({
        cat(paste('Cell type-specific Gene Expression Heatmaps'))
      })
      
      # GABAergic heatmap ----
      output$AIBS_figure_GABA <- renderPlot({
        ggplot(data = AIBS_GABA_heatmap_data, 
               mapping = aes(x = gene_symbol, y = layer, fill = expression)) +
          geom_tile() +
          scale_fill_distiller(palette = "RdYlBu", limits = c(-3, 3)) +
          scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
          labs(y = "Cortical Layer", x = "", title = "GABAergic expression",
               fill = "Z-Score\nNormalized\nExpression") +
          theme(axis.text.x = element_text(angle = 25, size = 13, hjust = 0.95), 
                axis.text.y = element_text(size = 13),
                axis.title.x = element_text(size = 17),
                axis.title.y = element_text(size = 17),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 11),
                plot.title = element_text(size=17))
      }, height = heatmapHeight)
      
      # Glutamatergic heatmap ----
      output$AIBS_figure_GLUT <- renderPlot({
        ggplot(data = AIBS_GLUT_heatmap_data, 
               mapping = aes(x = gene_symbol, y = layer, fill = expression)) +
          geom_tile() +
          scale_fill_distiller(palette = "RdYlBu", limits = c(-3, 3)) +
          scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
          labs(y = "Cortical Layer", x = "", title = "Glutamatergic expression",
               fill = "Z-Score\nNormalized\nExpression") +
          theme(axis.text.x = element_text(angle = 25, size = 13, hjust = 0.95), 
                axis.text.y = element_text(size = 13),
                axis.title.x = element_text(size = 17),
                axis.title.y = element_text(size = 17),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 11),
                plot.title = element_text(size=17))
      }, height = heatmapHeight)
      
      # Non-neuronal heatmap ----
      output$AIBS_figure_NONN <- renderPlot({
        ggplot(data = AIBS_NONN_heatmap_data, 
               mapping = aes(x = gene_symbol, y = layer, fill = expression)) +
          geom_tile() +
          scale_fill_distiller(palette = "RdYlBu", limits = c(-3, 3)) +
          scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
          labs(y = "Cortical Layer", x = "", title = "Non-neuronal expression",
               fill = "Z-Score\nNormalized\nExpression") +
          theme(axis.text.x = element_text(angle = 25, size = 13, hjust = 0.95), 
                axis.text.y = element_text(size = 13),
                axis.title.x = element_text(size = 17),
                axis.title.y = element_text(size = 17),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 11),
                plot.title = element_text(size=17))
      }, height = heatmapHeight)
      
      #AIBS heatmap caption
      output$AIBS_figure_caption <- renderPrint({
        cat(paste("Fig 3. Heatmaps were created using the snRNA-seq data from 
                   the Allen Institute of Brain Science (AIBS), specifically 
                   from the middle temporal gyrus (MTG). The raw RNA-seq data 
                   was pseudo-bulked at the layer level and normalized on using 
                   counts per million (CPM), log2-transformed and z-score 
                   normalized. The heatmaps represent the expression of the 
                   selected genes on a per cell-type basis, using the labels 
                   provided by AIBS."))
      })
      
    # AIBS scatterplots ----
      
    } else {
      
      # AIBS figure title
      output$AIBS_figure_title <- renderPrint({
        cat(paste('Cell type-specific Gene Expression Scatterplots'))
      })
      
      # GABAergic scatterplot ----
      output$AIBS_figure_GABA <- renderPlot({
        AIBS_GABA_heatmap_data %<>% 
          inner_join(AIBS_GABA_heatmap_data %>% group_by(layer) %>% 
                       summarize(median_rank = median(expression)), 
                                   by = "layer") %>%
          mutate(layer = factor(layer, levels = c("L1", "L2", "L3", "L4", "L5", 
                                                  "L6", "WM"))) %>%
          ggplot(aes(x = layer, y = expression, group = layer, 
                     names = gene_symbol, fill = layer)) +
          geom_errorbar(aes(ymax = median_rank, ymin = median_rank), 
                        colour = "black", linetype = 1) +
          geom_jitter(width = .05, alpha = 0.4) + ylim(c(-3,3)) +
          guides(fill = "none") +
          theme_bw() +
          labs(x = "\nCortical Layer", y = "Z-Score Normalized Expression", 
               title = "GABAergic Expression") +
          theme(axis.text.x = element_text(size = 13), 
                axis.text.y = element_text(size = 13),
                axis.title.x = element_text(size = 17),
                axis.title.y = element_text(size = 17),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 11),
                plot.title = element_text(size=17))
      }, height = heatmapHeight)
      
      # Glutamatergic scatterplot ----
      output$AIBS_figure_GLUT <- renderPlot({
        AIBS_GLUT_heatmap_data %<>% inner_join(AIBS_GLUT_heatmap_data %>% group_by(layer) %>% 
                                     summarize(median_rank = median(expression)), 
                                   by = "layer") %>%
          mutate(layer = factor(layer, levels = c("L1", "L2", "L3", "L4", "L5",
                                                  "L6", "WM"))) %>%
          ggplot(aes(x = layer, y = expression, group = layer, 
                     names = gene_symbol, fill = layer)) +
          geom_errorbar(aes(ymax = median_rank, ymin = median_rank), 
                        colour = "black", linetype = 1) +
          geom_jitter(width = .05, alpha = 0.4) + ylim(c(-3,3)) +
          guides(fill = "none") +
          theme_bw() +
          labs(x = "\nCortical Layer", y = "Z-Score Normalized Expression", 
               title = "Glutamatergic Expression") +
          theme(axis.text.x = element_text(size = 13), 
                axis.text.y = element_text(size = 13),
                axis.title.x = element_text(size = 17),
                axis.title.y = element_text(size = 17),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 11),
                plot.title = element_text(size=17))
      }, height = heatmapHeight)
      
      # Non-neuronal scatterplot ----
      output$AIBS_figure_NONN <- renderPlot({
        AIBS_NONN_heatmap_data %<>% inner_join(AIBS_NONN_heatmap_data 
                                                %>% group_by(layer) %>% 
                                    summarize(median_rank = median(expression)), 
                                  by = "layer") %>%
          mutate(layer = factor(layer, levels = c("L1", "L2", "L3", "L4", "L5", 
                                                  "L6", "WM"))) %>%
          ggplot(aes(x = layer, y = expression, group = layer, 
                     names = gene_symbol, fill = layer)) +
          geom_errorbar(aes(ymax = median_rank, ymin = median_rank), 
                        colour = "black", linetype = 1) +
          geom_jitter(width = .05, alpha = 0.4) + ylim(c(-3,3)) +
          guides(fill = "none") +
          theme_bw() +
          labs(x = "\nCortical Layer", y = "Z-Score Normalized Expression", 
               title = "Non-neuronal Expression") +
          theme(axis.text.x = element_text(size = 13), 
                axis.text.y = element_text(size = 13, hjust = 0.05),
                axis.title.x = element_text(size = 17),
                axis.title.y = element_text(size = 17),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 11),
                plot.title = element_text(size=17))
      }, height = heatmapHeight)
      
      output$AIBS_figure_caption <- renderPrint({
        cat(paste("Fig 3. The scatter plots were created using the scRNA-seq data 
                  from the Allen Institute of Brain Science which was created from 
                  live tissue, specifically from the middle temporal gyrus (MTG). 
                  The data was first log(x+1) normalized on a per-gene basis and 
                  the mean of all samples were taken on a per-gene and layer basis. 
                  The mean was then transformed using z-score normalization. As 
                  there were three cell types identified, the scatter plots represent 
                  the expression of the selected genes on a per-cell type basis. 
                  The line represents the median of the gene expression values 
                  for all genes in that layer."))
      })
    }
    
    # Summary textbox ----
    
    # Agreement Statistics textbox title
    output$summary_statistics_multiple_title <- renderPrint({
      cat(paste('Agreement Statistics'))
    })
    
    # Summary textbox
    output$summary_statistics_multiple <- renderPrint({
      cat(paste0(
        assayed_gene_string(selected_gene_list_multiple, He_DS1_logCPM_dataset,
                            Maynard_logCPM_dataset, "multiple"),
        #Compared genome-wide, the AUC value for the input genes is [?] (p = <as currently setup>).",
        stats_string(selected_gene_list_multiple, multi_gene_cor, 
                     p_value_multiple_gene, multi_gene_quantile, "multiple")
      ))
    })
    
    # Layer marker textbox title
    output$summary_markers_multiple_title <- renderPrint({
      cat(paste('Layer Markers'))
    })
    
    # Summary textbox
    output$summary_markers_multiple <- renderPrint({
      cat(paste0(
        layer_marker_preempt_string("He", "multiple", selected_gene_list_multiple),
        paste0(layer_marker_string(He_layer_marker, "multiple")),
        layer_marker_preempt_string("Maynard", "multiple", selected_gene_list_multiple),
        paste0(layer_marker_string(Maynard_layer_marker, "multiple"))
      ))
    })
  })
}




# Run the app ----
shinyApp(ui = ui, server = server)
