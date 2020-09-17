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


# Define UI ----
ui <- fluidPage(
  shinyjs::useShinyjs(),
  tags$head(includeHTML("google-analytics.html")),
  navbarPage(title = "Gene Expression Profile Comparison", 
             
    # Visualize gene expression across layers through heatmap or barplot    
    tabPanel(title = "Gene Visualization",
      sidebarLayout(
        sidebarPanel(
          # Input: Selector for which genes to visualize
          radioButtons(
            inputId = "selector", label =  "Single or multiple genes?",
            choices = c("Single", "Multiple")
          ),
          # UI - Single Input  ----
          conditionalPanel(
            condition = "input.selector == 'Single'",
            selectizeInput(
              inputId = "genelist", label = "Input gene:", choices = NULL,
              selected = NULL, multiple = FALSE, options = NULL),
            actionButton(inputId = "submit_barplot", label = "Submit")
          ),
          # UI - Multiple Input ----     
          conditionalPanel(
            condition = "input.selector == 'Multiple'",
            textAreaInput(
              inputId = "multiple_genelist", 
              label = "Input your gene list:", 
              placeholder = "GAD1, CCK, GRIN1"),
            selectInput(
              inputId = "region_label_choice",
              label = "Choose a region to examine with scRNA-seq data:",
              choices = c("A1C", "MTG", "V1C", "CgG", "M1lm", "M1ul", "S1lm", "S1ul")
            ),
            actionButton(inputId = "submit_heatmap", label = "Submit")
          ),
                        
        ),
                          
        mainPanel(
          tabsetPanel(type = "tabs", id = "tabset",
            tabPanel(title = "Dataset Overview", value = "overview",
              br(),
              h3("Welcome to the Gene Visualization project!"),
              br(),
              h4("This web application allows you to examine layer-specific gene expression
                across the cortex and determine layer annotations. "),
              br(),
              h4("The data for this application has been sourced from these following studies:"),
              br(),
              br(),
              h3(a("Zeng et al. (2012)", href = "https://pubmed.ncbi.nlm.nih.gov/22500809/", target = "_blank")),
              h4("This study from the Allen Brain Institute examined 46 neurotypical brains and assayed roughly 1000 genes 
                through in-situ hybridization."),
              br(),
              h3(a("He et al. (2017)", href = "https://pubmed.ncbi.nlm.nih.gov/28414332/", target = "_blank")),
              h4("This study assayed the whole genome using high-throughput RNA-seq."),
              br(),
              h3(a("Maynard et al. (2020)*", href = "https://www.biorxiv.org/content/10.1101/2020.02.28.969931v1", target = "_blank")),
              h4("*This study is currently a pre-print; it assayed the whole genome through the 10X Genomics Visium Platform.")
            ),
            tabPanel(title = "Gene Visualization", value = "visualization",
             # Multiple Input - Heatmaps and AUC ----
              br(),
              conditionalPanel(
                h4("Layer-specific Heatmaps"),
                condition = "input.selector == 'Multiple'",
                br(),
                plotOutput("He_heatmap", height = "auto") %>% withSpinner(),
                plotOutput("Maynard_heatmap", height = "auto") %>% withSpinner(),
                br(),
                h5(textOutput("heatmap_caption")),
                br(),
                h4(verbatimTextOutput("summary_multiple")),
                br(),
                br(),
                DT::dataTableOutput("table", width = "100%", height = "auto") %>% withSpinner(),
                br(),
                br(),
                br(),
                plotOutput("scRNA_barplot", height = "800px") %>% withSpinner(),
                h5(textOutput("scRNA_caption")),
                br(),
                br(),
              ),
              # Single Input - Barplot ----
              conditionalPanel(
                h3("Layer-specific gene expression"),
                condition = "input.selector == 'Single'",
                plotOutput("Barplot") %>% withSpinner(),
                br(),
                br(),
                br(),
                br(),
                h4(verbatimTextOutput("summary_single"))
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
  output$He_heatmap <- NULL
  output$Maynard_heatmap <- NULL 
  output$heatmap_caption <- NULL
  output$table <- DT::renderDT(NULL)
  output$scRNA_barplot <- NULL
  output$scRNA_caption <- NULL

  # Single gene input ----
  observeEvent(input$submit_barplot, {
    updateTabsetPanel(session, "tabset", selected = "visualization")

    # List of selected gene(s)
    selected_gene_list_single <- isolate(process_gene_input(input$genelist))
    
    # Process dataset to correct format for heatmap and barplot
    Barplot_data <- process_barplot_data(selected_gene_list_single, He_DS1_Human_averaged, Maynard_dataset_average)
    
    ## Single gene input:
    if (input$selector == "Single") {
      output$Barplot <- renderPlot({
        ggplot(data = Barplot_data, aes(x = Layer, y = Z_score, fill = Dataset, group = Dataset)) +
          geom_bar(stat = "identity", position = "dodge", width = 0.75) + theme(axis.text.x = element_text(angle = 45)) +
          geom_text(aes(label = layer_label, group = Dataset, 
                        vjust = ifelse(Z_score >= -0.1, 0, 2.5)), 
                    position = position_dodge(width = 0.75)) +
          geom_hline(yintercept = 1.681020) +
          geom_hline(yintercept = -1.506113) +
          labs(caption = "Barplot of z-scored gene expression levels.") 
      }) 
      
      # Filter for selected genes from table containing Zeng layer marker annotations
      layer_marker_table_single <- layer_marker_table %>%
        dplyr::filter(gene_symbol %in% selected_gene_list_single)
      
      layer_specific_gene_list_single <- separate_layers(layer_marker_table_single, selected_gene_list_single)
      names(layer_specific_gene_list_single) <- c("Layer 1", "Layer 2", "Layer 3", "Layer 4",
                                                  "Layer 5", "Layer 6", "White_matter")
      
      # Generate correlation value for single gene, the quantile and associated p-value
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
              "Specificaly in the Zeng dataset, ",
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
    }
  })
  
  #Multiple gene input ----
  observeEvent(input$submit_heatmap, {
    updateTabsetPanel(session, "tabset", selected = "visualization")

    # Input genes
    selected_gene_list_multiple <- isolate(process_gene_input(input$multiple_genelist))
    # Process He & Maynar data according to input genes
    He_heatmap_data <- process_heatmap_function(He_DS1_Human_averaged, selected_gene_list_multiple)
    Maynard_heatmap_data <- process_heatmap_function(Maynard_dataset_average, selected_gene_list_multiple)
    # Input area (for scRNA)
    selected_scRNA_region <- input$region_label_choice
    
    ## Multiple gene Input
    
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
    
    AUROC_table <- AUROC_function(He_DS1_Human_averaged, Maynard_dataset_average, selected_gene_list_multiple) 
    
    #Heatmaps ---- 
    heatmapHeight <- heatmap_height(selected_gene_list_multiple)
    output$He_heatmap <- renderPlot({
      ggplot(data = He_heatmap_data, mapping = aes(x = layer, y = gene_symbol, fill = Z_score)) +
        geom_tile() +
        scale_fill_distiller(palette = "RdYlBu", limits = c(-1,1)*max(abs(He_heatmap_data$Z_score))) +
        scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
        labs(y = "", x = "", title = "He et al Heatmap") +
        #Puts stars on layer marker annotations
        geom_text(aes(label = layer_label), size = 7, vjust = 1)
    }, height = heatmapHeight)
    
    output$Maynard_heatmap <- renderPlot({
      ggplot(data = Maynard_heatmap_data, mapping = aes(x = layer, y = gene_symbol, fill = Z_score)) +
        geom_tile() +
        scale_fill_distiller(palette = "RdYlBu", limits = c(-1,1)*max(abs(Maynard_heatmap_data$Z_score))) +
        scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
        labs(y = "", x = "", title = "Maynard et al Heatmap") +
        #Puts stars on layer marker annotations
        geom_text(aes(label = layer_label), size = 7, vjust = 1)
    }, height = heatmapHeight)
    
    output$heatmap_caption <- renderPrint({
      cat(paste("Fig 1. The heatmaps were created using the data from He et al and Maynard et al studies. The raw 
                RNA-seq data was normalized using z-score normalization. The stars (*) indicate if the
                source paper denoted the gene to be highly enriched within the specific cortical layer."))
    })
    
    output$scRNA_caption <- renderPrint({
      cat(paste("Fig 2. The barplots were created using the scRNA-seq data from the Allen Brain Institute. The 
                data is filtered by the chosen region and the data was normalized by on a per-gene basis. The median of 
                the normalized data was then calculated per gene."))
    })
    
    # AUROC table ----
    output$table <- DT::renderDataTable({
      AUROC_table
    }, escape = FALSE, )
    
    output$AUC_table_caption <- renderPrint({
      cat(paste("Table 1. The layers were ranked in each dataset with respect to the gene expression of the chosen genes, and
                #the AUC score was calculated per layer. P-values were calcuated using the Mann-Whitney U test."))
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
        "It was found that between the He and Maynard datasets, the genes had a mean Pearson correlation value of ",
        multi_gene_cor,
        " (p = ",
        p_value_multiple_gene,
        "), which ranks in the ",
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
    
    # scRNA barplot ----
    scRNA_matrix <- load_scRNA_region(selected_scRNA_region, selected_gene_list_multiple)
    output$scRNA_barplot <- renderPlot({
      ggplot(data = scRNA_matrix, mapping = aes(x = class_label, y = median_exp_value, fill = gene)) +
        geom_col(position = "dodge") +
        facet_wrap( ~ cortical_layer_label, ncol = 2)
      })
  })
}



# Run the app ----
shinyApp(ui = ui, server = server)
