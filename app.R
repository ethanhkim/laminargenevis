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
source("string_processing.R") 
source("data_processing.R")
load(here("data", "processed", "He_DS1_Human_averaged.Rdata"), verbose = TRUE)
load(here("data", "processed", "Maynard_dataset_average.Rdata"), verbose = TRUE)
load(here("data", "processed", "Zeng_dataset_updated.Rdata"), verbose = TRUE)
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
          # Input: Selector for which genes to visualize ----
          radioButtons(
            inputId = "selector", label =  "Single or multiple genes?",
            choices = c("Single", "Multiple")
          ),
          # Only show if Input is Single
          conditionalPanel(
            condition = "input.selector == 'Single'",
            selectizeInput(
              inputId = "genelist", label = "Input gene:", choices = NULL,
              selected = NULL, multiple = FALSE, options = NULL)
          ),
          # Only show if Input is Multiple          
          conditionalPanel(
            condition = "input.selector == 'Multiple'",
            textAreaInput(
              inputId = "multiple_genelist", 
              label = "Input your gene list:", 
              placeholder = "GAD1, CCK, GRIN1")
          ),
          # Submit button
          actionButton(inputId = "submit_heatmap", label = "Submit"),
          
          br(),
          br(),
                        
      ),
                          
        mainPanel(
          tabsetPanel(type = "tabs",
            tabPanel("Dataset Overview",
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
              h4("This study assayed the whole genome using high-throughout RNA-seq."),
              br(),
              h3(a("Maynard et al. (2020)*", href = "https://www.biorxiv.org/content/10.1101/2020.02.28.969931v1", target = "_blank")),
              h4("*This study is currently a pre-print; it assayed the whole genome through the 10X Genomics Visium Platform.")
            ),
            tabPanel("Gene Visualization",
             # Only show if Input from checkboxGroupInput is Multiple - only shows rendered heatmaps
              br(),
              conditionalPanel(
                h3("Layer-specific Heatmaps"),
                condition = "input.selector == 'Multiple'",
                plotlyOutput("He_heatmap"),
                br(),
                br(),
                plotlyOutput("Maynard_heatmap"),
                br(),
                h4(verbatimTextOutput("summary_multiple")),
                br(),
                dataTableOutput("table")
              ),
              # Only show if input from checkboxGroupInput is Single - only show layer-specific barplot
              conditionalPanel(
                h3("Layer-specific gene expression"),
                condition = "input.selector == 'Single'",
                plotlyOutput("Barplot"),
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

  observeEvent(input$submit_heatmap, {
    
    # List of selected gene(s)
    selected_gene_list_single <- isolate(process_gene_input(input$genelist))
    selected_gene_list_multiple <- isolate(process_gene_input(input$multiple_genelist))
    
    # Process dataset to correct format for heatmap and barplot
    He_heatmap_data <- process_heatmap_function(He_DS1_Human_averaged, selected_gene_list_multiple)
    Maynard_heatmap_data <- process_heatmap_function(Maynard_dataset_average, selected_gene_list_multiple)
    Barplot_data <- process_barplot_data(selected_gene_list_single, He_DS1_Human_averaged, Maynard_dataset_average)
    
    ## Single gene input:
    if (input$selector == "Single") {
      output$Barplot <- renderPlotly({
        p <- ggplot(data = Barplot_data, aes(x = Layer, y = Z_score, fill = Dataset, group = Dataset)) +
          geom_bar(stat = "identity", position = "dodge", width = 0.75) + theme(axis.text.x = element_text(angle = 45)) +
          geom_text(aes(label = layer_label, group = Dataset, 
                        vjust = ifelse(Z_score >= -0.1, 0, 2.5)), 
                    position = position_dodge(width = 0.75)) +
          geom_hline(yintercept = 1.681020) +
          geom_hline(yintercept = -1.506113) +
          labs(caption = "Barplot of z-scored gene expression levels.")
        p <- ggplotly(p, height = 400)
        p
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
      
      
      output$summary_single <- renderPrint({
        cat(paste0(
          selected_gene_list_single,
          " has a Pearson correlation value of ",
          single_gene_cor,
          " between the He and Maynard datasets, ranking in the ",
          single_gene_quantile,
          "th quantile (p = ",
          p_value_single_gene,
          ").\n\n",
          selected_gene_list_single,
          " was found to be:\n\n",
          if (sum(selected_gene_list_single %in% unique(He_DS1_Human_averaged$gene_symbol)) == 0) {
            "Not assayed by He et al.,\n"
          } else
            "Assayed by He et al.,\n",
          if (sum(selected_gene_list_single %in% unique(Maynard_dataset_average$gene_symbol)) == 0) {
            "Not assayed by Maynard et al.,\n"
          } else
            "Assayed by Maynard et al.,\n",
          if (sum(selected_gene_list_single %in% unique(Zeng_dataset_updated$gene_symbol)) == 0) {
            "Not assayed by Zeng et al.\n"
          } else {
            paste0(
              "Assayed by Zeng et al.\n\nIn the Zeng dataset, ",
              selected_gene_list_single,
              " was found to:\n\n",
              if (length(unlist(layer_specific_gene_list_single $`Layer 1`)) == 0 &
                  length(unlist(layer_specific_gene_list_single $`Layer 2`)) == 0 &
                  length(unlist(layer_specific_gene_list_single $`Layer 3`)) == 0 &
                  length(unlist(layer_specific_gene_list_single $`Layer 4`)) == 0 & 
                  length(unlist(layer_specific_gene_list_single $`Layer 5`)) == 0 &
                  length(unlist(layer_specific_gene_list_single $`Layer 6`)) == 0 &
                  length(unlist(layer_specific_gene_list_single $White_matter)) == 0) {
                "Not mark any layer."
              } else {
                paste0(
                  if (length(unlist(layer_specific_gene_list_single$`Layer 1`)) == 0) {
                    ""
                  } else "Mark layer 1.\n",
            
                  if (length(unlist(layer_specific_gene_list_single$`Layer 2`)) == 0) {
                    ""
                  } else "Mark layer 2.\n",
                  
                  if (length(unlist(layer_specific_gene_list_single$`Layer 3`)) == 0) {
                    ""
                  } else "Mark layer 3.\n",
                  
                  if (length(unlist(layer_specific_gene_list_single$`Layer 4`)) == 0) {
                    ""
                  } else "Mark layer 4.\n",
                  
                  if (length(unlist(layer_specific_gene_list_single$`Layer 5`)) == 0) {
                    ""
                  } else "Mark layer 5.\n",
                  
                  if (length(unlist(layer_specific_gene_list_single$`Layer 6`)) == 0) {
                    ""
                  } else "Mark layer 6.\n",
                  
                  if (length(unlist(layer_specific_gene_list_single$`WM`)) == 0) {
                    ""
                  } else "Mark white matter.",
                  sep = "\n"
                )
              }
            )
          }
        ))
      })
      
    } else {
      ## Multiple gene Input
      
      layer_marker_table_multiple <- layer_marker_table %>%
        dplyr::filter(gene_symbol %in% selected_gene_list_multiple)
      layer_specific_gene_list_multiple <- separate_layers(layer_marker_table_multiple, selected_gene_list_multiple)
      names(layer_specific_gene_list_multiple) <- c("Layer 1", "Layer 2", "Layer 3", "Layer 4",
                                                    "Layer 5", "Layer 6", "White_matter")
      
      
      # Generate correlation value for multiple gene and the quantile that value belongs in

      multi_gene_cor <- multi_gene_correlation(selected_gene_list_multiple, He_DS1_Human_averaged, Maynard_dataset_average)
      multi_gene_quantile <- quantile_distribution(He_Maynard_cor_diagonal, multi_gene_cor)
      p_value_multiple_gene <- wilcoxtest(selected_gene_list_multiple, He_DS1_Human_averaged, Maynard_dataset_average, He_Maynard_cor_diagonal)
      
      ## Code for AUROC analysis adapted from Derek Howard & Leon French
      
      He_df <- rank_dataset(He_DS1_Human_averaged)
      Maynard_df <- rank_dataset(Maynard_dataset_average)
      
      He_indices <- return_indices(He_df, selected_gene_list_multiple)
      Maynard_indices <- return_indices(Maynard_df, selected_gene_list_multiple)
      
      He_df %<>% select(-gene_symbol)
      Maynard_df %<>% select(-gene_symbol)
      
      AUROC_He <- map_df(He_df, auroc_analytic, He_indices)
      wilcox_AUROC_He <- map_df(He_df, apply_MWU, He_indices)
      AUROC_Maynard <- map_df(Maynard_df, auroc_analytic, Maynard_indices)
      wilcox_AUROC_Maynard <- map_df(Maynard_df, apply_MWU, Maynard_indices)
      
      AUROC_table <- bind_cols(gather(AUROC_He, key = Layers, value = AUROC_He),
                               gather(wilcox_AUROC_He, value = pValue),
                               gather(AUROC_Maynard, key = Layers, value = AUROC_Maynard),
                               gather(wilcox_AUROC_Maynard, value = pValue)) %>%
        rename(Layers = Layers...1) %>%
        rename(key = key...3) %>%
        rename(pValue_He = pValue...4) %>%
        rename(pValue_Maynard = pValue...8) %>%
        select(Layers, AUROC_He, pValue_He, AUROC_Maynard, pValue_Maynard)
      
      AUROC_table %<>% mutate(pValue_He = signif(pValue_He, digits = 3),
                              pValue_Maynard = signif(pValue_Maynard, digits = 3),
                              AUROC_He = signif(AUROC_He, digits = 3),
                              AUROC_Maynard = signif(AUROC_Maynard, digits = 3),
                              adjusted_P_He = signif(p.adjust(pValue_He), digits = 3),
                              adjusted_P_Maynard = signif(p.adjust(pValue_Maynard), digits = 3)) %<>%
        select(Layers, AUROC_He, pValue_He, adjusted_P_He, AUROC_Maynard, pValue_Maynard, adjusted_P_Maynard)
      
      # Heatmaps
      
      output$He_heatmap <- renderPlotly({
        p <- ggplot(data = He_heatmap_data, mapping = aes(x = layer, y = gene_symbol, fill = Z_score)) +
          geom_tile() +
          scale_fill_distiller(palette = "RdYlBu") + 
          scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
          labs(y = "", x = "", title = "He et al Heatmap") +
          labs(caption = "(based on data from He et al, 2017)") +
          geom_text(aes(label = layer_label))
        p <- ggplotly(p)      
        p
      })
      
      output$Maynard_heatmap <- renderPlotly({
        p <- ggplot(data = Maynard_heatmap_data, mapping = aes(x = layer, y = gene_symbol, fill = Z_score)) +
          geom_tile() +
          scale_fill_distiller(palette = "RdYlBu") +
          scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
          labs(y = "", x = "", title = "Maynard et al Heatmap",
               caption = "(based on data from Maynard et al, 2020)") +
          geom_text(aes(label = layer_label))
        p <- ggplotly(p)      
        p
      })
      
      # AUROC table
      
      output$table <- renderDataTable({
        AUROC_table
      }, escape = FALSE)
      
      # Summary textbox
      
      output$summary_multiple <- renderPrint({
        #count of intersection of submitted genes with total gene list
        cat(paste0(
          "Of the ",
          length(selected_gene_list_multiple),
          " input genes:\n\n",
          "The genes had a mean Pearson correlation value of ",
          multi_gene_cor,
          " (p = ",
          p_value_multiple_gene,
          "), which ranks in the ",
          multi_gene_quantile,
          "th quantile.\n\n",
          sum(selected_gene_list_multiple %in% unique(He_DS1_Human_averaged$gene_symbol)),
          " were assayed by He et al.,\n",
          sum(selected_gene_list_multiple %in% unique(Maynard_dataset_average$gene_symbol)),
          " were assayed by Maynard et al., \nand ",
          
          if (sum(selected_gene_list_multiple %in% unique(Zeng_dataset_updated$gene_symbol)) == 0) {
            "0 genes were assayed by Zeng et al."
          } else {
            paste0(
              sum(selected_gene_list_multiple %in% unique(Zeng_dataset_updated$gene_symbol)),
              " were assayed by Zeng et al.\n\nIn the Zeng dataset:\n\n",
              
              if (length(unlist(layer_specific_gene_list_multiple $`Layer 1`)) == 0 &
                  length(unlist(layer_specific_gene_list_multiple $`Layer 2`)) == 0 &
                  length(unlist(layer_specific_gene_list_multiple $`Layer 3`)) == 0 &
                  length(unlist(layer_specific_gene_list_multiple $`Layer 4`)) == 0 & 
                  length(unlist(layer_specific_gene_list_multiple $`Layer 5`)) == 0 &
                  length(unlist(layer_specific_gene_list_multiple $`Layer 6`)) == 0 &
                  length(unlist(layer_specific_gene_list_multiple $White_matter)) == 0) {
                "No genes were found to mark any layer."
              } else {
                paste0( if (length(unlist(layer_specific_gene_list_multiple$`Layer 1`)) == 0) {
                  ""
                } else paste0(length(unlist(layer_specific_gene_list_multiple $`Layer 1`))," marked layer 1 (", 
                              paste(layer_specific_gene_list_multiple $`Layer 1`, collapse = ", "), ").\n"),
                
                if (length(unlist(layer_specific_gene_list_multiple $`Layer 2`)) == 0) {
                  ""
                } else paste0(length(unlist(layer_specific_gene_list_multiple $`Layer 2`)), " marked layer 2 (", 
                              paste(layer_specific_gene_list_multiple $`Layer 2`, collapse = ", "), ").\n"),
                
                if (length(unlist(layer_specific_gene_list_multiple $`Layer 3`)) == 0) {
                  ""
                } else paste0(length(unlist(layer_specific_gene_list_multiple $`Layer 3`)), " marked layer 3 (",
                              paste(layer_specific_gene_list_multiple $`Layer 3`, collapse = ", "), ").\n"),
                
                if (length(unlist(layer_specific_gene_list_multiple $`Layer 4`)) == 0) {
                  ""
                } else paste0(length(unlist(layer_specific_gene_list_multiple $`Layer 4`)), " marked layer 4 (",
                              paste(layer_specific_gene_list_multiple $`Layer 4`, collapse = ", "), ").\n"),
                
                if (length(unlist(layer_specific_gene_list_multiple $`Layer 5`)) == 0) {
                  ""
                } else paste0(length(unlist(layer_specific_gene_list_multiple $`Layer 5`)), " marked layer 5 (",
                              paste(layer_specific_gene_list_multiple $`Layer 5`, collapse = ", "), ").\n"),
                
                if (length(unlist(layer_specific_gene_list_multiple $`Layer 6`)) == 0) {
                  ""
                } else paste0(length(unlist(layer_specific_gene_list_multiple $`Layer 6`)), " marked layer 6 (",
                              paste(layer_specific_gene_list_multiple $`Layer 6`, collapse = ", "), ").\n"),
                
                if (length(unlist(layer_specific_gene_list_multiple $White_matter)) == 0) {
                  ""
                } else paste0(length(unlist(layer_specific_gene_list_multiple $White_matter)), " marked white matter (",
                              paste(layer_specific_gene_list_multiple $White_matter, collapse = ", "), ")."),
                sep = "\n")
              }
            )
          }
        ))
      })
    }
  })
}


# Run the app ----
shinyApp(ui = ui, server = server)
