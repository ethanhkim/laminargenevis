## App script ##

library(shiny)
library(ggplot2)
library(magrittr)
library(plotly)
library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
library(reshape2)
source("./string_processing.R") 
load('./data/processed/He_DS1_Human_averaged.Rdata', verbose = TRUE)
load('./data/processed/Maynard_dataset_average.Rdata', verbose = TRUE)
load('./data/processed/Compared_Layer_markers.Rdata', verbose = TRUE)

# Define UI ----
ui <- fluidPage(
  navbarPage(title = "Gene Expression Profile Comparison",
             
      # Visualize gene expression across layers through heatmap or barplot    
      tabPanel(title = "Heatmaps",
          sidebarLayout(
            sidebarPanel(
            # Input: Selector for which genes to visualize ----
              radioButtons(
                inputId = "selector", label =  "Single or multiple gnees?",
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
                          
            ),
                        
              mainPanel(
                # Only show if Input from checkboxGroupInput is Multiple - only shows rendered heatmaps
                conditionalPanel(
                  h3("Layer-specific Heatmaps"),
                  condition = "input.selector == 'Multiple'",
                  plotlyOutput("He_heatmap"),
                  br(),
                  br(),
                  plotlyOutput("Maynard_heatmap")
                ),
                # Only show if input from checkboxGroupInput is Single - only show layer-specific barplot
                conditionalPanel(
                  h3("Layer-specific gene expression"),
                  condition = "input.selector == 'Single'",
                  plotlyOutput("Barplot")
                )
                          
              ),
            ),
        ),
        
        # Show if layer markers are consistent across various datasets
        tabPanel(title = "Layer Markers",
            sidebarLayout(
              sidebarPanel(
                # Input: Selector for which genes to visualize ----
                textAreaInput(
                  inputId = "layer_marker_genelist", 
                  label = "Input your gene list:", 
                  placeholder = "GAD1, CCK, GRIN1"
                ),
                          
                br(),
                br(),
                          
                # Submit button
                actionButton(inputId = "submit_layer_marker", label = "Submit")
              ),
                        
              mainPanel(h3("Layer Enrichment"),
                br(),
                br(),
    
                # Output: List of layer markers
                dataTableOutput("Layer_marker")
              )
            )
        )
  )
)

# Define server logic ----
server <- function(input, output, session) {
  
  common_genelist <- intersect(He_DS1_Human_averaged$gene_symbol, Maynard_dataset_average$gene_symbol) %>%
    sort()
  
  layer_marker_table <- Compared_Layer_Markers %>%
    replace(is.na(.), "N/A")
  
  updateSelectizeInput(session, inputId = "genelist", 
                       choices = common_genelist, server = TRUE)
  
  observeEvent(input$submit_heatmap, {
    
    selected_gene_list_barplot <- isolate(process_gene_input(input$genelist))
    selected_gene_list_heatmap <- isolate(process_gene_input(input$multiple_genelist))
    
    He_heatmap_data <- He_DS1_Human_averaged %>%
      dplyr::filter(gene_symbol %in% selected_gene_list_heatmap) %>%
      distinct(gene_symbol, .keep_all = TRUE)  %>%
      arrange(gene_symbol) %>%
      column_to_rownames(var = "gene_symbol") %>%
      t() %>%
      melt() %>%
      as_tibble() %>%
      dplyr::rename(Cuts = "Var1", Gene_symbol = "Var2", Z_score = "value")
    
    Maynard_heatmap_data <- Maynard_dataset_average %>%
      dplyr::filter(gene_symbol %in% selected_gene_list_heatmap) %>%
      distinct(gene_symbol, .keep_all = TRUE)  %>%
      arrange(gene_symbol) %>%
      column_to_rownames(var = "gene_symbol") %>%
      t() %>%
      melt() %>%
      as_tibble() %>%
      dplyr::rename(Cuts = "Var1", Gene_symbol = "Var2", Z_score = "value")
    
    He_barplot_data <- He_DS1_Human_averaged %>%
      dplyr::filter(gene_symbol %in% selected_gene_list_barplot)
    Maynard_barplot_data <- Maynard_dataset_average %>%
      dplyr::filter(gene_symbol %in% selected_gene_list_barplot)
    
    Barplot_data <- He_barplot_data %>%
      rename(He_Layer_1 = "Layer_1", He_Layer_2 = "Layer_2", He_Layer_3 = "Layer_3",
             He_Layer_6 = "Layer_6", He_Layer_5 = "Layer_5", He_Layer_4 = "Layer_4",
             He_WM = "WM"
      ) %>%
      add_column(Maynard_Layer_1 = Maynard_barplot_data$Layer_1, Maynard_Layer_2 = Maynard_barplot_data$Layer_2,
                 Maynard_Layer_3 = Maynard_barplot_data$Layer_3, Maynard_Layer_4 = Maynard_barplot_data$Layer_4,
                 Maynard_Layer_5 = Maynard_barplot_data$Layer_5, Maynard_Layer_6 = Maynard_barplot_data$Layer_6,
                 Maynard_WM = Maynard_barplot_data$WM
      ) %>%
      melt(id.vars = "gene_symbol") %>%
      dplyr::rename(Z_score = "value", Layers = "variable") %>%
      mutate(Dataset = ifelse(test = str_detect(Layers, "He"),
                              yes = "He",
                              no = "Maynard")) %>%
      mutate(Layer = case_when(
        str_detect(Layers, "1") ~ "1", str_detect(Layers, "2") ~ "2", str_detect(Layers, "3") ~ "3",
        str_detect(Layers, "4") ~ "4", str_detect(Layers, "5") ~ "5", str_detect(Layers, "6") ~ "6",
        str_detect(Layers, "WM") ~ "WM"
      ))
    
    if (input$selector == "Single") {
      
      output$Barplot <- renderPlotly({
        p <- ggplot(data = Barplot_data, aes(x = Layer, y = Z_score, fill = Dataset)) +
          geom_bar(stat = "identity", position = position_dodge(0.5), width = 0.5) +
          theme(axis.text.x = element_text(angle = 45))
        
        p <- ggplotly(p)
        p
      }) 
    } else {
      output$He_heatmap <- renderPlotly({
        p <- ggplot(data = He_heatmap_data, mapping = aes(x = Cuts, y = Gene_symbol, fill = Z_score)) +
          geom_tile() +
          scale_fill_distiller(palette = "RdYlBu") +
          labs(y = "", x = "", title = "He et al Heatmap") +
          labs(caption = "(based on data from He et al, 2017")
        p <- ggplotly(p)      
        p
      })
      
      
      
      output$Maynard_heatmap <- renderPlotly({
        p <- ggplot(data = Maynard_heatmap_data, mapping = aes(x = Cuts, y = Gene_symbol, fill = Z_score)) +
          geom_tile() +
          scale_fill_distiller(palette = "RdYlBu") +
          labs(y = "", x = "", title = "Maynard et al Heatmap",
               caption = "(based on data from Maynard et al, 2020)")
        p <- ggplotly(p)      
        p
      })
    }
    
  })
  
  # Take input for Layer Enrichment subpage and output data table with inputted genes
  observeEvent(input$submit_layer_marker, {
    
    selected_genelist_layer_marker <- isolate(process_gene_input(input$layer_marker_genelist))
    
    table <- layer_marker_table %>%
      dplyr::filter(gene_symbol %in% selected_genelist_layer_marker)
    
    output$Layer_marker <- renderDataTable({
      table
    }, escape = FALSE)
    
  })
}


# Run the app ----
shinyApp(ui = ui, server = server)
