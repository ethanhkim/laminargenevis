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

# Define UI ----
ui <- fluidPage(
  navbarPage(title = "Gene Expression Profile Comparison",
             
      tabPanel(title = "Heatmaps",
          sidebarLayout(
              sidebarPanel(
                
                # Input: Selector for which genes to visualize ----
                selectizeInput(
                  inputId = "genelist", label = "Gene", choices = NULL,
                  selected = NULL, multiple = TRUE, options = NULL
                  ),
                
                # Submit button
                actionButton(inputId = "submit", label = "Submit"),
                
                # Spaces
                br(),
                br(),
                br(),
                            
                ),
              
              mainPanel(h3("Layer-specific Heat Map"),
                  br(),
                  br(),
                  br(),
                  div(id = "main",
                      # Output: Heatmap of He et al gene expression
                      plotlyOutput("He_heatmap")),
                  
                  br(),
                  br(),
                  br(),

                  div(id = "second",
                      # Output: Heatmap of Maynard et al gene expression
                      plotlyOutput("Maynard_heatmap"))
                )
            )
        )
    )
)

# Define server logic ----
server <- function(input, output, session) {
  
  common_genelist <- intersect(He_DS1_Human_averaged$gene_symbol, Maynard_dataset_average$gene_symbol) %>%
    sort()
  
  updateSelectizeInput(session, inputId = "genelist", 
                       choices = common_genelist, server = TRUE)
  
  observeEvent(input$submit, {
    
    selected_gene_list <- isolate(process_gene_input(input$genelist))
    
    He_heatmap_data <- He_DS1_Human_averaged %>%
      dplyr::filter(gene_symbol %in% selected_gene_list) %>%
      distinct(gene_symbol, .keep_all = TRUE)  %>%
      arrange(gene_symbol) %>%
      column_to_rownames(var = "gene_symbol") %>%
      t() %>%
      melt() %>%
      as_tibble() %>%
      dplyr::rename(Cuts = "Var1", Gene_symbol = "Var2", Z_score = "value")
    
    output$He_heatmap <- renderPlotly({
      p <- ggplot(data = He_heatmap_data, mapping = aes(x = Cuts, y = Gene_symbol, fill = Z_score)) +
        geom_tile() +
        scale_fill_distiller(palette = "RdYlBu") +
        labs(y = "", x = "", title = "He et al Heatmap") +
        labs(caption = "(based on data from He et al, 2017")
      p <- ggplotly(p)      
      p
    })
    
    Maynard_heatmap_data <- Maynard_dataset_average %>%
      dplyr::filter(gene_symbol %in% selected_gene_list) %>%
      distinct(gene_symbol, .keep_all = TRUE)  %>%
      arrange(gene_symbol) %>%
      column_to_rownames(var = "gene_symbol") %>%
      t() %>%
      melt() %>%
      as_tibble() %>%
      dplyr::rename(Cuts = "Var1", Gene_symbol = "Var2", Z_score = "value")
    
    output$Maynard_heatmap <- renderPlotly({
      p <- ggplot(data = Maynard_heatmap_data, mapping = aes(x = Cuts, y = Gene_symbol, fill = Z_score)) +
        geom_tile() +
        scale_fill_distiller(palette = "RdYlBu") +
        labs(y = "", x = "", title = "Maynard et al Heatmap",
             caption = "(based on data from Maynard et al, 2020)")
      p <- ggplotly(p)      
      p
    })
  })
}
  

# Run the app ----
shinyApp(ui = ui, server = server)

