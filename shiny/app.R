library(rsconnect)
library(tidyverse)
library(plotly)
library(shiny)
library(shinydashboard)

#shinylive::export(appdir = "~/Library/CloudStorage/Box-Box/Robbi_PhD/TWIST_Tests/cost_benefit/shiny", destdir = "~/Library/CloudStorage/Box-Box/Robbi_PhD/TWIST_Tests/cost_benefit/docs")
#httpuv::runStaticServer("~/Library/CloudStorage/Box-Box/Robbi_PhD/TWIST_Tests/cost_benefit/docs/", port=8008)

#setwd("~/Library/CloudStorage/Box-Box/Robbi_PhD/TWIST_Tests/cost_benefit/shiny")

# Define UI
ui <- fluidPage(
  titlePanel("Cost Per SNP"),
  sidebarLayout(
    sidebarPanel(
      numericInput("cost_enrichment", "Cost of one Enrichment reaction:", value = 200.42),
      numericInput("total_sequencing_cost", "Total Cost of Sequencing Run:", value = 7284.38),
      numericInput("total_reads_sequenced", "Total Number of Reads Sequenced:", value = 1022807058)
    ),
    mainPanel(
     # plotly::plotlyOutput("scatter_plot"), plotly::plotlyOutput("plot9b")
      plotlyOutput("scatter_plot"), plotlyOutput("plot9b")
    )
  )
)


# Define server
server <- function(input, output) {
  #Define colours
  colours <- c("SG"="#ff0000","TW1"="#01cc54", "TW2"="#066601","1"="#fe7a00", "2"="#019dfe", "4"="#001bcc")

  # Plot 1
  output$scatter_plot <- renderPlotly({
    # Read data from CSV
    df <- read_tsv("TableS1.csv")
    
    # Subset data
    df_1 <- dplyr::select(df, Deidentified_IDs, n_rounds, method, n_pool, MappableBioEndoDNA_SG, SNPs_Covered, Nr._Sequenced_Read_Pairs)
    
    # Calculate other values based on user input
    CostOfEnrichment <- input$cost_enrichment
    TotalCostOfSequencingRun <- input$total_sequencing_cost
    TotalNumberOfReadsSequenced <- input$total_reads_sequenced
    
    df_2 <- df_1 %>% 
      mutate(CostPerSNP =  (((CostOfEnrichment*as.numeric(n_rounds))/as.numeric(n_pool)) + (TotalCostOfSequencingRun/TotalNumberOfReadsSequenced)*Nr._Sequenced_Read_Pairs )/as.numeric(SNPs_Covered) ) %>%
      mutate(CostPerSample =  ((CostOfEnrichment*as.numeric(n_rounds))/as.numeric(n_pool)) + (TotalCostOfSequencingRun/TotalNumberOfReadsSequenced)*Nr._Sequenced_Read_Pairs ) 
    
    # Create plot 1
    ggplot(df_2, aes(x = MappableBioEndoDNA_SG, y = CostPerSNP, fill = method, colour = method)) +
      geom_smooth(formula = y ~ log(x), method = "lm", alpha = 0.2) +
      geom_point(size = 4, alpha = 1, shape = 21, colour = "black") +
      scale_colour_manual(values = colours) + scale_fill_manual(values = colours) + 
      scale_shape_manual(values = as.character(letters), guide = "none") +
      scale_y_log10(labels = comma_format(), n.breaks = 12) +
      theme_bw() +
      theme(legend.position = c(0.8, 0.8), text = element_text(size = 12)) +
      labs(x = "Mappable Endo%, Screening SG", y = "Cost per SNP", title = "")
  })
  
  # Plot 2
  output$plot9b <- renderPlotly({
    # Calculate other values based on user input
    CostOfEnrichment <- input$cost_enrichment
    TotalCostOfSequencingRun <- input$total_sequencing_cost
    TotalNumberOfReadsSequenced <- input$total_reads_sequenced
    # Read data from CSV
    df <- read_tsv("TableS1.csv")
    
    # Subset data
    df_1 <- dplyr::select(df, Deidentified_IDs, n_rounds, method, n_pool, MappableBioEndoDNA_SG, SNPs_Covered, Nr._Sequenced_Read_Pairs)
    
    df_2 <- df_1 %>% 
      mutate(CostPerSNP =  (((CostOfEnrichment*as.numeric(n_rounds))/as.numeric(n_pool)) + (TotalCostOfSequencingRun/TotalNumberOfReadsSequenced)*Nr._Sequenced_Read_Pairs )/as.numeric(SNPs_Covered) ) %>%
      mutate(CostPerSample =  ((CostOfEnrichment*as.numeric(n_rounds))/as.numeric(n_pool)) + (TotalCostOfSequencingRun/TotalNumberOfReadsSequenced)*Nr._Sequenced_Read_Pairs ) 
    
    # Subset data
    df_3 <- dplyr::select(df_2, Deidentified_IDs, CostPerSNP, method, MappableBioEndoDNA_SG) 
df_3 <- pivot_wider(df_3, names_from = "method", values_from = "CostPerSNP", names_prefix = "costperSNP_")
df_3 <- pivot_longer(df_3, values_to = "cost_TW", names_to = "method", names_prefix = "costperSNP_", cols = c("costperSNP_TW1", "costperSNP_TW2")) %>%
  mutate(CostperSNP_TW = cost_TW / costperSNP_SG) %>%
  mutate(CostperSNP_TW_INV = 1 / (cost_TW / costperSNP_SG))

# Create plot 2
ggplot(df_3, aes(x = MappableBioEndoDNA_SG, y = CostperSNP_TW_INV, fill = method, colour = method)) +
  geom_hline(yintercept = 0) +
  geom_smooth(data = subset(df_3, method == "TW2"), formula = y ~ log(x), method = 'lm', alpha = 0.2) +
  geom_smooth(data = subset(df_3, method == "TW1"), formula = y ~ x, method = "lm", alpha = 0.2) +
  geom_point(size = 4, colour="black", alpha = 1, shape=21) +
  scale_colour_manual(values = colours) + scale_fill_manual(values = colours) + 
  scale_shape_manual(values = as.character(letters), guide = "none") +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "none") +
  labs(x = "Mappable Endo%, Screening SG", y = "Fold cost saving per SNP compared to SG", title = "")
  })
}

# Run the application
shinyApp(ui = ui, server = server)
