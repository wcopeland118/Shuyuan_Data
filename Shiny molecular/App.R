#Shiny Molecular
# a collaboration with Will to develop a web interface for Simon Griffin
# purpose is to build capacity in Shiny and test to see if it can help produce dendograms and 
# MDS plots of samples lookign at frequencies of bacteria phages analysed using Bray Curtis 
# dissimilarity.

library(shiny)
library(tidyverse)
library(lubridate)
library(plotly)
library(hms)
library(vegan)


# Define UI for ISF molecular bio app
ui <- fluidPage(
  
  # App title with ISF logo goes over the top of the webpage
  fluidRow(
    column(2, tags$img(height = 70, width = 70, src = "ISF_logo.png")),
    column(9, offset = 1, tags$h2("ISF Shuyuan Molecular Biology", style = "color:CornflowerBlue"))  # need official title from Simon
  ),
  
  
  # Sidebar layout with input controls
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      tabsetPanel(type = "tabs",
      
      # Input: Select the file for input
      tabPanel("Input", wellPanel(
        
        fileInput("file", label = tags$h4("Upload phage data", style = "color:CornflowerBlue"),
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        
        
        tags$hr(), #horizontal line
        
        # Input: Checkbox if file has header ----
        checkboxInput("header", "Header", TRUE),
        
        # Input: Select separator ----
        radioButtons("sep", "Separator",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = ","),
        
        # Input: Select quotes ----
        radioButtons("quote", "Quote",
                     choices = c(None = "",
                                 "Double Quote" = '"',
                                 "Single Quote" = "'"),
                     selected = '"'),
        
        # Horizontal line ----
        tags$hr(),
        
        # Input: Select number of rows to display ----
        radioButtons("disp", "Display",
                     choices = c(Head = "head",
                                 All = "all"),
                     selected = "head")



      
      )),
      
      tabPanel("Identify", wellPanel(  
      #input for column number of first phage
      selectInput("phage_col", label = tags$h4("Column # of first phage", style = "color:CornflowerBlue"), 
                  choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                  selected = 3),
      
      tags$hr(),
      fluidRow(column(6, verbatimTextOutput("choice")))
    ))
      

  
    )),
    
    # Main panel for displaying outputs, inc Bray Curtis table, MDS plot etc
    mainPanel(
      
      # Output: Tabset 
      tabsetPanel(type = "tabs",
                  
                  tabPanel("raw data", tableOutput("raw_data")), # show raw data file uploaded by user

                  tabPanel("BC Matrix", tableOutput("BC_table")), # show Bray Curtis matrix
                  
                  tabPanel("NMDS", plotlyOutput("NMDS")), # show plotly
                  
                  tabPanel("Cluster", plotOutput("dend")) # show cluster dendrogram 
      )
      
    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output) {
  
  # Defining & initializing reactiveValues objects
  #raw_data <- reactiveValues(value = matrix(nrow = 1, ncol = 1)) 
  
  # which column number for 1st phage from second tab in sidebar
  phage_col <- reactive({
    phage_col <- input$phage_col  
 
  })
  
  #for the chosen number display in second tab in sidebar
  output$choice <- renderText({phage_col()})
  
  
  #output for main
  
  # first pull in the data from the csv file into reactive expression df
  df <- reactive({
    req(input$file)
    
    df <- read.csv(input$file$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote)
    
  })
  
  # separate phage data from first few columns of df
  x <- reactive(df()[,phage_col():length(colnames(df()))])
  
  # output of rawdata runs when reactive expression df changes
  output$raw_data <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file)
    
    if(input$disp == "head") {
      return(head(df()))
    }
    else {
      return(df())
    }

  })
  
  #output for Bray Curtis table
  output$BC_table <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file)
    
    BC_matrix <- as.matrix(vegdist(x(), method="bray", binary=FALSE, diag=FALSE, upper=T,
                                   na.rm = FALSE, rownames = T))
    
    return(BC_matrix)
  
  
  }, rownames = T)
  
  #output for NMDS plot
  output$NMDS <- renderPlotly({
    
    # shows the MDS plot based on BC matrix
    NMDS <- metaMDS(x(), distance = "bray", k = 2)
    Habitat <- df()[,2]
    
    plot_df <- tibble(x = NMDS$points[,1], y = NMDS$points[,2], 
                      Habitat) %>% 
      rownames_to_column(var = "n")
    
    p <- plot_df %>% 
      ggplot() +
      geom_point(mapping = aes(x = x, y = y, colour = Habitat, 
                               label = n)) + # add label to aes to show in plotly tooltip
      labs(
        x = "Axis 1",
        y = "Axis 2"
      ) +
      theme(plot.title = element_text(color = "steelblue4", size = 14, face = "bold"))
    return(ggplotly(p, tooltip = "n"))
  })
  
  #output for cluster dendogram
  output$dend <- renderPlot({
    dend <- vegdist(x(), method="bray", binary=FALSE, diag=FALSE, upper=T,
                           na.rm = FALSE) %>% 
      hclust %>% 
      as.dendrogram
    
    
      c_plot <- dend %>% highlight_branches %>% plot(main = "Dendrogram")
      
      return(c_plot)
    
  })
  
}

# Create Shiny app ----
shinyApp(ui, server)