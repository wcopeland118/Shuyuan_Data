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
    column(9, offset = 1, tags$h2("ISF Shuyuan Molecular Biology"))  # need official title from Simon
  ),
  
  
  # Sidebar layout with input controls
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select the time series
      wellPanel(
        
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



      
      ),
    wellPanel(  
      #input for column number of first phage
      selectInput("phage_col", label = tags$h4("Column # of first phage", style = "color:CornflowerBlue"), 
                  choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                  selected = 1),
      
      tags$hr(),
      fluidRow(column(3, verbatimTextOutput("choice")))
    )
      

  
    ),
    
    # Main panel for displaying outputs, inc Bray Curtis table, MDS plot etc
    mainPanel(
      
      # Output: Tabset 
      tabsetPanel(type = "tabs",
                  
                  tabPanel("raw data", tableOutput("raw_data")),

                  tabPanel("BC Matrix", tableOutput("BC_table")) # show Bray Curtis matrix
      )
      
    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output) {
  
  # Defining & initializing reactiveValues objects
  #raw_data <- reactiveValues(value = matrix(nrow = 1, ncol = 1)) 
  
  # which column number for 1st phage 
  output$choice <- renderText({input$phage_col})
  
  #output for main
  
  # first pull in the data from the csv file into reactive expression df
  df <- reactive({
    req(input$file)
    
    df <- read.csv(input$file$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote)
    
  })
  
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
  
  output$BC_table <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file)
    
    NumCols <- length(colnames(df()))
    x <- df()[,3:11] 
    BC_matrix <- as.matrix(vegdist(x, method="bray", binary=FALSE, diag=FALSE, upper=T,
                                   na.rm = FALSE, rownames = T))
    
    return(BC_matrix)
  
  
  }, rownames = T)
}

# Create Shiny app ----
shinyApp(ui, server)