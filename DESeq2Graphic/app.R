#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

#### UI ####
ui <- shinyUI(fluidPage(
  
  #### FILE UPLOAD ####
  titlePanel("Uploading Files"),
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Choose CSV File',
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv')),
      tags$hr(),
      checkboxInput('header', 'Header', TRUE),
      radioButtons('sep', 'Separator',
                   c(Comma=',',
                     Semicolon=';',
                     Tab='\t'),
                   ','),
      radioButtons('quote', 'Quote',
                   c(None='',
                     'Double Quote'='"',
                     'Single Quote'="'"),
                   '"')
    ),
    
    
    #### MAIN PANEL ####
    mainPanel(
      tabsetPanel(
        tabPanel('Content',tableOutput('contents')),
        tabPanel('PCA', tableOutput('pca'))
      )
    )
  )
))




#### Server ####
server <- shinyServer(function(input, output) {
  data <- reactive({
    inFile <- input$file1
    if(is.null(inFile)){return(NULL)}
    data <- read.csv(inFile$datapath, header=input$header, sep=input$sep, 
                     quote=input$quote)
  })  
  
  
  output$contents <- renderTable({
      
      # input$file1 will be NULL initially. After the user selects
      # and uploads a file, it will be a data frame with 'name',
      # 'size', 'type', and 'datapath' columns. The 'datapath'
      # column will contain the local filenames where the data can
      # be found.
      
      if (is.null(data())){return(NULL)}
      head(data())
    })
    
    output$sum <- renderTable({
      if (is.null(data())){return(NULL)}
      summary(data())
    })
})

# Run the application 
shinyApp(ui = ui, server = server)

