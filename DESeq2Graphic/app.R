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
  titlePanel("DESeq2 Analysis"),
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
                   '\t'),
      radioButtons('quote', 'Quote',
                   c(None='',
                     'Double Quote'='"',
                     'Single Quote'="'"),
                   '"'),
      tags$hr(),
      textInput('colIDs', label = 'Column Description', placeholder = 'eg.: FT,FT,E,E - NO spaces, comma seperated'),
      textInput('norm', 'Normalization (Exp/Crtl)', placeholder = 'eg.: E/FT'),
      numericInput('signif', 'Significance level (p_adjusted)', value = 0.1),
      radioButtons('althyp', label = 'Alternative Hypothesis',
                   c(Greater = 'greater',
                     Less = 'less',
                     Both = ''),
                   ''
                   ),
      submitButton(text = 'Apply')
    ),
    
    
    #### MAIN PANEL ####
    mainPanel(
      tabsetPanel(
        tabPanel('Content',tableOutput('contents')),
        tabPanel('PCA', 
                 h4('PCA-Plot', align='center'),
                 br(),
                 plotOutput('pca')
         ),
        tabPanel('DESeq',
                 h4('MA-Plot', align='center'),
                 plotOutput('plotma')
         ),
        tabPanel('Result Table',
                 dataTableOutput('resultTable'))
      )
    )
  )
))




#### Server ####
library(DESeq2)
server <- shinyServer(function(input, output) {
  data <- reactive({
    inFile <- input$file1
    if(is.null(inFile)){return(NULL)}
    data <- as.matrix(read.csv(inFile$datapath, header=input$header, sep=input$sep, 
                     quote=input$quote, row.names = 1))
  })
  
  output$contents <- renderTable({
      
      # input$file1 will be NULL initially. After the user selects
      # and uploads a file, it will be a data frame with 'name',
      # 'size', 'type', and 'datapath' columns. The 'datapath'
      # column will contain the local filenames where the data can
      # be found.
      
      if (is.null(data())){return(NULL)}else{head(data(), 100)}
    })
    
  colData <- reactive({
    names <- input$colIDs
    names <- strsplit(names, split = ',')[[1]]
    colData <- data.frame('Samples' = colnames(data()), 'Conditions' = names)
  })
  
  dds <- reactive({
    dds <- DESeqDataSetFromMatrix(countData = data(), colData = colData(), design=~Conditions)
    dds <- DESeq(dds)
  })
  
  res <- reactive({
    res <- dds()
    padj <- input$signif
    norm <- strsplit(input$norm, split = '/')[[1]]
    norm <- c('Conditions', norm)
    if(input$althyp == ''){
      res <- results(res, contrast = norm, alpha = padj)
    }else{
      res <- results(res, contrast = norm, alpha = padj, altHypothesis = input$althyp)
    }
  })
  
  output$pca <- renderPlot({
    if(is.null(data()) && is.null(colData())){
      NULL
    }else{
      dds.df <- dds()
      rld <- rlog(dds.df)
      plotPCA(rld, intgroup = 'Samples')
    }
      
  })
    
  output$plotma <- renderPlot({
    if(length(colData())>1){
      res.df <- res()
      plotMA(res.df, main='MA-Plot')
    }else{
      NULL
    }
  })
  
  output$resultTable <- renderDataTable({
    res.df <- as.data.frame(res())
    res.df <- data.frame('GeneID'=row.names(res.df), res.df)
    res.df <- subset(res.df, res.df$padj < input$signif)
    res.df
  })
})

# Run the application 
shinyApp(ui = ui, server = server)

