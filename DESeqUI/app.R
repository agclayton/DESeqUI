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
      fileInput('file1', 'Choose Read-Count File',
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
      submitButton(text = 'Apply/Update')
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
        tabPanel('Classes',
                 plotOutput('classes'),
                 br(),
                 br(),
                 p('The gene-class ZUnknown is remove here, since it often represents the vast majority and therefore makes a evaluation of other classes hard')),
        tabPanel('Result Table',
                 downloadButton('dload', label='Download'),
                 dataTableOutput('resultTable')
                 )
      )
    )
  )
))




#### Server ####
require(DESeq2)
require(ggplot2)
server <- shinyServer(function(input, output) {
  gene.classes <- read.table('GeneClasses.txt', sep='\t', header=T)
  single.genes <- read.table('Tb_singleGenes.txt')
  colnames(single.genes) <- 'GeneID'
  single.genes.classes <- merge(single.genes, gene.classes)
  
  data <- reactive({
    inFile <- input$file1
    if(is.null(inFile)){return(NULL)}
    data <- read.csv(inFile$datapath, header=input$header, sep=input$sep, 
                     quote=input$quote)
    colnames(data)[1] <- 'GeneID'
    data <- merge(data, single.genes)
    data <- data.frame(row.names = data$GeneID, data[,-1])
  })
  
  # Preview of input data
  output$contents <- renderTable({
      
      # input$file1 will be NULL initially. After the user selects
      # and uploads a file, it will be a data frame with 'name',
      # 'size', 'type', and 'datapath' columns. The 'datapath'
      # column will contain the local filenames where the data can
      # be found.
      
      if (is.null(data())){return(NULL)}else{head(data(), 100)}
    })
  
  # Generating the Column data frame
  colData <- reactive({
    names <- input$colIDs
    names <- strsplit(names, split = ',')[[1]]
    colData <- data.frame('Samples' = colnames(data()), 'Conditions' = names)
  })
  
  # Doing the DESeq analysis
  dds <- reactive({
    dds <- DESeqDataSetFromMatrix(countData = data(), colData = colData(), design=~Conditions)
    dds <- DESeq(dds)
  })
  
  # Extracting the results of the DESeq analysis
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
  
  # Generating the PCA plot
  output$pca <- renderPlot({
    if(is.null(data()) && is.null(colData())){
      NULL
    }else{
      dds.df <- dds()
      rld <- rlog(dds.df)
      plotPCA(rld, intgroup = 'Samples')
    }
      
  })
  
  # Generating the MA-plot
  output$plotma <- renderPlot({
    if(length(colData())>1){
      res.df <- res()
      plotMA(res.df, main='MA-Plot')
    }else{
      NULL
    }
  })
  
  # Generating the results table
  output$resultTable <- renderDataTable({
    res.df <- classes.df()
    #res.df <- subset(res.df, res.df$padj < input$signif)
    res.df <- data.frame('GeneID' = res.df$GeneID, 'log2FoldChange' = res.df$log2FoldChange, 
                         'p_adj' = res.df$padj, 
                         'Class' = res.df$Class)
    res.df
  })
  
  output$dload <- downloadHandler(
    filename = function(){'Resutls.txt'},
    content = function(file){
      res.all <- res()
      res.df <- as.data.frame(res.all)
      res.df <- data.frame('GeneID' = row.names(res.df), res.df)
      res.df <- merge(res.df, gene.classes, all.x =T)
      #res.df <- subset(res.df, res.df$padj < input$signif)
      write.table(res.df, file, quote=F, row.names=F, sep='\t')
    }
  )
  
  classes.df <- reactive({
    df <- as.data.frame(res())
    df <- data.frame('GeneID' = row.names(df), df)
    df <- merge(df, gene.classes)
    df <- subset(df, df$padj < input$signif)
  })
  
  output$classes <- renderPlot({
    df <- classes.df()
    occurence <- NULL
    occurence.bg <- NULL
    for(class in levels(df$Class)){
      #print(subset(df, Class == class))
      #print(length(subset(df, Class == class)$Class))
      occurence <- c(occurence, length(subset(df, Class == class)$Class))
      occurence.bg <- c(occurence.bg, length(subset(single.genes.classes, Class == class)$Class))
    }
    class.occ <- data.frame('Class' = levels(df$Class), 'Occurence' = occurence, 'Sample' = 'Experiment')
    class.occ <- subset(class.occ, Occurence > 0)
    
    class.occ.bg <- data.frame('Class' = levels(df$Class), 'Occurence' = occurence.bg, 'Sample' = 'Background')
    class.occ.bg <- subset(class.occ.bg, Occurence > 0)
    
    class.occ.norm <- class.occ$Occurence / sum(class.occ$Occurence)
    class.occ <- data.frame(class.occ, 'Normalized' = class.occ.norm)
    class.occ.bg.norm <- class.occ.bg$Occurence / sum(class.occ.bg$Occurence)
    class.occ.bg <- data.frame(class.occ.bg, 'Normalized' = class.occ.bg.norm)
    
    class.occ.sum <- rbind(class.occ, class.occ.bg)
    class.occ.sum <- subset(class.occ.sum, Class != 'ZUnknown')
    ggplot(class.occ.sum, aes(x=Class, y=Normalized, fill=Sample)) + geom_bar(stat='identity', position='dodge') +
      theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1, size=11)) +
      ggtitle('Gene-Classes')
    
  })
  
# End of Server
})

# Run the application 
shinyApp(ui = ui, server = server)

