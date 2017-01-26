
library(shiny)
require(plotly)
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
      tags$h4('DESeq2 Parameters'),
      textInput('colIDs', label = 'Column Description', placeholder = 'eg.: FT,FT,E,E - NO spaces, comma seperated'),
      textInput('norm', 'Normalization (Exp/Crtl)', placeholder = 'eg.: E/FT'),
      numericInput('signif', 'Significance level (p_adjusted)', value = 0.1, step = 0.05),
      radioButtons('althyp', label = 'Alternative Hypothesis',
                   c(Greater = 'greater',
                     Less = 'less',
                     Both = ''),
                   ''
                   ),
      tags$hr(),
      tags$h4('Class enrichment Parameters'),
      numericInput('class.signif', 'Significance level', value=0.05),
      radioButtons('class.althyp', 'log2Fold change',
                   c(Greater = 'greater',
                     Less = 'less',
                     both = 'two.sided'),
                   'greater'
                   ),
      numericInput('MinLog2', 'Cutoff log2Fold-Change', value=0),
      submitButton(text = 'Apply/Update')
    ),
    
    
    #### MAIN PANEL ####
    mainPanel(
      #### Content panel ####
      tabsetPanel(
        tabPanel('Content',
                 downloadButton('rpm.dload',
                                label = 'RPM data download'),
                 tableOutput('contents'),
                 tableOutput('rpm.content')),
        
        #### PCA Panel ####
        tabPanel('PCA', 
                 h4('PCA-Plot', align='center'),
                 plotlyOutput('pca.groups'),
                 downloadButton('dload.pca.groups')
         ),
        
        #### DESeq Panel ####
        tabPanel('DESeq',
                 h4('MA-Plot', align='center'),
                 plotOutput('plotma'),
                 downloadButton('maplot')
         ),
        
        #### Classes Panel ####
        tabPanel('Classes',
                 plotOutput('classes'),
                 downloadButton('classesdload'),
                 br(),
                 p('The gene-class ZUnknown is remove here, since it often represents the vast majority and therefore makes a evaluation of other classes hard'),
                 br(),
                 br(),
                 plotOutput('boxclass'),
                 downloadButton('boxclass.dload'),
                 p('This plot shows the behaviour of each class, considering significant DE-genes.'),
                 br(),
                 br(),
                 plotOutput('c.cycle'),
                 downloadButton('c.cycle.dload'),
                 br(),
                 p('The plot shows the number of selected genes peaking at a specific time during the cell-cycle'),
                 br(),
                 plotlyOutput('len2log'),
                 downloadButton('len2log.dload')
                 
        ),
        
        ##### Results table panel ####
        tabPanel('Result Table',
                 downloadButton('dload', label='Download'),
                 dataTableOutput('resultTable')
                 )
      )
    )
  )
# End of UI
))




#### Server ####
require(DESeq2)
require(ggplot2)

server <- shinyServer(function(input, output) {
  single.genes <- read.table('unique_list_transcLength.txt', 
                             sep='\t', header=T, row.names = 1, quote='')
  cell.cycle <- read.table('Gene_CellCycle_Peak.txt',
                           sep = '\t', header=T,
                           row.names = 1)
  
  data <- reactive({
    inFile <- input$file1
    if(is.null(inFile)){return(NULL)}
    data <- read.csv(inFile$datapath, header=input$header, sep=input$sep, 
                     quote=input$quote, row.names = 1)
    
    data <- data[row.names(data) %in% row.names(single.genes), ]
    data <- data[complete.cases(data), ]
    return(data)
  })
  
  data.rpm <- reactive({
    df <- data()
    df.colsum <- colSums(df)
    
    for(colName in colnames(df)){
      df[[colName]] <- (df[[colName]] / df.colsum[[colName]]) *1000000
    }
    return(df)
  })
  
  output$rpm.content <- renderTable({
    if(is.null(data.rpm())){return(NULL)}else{head(data.rpm(), 30)}
  })
  
  
  output$rpm.dload <- downloadHandler(
    filename = 'RPM_data.txt',
    content = function(file){
      df <- data.rpm()
      df <- data.frame('GeneID' = row.names(df),
                       df)
      write.table(df, file=file, sep='\t', row.names = F, quote = F)
    }
  )
  
  # Preview of input data
  output$contents <- renderTable({
      
      # input$file1 will be NULL initially. After the user selects
      # and uploads a file, it will be a data frame with 'name',
      # 'size', 'type', and 'datapath' columns. The 'datapath'
      # column will contain the local filenames where the data can
      # be found.
      
      if (is.null(data())){return(NULL)}else{head(data(), 30)}
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
  
  rld <- reactive({
    dds.df <- dds()
    rld <- rlog(dds.df)
  })
  
  #### PCA plot ####
  
  # This plot should help to see
  # whether replicates of the same condition cluster
  # well or if they mix
  output$pca.groups <- renderPlotly({
    if(is.null(data()) && is.null(colData())){
      NULL
    }else{
      pca <- plotPCA(rld(), intgroup = 'Conditions', returnData=TRUE)
      pca.s <- plotPCA(rld(), intgroup = 'Samples', returnData=TRUE)
      percentVar <- round(100 * attr(pca, "percentVar"))
      
      pca.pl <- ggplot(pca, aes(PC1, PC2, color=Conditions, text = pca.s$Samples)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance"))
      
      ggplotly(pca.pl, tooltip = 'text')
    }
  })
  
  # Download the intgroup = samples plot
  output$pcaplot <- downloadHandler(
    filename = 'PCA-Plot.pdf',
    content = function(file) {
      pca.plot <- plotPCA(rld(), intgroup = 'Samples')
      ggsave(file, plot = pca.plot, device='pdf')
    }
  )
  
  output$dload.pca.groups <- downloadHandler(
    filename = 'PCA-Group-Plot.pdf',
    content = function(file) {
      pca.s <- plotPCA(rld(), intgroup = 'Samples', returnData=TRUE)
      
      percentVar <- round(100 * attr(pca.s, "percentVar"))
      
      pca.pl <- ggplot(pca.s, aes(PC1, PC2, color=Samples)) +
        geom_point(size=5, shape=15) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance"))
      
      ggsave(file, plot = pca.pl, device='pdf',
             height = 3,
             width = 10)
    }
  )
  
  #### MA-plot ####
  output$plotma <- renderPlot({
    if(length(colData())>1){
      res.df <- res()
      plotMA(res.df, main='MA-Plot')
    }else{
      NULL
    }
  })
  
  # Save the MA-plot
  output$maplot <- downloadHandler(
    filename = 'MA-Plot.pdf',
    content = function(file) {
      res.df <- res()
      pdf(file)
      plotMA(res.df, main='MA-Plot')
      dev.off()
    }
  )
  
  #### Results table ####
  # The following generates the result table
  output$resultTable <- renderDataTable({
    res.df <- classes.df()
    res.df <- data.frame('GeneID' = res.df$Row.names, 
                         'log2FoldChange' = res.df$log2FoldChange, 
                         'p_adj' = res.df$padj,
                         'Annotation' = res.df$Annotation,
                         'Class' = res.df$Class)
    res.df
  })
  
  output$dload <- downloadHandler(
    filename = function(){'Results.txt'},
    content = function(file){
      res.all <- res()
      res.df <- as.data.frame(res.all)
      res.df <- merge(res.df, single.genes, by='row.names')
      res.df <- data.frame('GeneID' = res.df$Row.names, res.df[,-1])
      #res.df <- subset(res.df, res.df$padj < input$signif)
      write.table(res.df, file, quote=F, row.names=F, sep='\t')
    }
  )
  
  #### Subset f. Class enr. ####
  classes.df <- reactive({
    df <- as.data.frame(res())
    df <- merge(df, single.genes, by='row.names') # pack the result table and the unique gene list together in one df
    df <- subset(df, df$padj < input$signif) # Only look at data which is below the significance value
    
    if(input$class.althyp == 'greater'){
      df <- subset(df, df$log2FoldChange > input$MinLog2) # Remove everything below log2FoldChange cutoff
      return(df)
    }
    
    if(input$class.althyp == 'less'){
      df <- subset(df, df$log2FoldChange < input$MinLog2) # Remove everything above log2FoldChange cutoff
      return(df)
    }
    
    if(input$class.althyp == ''){
      return(df) # return all significant genes, without cutoff
    }
  })
  
  #### Class Enrichment ####
  # The following function tries to find significant enriched gene classes
  class.occ.sum <- reactive({
    df <- classes.df() # Get the DESeq2 results data frame with Class and Annotation column
    occurence <- NULL # Initialize the counting vector for the sample
    occurence.bg <- NULL # Initialize the counting vector for the background
    # As background the whole unique gene list is used. This gives an idea about the distribution
    # of classes, if a random picking would be applied.
    
    # Loop over all classes and count their occurence in the sample (genes withing the previous specified cutoffs,
    # padj and log2FoldChange)
    for(class in levels(df$Class)){
      occurence <- c(occurence, length(subset(df, Class == class)$Class))
      occurence.bg <- c(occurence.bg, length(subset(single.genes[row.names(single.genes) %in% row.names(as.data.frame(res())), ], 
                                                    Class == class)$Class))
    }
    
    class.occ <- data.frame('Class' = levels(df$Class), 'Occurence' = occurence, 'Sample' = 'Experiment')
    
    class.occ.bg <- data.frame('Class' = levels(df$Class), 'Occurence' = occurence.bg, 'Sample' = 'Background')
    
    # Fishers exact test to find only significantly enriched classes
    fisher.pvalue <- NULL
    for(class in levels(class.occ$Class)){
      # The following generates the fisher.test matrix
      sam.class <- class.occ[class.occ$Class == class, ]$Occurence
      sam.class <- c(sam.class, sum(class.occ$Occurence) - sam.class)
      bg.class <- class.occ.bg[class.occ.bg$Class == class, ]$Occurence
      bg.class <- c(bg.class, sum(class.occ.bg$Occurence) - bg.class)
      
      f.pvalue <- fisher.test(
        matrix(c(sam.class, bg.class), ncol=2), 
        alternative = 'greater')$p.value # We are only interested in enrichment in the specified subset
      
      fisher.pvalue <- c(fisher.pvalue, f.pvalue) # save the fisher test pvalues
    }
    
    # Since we test multiple times for significance, the FDR must be calculated.
    # Here I use the p.adjust function and the method according to Benjamini-Hochberg
    fisher.pvalue <- p.adjust(fisher.pvalue, method = 'BH')
    
    # Normailize the appearance of each class with the number of genes being looked at.
    # This generates a nice plot, comparing the percentage in the unique-gene list and the 
    # data subset.
    class.occ.norm <- class.occ$Occurence / sum(class.occ$Occurence)
    class.occ <- data.frame(class.occ, 'Normalized' = class.occ.norm, 'Fisher.padj' = fisher.pvalue)
    
    class.occ.bg.norm <- class.occ.bg$Occurence / sum(class.occ.bg$Occurence)
    class.occ.bg <- data.frame(class.occ.bg, 'Normalized' = class.occ.bg.norm, 'Fisher.padj' = fisher.pvalue)
    
    # Putting background and data subset in one data frame for plotting
    class.occ.sum <- rbind(class.occ, class.occ.bg)
    return(class.occ.sum)
  })
  
  # Generate the significantly enriched (fisher exact test) classes barplot
  class.plot <- reactive({
    df <- class.occ.sum()
    df <- subset(df, Class != 'ZUnknown' & Fisher.padj < input$class.signif)
    print(df)
    c.plot <- ggplot(df, aes(x=Class, y=Normalized, fill=Sample)) + 
      geom_bar(stat='identity', position='dodge') +
      theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1, size=11)) +
      geom_text(aes(label=Occurence), position=position_dodge(width=0.9), vjust = -0.25) +
      geom_text(aes(label=round(Fisher.padj, digits = 3), y=0), vjust=1.25) +
      ggtitle('Gene-Classes')
    return(c.plot)
  })
  
  output$classes <- renderPlot({
    class.plot()
  })
  
  output$classesdload <- downloadHandler(
    filename = 'ClassEnrichment.pdf',
    content = function(file){
      p.class <- class.plot()
      ggsave(file, plot = p.class, device = 'pdf', width = 9, height = 5)
    }
  )
  
  # Generate the boxplot which observes all significant (according to DESeq2) genes' classes
  box.class <- reactive({
    res.all <- res()
    res.df <- as.data.frame(res.all)
    res.df <- merge(res.df, single.genes, by='row.names')
    res.df <- data.frame('GeneID' = res.df$Row.names, res.df[,-1])
    res.df <- subset(res.df, padj < input$signif)
    box.plot <- ggplot(res.df, aes(x=Class, y=log2FoldChange)) + 
      geom_boxplot() +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
      geom_hline(yintercept = input$MinLog2, color='red') +
      ggtitle('Subset Fraction: Gene Class distribution')
    
    return(box.plot)
  })
  
  output$boxclass <- renderPlot({
    box.class()
  })
  
  output$boxclass.dload <- downloadHandler(
    filename = 'ClassBoxplot.pdf',
    content = function(file){
      ggsave(file, plot = box.class(), device='pdf', width=12,
             height = 5, 
             units = 'in')
    }
  )
  
  #### Cell cycle plot ####
  c.cycle <- reactive({
    df <- classes.df()
    row.names(df) <- df$Row.names
    df <- df[,-1]
    df <- merge(df, cell.cycle, by='row.names')  # annotate each gene with its peak time of expression
    p.plot <- ggplot(df, aes(x=peak.time)) + geom_bar() +
      ggtitle('Subset Fraction: Peak of expression in cellcycle')
    return(p.plot)
  })
  
  output$c.cycle <- renderPlot({
    c.cycle()
  })
  
  output$c.cycle.dload <- downloadHandler(
    filename = 'CellCycle.pdf',
    content = function(file){
      ggsave(file, plot = c.cycle(), device = 'pdf', width = 8,
             height = 5, 
             units = 'in')
    }
  )
  
  len2log.plot <- reactive({
    df <- as.data.frame(res())
    df <- merge(df, single.genes, by='row.names') # pack the result table and the unique gene list together in one df
    #df <- subset(df, df$padj < input$signif)
    row.names(df) <- df$Row.names
    df <- df[,-1]

    df$TranscriptLength <- as.numeric(as.character(df$TranscriptLength))
    
    df.p <- ggplot(df, aes(x=TranscriptLength, y=log2FoldChange, text=Annotation,
                           color= padj < input$signif)) +
      scale_color_manual(name = paste('padj <', input$signif, sep = ' '), 
                         values = setNames(c('red', 'black'), c(T,F))) +
      geom_point(size = 0.5) +
      scale_x_continuous(trans = 'log2') +
      ylab('Log(2)-Fold Change') +
      xlab('Transcript Length [bp] log(2)-scale') +
      geom_hline(yintercept = input$MinLog2, color='gray')
    
    return(df.p)
  })
  
  output$len2log <- renderPlotly({
    # Renders the plot showing CDS length vs log2FoldChange
    ggplotly(len2log.plot())
  })
  
  output$len2log.dload <- downloadHandler(
    filename = 'Length_vs_log2FoldChange.pdf',
    content = function(file){
      ggsave(file, plot = len2log.plot(), device = 'pdf', width = 8,
             height = 5, 
             units = 'in')
    }
  )
  
# End of Server
})

# Run the application 
shinyApp(ui = ui, server = server)

