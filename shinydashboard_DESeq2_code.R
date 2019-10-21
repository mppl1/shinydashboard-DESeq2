# To access shiny-server at 0.0.0.0:3838, I will need to start another process of container for ex. sudo docker exec -ti -e ROOT=TRUE 0ca61e5e2337 shiny
# It is already running(if navigate to 0.0.0.0:3838) but I am unable to make changes without command line enquiries
# This script is to output basic Mov10_DESeq2 analysis into a shiny server for display purposes



# First stage is to just perform analysis and then to allow users to interact with the data after DESeq2 statistics

##### Load data #####
#
#metadata<-read.csv(url("https://raw.githubusercontent.com/hbc/NGS_Data_Analysis_Course/master/sessionIII/data/Mov10_full_meta.txt"), sep = "\t", row.names = 1)

if( !exists("countdata")){
  countdata<-read.csv(url("https://raw.githubusercontent.com/hbc/NGS_Data_Analysis_Course/master/sessionIII/data/Mov10_full_counts.txt"), sep = "\t", row.names = 1)
}

if( !exists("metadata")){
  metadata<-read.csv(url("https://raw.githubusercontent.com/hbc/NGS_Data_Analysis_Course/master/sessionIII/data/Mov10_full_meta.txt"), sep = "\t", row.names = 1)
}


# check that the sample names of count matrix match row names of phenodata. DESeq2 expects this 
all(rownames(metadata) %in% colnames(countdata))
all(rownames(metadata) == colnames(countdata))


##############Once algorithms are done it is best to save and reload data for display####################
# ## For EDA
# # create DESeq2 object. For design give the column with factor variable in metadata or phenodata, which indicates the grouping of samples
# dds <- DESeqDataSetFromMatrix(countData = countdata, colData = metadata, design = ~ sampletype)
# # add size factors to slot
# dds <- estimateSizeFactors(dds)
# # generate normalized counts (median of ratios method)
# normalized_counts <- counts(dds, normalized=TRUE)
# 
# #To dampen the variance mean correlation perform variance stabilization to normalized data
# #  This is more appropriate for unsupervised clustering etc.
# rld<-rlog(dds, blind = T)

###########################################



##########################################################################################
# Perform DESeq2 analysis
#dds<-DESeq(dds)


# ### Actual fold change test correcting for multiple testing
# contrast_oe <- c("sampletype", "MOV10_overexpression", "control")
# 
# res_tableOE_unshrunken <- results(dds, contrast=contrast_oe, alpha = 0.05)
# 
# #res_tableOE <- lfcShrink(dds, contrast=contrast_oe, res=res_tableOE_unshrunken) # Here with contrast default is type='normal'; type = c("normal", "apeglm", "ashr")
# # To get the coef names. Only available once "DESeq" is called
# colnames(coef(dds))
# res_tableOE_shrunk <- lfcShrink(dds, coef = "sampletype_MOV10_overexpression_vs_control", res= res_tableOE_unshrunken , type = 'apeglm')
# 


##################################################################################################################################

# These are overhead computations and to ensure that they are not done each time the app loads, will use a conditional for existence
# if exists conditional https://stackoverflow.com/questions/28218508/r-check-if-r-object-exists-before-creating-it

if( exists("res_table_shrunk")){
  print("res_table_shrunk was in global env. at start")
  rm(res_table_shrunk)
} else{
  print("res_table_shrunk was Not found so deletion worked.")
  }

# ## For EDA
# create DESeq2 object. For design give the column with factor variable in metadata or phenodata, which indicates the grouping of samples
if( !exists("ddsobj")){
  ddsobj <- DESeqDataSetFromMatrix(countData = countdata, colData = metadata, design = ~ sampletype)
}

if( !exists("ddsnorm")){
  ddsnorm <- estimateSizeFactors(ddsobj)
}
# add size factors to slot

if( !exists("ddss")){
  library(DESeq2)
  ddss <- DESeq(ddsnorm)
}

#To dampen the variance mean correlation perform variance stabilization to normalized data
#  This is more appropriate for unsupervised clustering etc.

if( !exists("rld")){
  rld<-rlog(ddsnorm, blind = T)
}


library(scales)

if( !exists("pca")){
  rld_mat <- assay(rld)
  pca <- prcomp(t(rld_mat)) # Single Value Decomposition
}


# screeplot(pca) # Not clear what the y-axis represents
# To plot screeplot with ggplot2
dfscreeplot<-data.frame(t(summary(pca)[6]$importance))
# This will be used for plotting PCAs
dfpca <- cbind(metadata, pca$x)


##################################################################################################################################



library(shiny)
library(shinydashboard)
library(shinyjs) # for dynamically displaying tab titles



ui <- dashboardPage(
 
  dashboardHeader(title = "Basic dashboard"),
  ## Sidebar content
  dashboardSidebar(
    sidebarMenu(id="tab", #style = "position: fixed; overflow: visible;", # To fix the sidebar so it doesn't scroll with body. but title still moves etc.
      menuItem("EDA", tabName = "EDA", icon = icon("dashboard")),
      menuItem("Actual DESeq2 analysis", tabName = "Actual_DESeq2_analysis", icon = icon("th")),
      menuItem("Third_selection", tabName = "thirdselection", icon = icon("bar-chart-o")) #,
      #menuItem("table_tab", tabName = "table_tab", icon = icon("user")) # this part is causing problems
    ),
    uiOutput("out1")
  ),
  

  
  ## Body content
  dashboardBody(
    # This is so that headers can change depending on tab selected
    tags$head(tags$style(HTML(
      '.myClass { 
            font-size: 20px;
            line-height: 50px;
            text-align: left;
            font-family: "Helvetica Neue",Helvetica,Arial,sans-serif;
            padding: 0 15px;
            overflow: hidden;
            color: white;
            }
            '))),
    tags$script(HTML('
                           $(document).ready(function() {
                           $("header").find("nav").append(\'<div id="pageHeader" class="myClass"></div>\');
                           })
                           ')),
  
  
    tabItems(
      # First tab content
      tabItem(tabName = "EDA", #header="EDA", #h2("Exploratory Data Analysis"), # this does not display it in header
              # Boxes need to be put in a row (or column)
              fluidRow(
                box(plotOutput("plot1", height = 400), width = 12)#,
                
                # box(
                #   title = "Controls",
                #   sliderInput("sliderbin", "Chose Number of Bins:", min = 50, max = 500, step=20, value=200),
                #   selectInput("selectsample", "Select Sample for histogram", choices = colnames(countdata), selected= colnames(countdata)[3] , multiple = F ),
                #   # sliderInput("sliderbin", "Number of observations:", 1, 100, 50)
                #   # selectInput("selectpcx", "Select PCs for x-axis", choices = rownames(dfscreeplot), selected= c(rownames(dfscreeplot)[1:2]) , multiple = F ),
                #   # selectInput("selectpcy", "Select PCs for y-axis", choices = rownames(dfscreeplot), selected= c(rownames(dfscreeplot)[1:2]) , multiple = F )
                #   br("*******************************************"),
                #   selectInput("selectpcx", "Select PCs for x-axis", choices = colnames(dfpca)[3:10], selected= c(colnames(dfpca)[3:10][1]) , multiple = F ),
                #   selectInput("selectpcy", "Select PCs for y-axis", choices = colnames(dfpca)[3:10], selected= c(colnames(dfpca)[3:10][2]) , multiple = F )
                #   
                # )
              ),
              
              br("The plots below are of dataset after it is normalized (sizefactor estimated using median of ratios method) \n and rlog transformed for variance stabilization "),
              
              fluidRow(
                
                box(plotOutput("plot2", height = 400)),
                box(plotOutput("plot3"))
              ),
              
              fluidRow(
                box(plotOutput("plot4", height=400))
              )
      ),
      
      # Second tab content
      tabItem(tabName = "Actual_DESeq2_analysis",
              fluidRow(
                # box(h2("Actual DESeq2 analysis tab content"),
                #     selectInput(inputId= "selectcontrast", "Select the design column[Condition column of metadata or phenodata]", choices=colnames(metadata), selected=colnames(metadata)[1]) ,
                #     uiOutput("out2"),
                #     textOutput(outputId = "contrastchoice"),
                #     # tags$head( 
                #     #   tags$style(HTML('#run{background-color:orange}'))
                #     # ),
                #     # actionButton("res_param_submit", "Submit parameters for calculations") #, icon = "table" The parameters for DESeq2 results and lfcshrink functions
                #     # The above was an alternative
                #     actionButton("res_param_submit", "Submit parameters for calculations", icon("paper-plane"), 
                #                  style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                #     ),
                # box(plotOutput("mtcars_table"))
                # box(textOutput("mtcars_table"))
                box(
                  dataTableOutput("mtcars_table")
                  # verbatimTextOutput("mtcars_table")
                  ),
                box(
                  dataTableOutput("pvaluecut_table")
                )
                
              ),
              
              fluidRow(
                # plotlyOutput("mtcars_plotMA")
                plotOutput("mtcars_plotMA"), height=600
              ),
              
              # fluidRow(
              #   box(h4("Analysing results of DESeq2 analysis"),
              #       textInput(inputId= "selectL2FC", "Select the log2FC threshold(enter positive number)", value = 1, placeholder = "Positive numbers from 0 to any foldexpression") ,
              #       textInput("selectpadj", label="Type in alpha or P-value threshold", value = 0.05, placeholder = "numbers only between 0 and 1 for example 0.05")
              #       # uiOutput("out2"), 
              #       # textOutput(outputId = "contrastchoice"),
              #       # tags$head( 
              #       #   tags$style(HTML('#run{background-color:orange}'))
              #       # ),
              #       # actionButton("res_param_submit", "Submit parameters for calculations") #, icon = "table" The parameters for DESeq2 results and lfcshrink functions
              #       # The above was an alternative
              #       # actionButton("res_param_submit", "Submit parameters for calculations", icon("paper-plane"), 
              #       #              style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
              #   )#,
              #   # box(plotOutput("mtcars_table"))
              #   # box(textOutput("mtcars_table"))
              #   # box(dataTableOutput("siggenes_table"))
              #   
              # ),
              fluidRow(
                box(
                dataTableOutput("siggenes_table"), width = 12
                  )
              ),
              
              fluidRow(
               
                  plotOutput("plot5", height = 400)
               
              ),
              fluidRow(
                
                  plotOutput("plot6", height = 500)
                
              ),
              fluidRow(plotOutput("plot7", height = 600))
              
      ),
      
      # Third tab content
      tabItem(tabName = "thirdselection",
              h2("Third selection content")
      #         ),
      # tabItem(tabName = "table_tab",
      #         # fluidRow(
      #         #   valueBox("Unread Mail", 144, icon("mail"), color = "blue", width = 6, size = "small"),
      #         #   valueBox("Spam", 20, icon("mail"), color = "red", width = 5, size = "small"),
      #         #   valueBox("Readed Mail", 666, icon("mail"), color = "green", width = 5, size = "small")
      #         # ),
      #         fluidRow(
      #           box(title = "Classic box", color = "blue", ribbon = FALSE,
      #               title_side = "top left", width = 14,
      #               dataTableOutput("mtcars_table")
      #               # tags$div(
      #               #   dataTableOutput("mtcars_table")
      #               #   , style = paste0("color:", semantic_palette[["blue"]], ";"))
      # 
      #           ))),
    ))
  ),
  useShinyjs() # this required installation of shinyjs package. Moving it here (outside scope of "dashboardbody" but inside "dashboardpage") also allowed for proper function. This is for dynamic header title change
  
  
  
  
)







server <- function(input, output, session) {
  
  output$plot1 <- renderPlot({
    # data <- histdata[seq_len(input$slider)]
    # hist(data)
    # This is to plot histogram of raw count data
    ggplot(countdata)+geom_histogram(aes(x=countdata[, input$selectsample]), stat = "bin", bins=input$sliderbin)+ xlab("Raw expression counts") +
      ylab("Number of genes")+ ggtitle("Histogram of Raw Expression Values")
  })
  # Screeplot of the PCs
  output$plot2 <- renderPlot({
    ggplot(dfscreeplot, aes(x=rownames(dfscreeplot), y=dfscreeplot$Proportion.of.Variance, group=1))+
      geom_bar(stat = "identity", fill="steelblue", position = "dodge")+ # stat = "count" for geom_col
      ggtitle("Screeplot of Principle Components") +
      geom_text(aes(x = rownames(dfscreeplot), y = dfscreeplot$Proportion.of.Variance, label = round(dfscreeplot$Proportion.of.Variance, 2)), vjust=-0.3)
  })
  
  output$plot3 <- renderPlot({
    # This is to plot a 2D plot of the user chosen/selected PCs
    # dfpca <- cbind(metadata, pca$x)
    ggplot(dfpca) + geom_point(aes(x=dfpca[,input$selectpcx], y=dfpca[,input$selectpcy], color = sampletype), size=2)+
      labs(title=paste("Plot of", input$selectpcx, "vs", input$selectpcy), x=paste(input$selectpcx, " ",100*dfscreeplot$Proportion.of.Variance[rownames(dfscreeplot)==input$selectpcx], "%"), y = paste(input$selectpcy," ", 100*(dfscreeplot$Proportion.of.Variance[rownames(dfscreeplot)==input$selectpcy] ), "%")) +
      # labs(title="Plot of PC1 vs PC2", x=paste("PC1", " ",100*dfscreeplot$Proportion.of.Variance[[1]], "%"), y = paste("PC2"," ", 100*(dfscreeplot$Proportion.of.Variance[[2]] ), "%")) +
      theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1))
    
    #dfscreeplot$Proportion.of.Variance[rowname(dfscreeplot)==input$selectpcx]
    
  })
  
  output$plot4 <- renderPlot({
    rld_mat <- assay(rld) # although it was already available in screeplot
    rld_cor <- cor(rld_mat) # Using the rld_mat object created above
    library(pheatmap)
    library(RColorBrewer)

    # heat.colors <- brewer.pal(6, "Blues") #(6, "Blues")
    heat.colors <- brewer.pal(6, "RdYlGn") #(6, "Blues")
    pheatmap(rld_cor,  main = "Heatmap of correlation matrix" , color = heat.colors, border_color=NA, fontsize = 10,
             fontsize_row = 10) # color = heat.colors, , height=20
  })
  
  # This is to change subitems of tabs
  output$out1 <- renderUI({
    
    if (input$tab == "EDA") {
      
      dyn_ui <- tabsetPanel(id = "tabset_id", selected = "EDA_Histogram", 
                            tabPanel("EDA_Histogram", value = "EDA_Histogram",
                                     list(selectInput("selectsample", "Select Sample for histogram", choices = colnames(countdata), selected= colnames(countdata)[3] , multiple = F ),
                                          sliderInput("sliderbin", "Chose Number of Bins:", min = 50, max = 500, step=20, value=200)
                                     )
                                     # The above is aesthetically cleaner than a box
                                     # box(
                                     #   title = "Histogram Controls: \n choose sample and bins",
                                     #   
                                     #   selectInput("selectsample", "Select Sample for histogram", choices = colnames(countdata), selected= colnames(countdata)[3] , multiple = F ),
                                     #   sliderInput("sliderbin", "Chose Number of Bins:", min = 50, max = 500, step=20, value=200)
                                     #   # sliderInput("sliderbin", "Number of observations:", 1, 100, 50)
                                     #   # selectInput("selectpcx", "Select PCs for x-axis", choices = rownames(dfscreeplot), selected= c(rownames(dfscreeplot)[1:2]) , multiple = F ),
                                     #   # selectInput("selectpcy", "Select PCs for y-axis", choices = rownames(dfscreeplot), selected= c(rownames(dfscreeplot)[1:2]) , multiple = F )
                                     #   # br("*******************************************"),
                                     # , width = 12)
                                     ),
                            tabPanel("EDA_PCA", value = "EDA_PCA", 
                                     list(selectInput("selectpcx", "Select PCs for x-axis", choices = colnames(dfpca)[3:10], selected= c(colnames(dfpca)[3:10][1]) , multiple = F ),
                                          selectInput("selectpcy", "Select PCs for y-axis", choices = colnames(dfpca)[3:10], selected= c(colnames(dfpca)[3:10][2]) , multiple = F ))
                                     # box(
                                     #   title = "PCA Controls: \n choose PCs",
                                     #   # selectInput("selectpcx", "Select PCs for x-axis", choices = rownames(dfscreeplot), selected= c(rownames(dfscreeplot)[1:2]) , multiple = F ),
                                     #   # selectInput("selectpcy", "Select PCs for y-axis", choices = rownames(dfscreeplot), selected= c(rownames(dfscreeplot)[1:2]) , multiple = F )
                                     #   # br("*******************************************"),
                                     #   selectInput("selectpcx", "Select PCs for x-axis", choices = colnames(dfpca)[3:10], selected= c(colnames(dfpca)[3:10][1]) , multiple = F ),
                                     #   selectInput("selectpcy", "Select PCs for y-axis", choices = colnames(dfpca)[3:10], selected= c(colnames(dfpca)[3:10][2]) , multiple = F )
                                     #   
                                     #   , width = 12)
                                     ))
      
    }
    if (input$tab == "Actual_DESeq2_analysis") {
      
      dyn_ui <-
      # (tags$style(HTML("
      #               .tabbable > .nav > li > a                  {background-color: aqua;  color:black}
      #               .tabbable > .nav > li > a[data-value='Waldtest'] {background-color: red;   color:white}
      #               .tabbable > .nav > li > a[data-value='t2'] {background-color: blue;  color:white}
      #               .tabbable > .nav > li > a[data-value='t3'] {background-color: green; color:white}
      #               .tabbable > .nav > li[class=active]    > a {background-color: black; color:white}
      #             ")),
      #   
          # tabPanel("t2",h2("blue tab")), 
          tabsetPanel(id = "tabset_id", selected = "results_waldtest", 
                            tabPanel("Waldtest", value = "sumtable", 
                                     list(selectInput(inputId= "selectcontrast", "Select the design column[Condition column of metadata or phenodata]", choices=colnames(metadata), selected=colnames(metadata)[1]) ,
                                          uiOutput("out2"),
                                          textOutput(outputId = "contrastchoice"),
                                          # tags$head( 
                                          #   tags$style(HTML('#run{background-color:orange}'))
                                          # ),
                                          # actionButton("res_param_submit", "Submit parameters for calculations") #, icon = "table" The parameters for DESeq2 results and lfcshrink functions
                                          # The above was an alternative
                                          actionButton("res_param_submit", "Submit parameters for calculations", icon("paper-plane"), 
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                     )
                                     # The above is aesthetically cleaner than a box
                                     # box(
                                     #   title = "Histogram Controls: \n choose sample and bins",
                                     #   
                                     #   selectInput("selectsample", "Select Sample for histogram", choices = colnames(countdata), selected= colnames(countdata)[3] , multiple = F ),
                                     #   sliderInput("sliderbin", "Chose Number of Bins:", min = 50, max = 500, step=20, value=200)
                                     #   # sliderInput("sliderbin", "Number of observations:", 1, 100, 50)
                                     #   # selectInput("selectpcx", "Select PCs for x-axis", choices = rownames(dfscreeplot), selected= c(rownames(dfscreeplot)[1:2]) , multiple = F ),
                                     #   # selectInput("selectpcy", "Select PCs for y-axis", choices = rownames(dfscreeplot), selected= c(rownames(dfscreeplot)[1:2]) , multiple = F )
                                     #   # br("*******************************************"),
                                     # , width = 12)
                            ),
                            tabPanel("narrow_genes", value = "Narrow down genes", 
                                     list(textInput(inputId= "selectL2FC", "Select the log2FC threshold(enter positive number)", value = 1, placeholder = "Positive numbers from 0 to any foldexpression") ,
                                          textInput("selectpadj", label="Type in alpha or P-value threshold", value = 0.05, placeholder = "numbers only between 0 and 1 for example 0.05")
                                          )
                                     # box(
                                     #   title = "PCA Controls: \n choose PCs",
                                     #   # selectInput("selectpcx", "Select PCs for x-axis", choices = rownames(dfscreeplot), selected= c(rownames(dfscreeplot)[1:2]) , multiple = F ),
                                     #   # selectInput("selectpcy", "Select PCs for y-axis", choices = rownames(dfscreeplot), selected= c(rownames(dfscreeplot)[1:2]) , multiple = F )
                                     #   # br("*******************************************"),
                                     #   selectInput("selectpcx", "Select PCs for x-axis", choices = colnames(dfpca)[3:10], selected= c(colnames(dfpca)[3:10][1]) , multiple = F ),
                                     #   selectInput("selectpcy", "Select PCs for y-axis", choices = colnames(dfpca)[3:10], selected= c(colnames(dfpca)[3:10][2]) , multiple = F )
                                     #   
                                     #   , width = 12)
                            ))
      # )
      
      
    }
    # if (input$tab == "Actual_DESeq2_analysis") {
    #   
    #   dyn_ui <- list(selectInput("s1", label = "Select", choices = letters[1:3]), #choices = letters[1:3]
    #                  selectInput("s2", label = "Select2", choices = letters[4:6]))  #choices = letters[4:6]
    # }
    if (input$tab == "thirdselection") {
      
      dyn_ui <- list(selectInput("z1", label = "Select", choices = letters[1:3]), #choices = letters[1:3]
                     selectInput("z2", label = "Select2", choices = letters[4:6]))  #choices = letters[4:6]
    }
    return(dyn_ui)
  })
  
  # This is for changing the header as tabs change. Typically without this a dropdown menu is expected in order to display
  observeEvent(input$tab, {
    header <- switch(input$tab,
                     EDA = "EDA",
                     Actual_DESeq2_analysis = "Actual DESeq2 analysis",
                     thirdselection = "Tab 3"#,
                     #table_tab="Table"
                     )
    
    # you can use any other dynamic content you like
    shinyjs::html("pageHeader", header)
  })
  
  # This is no longer necessary as I have incorporated it into uiOutput("out2") below and in UI
  # output$value <- renderText(paste0("A + B = ", input$A + input$B, input$selectcontrast))
  # output$slider <- renderUI({
  #   sliderInput(inputId = "B", label = "B", min = 0, max = 2*input$A, value = 5)
  # })
  
  # To dynamically populate the 
  output$out2 <- renderUI({
    xvar <- input$selectcontrast
    # dyn_ui2 <- selectInput("condition1", label="Select3", choices=metadata$input$selectcontrast, selected = "testing")
    dyn_ui2 <- list(selectInput("condition_comp", label="select cond'n to compare", choices=metadata[,(as.character(xvar))], selected = "testing"),
                    selectInput("condition_ref", label="select Reference (ex. control)", choices=metadata[,(as.character(xvar))], selected = "testing"),
                    textInput("alpha_threshold", label="Type in alpha or P-value threshold  for wald test", value = 0.05, placeholder = "numbers only between 0 and 1 for example 0.05") #,
                    #selectInput("alpha_threshold", "select alpha threshold", choices = seq(0,1, 0.025))
                    )
 
    return(dyn_ui2)
  })
  
  # output$contrastchoice <- renderText(paste0("A + B = ", input$selectcontrast, input$condition_comp,  input$condition_ref, input$alpha_threshold))
  output$contrastchoice <- renderText(paste("The conditions chosen for comparison are:  ", input$condition_comp, " And "  , input$condition_ref, "\n", "The chosen alpha is: ", input$alpha_threshold))
  
  
  observeEvent(input$res_param_submit, ignoreInit=TRUE, {
    contrast_exp <- c(input$selectcontrast, input$condition_comp, input$condition_ref)
    # paste(input$selectcontrast, input$condition_comp, "vs", input$condition_ref, collapse = "_")->tempvar
    # print(tempvar)
    # print(paste(contrast_exp, collapse = "_"))
    if (input$condition_comp != input$condition_ref){
      library(DESeq2)
      res_table_unshrunken <- results(ddss, contrast=contrast_exp, alpha = as.character(input$alpha_threshold))
      res_table_shrunk <<- lfcShrink(ddss, coef = paste(c(input$selectcontrast, input$condition_comp, "vs", input$condition_ref), collapse = "_"), res= res_table_unshrunken , type = 'apeglm')
      
      capture.output(summary(res_table_shrunk))->tempsumtab
      tempsumtab[3:length(tempsumtab)-1] -> tempsumtab
      data.frame( summary_of_DESeq_stat_results=tempsumtab)->tempsumtabdf
      # as.data.frame(index= c(1:length(tempsumtab)), summary_of_results=tempsumtab)->tempsumtabdf
      # output$mtcars_table <- renderPlot(tempsumtab)
      # output$mtcars_table <- renderText(tempsumtab)
      output$mtcars_table <- renderDataTable(tempsumtabdf)
      # output$mtcars_table <- renderText(tempsumtab)
      # output$mtcars_table <- renderText(tempsumtabdf)
      
      # Another table to show how the number of genes passing threshold change as p-value threshold changes
      output$pvaluecut_table <- renderDataTable({
        ## Split features by different p-value cutoffs
        # removes the NA values from calculations
        data.frame(res_table_shrunk) -> res_shrunk_df
        pval_table <- lapply(c(1e-04, 0.001, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
                               0.6, 0.7, 0.8, 0.9, 1), function(x) {
                                 data.frame('P-value' = x, 'Num_of_Genes' = sum(res_shrunk_df$pvalue <= x, na.rm = TRUE))
                               })
        pval_table <- do.call(rbind, pval_table)
        
        pval_table
        
        
      })
      
      
      
      # output$mtcars_plotMA <- renderPlotly({
      output$mtcars_plotMA <- renderPlot({
        # here is the MA plot of just with user input alpha=0.05
        library(ggplot2)
        # library(plotly)
        # library(cowplot)
        
        data.frame(res_table_shrunk) -> res_shrunk_df
        
        #plotMA(res_table_shrunk)
        
        res_shrunk_tempdf<-data.frame(res_table_shrunk)
        res_shrunkdf <- data.frame(res_shrunk_tempdf[,c(1,2,3,4)], padj=(res_shrunk_tempdf$padj < as.double(input$alpha_threshold))) #0.05
        rm(res_shrunk_tempdf)
        res_shrunksubdf <- res_shrunkdf[!is.na(res_shrunkdf$baseMean) & !is.na(res_shrunkdf$log2FoldChange), ]
        # res_shrunksubdf2 <- res_shrunksubdf
        res_shrunksubdf[is.na(res_shrunksubdf$padj),"padj"]<-"FALSE" # to resolve issue of padj=NA being dropped
        res_shrunksubdf[(abs(res_shrunksubdf$log2FoldChange) > 2),"log2FoldChange"]<-2 # to resolve issue of values beyond limit being chopped off
        
        paste(tail(strsplit(res_table_shrunk@elementMetadata$description[[2]], split=" ")[[1]], -5), collapse = " ") -> conditionstatstest
        
        MAplot<-ggplot(res_shrunksubdf, aes((baseMean), log2FoldChange, colour = padj, shape=(abs(res_shrunksubdf$log2FoldChange)==2))) + # as.double(input$alpha_threshold) for plotly, qval=qval, geneNames=geneNames, id=id
          scale_color_manual(values=c("#999999", "#FF0000")) +
          geom_point(size=1) + # size=0.5
          labs(shape = "L2FC>2")+
          guides(shape = guide_legend(override.aes = list(size=3)), color = guide_legend(override.aes = list(size=3)) )+ # THIS WORKED BEST
          # The call on text size worked but not for symbols themselves
          theme(legend.title = element_text(size = 12), legend.key.size = unit(0.2, "cm"), legend.text = element_text(color = "red", size = 12))+ # , legend.key.width = unit(1.5,"cm") note useful
          labs(color = paste(c("padj", "<", input$alpha_threshold ), collapse = " "))+ #as.double(input$alpha_threshold)
          ylim(-2, 2)+
          xlab(label="log10(Mean of normalized counts)")+
          # scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
          #               labels = trans_format("log10", math_format(10^.x)))+
          scale_x_log10(labels=scales::scientific)+
          geom_hline(yintercept=0, size=1)+
          ggtitle(paste("MA plot of mean of normalized counts vs log2 Fold Change after \n lfcshrinkage (apeglm)", conditionstatstest) )+#+ #  ggtitle(paste("MA plot of mean of normalized counts vs log2 Fold Change", input$alpha_threshold))+
          theme (panel.border = element_blank(),
                 axis.line    = element_line(color='black'),
                 panel.background = element_rect(fill="transparent"),
                 axis.text=element_text(size=12, colour = "black")
          )
        
        # ggplotly(MAplot) # , tooltip = c("text", "size")
        MAplot
        
        
        })
      # output$mtcars_table <- renderDataTable(summary(res_tableOE_shrunk)) # , options = list(dom = 't')
    } else {
      # output$mtcars_table <- renderText("Error! You can't compare the condition to itself")
      output$mtcars_table <- DT::renderDataTable({
        library(DT)

        daterror<-data.frame(Error=c("Error! You didn't make correct comparison", "both conditions CAN'T be the same"))
        daterror2<-DT::datatable(daterror) %>% #,
          formatStyle('Error',  color = 'red', backgroundColor = 'yellow', fontWeight = 'bold') #'orange'
        # span(style = "color:green;font-weight:bold;", textOutput("right")) %>% formatStyle(backgroundColor=styleEqual(1, "red"))
        return(daterror2)
        })

      # output$Table<- DT::renderDataTable({ 
      #   dat <- datatable(iris, options = list(paging=FALSE)) %>%
      #     formatStyle('Sepal.Length',  color = 'red', backgroundColor = 'orange', fontWeight = 'bold')
      #   return(dat)
      # })
      
      } #"test1",
      # print("Error")
    
    
  })

  toListen <- reactive({
    list(input$selectL2FC, input$selectpadj, input$res_param_submit) #, input$condition_comp, input$condition_ref
  })
  
  # observeEvent(input$selectL2FC & input$selectpadj, { # ignoreInit=TRUE,
  # by wrapping the events in the reactive call, this event is always dependent on user inputs
  observeEvent(toListen(), ignoreInit=TRUE, { # ignoreInit=TRUE,
    
    output$siggenes_table <- DT::renderDataTable({
    
    # text input ids: selectL2FC  selectpadj
    
    DEGdf<-res_table_shrunk %>% data.frame() #%>% data.frame
    #dim(DEG__padj5<-DEG__padj5 %>% base::filter(padj < padj.cutoff))
    subset.data.frame(DEGdf, padj < as.double(input$selectpadj))->DEG_padj_df
    DEG_padj_lfc_df<-subset.data.frame(DEG_padj_df, padj < as.double(input$selectpadj) & abs(log2FoldChange)>as.double(input$selectL2FC))
   
    # reorder the table with respect to padj.
    # But it turns out user can interact with the table and order the columns accordingly
    # DEG_KD_padj5[head(order(DEG_padj_lfc_df$padj), 20),]
    
    
    library(DT)
    DEG_padj_lfc_dt<-DT::datatable( DEG_padj_lfc_df)
    return( DEG_padj_lfc_dt)
  
    })
    output$plot5 <- renderPlot({
      
      
      DEG_df<-res_table_shrunk %>% data.frame() #%>% data.frame
      top20_sigOE_genes<-DEG_df[head(order(DEG_df$padj), 20),]
      # top20_sigOE_genes<-DEG_padj_lfc_df[head(order(DEG_padj_lfc_df$padj), 20),]
      # assumes the normalized count matrix was already generated in the beginning of session as global variable. See top
      top20_sigOE_normdf<-data.frame(normalized_counts[(rownames(normalized_counts) %in% rownames(top20_sigOE_genes)),])
      
      # get top20 genes normalized counts to plot with ggplot2
      library(tidyr)
      rownames(top20_sigOE_normdf)->top20_sigOE_normdf$gene
      top20_sigOE_normdf<-top20_sigOE_normdf[,c(ncol(top20_sigOE_normdf),(1:ncol(top20_sigOE_normdf)-1))] # reorder
      top20_sigOE_normdfl<-top20_sigOE_normdf %>% gather("Sample", "Normalized_Counts", (colnames(top20_sigOE_normdf)[-1]))
      
      # Since it will be grouped based on factors it is important to let ggplot2 know that different replicates of same group are the same, for this . Otherwise each replicate will be a separate group
      metadata->mov10_meta
      rownames(metadata)->mov10_meta$Sample
      library(tidyr)
      top20_sigOE_normdfl <- inner_join(mov10_meta, top20_sigOE_normdfl)
      
      ## plot using ggplot2
      ggplot(top20_sigOE_normdfl) +
        geom_point(aes(x = gene, y = Normalized_Counts, color = sampletype)) +
        scale_y_log10() +
        xlab("Genes") +
        ylab("log10 Normalized Counts") +
        ggtitle("Top 20 Significant DE Genes by padj value") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5))
      
    })
    output$plot6 <- renderPlot({
      
      library(pheatmap)
      
      ### Set a color palette
      library(RColorBrewer)
      heat_colors <- brewer.pal(6, "YlOrRd")
      
      # sig_normdf<-data.frame(normalized_counts[(rownames(normalized_counts) %in% rownames(DEG_OE_padj5_lfc58)),c(3:8)]) 
      # these are columns in normalized counts and exclude KD samples but it should depend on user input in results calculation
      # take advantage of res_table_shrunk being in global environement and that this renderPlot is within scope of 
      # observeEvent call, which means it is responsive to the input value as
      
      DEG_df<-res_table_shrunk %>% data.frame() #%>% data.frame
      top100_sigOE_genes<-DEG_df[head(order(DEG_df$padj), 100),]
      sig_normdf<-data.frame(normalized_counts[(rownames(normalized_counts) %in% rownames(top100_sigOE_genes)), ]) # c(3:8)
      
      # To get condition compared by user, the "res_table_shrunk" datastructure(which is in global env.) has it and can be accessed as follows
      paste(tail(strsplit(res_table_shrunk@elementMetadata$description[[2]], split=" ")[[1]], -5), collapse = " ") -> conditionstatstest
      ###Annotation data####
      annotation= data.frame(mov10_meta$sampletype, row.names = rownames(mov10_meta))
      print(conditionstatstest)
      # print( c("Hello to see if input$condition_com value is still available", input$condition_comp, input$condition_ref ))
      ### Run pheatmap
      library(pheatmap)
      pheatmap(sig_normdf, #sigOE_normdf[, c(2:ncol(sigOE_normdf))]
               color = heat_colors, 
               cluster_rows = T, 
               show_rownames = F,
               annotation = annotation, 
               border_color = NA, 
               fontsize = 10, 
               scale = "row", # this avoids manual z-score for normalized counts
               fontsize_row = 10, 
               height = 20,
               main = paste("heatmap clustering of samples by \n normalized countdata of \n top100 genes with small padj value based on \n ",   conditionstatstest) #control vs Mov10 O/E"
      )
      
      
    })
    
    output$plot7 <- renderPlot({
      
      # This a volcano plot to give another perspective on how data is visualized using 
      # padj and chosen log2fold change (on top of top20 genes and sample based on 
      # log2foldchange etc. clustering)
      library("ggplot2") #Best plots
      library("ggrepel") #Avoid overlapping labels
      
      # user input
      # input$selectL2FC, input$selectpadj
      
      DEG_df<-res_table_shrunk %>% data.frame()
      
      DEG_df <- cbind(gene=rownames(DEG_df), DEG_df)
      
      DEG_df <- DEG_df[!is.na(DEG_df$log2FoldChange) & !is.na(DEG_df$pvalue),]
      
      # may have to coerce the log2FC, pvalue, and padj into double
      print(head(DEG_df))
      l2FCthreshold <- as.double(input$selectL2FC)
      pvaluethreshold <- as.double(input$selectpadj)
      
      log10_pval<-round(-log10(pvaluethreshold), 2)
      
      mutateddf <- mutate(DEG_df, sig=ifelse(is.na(DEG_df$padj),"NA_", ifelse(abs(DEG_df$log2FoldChange)<l2FCthreshold & DEG_df$padj>pvaluethreshold, "Not_sign", ifelse(abs(DEG_df$log2FoldChange)<l2FCthreshold & DEG_df$padj<pvaluethreshold, "P-value", ifelse(abs(DEG_df$log2FoldChange)>l2FCthreshold & DEG_df$padj<pvaluethreshold, "P-value & log2FC", "log2FC")))))
      # Adjust the factor levels so that get consistent coloring
      mutateddf$sig<-factor(mutateddf$sig, levels = c("NA_", "Not_sign","P-value", "P-value & log2FC", "log2FC")) # "padj<0.05", ">2X but \n NotSig"))
      
      volc <- ggplot(mutateddf, aes(log2FoldChange, -log10(pvalue), gene=gene)) + # , na.rm=T volcanoplot with log2Foldchange versus pvalue
        geom_point(aes(col=sig), alpha=0.5) + #add points colored by significance  aes(col=sig:FCg)
        scale_color_manual(values=c("black", "grey", "blue", "green", "red")) + #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
        # scale_color_manual(values = c('Red', 'Pink', 'Blue', 'LightBlue', 'Green', 'LightGreen'))+ #, 'Black', 'Grey'
        theme_bw()+
        labs(title=paste("Volcano Plot of \n ",input$condition_comp, " vs ", input$condition_ref ))+ #input$condition_comp, input$condition_ref
        theme(plot.title = element_text(hjust = 0.5))+
        # ggtitle(paste("Volcano plot", )) + #e.g. 'Volcanoplot DESeq2'
        # scale_y_log10(labels=scales::scientific)+
        geom_hline(yintercept=-log10(pvaluethreshold), linetype="dashed", color = "green", size=1.5)+
        annotate("text", x=(max(mutateddf$log2FoldChange)/2), y= -log10(pvaluethreshold), label=paste("log10(1/Pvalue)", "threshold=", log10_pval), angle=90, hjust=0, size=6)+ # "(-)log10(Pvalue) of threshold is:"
        geom_vline(xintercept = c(-(l2FCthreshold),l2FCthreshold), linetype="dotted", 
                   color = "green", size=1.5)+
        annotate("text", x = c(-(l2FCthreshold),l2FCthreshold), y=quantile(-log10(mutateddf$pvalue), 0.9999, na.rm=T), label=paste("log2FC-threshold is: ", l2FCthreshold), angle=45, size=6)+ #y = max(-log10(mutateddf$pvalue[!is.na(mutateddf$padj)]))/2
        #theme(axis.ticks = element_line(size = 2))+ # This is just for axis ticks demarcking numbers
        theme(axis.line = element_line(size = 3, colour = "grey80"))+
        #theme(panel.border = )
        theme(panel.border = element_rect(linetype = "blank", fill = NA))+ #linetype = "dashed", fill = NA)
        theme(panel.grid.major = element_line(colour = "grey90", size=1))+
        theme(panel.grid.minor = element_line(colour = "grey90", size=0.8))+
        theme(axis.text.x = element_text(family="Arial", color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"))+ # size=rel(0.5)
        theme(axis.text.y = element_text(family="Arial", color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"))+
        theme(axis.title.x = element_text(family="Arial", face="plain", size=20))+ #, colour,face="bold"
        theme(axis.title.y = element_text(family="Arial", face="plain", size=20)) #, colour
      
      #volc+geom_text_repel(data=filter(mutateddf, padj<pvaluethreshold && (abs(log2FoldChange) > l2FCthreshold)), aes(label=gene))+geom_rug() #adding text for the top 20 genes
      
     volcs<- volc + geom_text_repel(
        data = subset(mutateddf, padj < pvaluethreshold & (abs(log2FoldChange) > l2FCthreshold)),
        aes(label = gene),
        size = 5,
        alpha = 0.9,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      ) + geom_rug()
      
     
     # position = dodge,
     # size = 4,
     # alpha = 0.9,
     # segment.size = .25,
     # segment.alpha = .8,
     # force = 1
     
      # geom_text_repel(
      #   data = subset(genes, padj < 0.05),
      #   aes(label = Gene),
      #   size = 5,
      #   box.padding = unit(0.35, "lines"),
      #   point.padding = unit(0.3, "lines")
      # )
      
      volcs
      
    })
    
  }
  )
  
  
  
  # output$plot5 <- renderPlot({
  #   
  #   library(ggplot2)
  #   library(cowplot) # for publication theme
  #   observeEvent(input$res_param_submit, ignoreInit=TRUE, {
  #     contrast_exp <- c(input$selectcontrast, input$condition_comp, input$condition_ref)
  #     # paste(input$selectcontrast, input$condition_comp, "vs", input$condition_ref, collapse = "_")->tempvar
  #     # print(tempvar)
  #     # print(paste(contrast_exp, collapse = "_"))
  #     if (input$condition_comp != input$condition_ref){
  #       library(DESeq2)
  #       res_table_unshrunken <- results(dds, contrast=contrast_exp, alpha = as.character(input$alpha_threshold))
  #       res_table_shrunk <- lfcShrink(dds, coef = paste(c(input$selectcontrast, input$condition_comp, "vs", input$condition_ref), collapse = "_"), res= res_table_unshrunken , type = 'apeglm')
  #       
  #     } else  {output$mtcars_table <- renderDataTable(data.frame(Error=c("Error! You didn't make correct comparison", "both conditions CAN'T be the same")))}
  #     # print("Error")
  #     
  #     # Turn the results into a dataframe
  #     DEG_OE_padj5<-res_tableOE_shrunk %>% data.frame()
  #     subset.data.frame(DEG_OE_padj5, padj < input$alpha_threshold)->DEG_OE_padj5
  #     # incorporate the Log2FC user supplied
  #     DEG_OE_padj5_lfc58<-subset.data.frame(DEG_OE_padj5, padj < padj.cutoff & abs(log2FoldChange)>lfc.cutoff)
  #     
  #     return()
  #     
  #   })
  # 
  # })
  # # # contrast_oe <- c("sampletype", "MOV10_overexpression", "control")
  # contrast_exp <- c(input$selectcontrast, input$selectcondition_comp, input$selectcondition_ref)
  # # #
  # res_table_unshrunken <- results(dds, contrast=contrast_exp, alpha = input$alpha_threshold)
  # # res_table_unshrunken <- results(dds, contrast=contrast_exp, alpha = 0.05)
  # # #
  # # # #res_tableOE <- lfcShrink(dds, contrast=contrast_oe, res=res_tableOE_unshrunken) # Here with contrast default is type='normal'; type = c("normal", "apeglm", "ashr")
  # # # # To get the coef names. Only available once "DESeq" is called
  # # # colnames(coef(dds))
  # res_table_shrunk <- lfcShrink(dds, coef = paste(input$selectcontrast, input$selectcontrastcomp, input$selectcontrastref, collapse = "_"), res= res_table_unshrunken , type = 'apeglm')
  # # data("mtcars")
  # # output$mtcars_table <- renderDataTable(summary(mtcars), options = list(dom = 't'))
  # #
  # # #output$mtcars_table <- renderDataTable(summary(res_table_shrunk), options = list(dom = 't'))
  rm(res_table_shrunk) # needs to be within the "server" function. otherwise each new session some functions that are dependent won't wait for other events
}


shinyApp(ui, server)



#################### References ###############

# https://shiny.rstudio.com/images/shiny-cheatsheet.pdf

### How to change dashboard body Tabitem title
# https://stackoverflow.com/questions/45176030/add-text-on-right-of-shinydashboard-header
### Not used in the end since I was able to just call another renderPlot call
# https://gist.github.com/wch/5436415/
### How to dynamically(really fixed and changes with respect to tab) generate options based on tab selected
# https://stackoverflow.com/questions/38979146/shiny-dashboardsidebar-change-when-different-tabitem-selected
# https://shiny.rstudio.com/reference/shiny/1.2.0/updateSelectInput.html
# https://stackoverflow.com/questions/48519138/selectinput-value-update-based-on-previous-selectinput-in-r-shiny # This is probably what works
# How to fix sidebar
# https://stackoverflow.com/questions/31013769/locking-r-shiny-dashboard-sidebar-shinydashboard


### How to capture output to a variable or text file
# https://www.r-bloggers.com/export-r-output-to-a-file/

## Rendering the number of genes that pass p-value threshold as threshold changes
# http://research.libd.org/recountWorkshop/compiled_vignette/SRP056604_main-results.html
