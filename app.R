
library(BiocManager)
options(repos = BiocManager::repositories())

library('ggseqlogo');library('ggplot2');library('Biostrings'); library('DescTools'); library('parallel'); library("foreach"); library('doParallel'); library('cowplot');library('tools')   
library('MASS'); library('ggpubr'); library('randomcoloR'); library('DescTools'); library('DT') ; library('htmlwidgets'); library('shinyFiles'); library('BiocManager');library('readtext'); library('greekLetters'); library('universalmotif')
library('shinyURL')


source('gibbs_sampling.R')
source('seqs_TableData.R')
source('summary_motif_data.R')
source('freq_matrices.R')
source('calculate_MaxScoreInd_Of_SelectedSeq.R')
source('UVmotif_Obj.R')
source('mainCode.R')
source('scan_sequence.R')
source('mapToAlphabets.R')
source('alternate_string.R')
source('scan_Max_sequence.R')
source('CompileResults.R')
source('HSMotifDiscover.R')




ui <- fluidPage(
  
  titlePanel("HSMotifDiscover"),
  sidebarLayout(
    sidebarPanel(
      
      
        
        img(src = "main.png", height = 100, width = 280),
        #column(width=6, a(img(src="//rf.revolvermaps.com/h/m/a/0/ff0000/128/0/5o6vm4xrk24.png"), href="https://www.revolvermaps.com/livestats/5o6vm4xrk24/"))
        

    ),
    mainPanel(
      
        tabPanel("Global Visitors", a(img(src="http://rf.revolvermaps.com/h/m/a/0/ff0000/128/40/57n4o8karx4.png"), href="https://www.revolvermaps.com/livestats/57n4o8karx4/")),


    )
  ),
  
  tags$hr(),
  
  helpText("The app is also available in the form of an r-Package at",
  a("HSMotifDiscover",target="_blank",href="https://github.com/bioinfoDZ/HSMotifDiscover"), "and can be used on local server by downloading a standalone copy avaiable at", a("HSMotifDiscoverShiny",target="_blank",href="https://github.com/bioinfoDZ/HSMotifDiscover_ShinyApp")),
  
  
  h2(" Input Parameters"),
  
  
  tags$hr(),
  
  
  fluidRow(
    
    
    column(width=3,
           fileInput("file1", h3("Heparan Sulfate sequence file*"))),

  
  column(width=3, 
         h4("File Discription: "),
         helpText("File provides sequence data,", 
                  "Where, file consists of header and sequences information ",
                  a("click here",target="_blank",href="Ex_simulated_HSseqs_new_M1c_and_M3c_S100E.txt"),
                  "for sample data",
                  br(),
                  "To get desired file mapping ",
                  a("click Here",target="_blank",href="https://hsmotifdiscover.shinyapps.io/SeqMap_ShinyApp/")
                  
                  )),
  

  #a("click here for sample data",target="_blank",href="Ex_simulated_HSseqs_new_M1c_and_M3c_S100E.txt")
  
  
  ),
  
  tags$hr(),
  
  
  fluidRow(
    column(width=3,sliderInput(inputId = "bins",label = h3("Motif Length Range*:"),min = 3,max = 30, value = c(5,7))),
    column(width=3,
           h4("Input Discription: "),
           helpText("Length of the motif to be discovered in the Input sequences"))
  ),
  
  tags$hr(),

  
  fluidRow(
    
    column(width=3,
           fileInput("file2", h3("Sequence Character file "))),
    
    column(width=3, 
           h4("File Discription: "),
           helpText("File provides Dimer or Trimer characters and their group,", 
                    a("click here for sample file",target="_blank",href="Chars.txt")
           )),
  ),
  
  tags$hr(),
  
  
  
  
  fluidRow(
    
    column(width=3,
           fileInput("file3", h3("Sequence Weight file"))),
    
    column(width=3, 
           h4("File Discription: "),
           helpText("Weight for the sequences.", "If input in not given, all sequences will have equal weights" ),
           a("click here for sample file",target="_blank",href="motif1c_Weight.txt"))
  ),
  
  tags$hr(),
  
 
 fluidRow(
   column(width=3, numericInput("num2", h3("Weight Threshold" ), value = 0)),
   
   column(width=3, 
          h4("Input Discription: "),
          helpText("Sequences with weight greater than the threshold are used to discover motif. (Applicable only if weight file is uploaded) " ))
 ),
 
 tags$hr(),
 
  fluidRow(
    column(width=3, numericInput("num3", h3("Itrations"), value = "40000")),
    column(width=3, 
      h4("Input Discription: "),
      helpText("Number of itrations for Gibbs sampling" ))
 
  ),   

 tags$hr(),
 
 fluidRow(
   column(width=3, checkboxInput("RC_choice", h3("Reverse Complementry"), value = FALSE)),
   column(width=3, 
          h4("Input Discription: "),
          helpText("Tick if reverse complementry strand has to be considered while motif discovery (Note:Only in case of DNA sequences)" ))
   
 ),   

tags$hr(),

fluidRow(
  column(3,
           actionButton("start", "Start", class = "btn-primary")
           #actionButton("stop", "Stop", class = "btn-danger")
         )
      ),

tags$hr(),

########
h2(" Output"),
h6("(Visible after the completion of the programme) "),
# 

uiOutput('tabs'),
tags$hr(),

fluidRow(
  tags$script(src = "myscript.js")
),


)


server <- function(input, output, session) {
  
  
  ntext1 <- eventReactive(req(input$start), ignoreNULL = T,{

    MotifLenVec=input$bins[1]:input$bins[2]
    MotifInd=seq(1,length(MotifLenVec))
    
    
    out=list()
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...Depending on the computer hardware, input data size and number of iterations', value = 0, {
                   for (i in MotifInd) {
                     
                     #arg1 <- list(input_HSseq_file=input$file1$datapath, motif_length=MotifLenVec[i], out_dir=out_dir,charMapFile=input$file2$datapath, seq_weight_file= input$file3$datapath, affinity_threshold=input$num2, randmiseWeight=FALSE )
                     arg1 <- list(input_HSseq_file=input$file1$datapath, motifLenVec=MotifLenVec[i], charGrpFile=input$file2$datapath, seq_weight_file= input$file3$datapath,  affinity_threshold=input$num2, itr=input$num3, if_revCompStrand=input$RC_choice)
                     ## Filter out NULL values
                     ## define Filter fun which returns TRUE if you consider the value valid
                     is_valid <- function(.) !is.null(.) && . != ""
                     ## Filter out non valid args
                     arg1 <- Filter(is_valid, arg1)
                     ## Use do.call to run your fun (non valid args are not set)
                     print( arg1)
    
                     out[[i]]=isolate(do.call(HSMotifDiscover,arg1))
              
                     #out[[i]]=isolate(readRDS('OutTest.RDS')[[i]])
                     #out[[i]][[1]]$arg=arg1
                     
                     incProgress(i/length(MotifInd))
                     
             
    
    }
    
    })
    
    #saveRDS(object = out, file = 'OutTest.RDS')
    
    shinyURL.server()
    
    return(out)
  
  })
  
  
  
  
  
  
  output$tabs <- renderUI({
    #req(input$slider)
    MotifLenVec=input$bins[1]:input$bins[2]
    
    Tabs <- lapply(1:length(MotifLenVec), function(i) {
    #Tabs <- lapply(1:1, function(i) {
      
      {
        
        output[[paste0('plot_',i)]] <- renderPlot({
          return(ntext1()[[i]][[1]]$MotifLogo)
        })
        
        
        output[[paste0('mytable1_',i)]] <- renderDT( ## <<
          ntext1()[[i]][[1]]$resultsTable_df, # data
          class = "display nowrap compact", # style
          filter = "top", # location of column filters
          options = list(  # options
            scrollX = TRUE # allow user to scroll wide tables horizontally
          ))
        
        #output$mytable1 <- DT::renderDataTable({ ## <<
        #  DT::datatable(ntext1()[[1]]$df, options = list(orderClasses = TRUE)) # data
        #})
        
        
        output[[paste0('downloadPlot_',i)]] <- downloadHandler(
          filename = function() { paste0('motif_',i,'.png') },
          content = function(file) {
            ggsave(file,ntext1()[[i]][[1]]$MotifLogo)
          }
        )
        
        
        output[[paste0('downloadDataM_',i)]] <- downloadHandler(
          filename = function() { paste0('Motif_PSSM_',i,'.tsv' )},
          content = function(file) {
            write.table(ntext1()[[i]][[1]]$PSSM, file, quote = FALSE, sep = '\t', row.names = TRUE, col.names = FALSE)
          }
        )
        
        output[[paste0('downloadData_',i)]] <- downloadHandler(
          filename = function() { paste0('resultsMatrix_',i,'.tsv') },
          content = function(file) {
            write.table(ntext1()[[i]][[1]]$resultsTable_df, file, quote = FALSE, sep = '\t', row.names = FALSE)
          }
        )
        
      }
      
      output[[paste0('text_',i)]] <- renderText({
        str1=paste0("Information content: ",ntext1()[[i]][[1]]$IC)
        str2=paste0('KL Divergerce between motif and background: ', signif(ntext1()[[i]][[1]]$MAX_motifScore,3))
        str3=paste0('P-Value of the Motif: ', signif(mean(ntext1()[[i]][[1]]$resultsTable_df$P_value,5)))
        HTML(paste(str1, str2, str3, sep = '<br/>'))
      })
      
      
      
      tabPanel(title = paste0("Motif Length", MotifLenVec[i]), # <<
               tabsetPanel(
                 
                 tabPanel("Discovered Motif", plotOutput(paste0("plot_",i)),
                          downloadButton(paste0('downloadPlot_',i), 'Download Motif'),
                          downloadButton(paste0('downloadDataM_',i), 'Download Motif PSSM'),
                 ),
                 #tabPanel("Results Table", DT::dataTableOutput("mytable1"))
                 tabPanel("Results Table", div(dataTableOutput(paste0("mytable1_",i))),
                          downloadButton(paste0('downloadData_',i), 'Download Data'),
                 ),
                 tabPanel("Motif Scores", htmlOutput(paste0('text_',i)))
                 
  
               )  ## >>
               
      )
      
    })
    do.call(tabsetPanel, Tabs)
  })
  
  
  shinyURL.ui()
  
}

shinyApp(ui, server)
