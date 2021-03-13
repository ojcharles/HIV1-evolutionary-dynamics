library(shiny)
library(data.table)

ui <- fluidPage(
  titlePanel("Multiple file uploads"),
  sidebarLayout(
    sidebarPanel(
      fileInput("vcf_in",
                label="Upload VCFs here",
                multiple = TRUE)
    ),
    mainPanel(
      textOutput("count")
    )
  )
)

server <- function(input, output) {
  
  ### load all functions into environment
  funfiles = list.files("R",full.names = T)
  sapply(funfiles, function(x) source(x))
  
  
  
  mycsvs<-reactive({
    
    
    for(i in 1:length(input$files[,1])){
      lst[[i]] <- read.csv(input$files[[i, 'datapath']])
    }
    
    
    rbindlist(lapply(input$csvs$datapath, fread),
              use.names = TRUE, fill = TRUE)
  })
}

shinyApp(ui = ui, server = server)