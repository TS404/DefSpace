shinyUI(
  fluidPage(
    titlePanel("DefSpace"),

    fluidRow(
      column(8,uiOutput("CisPlot"))
    ),

    fluidRow(
      column(8,uiOutput("TraPlot"))
    )


  )
)
