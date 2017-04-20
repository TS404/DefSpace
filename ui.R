shinyUI(
  verticalLayout(
    fluidPage(
      titlePanel("DefSpace: a tool for doing something cool"),

      wellPanel(
        div(class="header",
            p("brief description of what this tool is and what it can do for the user. What inputs does it require and what outputs does it produce")
        )
      ),

      radioButtons("alignment_type","Perform superfam alignment:",c("CIS only"="align.cis","Trans only"="align.tra","Find best"="align.best")),

      fluidRow(
        column(8,uiOutput("CisPlot"))
      ),

      fluidRow(
        column(8,uiOutput("TraPlot"))
      )

    )
  )
)
