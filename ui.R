shinyUI(
  verticalLayout(
    fluidPage(
      titlePanel("DefSpace: a tool for doing something cool"),

      wellPanel(
        div(class="header",
            p("brief description of what this tool is and what it can do for the user. What inputs does it require and what outputs does it produce")
        )
      ),

      textAreaInput("user_sequence", "Your Sequence", "RTCESQSHRFKGPCSRDSNCATVCLTEGFSGGRCPWIPPRCFFFCTSPC", width = "600px"),

     radioButtons("alignment_type","Perform superfam alignment:",c("Find best"="align.best","CIS only"="align.cis","Trans only"="align.tra")),

      uiOutput("PCAScatterPlot")

    )
  )
)
