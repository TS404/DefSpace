
#
#
# # Activate packages
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("Biostrings",
#            "DECIPHER"))
#
suppressMessages(library(quiet=TRUE,"mclust"))
library(quiet=TRUE,"ggplot2")
library(quiet=TRUE,"Biostrings")
library(quiet=TRUE,"DECIPHER")
library(gplots)
if(!require("devtools")) install.packages("devtools")
devtools::install_github("bwlewis/rthreejs")
library(quiet=TRUE,"threejs")


# library(shinyRGL)
# Load local data

# For the server we may need to automatically discover this
data_dir="data"

SAPCA.cis <- readRDS(file.path(data_dir,"CisDef.reference.PCA.RDS"))
# SAPCA.tra <- readRDS(file.path(data_dir,"TransDef.reference.PCA.RDS"))
#
# view.cis <- readRDS(file.path(data_dir,"CisDef.viewangle.RDS"))
# view.tra <- readRDS(file.path(data_dir,"TransDef.viewangle.RDS"))
#
# BLOSUM40 <- readRDS(file.path(data_dir,"BLOSUM.RDS"))

# User-inputted data is a single sequenece
# example below
sequence <- "RTCESQSHRFKGPCSRDSNCATVCLTEGFSGGRCPWIPPRCFFFCTSPC" # cis defensin example

source("R/DefSpace functions.R")
source("R/Plots.R")
#plot_js_3Dclusters(SAPCA.cis, plotPCs = 1:3,col=named_cols, axeslabels = c("","",""))
# Plot cis (before use inputs any data) -----------------------------------------------

cis.colours <- factor(c("blue",            #1 Plant extreme
  "darkolivegreen4", #2 Plant main
  "grey",            #3 Intermed
  "purple1",         #4 Plant sex
  "orange",          #5 Plant his
  "maroon",          #6 Invert
  "red"))

tra.colours <- factor(c("blue",     #1 Theta
                        "red",      #2 Aalpha
                        "orange",   #3 Beta
                        "purple"))  #4

shinyServer(function(input, output) {
  output$cisplot <- renderScatterplotThree({
    plot_js_3Dclusters(SAPCA.cis, plotPCs = 1:3,col=cis.colours)
  })

  output$traplot <- renderScatterplotThree({
    plot_js_3Dclusters(SAPCA.tra, plotPCs = c(1,2,4),col=tra.colours)
  })

  output$CisPlot <- renderUI({
    scatterplotThreeOutput("cisplot")
   })


  output$TraPlot <- renderUI({
    scatterplotThreeOutput("traplot")
  })

})
