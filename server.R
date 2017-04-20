source("functions.R")
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
# if(!require("devtools")) install.packages("devtools")
# devtools::install_github("bwlewis/rthreejs")
library(quiet=TRUE,"threejs")


plot_js_3Dclusters <- function(SAPCA,
                               plotPCs    = 1:3,
                               col        = "cluster",
                               radius     = 1,
                               labels     = NULL,
                               axeslabels = "default"){


  data <- SAPCA$seq.space.PCA$coordinates


  # Calculate radius size
  rad <- (range(SAPCA$seq.space.PCA$coordinates[,plotPCs])[2]-range(SAPCA$seq.space.PCA$coordinates[,plotPCs])[1])/100
  rad <- rad*radius

  if (all(axeslabels=="default")){
    axes <- paste("PC",plotPCs,sep="")
  }else{
    axes <- axeslabels
  }
  if (is.null(axeslabels)){
    axes <- c("","","")
  }

  print(paste(c("Making 3d scatterplot with data of size ",dim(data[,plotPCs]))))

  sp <- scatterplot3js(data[,plotPCs],
                 color = col2hex(col[SAPCA$seq.space.clusters$classification]),
                 renderer = "auto"
  )

  print("Done making 3d scatterplot")
  sp
}

plot_pca <- function(refdata,sequence,align_type,BLOSUM40){

  newseq.cis <- seq.MSA.add(refdata$cis,sequence,"cis-Defensins",BLOSUM40)
  newseq.tra <- seq.MSA.add(refdata$tra,sequence,"trans-Defensins",BLOSUM40)


  if ( align_type == "align.best" ){
    align_cis <- (newseq.cis$aln.hit.score>=newseq.tra$aln.hit.score)
  } else {
    align_cis <- (align_type == "align.cis")
  }

  if (align_cis){
      pca_data <- refdata$cis
      newseq.match <- newseq.cis
      plotPCs <- c(1,2,3)
      plot.colours <- cis.colours
  } else {
      pca_data <- refdata$tra
      newseq.match <- newseq.tra
      plotPCs <- c(1,2,4)
      plot.colours <- tra.colours
  }

  print(paste(c("Plotting cis ",align_cis)))

  newseq.r  <- seq.rotate    (pca_data, newseq.match)
  newseq.c  <- seq.clust.add (pca_data, newseq.r)
  SAPCA.add <- seq.SAPCA.add (pca_data, newseq.r, newseq.c)

  r = c(2,rep(0.3,nrow(SAPCA.add$numerical.alignment$MSA)-1))

  plot_js_3Dclusters(pca_data, plotPCs = 1:3,col=plot.colours, radius = r)

  # if(is.cis){
  #   SAPCA.match  <- SAPCA.cis
  #   newseq.match <- seq.MSA.add(SAPCA.cis,sequence,"cis-Defensins")
  #   match        <- "cis-Defensin"
  # }else if(is.tra){
  #   SAPCA.match  <- SAPCA.tra
  #   newseq.match <- seq.MSA.add(SAPCA.tra,sequence,"trans-Defensins")
  #   match        <- "trans-Defensin"
  # }else{
  #   newseq.cis <- seq.MSA.add(SAPCA.cis,sequence,"cis-Defensins")
  #   newseq.tra <- seq.MSA.add(SAPCA.tra,sequence,"trans-Defensins")
  #   # Which is the better match
  #   if(newseq.cis$aln.hit.score>=newseq.tra$aln.hit.score){
  #     SAPCA.match  <- SAPCA.cis
  #     newseq.match <- newseq.cis
  #     view         <- view.cis
  #     plotPCs      <- c(1,2,3)
  #     match        <- "cis-Defensin"
  #   }else{
  #     SAPCA.match  <- SAPCA.tra
  #     newseq.match <- newseq.tra
  #     view         <- view.tra
  #     plotPCs      <- c(1,2,4)
  #     match        <- "trans-Defensin"
  #   }
  # }
}

# Load local data

# For the server we may need to automatically discover this
data_dir="data"

#SAPCA.cis <- readRDS(file.path(data_dir,"CisDef.reference.PCA.RDS"))
# SAPCA.tra <- readRDS(file.path(data_dir,"TransDef.reference.PCA.RDS"))
#
# view.cis <- readRDS(file.path(data_dir,"CisDef.viewangle.RDS"))
# view.tra <- readRDS(file.path(data_dir,"TransDef.viewangle.RDS"))
#
# BLOSUM40 <- readRDS(file.path(data_dir,"BLOSUM.RDS"))

# User-inputted data is a single sequenece
# example below
# sequence <- "RTCESQSHRFKGPCSRDSNCATVCLTEGFSGGRCPWIPPRCFFFCTSPC" # cis defensin example

# source("R/DefSpace functions.R")
# source("R/Plots.R")
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

data_dir <- "data"


shinyServer(function(input, output) {


  #BLOSUM40
  blosum <- reactive({
    cat("Reading BLOSUM")
    readRDS(file.path(data_dir,"BLOSUM.RDS"))
  })


  reference_data <- reactive({
    cat("Reading reference data")
    SAPCA.cis <- readRDS(file.path(data_dir,"CisDef.reference.PCA.RDS"))
    SAPCA.tra <- readRDS(file.path(data_dir,"TransDef.reference.PCA.RDS"))
    SAPCA <- list(cis = SAPCA.cis, tra = SAPCA.tra)
    SAPCA
  })

  output$pcaplot <- renderScatterplotThree({
    plot_pca(reference_data(),input$user_sequence,input$alignment_type,blosum())
  })

  output$PCAScatterPlot <- renderUI({
    scatterplotThreeOutput("pcaplot")
   })

})
