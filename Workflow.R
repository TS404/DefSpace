
# Setup ----------------------------------------------------------------------

# Activate packages
source("https://bioconductor.org/biocLite.R")
biocLite(c("Biostrings",
           "DECIPHER"))

library(quiet=TRUE,"rgl")
library(quiet=TRUE,"ape")
library(quiet=TRUE,"cluster")
library(quiet=TRUE,"mclust")
library(quiet=TRUE,"gplots")
library(quiet=TRUE,"ggplot2")
library(quiet=TRUE,"Biostrings")
library(quiet=TRUE,"DECIPHER")




# Import necessary files ----------------------------------------------------------------------

folder    <- "C:\\Users\\T\\OneDrive\\1-Scripts\\GitHub\\DefSpace"
setwd(folder)

SAPCA.cis <- readRDS("data\\CisDef.reference.PCA.RDS")
SAPCA.tra <- readRDS("data\\TransDef.reference.PCA.RDS")

view.cis <- readRDS("data\\CisDef.viewangle.RDS")
view.tra <- readRDS("data\\TransDef.viewangle.RDS")

BLOSUM40 <- readRDS("data\\BLOSUM.RDS")

# User-inputted data is a single sequenece
# example below
sequence <- "RTCESQSHRFKGPCSRDSNCATVCLTEGFSGGRCPWIPPRCFFFCTSPC" # cis defensin example

# Plot cis (before use inputs any data) -----------------------------------------------

colours<-palette(c("blue",            #1 Plant extreme
                   "darkolivegreen4", #2 Plant main
                   "grey",            #3 Intermed
                   "purple1",         #4 Plant sex
                   "orange",          #5 Plant his
                   "maroon",          #6 Invert
                   "red"))            #7 Tox
plot_3Dclusters(SAPCA.cis, plotPCs = 1:3, axeslabels = c("","",""))
rgl::par3d(view.cis)



# Plot trans (before user inputs any data)  ------------------------------------------------

colours<-palette(c("blue",     #1 Theta
                   "red",      #2 Aalpha
                   "orange",   #3 Beta
                   "purple"))  #4
plot_3Dclusters(SAPCA.tra, plotPCs = c(1,2,4), axeslabels = c("","",""))
rgl::par3d(view.tra)


# Assign sequence class  ----------------------------------------------------------------------

# Checkbox for whether to check alignments to both superfams
if(is.cis){
  SAPCA.match  <- SAPCA.cis
  newseq.match <- seq.MSA.add(SAPCA.cis,sequence,"cis-Defensins")
  match        <- "cis-Defensin"
}else if(is.tra){
  SAPCA.match  <- SAPCA.tra
  newseq.match <- seq.MSA.add(SAPCA.tra,sequence,"trans-Defensins")
  match        <- "trans-Defensin"
}else{
  newseq.cis <- seq.MSA.add(SAPCA.cis,sequence,"cis-Defensins")
  newseq.tra <- seq.MSA.add(SAPCA.tra,sequence,"trans-Defensins")
  # Which is the better match
  if(newseq.cis$aln.hit.score>=newseq.tra$aln.hit.score){
    SAPCA.match  <- SAPCA.cis
    newseq.match <- newseq.cis
    view         <- view.cis
    plotPCs      <- c(1,2,3)
    match        <- "cis-Defensin"
  }else{
    SAPCA.match  <- SAPCA.tra
    newseq.match <- newseq.tra
    view         <- view.tra
    plotPCs      <- c(1,2,4)
    match        <- "trans-Defensin"
  }
}

# Rotation, cluster assignment and completion for whichever SAPCA matches sequence
newseq.r  <- seq.rotate    (SAPCA.match, newseq.match)
newseq.c  <- seq.clust.add (SAPCA.match, newseq.r)
SAPCA.add <- seq.SAPCA.add (SAPCA.match, newseq.r, newseq.c)


# Plot combined SAPCA -------------------------------------------------------------------
if(match=="cis-Defensin"){
  colours<-palette(c("blue",            #1 Plant extreme
                     "darkolivegreen4", #2 Plant main
                     "grey",            #3 Intermed
                     "purple1",         #4 Plant sex
                     "orange",          #5 Plant his
                     "maroon",          #6 Invert
                     "red"))            #7 Tox
}
if(match=="trans-Defensin"){
  colours<-palette(c("blue",     #1 Theta
                     "red",      #2 Aalpha
                     "orange",   #3 Beta
                     "purple"))  #4
}

plot_3Dclusters(SAPCA.add,
                plotPCs = plotPCs,
                labels = "query",
                radius = c(2,rep(0.3,nrow(SAPCA.add$numerical.alignment$MSA)-1)))
rgl::par3d(view)

# Info
as.fasta(SAPCA.add$numerical.alignment$MSA[rownames(closest(SAPCA.add,"query")),],
         decolgap = TRUE)
paste(sep="","Sequence appears to be a ",match," (max similarity = ",percent(newseq.match$aln.hit.score),")")


# Plot comparison histograms ------------------------------------------------------------------
aln.cistra <- c(newseq.tra$aln.all.score,newseq.cis$aln.all.score)
hist(newseq.cis$aln.all.score,
     xlim   = c(min(aln.cistra),1),
     ylim   = c(0,10),
     breaks = seq(min(range(aln.cistra)-0.5),
                  max(range(aln.cistra)+0.5),
                  0.02),
     col    = rgb(0,0.5,0,0.5), # green
     freq   = FALSE,
     xlab   = "similarity",
     main   = "")
hist(newseq.tra$aln.all.score,
     breaks = seq(min(range(aln.cistra)-0.5),
                  max(range(aln.cistra)+0.5),
                  0.02),
     col    = rgb(0,0,0.8,0.5), # blue
     freq   = FALSE,
     add    = TRUE)
box()



# End  ----------------------------------------------------------------------




