plot_3Dclusters(SAPCA.cis, plotPCs = 1:3, axeslabels = c("","",""))
rgl::par3d(view.cis)



# Plot trans (before user inputs any data)  ------------------------------------------------


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
     xlab   = "Similarity",
     main   = "")
hist(newseq.tra$aln.all.score,
     breaks = seq(min(range(aln.cistra)-0.5),
                  max(range(aln.cistra)+0.5),
                  0.02),
     col    = rgb(0,0,0.8,0.5), # blue
     freq   = FALSE,
     add    = TRUE)
box()



# Info  ----------------------------------------------------------------------

as.fasta(SAPCA.add$numerical.alignment$MSA[rownames(closest(SAPCA.add,"query")),],
         decolgap = TRUE)

paste(sep="","Sequence appears to be a ",match," (max similarity = ",percent(newseq.match$aln.hit.score),")")




