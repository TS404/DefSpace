
numericise_MSA <- function(MSA,
                           res.prop,
                           cys){
  seq.names <- rownames(MSA)
  aln.len   <- ncol(MSA)
  res.props <- colnames(res.prop)
  res.avail <- row.names(res.prop)
  # Numericise MSA based on res.prop
  MSA.num.tall <- res.prop[t(MSA),]
  # Name data types
  rownames(MSA.num.tall) <- NULL
  sequence     <- rep(x = seq.names, each  = aln.len)
  residue      <- rep(x = 1:aln.len, times = length(seq.names))
  MSA.num.tall <- cbind(sequence, residue, MSA.num.tall)
  # Stack data into list of matrices
  MSA.num.stack <- NULL
  for (x in 1:length(res.props)) {
    col.names <- paste(1:aln.len,
                       rep(res.props[x],aln.len),
                       sep = ".")
    MSA.num.stack[[res.props[x]]] <- matrix(MSA.num.tall[,x+2],
                                            ncol     = aln.len,
                                            byrow    = TRUE,
                                            dimnames = list(seq.names,
                                                            col.names))
  }
  # Also reflow into single wide matrix
  MSA.num.wide <- MSA.num.stack[[1]]
  for (x in 2:length(res.props)) {
    MSA.num.wide <- cbind(MSA.num.wide, MSA.num.stack[[res.props[x]]])
  }
  
  
  ############################.
  # Scaling by property type #
  ############################.
  
  # Take means and variances of each property type
  prop.means <- NULL
  prop.vars  <- NULL
  for (x in 1:length(res.props)) {
    prop.means[x] <- mean(MSA.num.stack[[x]],na.rm=1)
    prop.vars[x]  <- var(tidyr::gather(data.frame(MSA.num.stack[[x]]))[2],na.rm=1)
  }
  names(prop.means) <- res.props
  names(prop.vars)  <- res.props
  
  # Scale numericised MSA to prop.means and prop.vars
  MSA.scale.stack <- NULL
  for (x in 1:length(res.props)) {
    MSA.scale.stack[[res.props[x]]] <- (MSA.num.stack[[res.props[x]]]- prop.means[x]) /
      sqrt(prop.vars[x])
  }
  
  # Replace gaps (currently "NA") with column average
  # Create na.colmean function
  na.colmean<-function(x){
    x[is.na(x)] <- mean(as.matrix(x),na.rm = 1)
    x
  }
  # For each property of MSA.num.stack, apply na.colmean function to each matrix comlumn
  for (x in 1:length(res.props)) {
    MSA.scale.stack[[x]] <- apply(MSA.scale.stack[[x]],2,na.colmean)
  }
  
  # Also reflow into singe wide matrix for PCA
  MSA.scale.wide <- MSA.scale.stack[[1]]
  for (x in 2:length(res.props)) {
    MSA.scale.wide <- cbind(MSA.scale.wide, MSA.scale.stack[[x]])
  }
  
  
  ##################.
  # Alignment list #
  ##################.
  
  numerical.alignment <- list(MSA             = MSA,
                              res.prop        = res.prop,
                              MSA.num.stack   = MSA.num.stack,
                              MSA.num.wide    = MSA.num.wide,
                              MSA.scale.stack = MSA.scale.stack,
                              MSA.scale.wide  = MSA.scale.wide,
                              prop.means      = prop.means,
                              prop.vars       = prop.vars,
                              seq.names       = seq.names,
                              aln.len         = aln.len)
  numerical.alignment
}



closest <- function (SAPCA,
                     sequence,
                     PC = 1:3,
                     n  = 10){
  
  coords     <- SAPCA$seq.space.PCA$coordinates
  centre     <- coords[sequence,PC]
  centre.m   <- matrix(rep(centre,nrow(coords)),
                       nrow  = nrow(coords),
                       byrow = TRUE)
  distances  <- SAPCA$seq.space.PCA$coordinates[,PC]-centre.m
  rootsquare <- sqrt(rowSums(distances^2))
  
  sorted     <- as.matrix(rootsquare[order(rootsquare)])
  colnames(sorted) <- "distance"
  
  head(sorted,n)
}



read.MSA <- function(MSA){
  # Load sequence MSA
  # if a matrix, can be used straight away
  # if raw fasta file, use seqinr to convert to data frame
  if (!is.matrix(MSA)){
    MSA       <- data.frame(seqinr::read.fasta(MSA,set.attributes=FALSE))
  }
  # if a data frame, convert to matrix
  if (is.data.frame(MSA)){
    MSA       <- as.matrix(t(toupper(as.matrix(MSA))))
  }
  
  MSA
}



as.fasta <- function(matrix,degap=FALSE,decolgap=FALSE,write=FALSE,print=FALSE,name=NULL){
  
  # Remove empty columns
  if(decolgap){
    matrix<-matrix[,colMeans(matrix=="-")!=1]
  }
  
  # Convert alignment matrix to list of strings
  names <- paste(">",row.names(matrix),sep="")
  seqs  <- do.call("paste",c(data.frame(matrix),sep=""))
  
  # If just one sequence, this is how to name it
  if(is.null(dim(matrix))){
    names <- ">sequence"
    if(!is.null(name)){
      names <- paste(">",name,sep="")
    }
    seqs  <- paste(matrix,collapse="")
  }
  
  # Degap sequences
  if (degap){
    seqs <- gsub("-","",seqs)
  }
  
  # Interleave names and sequences
  ord1 <- 2*(1:length(names))-1
  ord2 <- 2*(1:length(seqs))
  
  # Output
  if (print==TRUE){
    paste0(c(names,seqs)[order(c(ord1,ord2))], sep = "<br/>")
  }else{
    if (write==FALSE){
      cat(c(names,seqs)[order(c(ord1,ord2))], sep = "\n")
    }
    else{
      if (!grepl(".fa",write,ignore.case=TRUE)){
        write<-paste(write,".fa",sep="")
      }
      cat(c(names,seqs)[order(c(ord1,ord2))], sep = "\n", file = write)
    }
  }
}



as.AAstring<-function(string, degap=FALSE){
  string <- paste(string,collapse="")
  if(degap==TRUE){
    string<-gsub("-","",string)
  }
  output <- Biostrings::AAString(string)
  output
}


as.AAstringSet<-function(MSA, degap=FALSE){
  MSA <- apply(MSA,1,paste,collapse="")
  if(degap==TRUE){
    MSA<-gsub("-","",MSA)
  }
  output <- Biostrings::AAStringSet(MSA)
  output
}



seq.MSA.add <- function(SAPCA,sequence,SAPCAname=NULL,smatrix=BLOSUM40){
  sequence <- casefold(sequence,upper=TRUE)
  MSA   <- SAPCA$numerical.alignment$MSA
  MSA2  <- as.AAstringSet(MSA,degap = TRUE)
  seqs  <- nrow(MSA)
  seq   <- as.AAstring(sequence, degap=FALSE)
  seq.d <- as.AAstring(sequence, degap=TRUE)
  BLOSUM40 <- smatrix
  
  aln.all <- Biostrings::pairwiseAlignment(MSA2,
                                           seq.d,
                                           substitutionMatrix = BLOSUM40,
                                           gapOpening   = 0,
                                           gapExtension = 4,
                                           scoreOnly    = TRUE)
  
  # Max possible similarity score
  aln.limit <- Biostrings::pairwiseAlignment(seq.d,
                                             seq.d,
                                             substitutionMatrix = BLOSUM40,
                                             gapOpening   = 0,
                                             gapExtension = 4,
                                             scoreOnly    = TRUE)
  
  # Similarity score as percentage of max
  aln.hit.score <- max(aln.all)/aln.limit
  
  # The sequence of the best matching member of the database
  aln.hit.num  <- which(aln.all==max(aln.all))[1]
  aln.hit.name <- SAPCA$numerical.alignment$seq.names[aln.hit.num]
  aln.hit.seq  <- paste(as.AAstring(MSA[aln.hit.num,],degap = 1))
  
  # Use "*" to indicate gaps in the best reference sequence (aln.hit)
  aln.hit <- gsub("-","*",as.AAstring(MSA[aln.hit.num,]))
  
  # Use "#" to anchor ends of sequences so that they are not trimmed
  aln.hit.anchor <- paste0("########",aln.hit,"########")
  seq.d.anchor   <- paste0("########",seq.d,  "########")
  
  aln.add <- Biostrings::pairwiseAlignment(aln.hit.anchor,
                                           seq.d.anchor,
                                           substitutionMatrix = BLOSUM40,
                                           gapOpening   = 0,
                                           gapExtension = 4)
  
  # Has the new sequence introduced exrta gaps into the hit sequence alignement?
  aln.hit.orig         <- as.AAstring(MSA[aln.hit.num,])
  aln.hit.new          <- gsub("#","",Biostrings::pattern(aln.add))
  aln.hit.seq.aln.orig <- unlist(strsplit(as.character(aln.hit.orig),""))
  aln.hit.seq.aln.new  <- unlist(strsplit(as.character(aln.hit.new),""))
  
  if(as.character(aln.hit.orig)!=as.character(aln.hit.new)){
    print(paste(sum(aln.hit.seq.aln.new=="-"),
                "residues of the new sequence were not alignable to the",
                SAPCAname,
                "reference MSA so were ignored"))
  }
  
  # Alignment addition as matrix
  aln.add.mat <- rbind(unlist(strsplit(as.character(gsub("#","",Biostrings::pattern(aln.add))),"")),
                       unlist(strsplit(as.character(gsub("#","",Biostrings::subject(aln.add))),"")))
  
  # Unmathcable resiues removed from aligned sequence
  aln.add2 <- aln.add.mat[2,][aln.add.mat[1,]!="-"]
  aln.add3 <- paste(as.AAstring(aln.add2))

  # Unalignable residues ignored from the middle of the sequence
  
  seq.alignable <- Biostrings::pairwiseAlignment(seq.d,
                                                 as.AAstring(aln.add3, degap=TRUE),
                                                 substitutionMatrix = BLOSUM40,
                                                 gapOpening   = 0,
                                                 gapExtension = 4)
  # res.mid.removed <- sum(unlist(strsplit(as.character(subject(seq.alignable)),""))=="-")
  # 
  # # Reidues removed from the start or finish template (aln.hit) during alignment
  # res.terminal.discrep <- nchar(gsub("[*]","",aln.hit)) - nchar(gsub("[*]","",aln.hit.new)) + res.mid.removed
  # 
  # aln.hit.d <- unlist(strsplit(as.character(gsub("[*]","",aln.hit)),""))
  # aln.hit.new.d <- unlist(strsplit(as.character(gsub("[*]","",aln.hit.new)),""))
  # 
  # aln.hit.d[x + 1:length(aln.hit.new.d)]
  # 
  # displacement <- NULL
  # for(x in 0:(res.terminal.discrep)){
  #   displacement <- append(displacement,mean(aln.hit.d[x + 1:length(aln.hit.new.d)]==aln.hit.new.d))
  # }
  # res.lead.missing <- which.max(displacement)-1
  # res.tail.missing <- res.terminal.discrep-which.max(displacement)+1
  # 
  # # Gaps in the hit sequence (original and newly aligned)
  # gaps.orig       <- unlist(strsplit(as.character(aln.hit.orig),"[A-Z]"))
  # gaps.count.orig <- nchar(gaps.orig)
  # gaps.lead.orig  <- gaps.count.orig[1]
  # gaps.tail.orig  <- gaps.count.orig[length(gaps.count.orig)]
  # 
  # gaps.new        <- unlist(strsplit(as.character(gsub("-","",aln.hit.new)),"[A-Z]"))
  # gaps.count.new  <- nchar(gaps.new)
  # gaps.lead.new   <- gaps.count.new[1]
  # gaps.tail.new   <- gaps.count.new[length(gaps.count.new)]
  # 
  # if(length(gaps.count.orig)>length(gaps.count.new)){
  #   gaps.count.new <- append(gaps.count.new,0)
  # }
  # 
  # gaps.discrep     <- rbind(gaps.count.orig,
  #                            c(rep(0,res.lead.missing),gaps.count.new))
  # gaps.discrep.num <- gaps.discrep[1,]-gaps.discrep[2,]
  # gaps.lead.add    <- gaps.discrep.num[1] + res.lead.missing
  # gaps.tail.add    <- gaps.discrep.num[length(gaps.discrep.num)] + res.tail.missing

  # # Final aligned sequence to add (with gps at beginning and end to fit)
  # query <- c(rep("-",gaps.lead.add),
  #            aln.add2,
  #            rep("-",gaps.tail.add))
  
  query <- aln.add2
  
  aln.final <- rbind(query,MSA)
  output    <- list(MSA             = aln.final,
                    aln.hit.name    = aln.hit.name,
                    aln.hit.seq     = aln.hit.seq,
                    aln.hit.score   = aln.hit.score,
                    aln.all.score   = aln.all/aln.limit,
                    seq.unalignable = seq.alignable)
  output
}



seq.rotate <- function(SAPCA,newseq){
  
  res.props  <- colnames(SAPCA$numerical.alignment$res.prop)
  prop.means <- SAPCA$numerical.alignment$prop.means
  prop.vars  <- SAPCA$numerical.alignment$prop.vars
  # Align new sequence with MSA
  
  # Numericise new sequence
  newseq.num <- numericise_MSA(MSA      = newseq$MSA[c("query",newseq$aln.hit.name),],
                               res.prop = SAPCA$numerical.alignment$res.prop)
  
  # Scale new sequnce using same scaling as SAPCA (gaps as "NA")
  newseq.scale.stack <- NULL
  for (x in 1:length(res.props)) {
    newseq.scale.stack[[res.props[x]]] <- (newseq.num$MSA.num.stack[[res.props[x]]]- prop.means[x]) /
      sqrt(prop.vars[x])
  }
  
  # Reflow into single wide matrix
  newseq.scale.wide <- newseq.scale.stack[[1]]
  for (x in 2:length(res.props)) {
    newseq.scale.wide <- cbind(newseq.scale.wide, newseq.scale.stack[[res.props[x]]])
  }
  
  # Replace gaps (currently "NA") with column average of the scaled SAPCA
  # Create na.colmean function
  gapvalues  <- colMeans(SAPCA$numerical.alignment$MSA.scale.wide)
  na.replace <- function(x,y){
    x[is.na(x)] <- y
    x
  }
  newseq.scale.wide.g <- NULL
  for(i in 1:ncol(newseq.scale.wide)){
    newseq.scale.wide.g <- cbind(newseq.scale.wide.g,
                                 na.replace(newseq.scale.wide[,i],gapvalues[i]))
  }
  
  # Rotate scaled sequence into same space as SAPCA
  newseq.rot <- scale(newseq.scale.wide.g,
                      SAPCA$seq.space.PCA$centre,
                      SAPCA$seq.space.PCA$scale) %*% SAPCA$seq.space.PCA$loadings
  
  # Output
  output <- list(seq       = newseq,
                 seq.num   = newseq.num$MSA.num.wide,
                 seq.scale = newseq.scale.wide.g,
                 seq.rot   = newseq.rot)
  output
}



seq.clust.add <- function(SAPCA,newseq.r){
  
  SAPCA.c  <- mclustrev(SAPCA)
  newseq.c <- mclust::predict.Mclust(SAPCA.c,newseq.r$seq.rot[,SAPCA$call$clusterPCs])
  newseq.c
}



mclustrev <- function(SAPCA){
  output <- SAPCA$seq.space.clusters$other
  
  output$classification <- SAPCA$seq.space.clusters$classification
  output$G              <- SAPCA$seq.space.clusters$optimal
  output$z              <- SAPCA$seq.space.clusters$likelihoods
  
  class(output) <- "Mclust"
  output
}



seq.SAPCA.add <- function (SAPCA,newseq.r,newseq.c){
  output <- SAPCA
  output$numerical.alignment$seq.names      <- rownames(newseq.r$seq$MSA)
  output$numerical.alignment$MSA            <- newseq.r$seq$MSA
  output$numerical.alignment$MSA.num.wide   <- rbind(newseq.r$seq.num[1,],
                                                     SAPCA$numerical.alignment$MSA.num.wide)
  output$numerical.alignment$MSA.scale.wide <- rbind(newseq.r$seq.scale[1,],
                                                     SAPCA$numerical.alignment$MSA.scale.wide)
  output$numerical.alignment$MSA.num.stack  <- NULL
  output$numerical.alignment$MSA.scale.stack<- NULL
  
  output$seq.space.PCA$coordinates          <- rbind(newseq.r$seq.rot[1,],
                                                     SAPCA$seq.space.PCA$coordinates)
  output$seq.space.clusters$likelihoods     <- rbind(newseq.c$z[1,],
                                                     SAPCA$seq.space.clusters$likelihoods)
  output$seq.space.clusters$classification  <- c(newseq.c$classification[1],
                                                 SAPCA$seq.space.clusters$classification)
  
  rownames(output$numerical.alignment$MSA.num.wide)[1]   <- "query"
  rownames(output$numerical.alignment$MSA.scale.wide)[1] <- "query"
  rownames(output$seq.space.PCA$coordinates)[1]          <- "query"
  rownames(output$seq.space.clusters$likelihoods)[1]     <- "query"
  
  output
}



seq.add.full <- function (SAPCA,sequence,SAPCAname=NULL){
  
  newseq   <- seq.MSA.add(SAPCA,sequence,SAPCAname)
  newseq.r <- seq.rotate(SAPCA,newseq)
  newseq.c <- seq.clust.add(SAPCA,newseq.r)
  SAPCA2   <- seq.SAPCA.add(SAPCA,newseq.r,newseq.c)
  
  SAPCA2
}



percent <- function(x, digits = 1, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}



plot_3Dclusters <- function(SAPCA,
                            plotPCs    = 1:3,
                            col        = "cluster",
                            radius     = 1,
                            labels     = NULL,
                            write      = FALSE,
                            axeslabels = "default"){
  if (!is.null(SAPCA$seq.space.PCA$coordinates)){
    data <- SAPCA$seq.space.PCA$coordinates
  }else{
    data <- SAPCA    
  }
  
  if (all(col=="cluster")){
    colour <- SAPCA$seq.space.clusters$classification
  }else{
    colour <- col
  }
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
  
  # Plot model-based clusters in 3D
  rgl::plot3d(data[,plotPCs],
              col      = colour,      # colour by clusters
              specular = "black",     # matte lighting
              type     = "s",         # "p" is points, "s" is spheres
              radius   = rad,         # sphere radius if using spheres
              size     = 4,           # point size
              axes     = FALSE,       # draw axes separately
              xlab     = axes[1],
              ylab     = axes[2],
              zlab     = axes[3])       
  # Draw axes
  if (write!=FALSE){
    rgl::axes3d(color = "black", labels = FALSE)                       
  }else{
    rgl::axes3d(color = "black", alpha=0.5, labels = FALSE) 
  }
  
  for (NAME in labels){
    SUB = row.names(SAPCA$seq.space.PCA$coordinates)==NAME      # Label based on its row.name
    rgl::text3d(subset(SAPCA$seq.space.PCA$coordinates[,plotPCs],subset=SUB), 
                text      = paste('---',NAME),   # data label text
                font      = 2,                   # bold
                color     = "black",             # colour
                adj       = -rad/2)              # offset
  }
  
  # Write html for interactive data
  if (write!=FALSE){
    rglwidget::.writeWebGL(write)                          
  }
}



plot_overlay_3Dlabel.A <- function(SAPCA,
                                   plotPCs = 1:3){
  selected     <- rgl::select3d()
  selected.set <- selected(SAPCA$seq.space.PCA$coordinates[,plotPCs])
 
  list(selected     = selected,
       selected.set = selected.set)
}


plot_overlay_3Dlabel.B <- function(selection,
                                   SAPCA,
                                   plotPCs = 1:3){
  selected     <- selection$selected
  selected.set <- selection$selected.set
  
  if(sum(selected.set)!=0){
    as.fasta(SAPCA$numerical.alignment$MSA[selected.set,],
             decolgap=TRUE)
    
    
    for (NAME in SAPCA$numerical.alignment$seq.names[selected.set]){
      # Which point will be labelled
      SUB <- NAME   
      # What is the label text
      TEXT <- NAME
      
      rad <- (range(SAPCA$seq.space.PCA$coordinates[,plotPCs])[2] -
                range(SAPCA$seq.space.PCA$coordinates[,plotPCs])[1]) /
        100
      
      rgl::text3d(SAPCA$seq.space.PCA$coordinates[SUB,plotPCs], 
                  text      = paste('---',TEXT),   # data label text
                  font      = 2,                   # bold
                  color     = "black",             # colour
                  adj       = -rad/2)              # offset
    }
  }
}
