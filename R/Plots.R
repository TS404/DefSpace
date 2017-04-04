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


plot_3Dclosest <- function(SAPCA,
                           sequence,
                           plotPCs    = 1:3,
                           measurePCs = 1:3,
                           radius     = 1,
                           n          = 10,
                           write      = FALSE){

  temp <- SAPCA$numerical.alignment$MSA[,1]
  temp[1:length(temp)]                    <- "white"
  temp[rownames(closest(SAPCA,sequence,PC = measurePCs, n = n))] <- "red"
  temp[sequence]                          <- "black"

  plot_3Dclusters(SAPCA,
                  radius  = radius,
                  plotPCs = plotPCs,
                  labels  = sequence,
                  col     = temp,
                  write   = write)
  
}
