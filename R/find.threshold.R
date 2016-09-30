create.net <- function(){
  
  Ad <- matrix(0,ncol = nrow(simil),nrow = nrow(simil))
  
  # Transforms to adjacency matrix
  
  for(i in 1:nrow(simil)){
    Ad[which(simil[,i]>=value),i]<-1
    Ad[which(simil[,i]<value),i]<-0
  }
  
  # Changes the names of adjacency matrix to the genes names
  
  colnames(Ad)<-rownames(Ad)<-rownames(simil)
  
  # Diagonal equal to zero
  
  diag(Ad)<-0
  
  # Creates the network from the adjacency matrix
  
  Gr=graph.adjacency(Ad,mode="undirected",add.colnames=NULL,diag=FALSE)
  
  # Shows the final values of the co-expression network
  
  cat("Final Network",sep = "\n")
  cat(paste("p-value =",fit,sep = " "),sep = "\n")
  cat(paste("threshold =",value,sep = " "),sep = "\n")
  
  # Create a plot with local maximum.
  
  plot(x = pcv,y = abs(Cis-C0s),t="l",xlab="Threshold",ylab="|C-C0|")
  
  abline(v=value,col="red")
  
  # Returns the network as an igraph object
  
  return(Gr)
}