#' @export create.net
#' @author Juan David Henao Sanchez <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia


create.net <- function(difexp,method, threshold){
  
  simil <- .correlation.matrix(difexp,method)
  
  Ad <- matrix(0,ncol = nrow(simil),nrow = nrow(simil))
  
  # Transforms to adjacency matrix
  
  for(i in 1:nrow(simil)){
    Ad[which(simil[,i]>=threshold),i]<-1
    Ad[which(simil[,i]<threshold),i]<-0
  }
  
  # Changes the names of adjacency matrix to the genes names
  
  colnames(Ad)<-rownames(Ad)<-rownames(simil)
  
  # Diagonal equal to zero
  
  diag(Ad)<-0
  
  # Creates the network from the adjacency matrix
  
  Gr=graph.adjacency(Ad,mode="undirected",add.colnames=NULL,diag=FALSE)
  
  de <- degree(Gr,loops = F)
  
  Ad <- Ad[which(de > 0), which(de > 0)]
  
  net <- graph.adjacency(Ad,mode="undirected",add.colnames=NULL,diag=FALSE)
  
  write.graph(graph = net,file = "~/testEN2.net","ncol")
  
  return(net)
  
}