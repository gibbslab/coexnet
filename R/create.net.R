#' @export create.net
#' @author Juan David Henao Sanchez <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Create a co-expression matrix using correlation matrix
#' @description Create a co-expression network from the expressed diferential genes matrix and a threshold value using a 
#' correlation matrix created from two different methods. The first one is the Pearson Correlation Coefficient and the second one is 
#' the Mutual Information test.
#' @param difexp A data.frame with the expressed diferential genes matrix.
#' @param method The method to create the correlation matrix, can be "correlation" to use the Pearson Correlation Coefficient or
#' "mutual information" to use the Mutual Information test.
#' @param threshold The threshold value to establish the edges in the co-expression network.
#' @return A igraph objecto with the co-expression network.
#' @seealso \code{\link{find.threshold}} to find a threshold value.

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
  
  return(Gr)
  
}