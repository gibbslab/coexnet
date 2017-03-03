#' @export create.net
#' @author Juan David Henao Sanchez <judhenaosa@unal.edu.co>
#' @author Liliana Lopez Kleine <llopezk@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Creating a co-expression network from expression matrix.
#' @description From an expression matrix, this function creates a co-expression network like a graph object using a threshold value and one similarity function.
#' @param difexp A whole expression matrix or differentially expressed genes matrix.
#' @param method A function to calculate the similarity matrix between genes. It can be "correlation" to use Pearson function or "mutual information" to use a based on entropy information function.
#' @param threshold A value between 0 and 1 to filter the similarity matrix and create the co-expression network.
#' @return An undirected co-expression network as igraph object.
#' @seealso \code{\link{find.threshold}} to obtain a threshold value based on biology network assumptions.
#' @examples 
#' 
#' # Loading data
#' 
#' pathfile <- system.file("extdata","expression_example.txt",package = "coexnet")
#' data <- read.table(pathfile,stringsAsFactors = FALSE)
#' 
#' # Building the network
#' 
#' cor_pearson <- create.net(difexp = data,threshold = 0.7,method = "correlation")
#' cor_pearson
#' 
#' mut_inf <- create.net(difexp = data,threshold = 0.5,method = "mutual information")
#' mut_inf

create.net <- function(difexp,method, threshold){
  
  # Obtains the similarity values
  
  simil <- .correlation.matrix(difexp,method)
  
  # Creates an empty matrix
  
  Ad <- matrix(0,ncol = nrow(simil),nrow = nrow(simil))
  
  # Transforms to adjacency matrix
  
  for(i in 1:nrow(simil)){
    Ad[which(simil[,i]>=threshold),i]<-1
    Ad[which(simil[,i]<threshold),i]<-0
  }
  
  # Changes the names of adjacency matrix to the names of the genes
  
  colnames(Ad)<-rownames(Ad)<-rownames(simil)
  
  # Diagonal equal to zero
  
  diag(Ad)<-0
  
  # Creates the network from the adjacency matrix
  
  Gr=graph.adjacency(Ad,mode="undirected",add.colnames=NULL,diag=FALSE)
  
  # Obtains the degree value for each node in the network
  
  de <- degree(Gr,loops = F)
  
  # Deletes the nodes without any edge
  
  Ad <- Ad[which(de > 0), which(de > 0)]
  
  # Creates the final network from the adjacency matrix filtered
  
  net <- graph.adjacency(Ad,mode="undirected",add.colnames=NULL,diag=FALSE)
  
  # Returns the network like igraph object
  
  return(net)
  
}