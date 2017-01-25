#' @export create.net
#' @author Juan David Henao Sanchez <judhenaosa@unal.edu.co>
#' @author Liliana Lopez Kleine <llopezk@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Create a co-expression network from expression matrix.
#' @description From a expression matrix, this function create a co-expression network like a igrph object using a threshold value and
#' one similarity function.
#' @param difexp  A whole expression matrix or the expression matrix to differentially expressed genes.
#' @param method  The function to calculate the similarity among the genes. "correlation" to use Pearson function or 
#' "mutual information" to use a function based on entropy information.
#' @param threshold  A value (between 0.1 to 0.99) to cut off the adjacency matrix and create the co-expression network.
#' @return A undirected co-expression network like igraph object.
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
#' cor_pearson <- create.net(difexp = a,threshold = 0.7,method = "correlation")
#' cor_pearson
#' 
#' mut_inf <- create.net(difexp = a,threshold = 0.5,method = "mutual information")
#' mut_inf

create.net <- function(difexp,method, threshold){
  
  # Obtains the similarity values
  
  simil <- .correlation.matrix(difexp,method)
  
  # Creates a empty matrix
  
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
  
  # Obtains the degree value for each node in the network
  
  de <- degree(Gr,loops = F)
  
  # Deletes the node without any edge
  
  Ad <- Ad[which(de > 0), which(de > 0)]
  
  # Creates the final network from the adjacency matrix filtered
  
  net <- graph.adjacency(Ad,mode="undirected",add.colnames=NULL,diag=FALSE)
  
  # return the network like igraph object
  
  return(net)
  
}