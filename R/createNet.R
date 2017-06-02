# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @export createNet
#' @author Juan David Henao Sanchez <judhenaosa@unal.edu.co>
#' @author Liliana Lopez Kleine <llopezk@unal.edu.co>
#' @title Creating a co-expression network from expression matrix.
#' @description From an expression matrix, this function creates a co-expression network like a graph object using a threshold value and one similarity function.
#' @param expData A whole expression matrix or differentially expressed genes matrix, it may be stored in a SummarizedExperiment object.
#' @param method A function to calculate the similarity matrix between genes. It can be "correlation" to use Pearson function or "mutual information" to use a based on entropy information function.
#' @param threshold A value between 0 and 1 to filter the similarity matrix and create the co-expression network.
#' @return An undirected co-expression network as igraph object.
#' @seealso \code{\link{findThreshold}} to obtain a threshold value based on biology network assumptions.
#' @examples 
#' 
#' # Loading data
#' 
#' pathfile <- system.file("extdata","expression_example.txt",package = "coexnet")
#' data <- read.table(pathfile,stringsAsFactors = FALSE)
#' 
#' # Building the network
#' 
#' cor_pearson <- createNet(expData = data,threshold = 0.7,method = "correlation")
#' cor_pearson
#' 
#' mut_inf <- createNet(expData = data,threshold = 0.5,method = "mutual information")
#' mut_inf

createNet <- function(expData,method,threshold){
  
  # Identifing the SummarizedExperiment object
  
  if(is(expData,"SummarizedExperiment")){
    expData <- assay(expData)
  }
  
  # Obtaining the similarity values
  
  simil <- .correlation.matrix(expData,method)
  
  # Create an empty matrix
  
  Ad <- matrix(0,ncol = nrow(simil),nrow = nrow(simil))
  
  # Transform to adjacency matrix
  
  for(i in seq_len(nrow(simil))){
    Ad[which(simil[,i]>=threshold),i]<-1
    Ad[which(simil[,i]<threshold),i]<-0
  }
  
  # Change the names of adjacency matrix to the names of the genes
  
  colnames(Ad)<-rownames(Ad)<-rownames(simil)
  
  # Diagonal equal to zero
  
  diag(Ad)<-0
  
  # Creating the network from the adjacency matrix
  
  Gr=graph.adjacency(Ad,mode="undirected",add.colnames=NULL,diag=FALSE)
  
  # Obtaining the degree value for each node in the network
  
  de <- degree(Gr,loops = FALSE)
  
  # Delete the nodes without any edge
  
  Ad <- Ad[which(de > 0), which(de > 0)]
  
  # Create the final network from the adjacency matrix filtered
  
  net <- graph.adjacency(Ad,mode="undirected",add.colnames=NULL,diag=FALSE)
  
  # Return the network like igraph object
  
  return(net)
  
}
