# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @export ppiNet
#' @author Juan David Henao <judhenaosa@unal.edu.co>
#' @title Create a protein-protein interaction network
#' @description Creates a protein-protein interaction network using an edge list with the relations between proteins or a vector with the gene symbols or any other molecular identifier type, widely used in biological databases, to create a predictive PPI network using information of evidence in STRING database.
#' @param molecularIDs A vector of IDs recognized by STRING database to create a PPI network from them.
#' @param file A file with an edge list to charge the PPI network.
#' @param speciesID The numerical ID from STRING database to the species of interest, by defect, is "9006" corresponding to human species.
#' @param evidence A vector with the evidence to support the interactions between the proteins, by default is all the evidence given in STRING database.
#' @return An igraph object as the protein-protein interaction network where the nodes are the molecular identifiers given in the input.
#' @examples
#' \dontrun{
#' # Creating a vector with identifiers
#' 
#' ID <- c("FN1","HAMP","ILK","MIF","NME1","PROCR","RAC1","RBBP7",
#' "TMEM176A","TUBG1","UBC","VKORC1")
#' 
#' # Creating the PPI network
#' 
#' ppi <- ppiNet(molecularIDs = ID,evidence = c("neighborhood","coexpression","experiments"))
#' ppi
 #' }
#' 
#' # Creating a PPI network from external data
#' 
#' ppi <- ppiNet(file = system.file("extdata","ppi.txt",package = "coexnet"))
#' ppi

ppiNet <- function(molecularIDs = NULL,file = NULL,speciesID = 9606,evidence = c("neighborhood","neighborhood_transferred",
            "fusion","cooccurence","homology","coexpression","coexpression_transferred",
            "experiments","experiments_transferred","database","database_transferred","textmining",
            "textmining_transferred","combined_score")){
  
  # Detecting the input type
  if(is.null(file) && molecularIDs > 0){
    # Creating the vector to store the unique identifiers
    for_gen <- unlist(sapply(molecularIDs,function(i){
      if(grepl("-",i) > 0){
        return(unlist(strsplit(i,"-")))
      }else{
        return(i)
      }
    }))
    
    # Remove duplicated identifiers
    for_gen <- unique(sort(for_gen))
    # Transform the vector into data frame 
    new_genes <- as.data.frame(for_gen,stringsAsFactors = FALSE)
    # The column name must be gene
    names(new_genes) <- "gene"
    # Loading the STRING database
    database <- STRINGdb$new(version="10",species=speciesID,score_threshold=0,input_directory="")
    # Obtaining the STRING ID to each identifier
    mapped <- database$map(new_genes,"gene",removeUnmappedRows = TRUE)
    # Removing the STRING ID duplicated
    mapped <- mapped[!duplicated(mapped$STRING_id),]
    # Obtaining the interactions among STRING IDs according to different types of evidence
    interactions <- database$get_interactions(mapped$STRING_id)
    # Extract the relations deleting the evidence information columns
    graph_relations <- data.frame(interactions$from,interactions$to,stringsAsFactors = FALSE)
    # Read each evidence given
    for(i in evidence){
      # Take each column with evidence information
      for(j in names(interactions)){
        # If the evidence in the "interactions" variable corresponds with one of the request evidence
        if(i == j){
          # Appends the column with the requested evidence in the data frame with interactions among STRING IDs
          graph_relations <- cbind(graph_relations,interactions[j])
        }
      }
    }
    # This data frame will be fill up with interactions with any evidence value greater than zero
    graph_ppi <- graph_relations[rowSums(graph_relations[,seq(3,ncol(graph_relations))]) > 0,]
    # This loop replace the STRING IDs with the original identifiers in the first column 
    for(n in seq_len(nrow(graph_ppi))){
      graph_ppi$interactions.from[n] <- mapped[
        graph_ppi$interactions.from[n] == mapped$STRING_id,][1]
    }
    # This loop replace the STRING IDs with the original identifiers in the second column
    for(n in seq_len(nrow(graph_ppi))){
      graph_ppi$interactions.to[n] <- mapped[
        graph_ppi$interactions.to[n] == mapped$STRING_id,][1]
    }
    
    # Is necessary transform the data frame into matrix
    edge_list <- as.matrix(as.vector(graph_ppi[,1],mode = "character"))
    edge_list <- cbind(edge_list,as.vector(graph_ppi[,2],mode = "character"))
    # Create a network based on the interactions in the matrix
    final_graph <- graph.edgelist(edge_list,directed = FALSE)
  }else if(is.null(molecularIDs) && file.exists(file)){
    # Read and create the igraph object
    final_graph <- read.graph(file = file,format = "ncol")
  }else{
    stop("A valid input was not found: Please, be sure to use an edge list or a vector of IDs as input")
  }
  # Return the PPI network
  return(final_graph)
}