#' @export ppi.net
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Create a protein-protein interaction network from STRING database
#' @description From a molecular identifiers (gene symbol or protein symbol), create a protein-protein
#' network using the different kind of evidence from STRING database.
#' @param genes A vector with the molecular identifiers, they could be two IDs united for the "-"
#' character.
#' @param species_ID The numerical ID from STRING database to the specie of interest, by defect is
#' "9006" corresponding to human specie.
#' @param evidence A vector with the evicences to support the interactions bewtween the proteins
#' by default are all te evidences given in STRING database.
#' @return An igraph object as the protein-protein interaction network where the nodes are the
#' molecular identifiers given in the vector input.
#' @examples
#' \dontrun{
#' # Creating vector with identifiers
#' 
#' ID <- c("FN1","HAMP","ILK","MIF","NME1","PROCR","RAC1","RBBP7",
#' "TMEM176A","TUBG1","UBC","VKORC1")
#' 
#' # Creating the PPI network
#' 
#' ppi <- ppi.net(genes = ID,evidence = c("neighborhood","coexpression","experiments"))
#' ppi
#' }

ppi.net <- function(genes,species_ID = 9606,evidence = c("neighborhood","neighborhood_transferred",
            "fusion","cooccurence","homology","coexpression","coexpression_transferred",
            "experiments","experiments_transferred","database","database_transferred","textmining",
            "textmining_transferred","combined_score")){
  
  # Detecting the input type
  if(is.vector(input)){
    # Replacing the name of input
    genes <- input
    # Creating the vector to store the unique identifiers
    for_gen <- vector()
    # To each ID in the vector input
    for(i in genes){
      # Split each ID with two or more different identifiers
      for(j in strsplit(i,"-")){
        # Add each of the ID such as 
        for_gen <- append(for_gen,j)
      }
    }
    # Remove de dupicated identifiers
    for_gen <- unique(sort(for_gen))
    # Transform the vector into data-frame 
    new_genes <- as.data.frame(for_gen,stringsAsFactors = FALSE)
    # The name of column must be "gene"
    names(new_genes) <- "gene"
    # Charge the STRING database
    database <- STRINGdb$new(version="10",species=species_ID,score_threshold=0,input_directory="")
    # Obtain the STRING ID to each identifier
    mapped <- database$map(new_genes,"gene",removeUnmappedRows = TRUE)
    # Remove the STRING ID duplicated
    mapped <- mapped[!duplicated(mapped$STRING_id),]
    # Obtain the interactions among STRING IDs from diferent types of evidences
    interactions <- database$get_interactions(mapped$STRING_id)
    # Extract the relations deleting the evidences
    graph_relations <- data.frame(interactions$from,interactions$to,stringsAsFactors = FALSE)
    # Read each evidence given
    for(i in evidence){
      # Take each column of all the evidences
      for(j in names(interactions)){
        # If the evidence in the "interactions" variable corresponde with one of the request evidence
        if(i == j){
          # Append the column with the request evidence in the data.frame with interactions
          # among STRING_ID
          graph_relations <- cbind(graph_relations,interactions[j])
        }
      }
    }
    # This data.frame will be fill up with interactions with any interaction greater than 0
    graph_ppi <- data.frame()
    # This loop fill up the "graph_ppi" data.frame
    for(i in 1:nrow(graph_relations)){
      if(any(graph_relations[i,3:ncol(graph_relations)] > 0)){
        graph_ppi <- rbind.data.frame(graph_ppi,graph_relations[i,],
                                      stringsAsFactors = FALSE)
      }
    }
    # This loop replace the STRING IDs with the original identigiers in the first column 
    for(n in 1:nrow(graph_ppi)){
      graph_ppi$interactions.from[n] <- mapped[
        graph_ppi$interactions.from[n] == mapped$STRING_id,][1]
    }
    # This loop replace the STRING IDs with the original identigiers in the second column
    for(n in 1:nrow(graph_ppi)){
      graph_ppi$interactions.to[n] <- mapped[
        graph_ppi$interactions.to[n] == mapped$STRING_id,][1]
    }
    # Is necesary transform the data.frame into matrix
    edge_list <- matrix()
    # Filling up the matrix with the interactions
    edge_list <- cbind(as.vector(graph_ppi[,1],mode = "character"))
    edge_list <- cbind(edge_list,as.vector(graph_ppi[,2],mode = "character"))
    # Creating a network from interactions in the matrix
    final_graph <- graph.edgelist(edge_list,directed = FALSE)
  }else{
    # Reading and creating the igraph object
    final_graph <- read.graph(file = input,format = "ncol")
  }
  # Returns the PPI network
  return(final_graph)
}

