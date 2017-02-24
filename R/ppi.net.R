ppi.net <- function(genes,species_ID = 9606,evidence = c("neighborhood","neighborhood_transferred",
            "fusion","cooccurence","homology","coexpression","coexpression_transferred",
            "experiments","experiments_transferred","database","database_transferred","textmining",
            "textmining_transferred","combined_score")){
  
  for_gen <- vector()
  
  for(i in genes){
    for(j in strsplit(i,"-")){
      for_gen <- append(for_gen,j)
    }
  }
  
  for_gen <- unique(sort(for_gen))
  
  new_genes <- as.data.frame(for_gen,stringsAsFactors = FALSE)
  names(new_genes) <- "gene"
  
  database <- STRINGdb$new(version="10",species=species_ID,score_threshold=0,input_directory="")
  mapped <- database$map(new_genes,"gene",removeUnmappedRows = TRUE)
  
  mapped <- mapped[!duplicated(mapped$STRING_id),]
  
  interactions <- database$get_interactions(mapped$STRING_id)
  
  graph_relations <- data.frame(interactions$from,interactions$to,stringsAsFactors = FALSE)
  
  for(i in evidence){
    for(j in names(interactions)){
      if(i == j){
        graph_relations <- cbind(graph_relations,interactions[j])
      }
    }
  }
  
  graph_ppi <- data.frame()
  
  for(i in 1:nrow(graph_relations)){
    if(any(graph_relations[i,3:ncol(graph_relations)] > 0)){
      graph_ppi <- rbind.data.frame(graph_ppi,graph_relations[i,],
                                    stringsAsFactors = FALSE)
    }
  }
  
  for(n in 1:nrow(graph_ppi)){
    graph_ppi$interactions.from[n] <- mapped[
      graph_ppi$interactions.from[n] == mapped$STRING_id,][1]
  }
  
  for(n in 1:nrow(graph_ppi)){
    graph_ppi$interactions.to[n] <- mapped[
      graph_ppi$interactions.to[n] == mapped$STRING_id,][1]
  }
  
  edge_list <- matrix()
  edge_list <- cbind(as.vector(graph_ppi[,1],mode = "character"))
  edge_list <- cbind(edge_list,as.vector(graph_ppi[,2],mode = "character"))
  
  final_graph <- graph.edgelist(edge_list,directed = FALSE)
  
  return(final_graph)
}

