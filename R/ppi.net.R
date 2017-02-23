PRK <- read.table(file = "GSE9807_0,3.txt",stringsAsFactors = FALSE)
genes <- as.data.frame(row.names(PRK),stringsAsFactors = FALSE)
names(genes) <- "gene"



### Function

for(n in 1:nrow(genes)){
  if(grepl("-",genes[n,]) == TRUE){
    print("find")
  }
}

database <- STRINGdb$new(version="10",species=9606,score_threshold=0,input_directory="")
mapped <- database$map(genes,"gene",removeUnmappedRows = TRUE)

interactions <- database$get_interactions(mapped$STRING_id)

evidence <- c("neighborhood","coexpression","experiments")
graph_relations <- data.frame(interactions$from,interactions$to,stringsAsFactors = F)

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

c(t(as.character(graph_ppi$interactions.from)), t(as.character(graph_ppi$interactions.to)))

edge_list <- as.matrix(graph_ppi[,1:2])

final_graph <- graph.edgelist(edge_list,directed = FALSE)
