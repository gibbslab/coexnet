#' @export CCP
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Finding Common Connection Pattern between different networks.
#' @description From the intersection of two or more networks, it obtains the connected components 
#' by deleting the solitary nodes in the intersection network. They must be igraph objects.
#' @param ... The networks (igraph objects) to obtain the Common Connection Patterns.
#' @param by It can be, a degree value or a range of values (c(min, max)), to defined the Common Connection Pattern. If you pass one value all the nodes above this degree will be used. By default is NULL, it calculates the CCPs using all the network.
#' @return An igraph object with all the Common Connection Pattern in the same network.
#' @details 
#' The Common Connection Pattern (CCP), is a new methodological proposal to identify molecular components linked together and common in several biological networks. The principal assumption behind Common Connection Pattern is that the networks to be compared must have the same molecular information from, i.e., must compare one layer of molecular abstraction at the same time, for example, co-expression layer, protein-protein layer, the gene regulation layer, among others.
#' 
#' For this, the intersection between biological networks is calculated whose result are the sub-networks with diameter greater than zero, being each of them considered as a Common Connection Pattern.
#' @examples 
#' # Loading data
#' 
#' data("net1")
#' data("net2")
#' 
#' # Obtaining Common Connection Patterns
#' 
#' ccp <- CCP(net1,net2)
#' ccp

CCP <- function(...,by = NULL){
  
  if(is.null(by)){
    # Identifying the intersection network
    ccp <- graph.intersection(...,keep.all.vertices = FALSE)
  }else{
    
    decom <- list(...)
    
    pre_ccp <- NULL
    
    for(i in decom){
      # Filtering the network from a data frame
      df <- as_data_frame(i)
      
      # Identifies the nodes with the degree given by the user
      if(length(by) == 2){
        v <- names(which(degree(i) >= by[1] & degree(i) <= by[2]))
      }else{
        v <- names(which(degree(i) >= by[1]))
      }
      
      # Filters the network as data frame
      ft <- df[df$from %in% v,]
      ft <- rbind(ft,df[df$to %in% v,])
      
      # Creating the igraph object
      newG <- graph_from_data_frame(unique(ft),directed = FALSE)
      
      # Obtaining the intersection network
      if(is.null(pre_ccp)){
        pre_ccp <- newG
      }else{
        ccp <- graph.intersection(pre_ccp,newG,keep.all.vertices = FALSE)
        pre_ccp <- ccp
      }
    }
  }
  
  # Obtains the names of vertices with diameter equal to zero
  members <- which(degree(ccp) == 0)
  # Deletes the solitary vertices in the intersection network
  finalccp <- delete.vertices(graph = ccp,v = names(members))
  
  if(diameter(finalccp) == 0){
    # Showing a messege
    stop("No Common Connection Patterns found between networks")
  }else{
    # Return only the nodes in the intersection network with at least one edge
    return(finalccp)
  }
}