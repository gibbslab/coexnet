#' @export CCP
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Finding Common Connection Pattern between different networks.
#' @description From the intersection of two or more networks, it obtains the connected components 
#' by deleting the solitary nodes in the intersection network. They must be graph objects.
#' @param ... The networks (igraph objects) to obtain the Common Connection Patterns.
#' @return An igraph object with all the Common Connection Pattern in the same network.
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

CCP <- function(...){
  # Obtains the intersection set of the networks
  intersect <- graph.intersection(...,keep.all.vertices = FALSE)
  # Selects the intersection nodes without any value of degree
  members <- which(degree(intersect) == 0)
  # Deletes the solitary vertices in the intersection network
  CCP <- delete.vertices(graph = intersect,v = names(members))
  # Returns only the nodes in the intersection network with at least one edge
  return(CCP)
}
