#' @export CCP
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Finding Common Connection Pattern between different networks.
#' @description From the ntersection of two or more networks, obntain the connected components deleting
#' the alone nodes in the intersection network. The networks must be igraph objects.
#' @param ... The networks (igraph objects) to obtain the Common Connection Patterns.
#' @return An igraph object with all the Common Connection Pattern in the same network.
#' @examples 
#' # Charge data
#' 
#' data("net1")
#' data("net2")
#' 
#' # Obtain Common Connection Patterns
#' 
#' ccp <- CCP(net1,net2)
#' ccp

CCP <- function(...){
  # Obtain the intersection set of the networks
  intersect <- graph.intersection(...,keep.all.vertices = FALSE)
  # Select the intersected nodes without any value of dregree
  members <- which(degree(intersect) == 0)
  # Delete the aone vertices in the intersection network
  CCP <- delete.vertices(graph = intersect,v = names(members))
  # Return only the nodes in the intersection network with at least one edge
  return(CCP)
}
