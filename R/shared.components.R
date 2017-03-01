#' @export shared.components
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Finding shared components between diferent networks.
#' @description From the intersected network, obtain the nodes without any value of degree. There are
#' the molecular shared components without biological relation among them.
#' @param ... The networks (igraph objects) to obtain the shared components.
#' @return A vector with the names of the shared components.
#' @examples 
#' # Loading data
#' 
#' data("net1")
#' data("net2")
#' 
#' # Obtain shared components
#' 
#' share <- shared.components(net1,net2)
#' share


shared.components <- function(...){
  # Obtain the intersection set of the networks
  intersect <- graph.intersection(...,keep.all.vertices = FALSE)
  # Select the nodes in the intersected network with degree equal cero
  members <- which(degree(intersect) == 0)
  # Return the names of the alone nodes in the intersected network
  return(names(members))
}
