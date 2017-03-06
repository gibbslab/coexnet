#' @export shared.components
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Finding shared components between diferent networks.
#' @description From the intersection network, obtains the nodes without any degree value. These are the shared molecular components without any biological relation among them.
#' @param ... The networks (igraph objects) to obtain the shared components.
#' @return A vector with the names of the shared components.
#' @examples 
#' # Loading data
#' 
#' data("net1")
#' data("net2")
#' 
#' # Obtaining shared components
#' 
#' share <- shared.components(net1,net2)
#' share


shared.components <- function(...){
  # Obtains the intersection set of the networks
  intersect <- graph.intersection(...,keep.all.vertices = FALSE)
  # Selects the nodes in the intersected network with degree equal zero
  members <- which(degree(intersect) == 0)
  # Returns the names of the solitary nodes in the intersection network
  return(names(members))
}
