# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @export sharedComponents
#' @author Juan David Henao <judhenaosa@unal.edu.co>
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
#' share <- sharedComponents(net1,net2)
#' share


sharedComponents <- function(...){
  # Obtaining the intersection set of the networks
  intersect <- graph.intersection(...,keep.all.vertices = FALSE)
  # Select the nodes in the intersected network with degree equal zero
  membrs <- which(degree(intersect) == 0)
  # To empty vector
  if(length(membrs) == 0){
    # Showing a message
    stop("No shared components found between the networks")
  }else{
    # Return the names of the solitary nodes in the intersection network
    return(names(membrs))
  }
}
