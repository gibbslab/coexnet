#' @export shared.components
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title 

shared.components <- function(...){
  # Obtain the intersection set of the networks
  intersect <- graph.intersection(...,keep.all.vertices = FALSE)
  # Select the nodes in the intersected network with degree equal cero
  members <- which(degree(intersect) == 0)
  # Return the names of the alone nodes in the intersected network
  return(names(members))
}
