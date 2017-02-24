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
