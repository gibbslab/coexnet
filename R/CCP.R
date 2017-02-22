CCP <- function(...){
  
  intersect <- graph.intersection(...,keep.all.vertices = FALSE)
  
  members <- which(degree(intersect) == 0)
  CCP <- delete.vertices(graph = intersect,v = names(members))
  
  return(CCP)
}