
shared.components <- function(...){
  
  intersect <- graph.intersection(...,keep.all.vertices = FALSE)
  
  members <- which(degree(intersect) == 0)
  
  return(names(members))
}