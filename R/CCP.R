#' @export CCP
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Finding Common Connection Pattern between different networks.
#' @description From the intersection of two or more networks, it obtains the connected components 
#' by deleting the solitary nodes in the intersection network. They must be graph objects.
#' @param ... The networks (igraph objects) to obtain the Common Connection Patterns.
#' @return An igraph object with all the Common Connection Pattern in the same network.
#' @details 
#' The Common Connection Pattern (CCP), is a new methodological proposal to identify molecular components linked together and common in several biological networks. The principal assumption behind Common Connection Pattern is that the networks to be compared must have the same molecular information from, i.e., must compare one layer of molecular abstraction at the same time, for example, co-expression layer, protein-protein layer, the gene regulation layer, among others.
#' In general, the comparison of biological networks is made to determine common elements or biomarkers among several related phenotypes. In this case, the Common Connection Patterns aims to identify common molecular elements between these phenotypes that are associated also with each other in a specific way. For this, the intersection between biological networks is calculated whose result can have two kinds of elements. On the one hand, the shared nodes without any connection with other nodes in the intersection network. On the other hand, nodes connected to one or more nodes in the intersection network. Each connected component in the intersection will be considered as a Common Connection Pattern.
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
