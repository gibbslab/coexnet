neta <- read.graph(file = "Esclerosis.net",format = "ncol")
netb <- read.graph(file = "Parkinson.net",format = "ncol")

intersect <- graph.intersection(neta,netb,keep.all.vertices = FALSE)


#### Funcion

members <- which(graph.coreness(intersect) > 0)

CCP <- remove.edge.attribute(graph = intersect,name = members)

lapply(seq_along(graph$csize)[graph$csize > 1], function(x) 
  V(intersect)$name[graph$membership %in% x])


table(names(V(ESC)) %in% names(V(PRK)))

V(intersect)
