neta <- read.graph(file = "Esclerosis.net",format = "ncol")
netb <- read.graph(file = "Parkinson.net",format = "ncol")

intersect <- graph.intersection(neta,netb,keep.all.vertices = FALSE)


#### Funcion

members <- which(degree(intersect) == 0)

remove.vertex.attribute(graph = intersect,name = names(members)[1])
