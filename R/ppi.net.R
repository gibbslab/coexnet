PRK <- read.table(file = "GSE9807_0,3.txt")
genes <- as.data.frame(row.names(PRK))
names(genes) <- "gene"

database <- STRINGdb$new(version="10",species=9606,score_threshold=0,input_directory="")
mapped <- database$map(genes,"gene",removeUnmappedRows = TRUE)

proteins <- mapped$STRING_id

interactions <- database$get_interactions(names(genes))
