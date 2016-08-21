# gene.symbol
# Juan David Henao Sanchez
# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

gene.symbol <- function(GPL, d = "."){

  # Information extracted from the file .soft
  
  gpl <- getGEO(filename = paste0(GPL,".soft"))
  
  # Creates a table object with all the data in the .soft file
  
  sym <- Table(gpl)
  
  # Creates a data.frame object with the gene symbos to the probes associated
  
  ta <- data.frame(sym$ID, gsub(" /// ","-",sym$`Gene Symbol`), stringsAsFactors = F)
  
  # Names of each column
  
  names(ta) <- c("probe","gene")
  
  return(ta)
}