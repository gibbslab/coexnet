#' @export gene.symbol
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Create a table relating probesets with genes.
#' @description From the information in the .soft file, creates a data.frame with two columns, the first one with the probeset names
#' and the second one with the name of the corresponding genes. 
#' The genes are annotated using the gene symbol system.
#' @param GPL The GPL ID.
#' @param d The pathfile where the .soft file are, by default is the current pathfile.
#' @return A data.frame with two columns, in the first one the probesets and in the second one the corresponding gene to each probeset.

gene.symbol <- function(GPL, d = "."){
  
  options(warn=-1)

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