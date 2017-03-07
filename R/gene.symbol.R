#' @export gene.symbol
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Create a table relating probesets with genes
#' @description From the information in the .soft file, creates a data.frame with two columns. In the first one, the probeset names. In the second one, the name of the corresponding genes.
#' @param GPL The GPL ID.
#' @param d The path file where the .soft file is. By default is the current path file.
#' @return A data.frame with two columns, in the first one, the probesets and in the second one, the corresponding gene to each probeset.
#' @examples 
#' 
#' # Creating the table with probesets and genes/IDs
#' 
#' gene_table <- gene.symbol(GPL = "GPL2025",d = system.file("extdata",package = "coexnet"))
#' 
#' # Cleaning the NA information
#' 
#' gene_na <- na.omit(gene_table)
#' 
#' # Cleaning gene/ID information empty
#' 
#' final_table <- gene_na[gene_na$ID != "",]
#' 
#' head(final_table)


gene.symbol <- function(GPL, d = "."){
  
  options(warn=-1)
  
  # Move to path file with the .soft file
  
  setwd(d)

  # Information extracted from the file .soft
  
  gpl <- getGEO(filename = paste0(GPL,".soft"))
  
  
  # Create a table object with all the data in the .soft file
  
  sym <- Table(gpl)
  
  # Create a data.frame with the gene symbols for the probes associated
  
  ta <- data.frame(sym$ID, gsub(" /// ","-",sym$`Gene Symbol`), stringsAsFactors = F)
  
  # Names of each column
  
  names(ta) <- c("probe","ID")
  
  return(ta)
}