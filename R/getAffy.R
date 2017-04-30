# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @export getAffy
#' @author Juan david Henao <judhenaosa@unal.edu.co> 
#' @title Charge and create an AffyBatch object 
#' @description Charges the data from a file with the GSM samples compressed using the "filelist.txt" file to identify the GSM ID
#' and create the AffyBatch object finally.
#' @param GSE The name of the file with the compressed samples data.
#' @param directory The path file where the samples and the "filelist" file are to create the AffyBatch object, by default is the current directory.
#' @seealso \code{\link{getInfo}} to download expression data from GEO DataSets included the "filelist" file.
#' @return An AffyBatch Object.
#' @examples 
#' 
#' # Creating the AffyBatch object from raw data downloaded
#' # The data is from partial experiment results from Greene JG(2006), GEO accession: GSE4773
#' 
#' # Data without CDF environment information
#' 
#' affy <- getAffy(GSE = "GSE1234",directory = system.file("extdata",package = "coexnet"))

getAffy <- function(GSE,directory="."){
  
  # Read the filelist.txt file with the samples information
  
  raw <- read.table(file = paste0(directory,"/",GSE,"/","filelist.txt"),sep = "\t",
                    header = TRUE,comment.char = "", stringsAsFactors = FALSE)
  
  # Obtaining all the names of samples
  
  GSMs <- raw$Name[grep(".CEL",raw$Name,ignore.case = TRUE)]
  
  # Read and return the raw data from each sample
  
  affy <- ReadAffy(filenames = as.character(GSMs), compress = TRUE,
                   celfile.path = paste0(directory,"/",GSE))
  
  return(affy)
}