#' @export get.affy
#' @author Juan david Henao <judhenaosa@unal.edu.co>
#' 
#' # Bioinformatics and Systems Biology | Universidad Nacional de Colombia
#' 
#' @title Charge and create an AffyBatch object 
#' @description Charges the data from a file with the GSM samples compressed using the "filelist.txt" file to identify the GSM ID
#' and create the AffyBatch object finally.
#' @param GSE The name of the file with the compressed samples data.
#' @param dir The pathfile where the samples and the filelist file are to create the AffyBatch object, by default is the current directory.
#' @seealso \code{\link{get.info}} to download expression data from GEO DataSet included the filelist file.
#' @examples 
#' 
#' # Load the AffyBatch from downloaded raw data
#' # The data is from partial experiment results from Greene JG(2006), GEO accession: GSE4773
#' 
#' # Data whitout CDF enviroment information
#' 
#' affy <- get.affy(GSE = "GSE1234",dir = system.file("extdata",package = "coexnet"))

get.affy <- function(GSE,dir="."){
  
  # Change the directory 
  
  setwd(dir)
  
  # Reads the filelist.txt file with the samples information
  
  raw <- read.table(file = paste0(GSE,"/","filelist.txt"),sep = "\t",
                    header = T,comment.char = "", stringsAsFactors = F)
  
  # Obtains all the names of samples
  
  GSMs <- raw$Name[grep(".CEL",raw$Name,ignore.case = T)]
  
  # Reads and returns the raw data from each sample
  
  affy <- ReadAffy(filenames = as.character(GSMs), compress = T,
                   celfile.path = GSE)
  
  return(affy)
}