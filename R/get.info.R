#' @export get.info
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Download raw expression data from GEO DataSet
#' @description Using a GSE ID list, from the GEO DataSets FTP, it searches and downloads each raw expression data in a separate file with the GSE ID as the name. After, the function uncompresses the downloaded file to be ready to use. Additionally, searches and downloaded the information of the chip using the GPL ID to obtain the .soft file.
#' @param GSE The list of GSE ID.
#' @param GPL The GPL ID of the chip used in the microarray experiment.
#' @param dir The path file to download and uncompress the raw expression data, by default is the current path file.
#' @return A set of files with the uncompressed data and ready to use.
#' @seealso \code{\link{get.affy}} to charge the expression data.
#' @examples 
#' 
#' # Extract data from GEO DataSets
#' 
#' get.info(GSE = "GSE8216", GPL = "GPL2025",dir = tempdir())

get.info <- function(GSE,GPL,dir="."){
  
  # Move to the specified directory
  
  setwd(dir)
  
  # Download the raw data from GEO DataSets database
  
  sapply(as.vector(t(GSE)), getGEOSuppFiles)
  
  # Obtaining the name of the compressed data
  
  files <- dir(".")[grep("^GSE[0-9]",dir("."),ignore.case = TRUE)]
  
  # Uncompress all the raw data
  
  for(j in files){
    untar(paste0(j,"/",j,"_RAW.tar"), exdir = paste0(j,"/"))
  }
  
  # Download the .soft file of the microarray chip
  
  getGEOfile(GPL, destdir = ".")
}