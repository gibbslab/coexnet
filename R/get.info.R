#' @export get.info
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Download raw expression data from GEO DataSet
#' @description Using a GSE ID list, from the GEO DtaSet ftp, search and download each one of the raw expression data in a 
#' separate file with the GSE ID as the name, after that, this function uncompress the downloaded file to be ready to use.
#' Adittionally, this function search and downloaded the information of the chip using the GPL ID to obtain the .soft file.
#' @usage get.info(GSE, GPL, dir = ".")
#' @param GSE The list of GSE ID.
#' @param GPL The GPL ID of the chip used in the microarray experiment.
#' @param dir The pathfile to download and uncompress the raw wxpewssion data, by default is the current pathfile.
#' @return A set of files with the uncompressed data and ready to use.
#' @seealso \code{\link{get.affy}} to charge the expression data.
#' @examples 
#' 
#' # Extract data from GEO DataSet
#' 
#' get.info(GSE = "GSE8216", GPL = "GPL2025")

get.info <- function(GSE,GPL,dir="."){
  
  # Moves to the specified directory
  
  setwd(dir)
  
  # Downloads the raw data from GEO DataSets database
  
  sapply(as.vector(t(GSE)), getGEOSuppFiles)
  
  # Obtains the name of the compressed data
  
  files <- dir(".")[grep("^GSE[0-9]",dir("."),ignore.case = T)]
  
  # Uncompress all the raw data
  
  for(j in files){
    untar(paste0(j,"/",j,"_RAW.tar"), exdir = paste0(j,"/"))
  }
  
  # Downloads the .soft file of the microarray chip
  
  getGEOfile(GPL, destdir = ".")
}