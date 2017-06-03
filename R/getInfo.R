# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @export getInfo
#' @author Juan David Henao <judhenaosa@unal.edu.co>
#' @title Download raw expression data from GEO DataSet
#' @description Using a GSE ID list, from the GEO DataSets FTP, it searches and downloads each raw expression data in a separate file with the GSE ID as the name. After, the function uncompresses the downloaded file to be ready to use. Additionally, searches and downloaded the information of the chip using the GPL ID to obtain the .soft file.
#' @param GSE The list of GSE ID.
#' @param GPL The GPL ID of the chip used in the microarray experiment.
#' @param directory The path file to download and uncompress the raw expression data, by default is the current path file.
#' @return A set of files with the uncompressed data and ready to use.
#' @seealso \code{\link{getAffy}} to charge the expression data.
#' @examples 
#' \dontrun{
#' # Extract data from GEO DataSets (Takes time)
#' 
#' getInfo(GSE = "GSE8216", GPL = "GPL2025",directory = tempdir())
#' }

getInfo <- function(GSE,GPL,directory="."){
  
  # Download the raw data from GEO DataSets database
  
  sapply(as.vector(t(GSE)), function(n){
    getGEOSuppFiles(n,baseDir = directory)
  })
  
  # Obtaining the name of the compressed data
  
  files <- dir(directory)[grep("^GSE[0-9]",dir(directory),ignore.case = TRUE)]
  
  # Uncompress all the raw data
  
  for(j in GSE){
    untar(paste0(directory,"/",j,"/",j,"_RAW.tar"), exdir = paste0(directory,"/",j,"/"))
  }
  
  # Download the .soft file of the microarray chip
  
  getGEOfile(GPL, destdir = directory)
}