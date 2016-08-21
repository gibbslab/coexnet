# get.info
# Juan David Henao Sanchez
# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

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