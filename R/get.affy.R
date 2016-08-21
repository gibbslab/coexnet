# get.affy
# Juan David Henao Sanchez
# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

get.affy <- function(GSE){
  
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