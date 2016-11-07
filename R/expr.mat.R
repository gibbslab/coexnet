#' @export expr.mat
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title gggg

expr.mat <- function(affy,genes,NormalizeMethod,SummaryMethod){
  
  if(NormalizeMethod == "vsn"){
    
    # Normalizing with vsn method
    
    pvsn <- normalize.AffyBatch.vsn(affy)
   
    vsn <- computeExprSet(x = pvsn,pmcorrect.method = "pmonly",summary.method = "medianpolish")
    
    dates <- protocolData(affy)$ScanDate
    
    strdates <- strsplit(dates," ")
    
    batch.dates <- vector()
    
    for (i in 1:length(strdates)) {
      batch.dates[i]  <- strdates[[i]][1]
    }
    
    tab <-names(table(batch.dates))
    
    for (n in 1:length(tab)) {
      batch.dates[batch.dates == tab[n]] <- paste0("b", n)
    }
    
    batch <- removeBatchEffect(vsn,batch.dates)
    
    cat("Summarizing")
    
    if(SummaryMethod == "max"){
      
      # Summarizing using the high expression value
      
      eset <- .max.probe(batch,genes) 
      
    }else if(SummaryMethod == "median"){
      
      # Summarizing using the median expression value
      
      eset <- .median.probe(genes,batch)
    }
    
  }else if(NormalizeMethod == "rma"){
    
    # Normalizing using ram method
    
    rma <- rma(affy) 
    
    cat("Summarizing",sep = "\n")
    
    if(SummaryMethod == "max"){
      
      # Summarizing using the high expression value
      
      eset <- .max.probe(rma,genes) 
      
    }else if(SummaryMethod == "median"){
      
      # Summarizing using the median expression value
      
      eset <- .median.probe(genes,rma)
    }
  }
  
  cat("DONE")
  
  # Returns the gene expression matrix
  
  return(eset)
}