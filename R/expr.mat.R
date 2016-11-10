#' @export expr.mat
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title gggg

expr.mat <- function(affy,genes,NormalizeMethod,SummaryMethod,BatchCorrect = TRUE){
  
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
  
  if(NormalizeMethod == "vsn"){
    
    # Normalizing with vsn method
    
    pvsn <- as.matrix.ExpressionSet(affy)
   
    norm <- normalizeVSN(pvsn) 
    
    exprs(affy) <- norm
    
    vsn <- computeExprSet(x = affy,pmcorrect.method = "pmonly",summary.method = "avgdiff")
    
    if (BatchCorrect == TRUE) {
      batch <- removeBatchEffect(vsn,batch.dates)
    }else{
      batch <- as.matrix.ExpressionSet(vsn)
    }
    
    cat("Summarizing",sep = "\n")
    
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