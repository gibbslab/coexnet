# expr.mat
# Juan David Henao Sanchez
# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

expr.mat <- function(affy,genes,NormalizeMethod,SummaryMethod){
  
  if(NormalizeMethod == "vsn"){
    
    # Normalizing with vsn method
    
    pvsn <- justvsn(affy)
    
    vsn <- computeExprSet(vsn,"pmonly","avgdiff")
    
    cat("Summarizing")
    
    if(SummaryMethod == "max"){
      
      # Summarizing using the high expression value
      
      eset <- .max.probe(vsn,genes) 
      
    }else if(SummaryMethod == "median"){
      
      # Summarizing using the median expression value
      
      eset <- .median.probe(gene,vsn)
    }
    
  }else if(NormalizeMethod == "rma"){
    
    # Normalizing using ram method
    
    rma <- rma(affy) 
    
    cat("Summarizing")
    
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