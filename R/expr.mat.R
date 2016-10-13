#' @export expr.mat
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Create a expression matrix 

expr.mat <- function(affy,genes,NormalizeMethod,SummaryMethod){
  
  if(NormalizeMethod == "vsn"){
    
    # Normalizing with vsn method
    
    pvsn <- normalize.AffyBatch.vsn(affy)
   
    vsn <- computeExprSet(x = pvsn,pmcorrect.method = "pmonly",summary.method = "avgdiff")
    
    cat("Summarizing")
    
    if(SummaryMethod == "max"){
      
      # Summarizing using the high expression value
      
      eset <- .max.probe(vsn,genes) 
      
    }else if(SummaryMethod == "median"){
      
      # Summarizing using the median expression value
      
      eset <- .median.probe(genes,vsn)
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