# expr.mat
# Juan David Henao Sanchez
# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

expr.mat <- function(affy,genes,NormalizeMethod,SummaryMethod){
  
  if(NormalizeMethod == "vsn"){
    
    # Normalizing with vsn method
    
    vsn <- expresso(affy,pmcorrect.method = "pmonly", bg.correct = F,
                    normalize.method = "vsn", summary.method = "avgdiff")
    
    print("summarizing")
    
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
    
    print("summarizing")
    
    if(SummaryMethod == "max"){
      
      # Summarizing using the high expression value
      
      eset <- ProbeFilter(rma,genes) 
      
    }else if(SummaryMethod == "median"){
      
      # Summarizing using the median expression value
      
      eset <- medianProbe(genes,rma)
    }
  }
  
  # Returns the gene expression matrix
  
  return(eset)
}