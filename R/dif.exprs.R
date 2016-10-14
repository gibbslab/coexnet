#' @export dif.exprs
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Differential expression analysis using two different methods.
#' @description Using the expression matrix calculate the differential expressed
#' genes to two class analysis and fixing an expected FDR value. The methods
#' are SAM and ACDE.
#' @param eset A data.frame with the expression matrix.
#' @param treatment A vector with the ientifiers of the classes, 0 to control
#' and 1 to case.
#' @param fdr The expected FDR value.
#' @param DifferentialMethod The method to calculate the differential expressed
#' genes, can be "sam" or "acde"
#' @return A data.frame with the expression matrix to the expressed diferential 
#' genes only.
#' @seealso \code{\link{expr.mat}} to obtain the expression matrix.

dif.exprs <- function(eset,treatment,fdr,DifferentialMethod){
  
  if(DifferentialMethod == "sam"){
    
    # Differential analysis using sam method
    
    samr <- sam(data = eset,cl = treatment,B=100,rand=100)
    
    # Obtains the fdr to different thresholds
    
    tab <- as.data.frame(samr@mat.fdr)
    
    # Filters the fdr values from the expected value
    
    tab <- tab[tab$FDR >= fdr,]
    
    # Message to empty result
    
    if(nrow(tab) == 0){stop("No differentially expressed genes found")}
    
    # Obtains the threshold value
    
    value <- tab[nrow(tab),]
    
    # Shows the result of differential analysis
    
    plot(samr,value$Delta)
    
    # Summarizing the results of differential analysis
    
    sum <- summary(samr,value$Delta,entrez=F)
    
    # Obtains the names of genes differentially expressed
    
    dife <- sum@row.sig.genes
    
    # Filters the expression matrix with the genes differentially expressed
    
    genes <- eset[dife,]
    
    # Shows the achieved fdr
    
    cat(paste0("Achieved FDR: ",value$FDR))
    
  }else if(DifferentialMethod == "acde"){
    
    # Changes the value of cases and controls
    
    treatment[treatment == 0] <- 2
    
    # Differential analysis using acde method
    
    acde <- stp(eset,treatment,R = 100, PER = T,alpha = fdr)
    
    # Shows the result of differential analysis
    
    plot(acde)
    
    # Shows the achieved fdr
    
    cat(paste0("Achieved FDR: ",acde$astar))
    
    # Shows the threshold value
    
    cat(paste0("delta value: ",acde$tstar))
    
    # Creates a data.frame object with the result of differential analysis
    
    list <- data.frame(acde$gNames, acde$dgenes)
    
    # Obtains the name of genes differentially expressed
    
    diff <- list[list$acde.dgenes != "no-diff.",]
    
    # Filters the expression matrix with the genes differentially expressed
    
    genes <- eset[diff$acde.gNames,]
  }
  
  # Returns the matrix of the genes differentially expressed
  
   return(genes)
}
