# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @export difExprs
#' @author Juan David Henao <judhenaosa@unal.edu.co>
#' @title Differential expression analysis using two different methods.
#' @description Using the expression matrix calculate the differential expressed
#' genes to two class analysis and fixing an expected FDR value. The methods
#' are SAM and ACDE.
#' @param expData A matrix with the expression matrix, it may be stored in a SummarizedExperiment object.
#' @param treatment A vector with the ientifiers of the classes, 0 to control
#' and 1 to case.
#' @param fdr The expected FDR value.
#' @param DifferentialMethod The method to calculate the differential expressed
#' genes, can be "sam" or "acde"
#' @param plotting The option to show the result in a plot. By default FALSE.
#' @return A data.frame with the expression matrix to the expressed diferential 
#' genes only.
#' @seealso \code{\link{exprMat}} to obtain the expression matrix.
#' @references Tusher, V. G., Tibshirani, R., & Chu, G. (2001). Significance analysis of 
#' microarrays applied to the ionizing radiation response. Proceedings of the National Academy of Sciences, 98(9), 5116-5121.
#' @references Acosta J and Lopez-Kleine L (2015). acde: Artificial Components Detection of Differentially Expressed Genes. R package version 1.4.0.
#' @examples 
#' 
#' ## Creating the expression matrix
#' 
#' # The matrix have 200 genes and 20 samples
#' 
#' n <- 200
#' m <- 20
#' 
#' # The vector with treatment samples and control samples
#' 
#' treat <- c(rep(0,10),rep(1,10))
#' 
#' # Calculating the expression values normalized
#' 
#' mat <- as.matrix(rexp(n, rate = 1))
#' norm <- t(apply(mat, 1, function(nm) rnorm(m, mean=nm, sd=1)))
#' 
#' ## Running the function using the two approaches
#' 
#' sam <- difExprs(expData = norm,treatment = treat,fdr = 0.2,DifferentialMethod = "sam")
#' acde <- difExprs(expData = norm,treatment = treat,fdr = 0.2,DifferentialMethod = "acde")

difExprs <- function(expData,treatment,fdr,DifferentialMethod,plotting=FALSE){
  
  # Identifing the SummarizedExperiment object
  
  if(is(expData,"SummarizedExperiment")){
    expData <- assay(expData)
  }
  
  if(DifferentialMethod == "sam"){
    
    # Differential analysis using sam method
    
    samr <- sam(data = expData,cl = treatment,B=100,rand=100)
    
    # Obtaining the fdr to different thresholds
    
    tab <- as.data.frame(samr@mat.fdr)
    
    # Filter the fdr values from the expected value
    
    tab <- tab[tab$FDR >= fdr,]
    
    # Message to empty result
    
    if(nrow(tab) == 0){stop("No differentially expressed genes found")}
    
    # Obtaining the threshold value
    
    value <- tab[nrow(tab),]
    
    # Showing the result of differential analysis
    
    if(plotting == TRUE){
      plot(samr,value$Delta)
    }
    
    # Summarizing the results of differential analysis
    
    sum <- summary(samr,value$Delta,entrez=FALSE)
    
    # Obtaining the names of genes differentially expressed
    
    dife <- sum@row.sig.genes
    
    # Filter the expression matrix with the genes differentially expressed
    
    genes <- expData[dife,]
    
    # Showing the achieved fdr
    
    cat(paste0("Achieved FDR: ",value$FDR))
    
  }else if(DifferentialMethod == "acde"){
    
    # Change the value of cases and controls
    
    treatment[treatment == 0] <- 2
    
    # Differential analysis using acde method
    
    acde <- stp(expData,treatment,R = 100, PER = TRUE,alpha = fdr)
    
    # Showing the result of differential analysis
    
    if(plotting == TRUE){
      plot(acde) 
    }
    
    # Showing the achieved fdr
    
    cat(paste0("Achieved FDR: ",acde$astar))
    
    # Create a data.frame object with the result of differential analysis
    
    lis <- data.frame(acde$gNames, acde$dgenes)
    
    # Obtaining the name of genes differentially expressed
    
    diff <- lis[lis$acde.dgenes != "no-diff.",]
    
    # Filter the expression matrix with the genes differentially expressed
    
    genes <- expData[diff$acde.gNames,]
  }
  
  # Return the matrix of the genes differentially expressed
  
   return(genes)
}
