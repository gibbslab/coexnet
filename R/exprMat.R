# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @export exprMat
#' @author Juan David Henao <judhenaosa@unal.edu.co>
#' @title Calculate the expression matrix from the raw expression data.
#' @description This function use a affyBatch object with the raw expression data to normalize and transform the matrix from
#' probeset to gene considering the option to remove the batch effect in the long microarray data.
#' @param affy  A AffyBatch object with the raw expression data.
#' @param genes  A table with two columns, in the firt one the name of each probe 
#' in the microarray with header "probe" and in the second one the corresponding gene or ID with header "ID".
#' @param NormalizeMethod  The method to normalize the raw data. Can be "vsn" to apply 
#' Variance Stabilizing Normalization function or "rma" to apply Robust Multi-Array Average function.
#' @param SummaryMethod  The method to pass from probeset to genes or ID. Can be "max" to selecto the probeset with the most
#' average expression value or "median" to obtain the median of each sample of the set of probeset corresponding to particular
#' gen or ID.
#' @param BatchCorrect  The option to apply batch effect correction, by default TRUE.
#' @return A SummarizedExperiment object with the expression matrix.
#' @seealso \code{\link{getAffy}} to obtain the affyBatch object.
#' @seealso \code{\link{geneSymbol}} to obtain the data.frame with probeset and genes/ID from .SOFT file.
#' @references Huber, W., Von Heydebreck, A., Sultmann, H., Poustka, A., & Vingron, M. (2002). Variance stabilization applied to microarray data calibration and to the quantification of differential expression. Bioinformatics, 18(suppl 1), S96-S104.
#' @references Irizarry, R. A., Hobbs, B., Collin, F., Beazer Barclay, Y. D., Antonellis, K. J., Scherf, U., & Speed, T. P. (2003). 
#' Exploration, normalization, and summaries of high density oligonucleotide array probe level data. Biostatistics, 4(2), 249-264.
#' @examples 
#' \dontrun{
#' 
#' # Loading data
#' 
#' if (require(affydata)) {
#'  data(Dilution)
#' }
#' 
#' # Loading table with probeset and gene/ID information
#' 
#' data(info)
#' 
#' # Calculating the expression matrix
#' 
#' ## RMA
#' 
#' rma <- exprMat(affy = Dilution,genes = info,NormalizeMethod = "rma",
#' SummaryMethod = "median",BatchCorrect = FALSE)
#' head(rma)
#' 
#' ## VSN
#' 
#' vsn <- exprMat(affy = Dilution,genes = info,NormalizeMethod = "vsn",
#' SummaryMethod = "median",BatchCorrect = FALSE)
#' head(vsn)
#' 
#' }

exprMat <- function(affy,genes,NormalizeMethod,SummaryMethod,BatchCorrect = TRUE){
  
  if(NormalizeMethod == "vsn"){
    
    # Normalizing with vsn method
    
    cat("Normalizing", sep = "\n")
    
    # Extract the raw expression matrix from AffyBatch object
    
    pvsn <- as.matrix.ExpressionSet(affy)
    
    # Normalize using vsn
   
    norm <- normalizeVSN(pvsn, verbose = FALSE)
    
    # Replace the raw expression matrix from the affyBaych object with the normalized expression matrix
    
    exprs(affy) <- norm
    
    # Pass from probes to probeset
    
    vsn <- computeExprSet(x = affy,pmcorrect.method = "pmonly",summary.method = "avgdiff")
    
    # Conditional to apply batch effect correction
    
    if (BatchCorrect == TRUE) {
      # Extract the date of microarray scan
      
      dates <- protocolData(affy)$ScanDate
      
      # Spliting all the information related to date of scan
      
      strdates <- strsplit(dates," ")
      
      # List to vector
      
      batch.dates <- vapply(seq_len(length(strdates)),function(i){strdates[[i]][1]},1)
      
      # Obtaining the unique dates
      
      tab <-names(table(batch.dates))
      
      # Joining samples in batchs according the date of scan
      
      for (n in seq_len(length(tab))) {
        batch.dates[batch.dates == tab[n]] <- paste0("b", n)
      }
      
      batch <- removeBatchEffect(vsn,batch.dates)
    }else{
      batch <- as.matrix.ExpressionSet(vsn)
    }
    
    cat("Summarizing",sep = "\n")
    
    if(SummaryMethod == "max"){
      
      # Summarizing using the high expression value
      
      eset <- .max.probe(genes,batch) 
      
    }else if(SummaryMethod == "median"){
      
      # Summarizing using the median expression value
      
      eset <- .median.probe(genes,batch)
    }
    
  }else if(NormalizeMethod == "rma"){
    
    # Normalizing using ram method
    
    rma <- rma(affy)
    
    # Conditional to apply batch effect correction
    
    if (BatchCorrect == TRUE) {
      # Extract the date of microarray scan
      
      dates <- protocolData(affy)$ScanDate
      
      # Spliting all the information related to date of scan
      
      strdates <- strsplit(dates," ")
      
      # List to vector
      
      batch.dates <- vapply(seq_len(length(strdates)),function(i){strdates[[i]][1]},'character')
      
      # Obtaining the unique dates
      
      tab <-names(table(batch.dates))
      
      # Joining samples in batchs according the date of scan
      
      for (n in seq_len(length(tab))) {
        batch.dates[batch.dates == tab[n]] <- paste0("b", n)
      }
      
      batch <- removeBatchEffect(rma,batch.dates)
    }else{
      batch <- as.matrix.ExpressionSet(rma)
    }
    
    cat("Summarizing",sep = "\n")
    
    if(SummaryMethod == "max"){
      
      # Summarizing using the high expression value
      
      eset <- .max.probe(genes,batch) 
      
    }else if(SummaryMethod == "median"){
      
      # Summarizing using the median expression value
      
      eset <- .median.probe(genes,batch)
    }
  }
  
  cat("DONE")
  
  # Return the gene expression matrix
  
  return(eset)
}
