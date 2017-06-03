#' @importFrom GEOquery getGEO Table getGEOSuppFiles getGEOfile
#' @importFrom affy ReadAffy expresso rma computeExprSet
#' @importFrom siggenes sam plot summary
#' @importFrom acde stp
#' @importFrom SummarizedExperiment assay makeSummarizedExperimentFromExpressionSet
#' @importFrom igraph graph.intersection V graph.edgelist as_data_frame graph_from_data_frame degree delete.vertices diameter graph.adjacency transitivity fit_power_law read.graph
#' @importFrom minet build.mim
#' @importFrom limma removeBatchEffect as.matrix.ExpressionSet normalizeVSN
#' @importFrom graphics abline text
#' @importFrom stats cor median na.omit sd
#' @importFrom utils read.table untar
#' @importFrom vsn vsn2
#' @import Biobase
#' @import rmarkdown
#' @import STRINGdb

# internalfunctions
# Juan David Henao Sanchez
# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

.median.probe <- function(gene,array){
  
  colnames(array) <- gsub(".CEL.gz","",colnames(array),ignore.case = TRUE)
  list <- unique(gene$ID)
  newList <- list[grep(paste0("^","$"),list,invert = TRUE)]
  
  g <- matrix(0,length(newList),dim(array)[2])
  
  g <- t(vapply(seq_len(nrow(g)), function(n){
    a <- array[as.vector(gene[grep(paste0("^",newList[n],"$"),gene$ID),]$probe),]
    a <- na.omit(a)
    if(!is.null(dim(a))){
      g[n,] <- apply(a,2,median)
    }else{
      g[n,] <- a
    }
  },rep(0.0,ncol(array))))
  
  rownames(g) <- newList
  colnames(g) <- colnames(array)
  
  g <- na.omit(g)
  
  return(makeSummarizedExperimentFromExpressionSet(ExpressionSet(g)))
}

.max.probe <- function(gene,array){
  
  colnames(array) <- gsub(".CEL.gz","",colnames(array),ignore.case = TRUE)
  list <- unique(gene$ID)
  newList <- list[grep(paste0("^","$"),list,invert = TRUE)]
  
  g <- matrix(0,length(newList),dim(array)[2])
  
  g <- t(vapply(seq_len(nrow(g)), function(n){
    a <- array[as.vector(gene[grep(paste0("^",newList[n],"$"),gene$ID),]$probe),]
    a <- na.omit(a)
    
    if(!is.null(dim(a))){
      if(dim(a)[1] > 0){
        g[n,] <- array[names(sort(rowMeans(a),decreasing = TRUE))[1],] 
      }else{
        g[n,] <- rep(NA,ncol(a))
      }
    }else{
      g[n,] <- a
    }
  },rep(0.0,ncol(array))))
  
  rownames(g) <- newList
  colnames(g) <- colnames(array)
  
  g <- na.omit(g)
  
  return(makeSummarizedExperimentFromExpressionSet(ExpressionSet(g)))
}

.correlation.matrix <- function(difexp,method){
  
  if(method == "correlation"){
    
    simil <- abs(cor(t(difexp),use =  "pairwise.complete.obs"))
    
  }else if(method == "mutual information"){
    
    presimil <- build.mim(t(difexp), estimator = "mi.shrink", disc = "globalequalwidth")
    
    simil<-sqrt(1-exp(-2*presimil))
    
    simil[which(is.na(simil))]<-0
  }
  
  return(simil)
  
}