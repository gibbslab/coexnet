#' @importFrom GEOquery getGEO Table getGEOSuppFiles getGEOfile
#' @importFrom affy ReadAffy expresso rma computeExprSet
#' @import siggenes
#' @importFrom acde stp
#' @import igraph
#' @import STRINGdb
#' @importFrom minet build.mim
#' @importFrom limma removeBatchEffect as.matrix.ExpressionSet normalizeVSN
#' @import Biobase
#' @importFrom graphics abline text
#' @importFrom stats cor median na.omit sd
#' @importFrom utils read.table untar

# internalfunctions
# Juan David Henao Sanchez
# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

.median.probe <- function(gene,array){
  
  colnames(array) <- gsub(".CEL.gz","",colnames(array),ignore.case = TRUE)
  list <- unique(gene$ID)
  newList <- list[grep(paste0("^","$"),list,invert = TRUE)]
  g <- matrix(0,length(newList),dim(array)[2])
  
  for(n in 1:nrow(g)){
    a <- array[as.vector(gene[grep(paste0("^",newList[n],"$"),gene$ID),]$probe),]
    a <- na.omit(a)
    if(!is.null(dim(a))){
      g[n,] <- apply(a,2,median)
    }else{
      g[n,] <- a
    }
  }
  rownames(g) <- newList
  colnames(g) <- colnames(array)
  
  return(g)
}

.max.probe <- function(gene,array){
  
  colnames(array) <- gsub(".CEL.gz","",colnames(array),ignore.case = TRUE)
  list <- unique(gene$ID)
  newList <- list[grep(paste0("^","$"),list,invert = TRUE)]
  g <- matrix(0,length(newList),dim(array)[2])
  
  for(n in 1:nrow(g)){
    a <- array[as.vector(gene[grep(paste0("^",newList[n],"$"),gene$ID),]$probe),]
    a <- na.omit(a)
    
    if(!is.null(dim(a))){
      g[n,] <- array[names(sort(rowMeans(a),decreasing = TRUE))[1],]
    }else{
      g[n,] <- a
    }
  }
  
  rownames(g) <- newList
  colnames(g) <- colnames(array)
  
  return(g)
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