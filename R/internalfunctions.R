#' @importFrom GEOquery getGEO Table getGEOSuppFiles getGEOfile
#' @importFrom vsn normalize.AffyBatch.vsn exprs
#' @importFrom affy ReadAffy expresso rma computeExprSet
#' @import siggenes
#' @importFrom acde stp
#' @import igraph
#' @importFrom minet build.mim
#' @importFrom limma removeBatchEffect
#' @import Biobase

# internalfunctions
# Juan David Henao Sanchez
# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

.median.probe <- function(gene,array){
  
  colnames(array) <- gsub(".CEL.gz","",colnames(array),ignore.case = T)
  list <- unique(gene$gene)
  newList <- list[grep(paste0("^","$"),list,invert = T)]
  g <- matrix(0,length(newList),dim(array)[2])
  
  for(n in 1:nrow(g)){
    a <- batch[as.vector(gene[grep(paste0("^",newList[n],"$"),gene$gene),]$probe),]
    a <- na.omit(a)
    if(is.null(dim(a))){
      g[n,] <- a
    }else{
      g[n,] <- apply(a,2,median)
    }
    print(n)
  }
  rownames(g) <- newList
  colnames(g) <- colnames(array)
  
  return(g)
}

.max.probe <- function(gene,array){
  
  eset <- exprs(array)
  
  #eset <- array
  
  rows <- rowMeans(eset,na.rm = T)
  
  da <- data.frame(gene,stringsAsFactors = F)
  
  names(da) <- c("a","b")
  
  probemean <- data.frame(names(rows),rows, stringsAsFactors = F)
  
  names(probemean) <- c("a","b")
  
  merge <- merge.data.frame(probemean, da, by.x = "a", by.y = "a")
  
  dat <- data.frame(merge$b.y,merge$b.x, row.names = merge$a,stringsAsFactors = F)
  
  total <- dat[-c(grep(paste0("^","$"),dat[,1])),]
  
  onlygenes <- unique(total[,1])
  
  result <- data.frame()
  
  for(x in onlygenes){
    
    genes  <- total[grep(x,total[,1],fixed = T),]
    
    max <- genes[grep(max(total[grep(paste0("^",x,"$"),total[,1]),2]),genes[,2]),]
    
    if(length(result) == 0){
      result <- rbind(max)
    }else{
      result <- rbind(result,max)
    }
  }
  
  eset2 <- as.data.frame(eset)
  
  final <- eset2[row.names(result),]
  
  row.names(final) <- result$merge.b.y
  
  return(final)
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