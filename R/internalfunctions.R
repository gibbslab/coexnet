# internalfunctions
# Juan David Henao Sanchez
# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

.median.probe <- function(gene,array){
  
  marray <- as.data.frame(exprs(array))
  
  names(marray) <- gsub(".CEL.gz","",names(marray),ignore.case = T)
  
  uni <- data.frame(gene$gene,marray)
  
  wowithw <- uni[grep(paste0("^","$"),uni$gene.gene,ignore.case = T,invert = T),]
  
  g <- data.frame()
  
  for(i in unique(na.omit(wowithw$gene.gene))){
    
    a <- as.data.frame(t(sapply(wowithw[grep(paste0("^",i,"$"),wowithw$gene.gene),
                                        2:ncol(wowithw)],median)),row.names = i)
    g <- rbind.data.frame(g,a)
  }
  
  return(g)
}

.max.probe <- function(array,gene){
  
  eset <- exprs(array)
  
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