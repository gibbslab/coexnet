# create.net
# Juan David Henao Sanchez
# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

find.threshold <- function(difexp, method){
  
  simil <- .correlation.matrix(method)
  
  # Creates a sequence of threshold values
  
  pcv <- seq(0.01,0.99,by = 0.01)
  
  # Creates two empty vectors
  
  Cis <- vector()
  C0s <- vector()
  
  cat("First filter",sep = "\n")
  
  # Initial counter
  
  count <- 1
  
  for (val in pcv) {
    
    # Creates an empty matrix
    
    ady <- matrix(0,ncol = ncol(simil), nrow = nrow(simil))
    
    # Transforms to adjacency matrix
    
    for(i in 1:nrow(simil)){
      ady[which(simil[,i]>=val),i]<-1
      ady[which(simil[,i]<val),i]<-0
    }
    
    # Diagonal equal to zero
    
    diag(ady)<-0
    
    # Creates the network from the adjacency matrix
    
    G = graph.adjacency(ady,mode="undirected",diag=FALSE)
    
    # Obtains the clustering coefficient value
    
    Ci <- transitivity(G,type = "globalundirected")
    
    # If clustering coefficient cannot calculate, then equal to cero
    
    if(is.nan(Ci)){Ci <- 0}
    
    # Obtains values to calculate an artificial clustering coefficient
    
    K1 <- sum(degree(G,loops = F))
    K2 <- sum(degree(G,loops = F)^2)
    k1 <- (1/length(V(G)))*K1
    k2 <- (1/length(V(G)))*K2
    
    # Obtains a value to simulate clustering coefficient of a network
    # randomly created
    
    C0 <- ((k2-k1)^2)/(k1^3*length(V(G)))
    
    # If the result of simulate clustering coefficient is na, tranform to cero
    
    if(is.nan(C0)){C0 <- 0}
    
    # Adds the clustering coefficient to threshold values in a vector
    
    Cis[count] <- Ci
    
    # Adds the artificial clustering coefficient in a vector
    
    C0s[count] <- C0
    
    # Increases the counter
    
    count <- count + 1
    
  }
  
  thr <- abs(C0s-Cis)
  
  othr <- sort(thr,decreasing = T)
  
  dthr <- data.frame(pcv,thr)
  
  pvalue <- 0
  
  cntr <- 1
  
  while (pvalue <= 0.05) {
    
    mthr <- dthr[which(dthr$thr == othr[cntr]),][1]
    
    # Creates an empty matrix
    
    ad <- matrix(0,ncol = nrow(simil),nrow = nrow(simil))
    
    # Transforms to adjacency matrix
    
    for(i in 1:nrow(simil)){
      ad[which(simil[,i]>=mthr),i]<-1
      ad[which(simil[,i]<mthr),i]<-0
    }
    
    # Diagonal equal to zero
    
    diag(ad)<-0
    
    # Creates the network from the adjacency matrix
    
    gr=graph.adjacency(ad,mode="undirected",diag=FALSE)
    
    # Uses the function fit_power_law
    
    fit <- fit_power_law(degree(gr))
    
    pvalue <- fit$KS.p
    
    cntr <- cntr + 1
  }
  
  # Compares the clustering coefficients
  
  plot(pcv,abs(C0s-Cis),t="l",xlab = "Threshold",ylab = "| C0-Ci |")
  
  abline(v=mthr, col="red")
  
  text(x=0.1,y=0.25,paste0("Threshold = ", mthr))

}