# create.net
# Juan David Henao Sanchez
# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

create.net <- function(difexp, method){
  
  if(method == "corelation"){
    
    # Creates the Pearson correlation matrix
    
    simil <- abs(cor(t(difexp),use =  "pairwise.complete.obs"))
    
  }else if(method == "mutual information"){
    
    # Creates the mutual information matrix
    
    presimil <- build.mim(t(difexp), estimator = "mi.shrink", disc = "globalequalwidth")
    
    # Transforms data [0,inf) to [0,1] data
    
    simil<-sqrt(1-exp(-2*presimil))
    
    # Transforms the na data in 0
    
    simil[which(is.na(simil))]<-0
  }
  
  # Creates a sequence of threshold values
  
  pcv <- seq(0.01,0.99,by = 0.01)
  
  # Creates two empty vectors
  
  Cis <- vector()
  C0s <- vector()
  
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
    
    # Shows the progress of the procedure
    
    print(paste0(val*100,"%"))
  }
  
  # Creates an empty vector
  
  C <- vector()
  
  # Compares the clustering coefficients
  
  for (counter in 1:(length(pcv)-1)) {
    if((Cis[counter] - C0s[counter]) > (Cis[counter+1] - C0s[counter+1])){
      
      # Adds only the thresholds that pass the condition
      
      C[counter] <- pcv[counter]
    } 
  }
  
  # Omits the possible na values
  
  C <- na.omit(C)
  
  # p-value to K-S test
  
  fit <-0.05
  
  for (z in as.vector(C)) {
    
    # Creates an empty matrix
    
    ad <- matrix(0,ncol = nrow(simil),nrow = nrow(simil))
    
    # Transforms to adjacency matrix
    
    for(i in 1:nrow(simil)){
      ad[which(simil[,i]>=z),i]<-1
      ad[which(simil[,i]<z),i]<-0
    }
    
    # Diagonal equal to zero
    
    diag(ad)<-0
    
    # Creates the network from the adjacency matrix
    
    gr=graph.adjacency(ad,mode="undirected",diag=FALSE)
    
    # Uses the function fit_power_law
    
    pvalue <- fit_power_law(degree(gr))
    
    # Shows the threshold value
    
    print(paste("Threshold:",z))
    
    # Shows the p-value to K-S test in the fit_power_law function
    
    print(paste("p-value:",pvalue$KS.p))
    
    # Obtains the higher p-value and the corresponding threshold value
    
    if(pvalue$KS.p >= fit ){
      fit <- pvalue$KS.p
      value <- z                     
    }
  }
  
  # Creates an empty matrix
  
  Ad <- matrix(0,ncol = nrow(simil),nrow = nrow(simil))
  
  # Transforms to adjacency matrix
  
  for(i in 1:nrow(simil)){
    Ad[which(simil[,i]>=value),i]<-1
    Ad[which(simil[,i]<value),i]<-0
  }
  
  # Changes the names of adjacency matrix to the genes names
  
  colnames(Ad)<-rownames(Ad)<-rownames(simil)
  
  # Diagonal equal to zero
  
  diag(Ad)<-0
  
  # Creates the network from the adjacency matrix
  
  Gr=graph.adjacency(Ad,mode="undirected",add.colnames=NULL,diag=FALSE)
  
  # Shows the final values of the co-expression network
  
  print("")
  print("Final Network")
  print("")
  print(paste("p-value =",fit,sep = " "))
  print(paste("threshold =",value,sep = " "))
  
  # Returns the network as an igraph object
  
  return(Gr) 
}