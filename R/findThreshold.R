# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @export findThreshold
#' @author Juan David Henao Sanchez <judhenaosa@unal.edu.co>
#' @author Liliana Lopez Kleine <llopezk@unal.edu.co>
#' @title Find the threshold value to create a co-expression network
#' @description Finds the threshold value to establish the cutoff in the process to define the edges in the co-expression network final from two steps. In the first one, obtains the subtraction from clustering coefficient values of the real and random networks created from the possible threshold values in the correlation matrix. In the second one, a Kolmogorov-Smirnov test is made to evaluate the degree distribution respect normality.
#' @param difexp A whole expression matrix or the expression matrix to differentially expressed genes.
#' @param method The method name to create the correlation matrix, this can be "correlation" to obtain the Pearson Correlation Coefficient. On the other hand, can be "mutual information" to obtain the correlation values from an entropy-based method.
#' @return The best threshold value found using the two criteria and a plot showing the result.
#' @seealso \code{\link{difExprs}} to find the differentially expressed genes matrix.
#' @references Elo, L. L., Jarvenpaa, H., Oresic, M., Lahesmaa, R., & Aittokallio, T. (2007). Systematic construction of gene coexpression networks with applications to human T helper cell differentiation process. Bioinformatics, 23(16), 2096-2103.
#' @examples 
#' 
#' # Loading data
#' 
#' pathfile <- system.file("extdata","expression_example.txt",package = "coexnet")
#' data <- read.table(pathfile,stringsAsFactors = FALSE)
#' 
#' # Finding threshold value
#' 
#' cor_pearson <- findThreshold(difexp = data,method = "correlation")
#' cor_pearson
 
findThreshold <- function(difexp, method){
  
  # Obtaining the similarity values
  
  simil <- .correlation.matrix(difexp,method)
  
  # Create a sequence of threshold values
  
  pcv <- seq(0.01,0.99,by = 0.01)
  
  # Create two empty vectors
  
  Cis <- vector()
  C0s <- vector()
  
  # Initial counter
  
  count <- 1
  
  for (val in pcv) {
    
    # Create an empty matrix
    
    ady <- matrix(0,ncol = ncol(simil), nrow = nrow(simil))
    
    # Transform to adjacency matrix
    
    for(i in 1:nrow(simil)){
      ady[which(simil[,i]>=val),i]<-1
      ady[which(simil[,i]<val),i]<-0
    }
    
    # Diagonal equal to zero
    
    diag(ady)<-0
    
    # Create the network from the adjacency matrix
    
    G = graph.adjacency(ady,mode="undirected",diag=FALSE)
    
    # Obtaining the clustering coefficient value
    
    Ci <- transitivity(G,type = "globalundirected")
    
    # If clustering coefficient cannot calculate, then equal to zero
    
    if(is.nan(Ci)){Ci <- 0}
    
    # Obtaining values to calculate an artificial clustering coefficient
    
    K1 <- sum(degree(G,loops = FALSE))
    K2 <- sum(degree(G,loops = FALSE)^2)
    k1 <- (1/length(V(G)))*K1
    k2 <- (1/length(V(G)))*K2
    
    # Obtaining a value to simulate clustering coefficient of a network randomly created
    
    C0 <- ((k2-k1)^2)/(k1^3*length(V(G)))
    
    # If the result of simulate clustering coefficient is NA, then transform to zero
    
    if(is.nan(C0)){C0 <- 0}
    
    # Add the clustering coefficient for threshold values in a vector
    
    Cis[count] <- Ci
    
    # Add the artificial clustering coefficient in a vector
    
    C0s[count] <- C0
    
    # Increase the counter
    
    count <- count + 1
    
  }
  
  thr <- vector()
  pass <- vector()
  
  # Finding the subtraction between the clustering coefficient of random network and the real network
  
  for(i in 1:(length(pcv)-1)){
    if(Cis[i]-C0s[i] > Cis[i+1]-C0s[i+1]){
     thr[i] <- pcv[i]
     pass[i] <- Cis[i]*100-C0s[i]*100
    }
  }
  
  # Delete NA values
  
  thr <- na.omit(thr)
  
  # Round the result of the difference between the clustering coefficients
  
  pass <- round(na.omit(pass))
  
  thre <- vector()
  
  # Delete the minimum values of the difference between the clustering coefficients
  
  for(j in 1:length(thr)){
    if(pass[j] > min(pass)){
      thre[j] <- thr[j]
    }
  }
  
  # Delete the NA values
  
  thre <- na.omit(thre)
  
  mtr <- vector()
  
  for(n in 1:length(thre)){
    
    ad <- matrix(0,ncol = nrow(simil),nrow = nrow(simil))
    
    # Transforming into adjacency matrix
    
    for(i in 1:nrow(simil)){
      ad[which(simil[,i]>=thre[n]),i]<-1
      ad[which(simil[,i]<thre[n]),i]<-0
    }
    
    # Diagonal equal to zero
    
    diag(ad)<-0
    
    # Create the network from the adjacency matrix
    
    gr=graph.adjacency(ad,mode="undirected",diag=FALSE)
    
    # Use the function fit_power_law
    
    fit <- fit_power_law(degree(gr))
    
    # Obtaining the p-value
    
    pvalue <- fit$KS.p
    
    # Add the threshold value if p-value > 0.05
    
    if(pvalue > 0.05){
      mtr[n] <- thre[n]
    }
  }
  
  # Delete NA values
  
  mtr <- na.omit(mtr)
  
  # Delete the minimum values in the subtraction of clustering coefficient values
  
  if(length(mtr) == 0){
    new_pass <- vector()
    Cls <- round(Cis-C0s,digits = 3)
    for(j in 1:length(pcv)){
      if(Cls[j] > min(Cls)){
        mtr[j] <- pcv[j]
     }
    }
  }
  
  # Delete NA values
  
  mtr <- na.omit(mtr)
  
  # Create a plot comparing the clustering coefficients 
  
  plot(pcv,abs(Cis-C0s),t="l",xlab = "Threshold",ylab = "| Ci - C0 |")
  
  # Create the line to show the final threshold value
  
  abline(v=mtr[1], col="red")
  
  # Text to complete the plot
  
  text(min(pcv)+0.1,max(abs(C0s-Cis))-0.1,paste0("Threshold = ", mtr[1]))
  
  # return the value corresponding to the final threshold value
  
  return(mtr[1])
  
}
