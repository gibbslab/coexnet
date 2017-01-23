#' @export find.threshold
#' @author Juan David Henao Sanchez <judhenaosa@unal.edu.co>
#' @author Liliana Lopez Kleine <llopezk@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Find the threshold value to create a co-expression network
#' @description Find the threshold value to establish the edges in a co-expression network based on correlation matrix
#' using two different criteria. First, use the differential value between the clustering coefficient for the real network and
#' random network build and second test the network using a Kolgomorov-Smirnov test to reject normality in the degree distribution of 
#' the edges in the network.
#' @param difexp  A whole expression matrix or the expression matrix to differentially expressed genes.
#' @param method  The method to create the correlation matrix, can be "correlation" to Pearson test or "mutual information" to test 
#' based on entropy information.
#' @return The best threshold value find using the two criteria and a plot with the result.
#' @seealso \code{\link{dif.exprs}} to find the differential expressed genes matrix.
#' @references Elo, L. L., Järvenpää, H., Orešič, M., Lahesmaa, R., & Aittokallio, T. (2007). Systematic construction of gene coexpression networks with applications to human T helper cell differentiation process. Bioinformatics, 23(16), 2096-2103.
#' @examples 
#' 
#' # Loading data
#' 
#' data <- read.table(system.file("extdata","expression_example.txt",package = "coexnet"),stringsAsFactors = F)
#' 
#' # Find threshold value
#' 
#' cor_pearson <- find.threshold(difexp = data,method = "correlation")
#' cor_pearson
#' 
#' mut_inf <- find.threshold(difexp = data,method = "mutual information")
#' mut_info
 
find.threshold <- function(difexp, method){
  
  # Obtains the similarity values
  
  simil <- .correlation.matrix(difexp,method)
  
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
    
  }
  
  thr <- vector()
  pass <- vector()
  
  # Find the differencies between the clustering coefficient to random network anthe real network
  
  for(i in 1:(length(pcv)-1)){
    if(Cis[i]-C0s[i] > Cis[i+1]-C0s[i+1]){
     thr[i] <- pcv[i]
     pass[i] <- Cis[i]*100-C0s[i]*100
    }
  }
  
  # Deletes NA values
  
  thr <- na.omit(thr)
  
  # Rounds the result of the difference between the clustering coefficients.
  
  pass <- round(na.omit(pass))
  
  thre <- vector()
  
  # Deletes the minimun values of the difference between the clustering coefficients.
  
  for(j in 1:length(thr)){
    if(pass[j] > min(pass)){
      thre[j] <- thr[j]
    }
  }
  
  # Deletes the NA values
  
  thre <- na.omit(thre)
  
  mtr <- vector()
  
  for(n in 1:length(thre)){
    
    ad <- matrix(0,ncol = nrow(simil),nrow = nrow(simil))
    
    # Transforms to adjacency matrix
    
    for(i in 1:nrow(simil)){
      ad[which(simil[,i]>=thre[n]),i]<-1
      ad[which(simil[,i]<thre[n]),i]<-0
    }
    
    # Diagonal equal to zero
    
    diag(ad)<-0
    
    # Creates the network from the adjacency matrix
    
    gr=graph.adjacency(ad,mode="undirected",diag=FALSE)
    
    # Uses the function fit_power_law
    
    fit <- fit_power_law(degree(gr))
    
    # Obtaines the p-value
    
    pvalue <- fit$KS.p
    
    # Adds the threshold value if p-value > 0.05
    
    if(pvalue > 0.05){
      mtr[n] <- thre[n]
    }
  }
  
  # Deletes NA values
  
  mtr <- na.omit(mtr)
  
  # Deletes the minimun values of the difference between the clustering coefficients.
  
  if(length(mtr) == 0){
    new_pass <- vector()
    Cls <- round(Cis-C0s,digits = 3)
    for(j in 1:length(pcv)){
      if(Cls[j] > min(Cls)){
        mtr[j] <- pcv[j]
     }
    }
  }
  
  # Deletes NA values
  
  mtr <- na.omit(mtr)
  
  # Creates a plot comparing the clustering coefficients 
  
  plot(pcv,abs(Cis-C0s),t="l",xlab = "Threshold",ylab = "| Ci - C0 |")
  
  # Creates the line to show the final threshold value
  
  abline(v=mtr[1], col="red")
  
  # Text to complete the plot
  
  text(min(pcv)+0.1,max(abs(C0s-Cis))-0.1,paste0("Threshold = ", mtr[1]))
  
  # return the value corrisponding to the final threshold value
  
  return(mtr[1])
  
}
