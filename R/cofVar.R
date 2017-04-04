# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @export cofVar
#' @author Juan David Henao <judhenaosa@unal.edu.co>
#' @title Calculating the coefficient of variation for expression matrix.
#' @description This function calculates the mean and the coefficient of variation to each row 
#' (genes or probesets) in an expression matrix in two ways: i) in the whole matrix ii) for the specific phenotype (case or control).
#' @param data The whole normalized expression matrix, rows: genes or probeset, columns: samples.
#' @param complete Boolean to define if the function uses the whole expression matrix, by default TRUE.
#' @param treatment A numeric vector with 0s and 1s for each sample in the expression matrix, the 0 expresses the control samples and 1 expresses the case samples, by default is NULL.
#' @param type It can be "case" to calculate the mean and the coefficient of variation for the case samples or, otherwise, "control" to obtain these two values for the control samples.
#' @return The expression matrix with two new columns, the first one with the averages and the other one with the coefficient of variation values.
#' @examples 
#' 
#' ## Creating the expression matrix
#' 
#' # The matrix have 200 genes and 20 samples
#' 
#' n <- 200
#' m <- 20
#' 
#' # The vector with treatment and control samples
#' 
#' t <- c(rep(0,10),rep(1,10))
#' 
#' # Calculating the expression values normalized
#' 
#' mat <- as.matrix(rexp(n, rate = 1))
#' norm <- t(apply(mat, 1, function(nm) rnorm(m, mean=nm, sd=1)))
#' 
#' ## Calculating the mean and the coefficient of variation
#' 
#' # For the whole expression matrix
#' 
#' complete <- cofVar(norm)
#' head(complete)
#' 
#' # For the case samples
#' 
#' case <- cofVar(data = norm,complete = FALSE,treatment = t,type = "case")
#' head(case)


cofVar <- function(data,complete=TRUE,treatment=NULL,type=NULL){
  
  # Replacing the sample name for case/control ID
  
  if (complete == FALSE) {
    
    data <- as.data.frame(data)
    
    names(data) <- treatment
    
    if(type == "control"){
      
      type <- "0"
    }else if(type == "case"){
      
      type <- "1"
    }
    tdata <- data[names(data) == type]
  }else{
    tdata = as.data.frame(data)
  }
  
  # Obtaining the mean of the expression values
  
  tdata$mean <- rowMeans(tdata, na.rm = TRUE)
  
  # Obtaining the coefficient of variation to each row
  
  CV <- function(x){sd(x,na.rm = TRUE)/mean(x,na.rm = TRUE)}
  
  # Apply the function to the expression matrix
  
  tdata$cv <- apply(tdata[,1:(ncol(tdata)-1)],1,CV)
  
  # Return the expression matrix with the mean and the CV
  
  return(tdata)
}