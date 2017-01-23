#' @export cof.var
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Calculate the coefficient of variation to expression matrix.
#' @description This function calculate the mean and the coefficient of variation to each one of the row (gene or probeset)
#' in a expression matrix in two ways: i) in the whole matrix ii) to specific phenotype (case or control).
#' @param data  The whole normalize expression matrix, rows: genes or probeset, columns: samples.
#' @param complete  Boolean to define if the function use the whole expression matrix, by default TRUE.
#' @param A vector with 0 and 1 to each one of the samples in the expression matrix, the 0 express 
#' the control samples and 1 express the case samples, by default NULL.
#' @param type  Can be "case" to obtain the mean and the coefficient of variation to the case samples or, otherwise, "control" to 
#' obtain these two values to the control samples.
#' @return The expression matrix with two new columns, the first one with the averages and 
#' the another one with the coefficient of variation values.
#' @examples 
#' 
#' ## Creating the expression matrix
#' 
#' # The matrix have 200 genes and 20 samples
#' 
#' n <- 200
#' m <- 20
#' 
#' # The vector with treatment samples and control samples
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
#' # In whole expression matrix
#' 
#' complete <- cof.var(norm)
#' head(complete)
#' 
#' # In case samples
#' 
#' case <- cof.var(data = norm,complete = FALSE,treatment = t,type = "case")
#' head(case)


cof.var <- function(data,complete=TRUE,treatment=NULL,type=NULL){
  
  # Replace the sample name for case/control ID
  
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
  
  # Obtains the mean of the expression values
  
  tdata$mean <- rowMeans(tdata, na.rm = T)
  
  # Obtains the coefficient of variation to each row
  
  CV <- function(x){sd(x,na.rm = T)/mean(x,na.rm = T)}
  
  # Applies the function at the expression matrix
  
  tdata$cv <- apply(tdata[,1:(ncol(tdata)-1)],1,CV)
  
  # Returns the expression matrix with the mean and the CV
  
  return(tdata)
}