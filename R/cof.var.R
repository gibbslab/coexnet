#'@export cof.var
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Calculate the mean and coefficient of variation of the expression matrix to cases, controls or all the samples
#' @description Using a expression matrix, calculates the mean and the coefficient of variation to the probeset or genes using
#' only the cases samples or control samples or using all the samples.
#' @usage cof.var(data, type = c("control", "case", "complete"), treatment)
#' @param data A data frame with the expression matrix.
#' @param type The samples to be used, can be "control", "case" or "complete".
#' @param treatment A vector of 0 and 1 values, 0 are the control samples and 1 are the cases samples.
#' @return A data frame with the expression matrix and two new colums. one to the mean and the another one to the 
#' coefficient of variation values.

cof.var <- function(data,type,treatment){
  
  # Changes the names of samples to "0" or "1"
  
  names(data) <- treatment
  
  # Replaces the value of "type" variable
  
  if(type == "control"){
    
    type <- "0"
  }else if(type == "case"){
    
    type <- "1"
  }
  
  if(type == "complete"){
    
    # Uses all the samples
    
    tdata = data
    
  }else{
    
    # Uses the control samples or the case samples
    
    tdata <- data[names(data) == type]
  
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