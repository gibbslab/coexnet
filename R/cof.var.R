# cof.var
# Juan David Henao Sanchez
# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

cof.var <- function(data,type,treatment,complete = FALSE){
  
  # Changes the names of samples to "0" or "1"
  
  names(data) <- treatment
  
  # Replaces the value of "type" variable
  
  if(type == "control"){
    
    type <- "0"
  }else if(type == "case"){
    
    type <- "1"
  }
  
  if(complete == FALSE){
    
    # Uses the control samples or the case samples
    
    tdata <- data[names(data) == type]
    
  }else{
    
    # Uses all the samples
    
    tdata = data
  
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