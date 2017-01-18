#' @export cof.var
#' @author Juan David Henao <judhenaosa@unal.edu.co>

# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

#' @title Calculate the coefficient of variation to expression matrix 
#' @description 


cof.var <- function(data,complete=TRUE,treatment=NULL,type=NULL){
  
  # Replace the sample name for case/control ID
  
  if (complete == FALSE) {
    
    names(data) <- treatment
    
    if(type == "control"){
      
      type <- "0"
    }else if(type == "case"){
      
      type <- "1"
    }
    tdata <- data[names(data) == type]
  }else{
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