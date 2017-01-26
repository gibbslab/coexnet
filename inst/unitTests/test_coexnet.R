## test for cor.var function ##

test_cof.var <- function(){
  
  ## Correct cases
  
  t <- c(0,1,0,1)
  
  a <- seq(0.1,0.5,length.out = 20)
  norm <- data.frame(a,a+0.2,a+0.1,a+0.3)
  
  # complete
  
  result <- cof.var(norm)
  checkEqualsNumeric(round(result$cv[1],7),0.5163978)
  
  # case
  
  result2 <- cof.var(norm,complete = FALSE,treatment = t,type = "case")
  checkEqualsNumeric(round(result2$cv[1],7),0.2020305)
  
  # control
  
  result3 <- cof.var(norm,complete = FALSE,treatment = t,type = "control")
  checkEqualsNumeric(round(result3$cv[1],7),0.4714045)
  
  ## Errors
  
  checkException(cof.var(),silent = TRUE)
  checkException(cof.var(norm[1,3]),silent = TRUE)
  checkException(cof.var(norm,complete = FALSE,treatment = t,type = "tratment"),silent = TRUE)
  checkException(cof.var(norm,complete = FALSE,treatment = c(0,1),type = "control"),silent = TRUE)
}