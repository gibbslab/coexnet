## test for cor.var function ##

test_cof.var <- function(){
  
  ## correct case
  
  t <- c(rep(0,10),rep(1,10))
  
  mat <- as.matrix(rexp(200, rate = 1))
  norm <- t(apply(mat, 1, function(nm) rnorm(20, mean=nm, sd=1)))
  
  # complete
  
  result <- cof.var(norm)
  checkEqualsNumeric(round(result[1,1],6),2.675794)
  
  # case
  
  result2 <- cof.var(norm,complete = FALSE,treatment = t,type = "case")
  checkEqualsNumeric(round(result2[1,1],6),2.137933)
  
  # control
  
  result3 <- cof.var(norm,complete = FALSE,treatment = t,type = "control")
  checkEqualsNumeric(round(result3[1,1],6),2.675794)
  
  ## error
  
  checkException(cof.var(),silent = TRUE)
  checkException(cof.var(norm[1,3]),silent = TRUE)
  checkException(cof.var(norm,complete = FALSE,treatment = t,type = "tratment"),silent = TRUE)
  checkException(cof.var(norm,complete = FALSE,treatment = c(0,1),type = "control"),silent = TRUE)
}