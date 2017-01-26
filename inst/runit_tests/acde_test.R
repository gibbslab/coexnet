## test for cor.var function ##

test_cof.var <- function(){
  
  ## correct case
  
  t <- c(rep(0,10),rep(1,10))
  
  mat <- as.matrix(rexp(200, rate = 1))
  norm <- t(apply(mat, 1, function(nm) rnorm(20, mean=nm, sd=1)))
  
  # complete
  
  result <- cof.var(norm)
  checkEqualsNumeric(round(result[1,1],6),2.675794)
  
  
}