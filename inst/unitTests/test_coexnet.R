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

## test for create.net function ##

test_create.net <- function(){
  
  ## Correct cases
  
  pathfile <- system.file("extdata","expression_example.txt",package = "coexnet")
  data <- read.table(pathfile,stringsAsFactors = FALSE)
  
  # correlation
  
  cor_pearson <- create.net(difexp = data,threshold = 0.7,method = "correlation")

  checkEquals(names(cor_pearson["RFC2"][4]),"PTPN21")
  checkEquals(names(cor_pearson["RFC2"][137]),"IL28A")
  checkEqualsNumeric(cor_pearson["RFC2"][137],1)
  
  # mutual information
  
  mut_inf <- create.net(difexp = data,threshold = 0.5,method = "mutual information")
  
  checkEquals(names(mut_inf["PXK"][40]),"NEDD1")
  checkEqualsNumeric(mut_inf["PXK"][40],0)
  checkEqualsNumeric(table(mut_inf["RFC2"] == 1)[2],3)
  
  ## Errors
  
  checkException(create.net(),silent = TRUE)
  checkException(create.net(difexp = data,threshold = 0.3),silent = TRUE)
  checkException(create.net(difexp = data[1,5],method = "mutual information",threshold = 0.5),silent = TRUE)
}