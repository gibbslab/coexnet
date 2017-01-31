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

## test for find.threshold ##

test_find.threshold <- function(){
  
  ## Correct cases
  
  pathfile <- system.file("extdata","expression_example.txt",package = "coexnet")
  data <- read.table(pathfile,stringsAsFactors = FALSE)
  
  # correlation
  
  checkEqualsNumeric(find.threshold(difexp = data,method = "correlation"),0.72)
  checkEqualsNumeric(find.threshold(difexp = data[1:5,],method = "correlation"),0.01)
  
  # mutual information
  
  checkEqualsNumeric(find.threshold(difexp = data,method = "mutual information"),0.03)
  checkEqualsNumeric(find.threshold(difexp = data[1:100,],method = "mutual information"),0.03)
  
  ## Errors
  
  checkException(find.threshold(),silent = TRUE)
  checkException(find.threshold(difexp = data),silent = TRUE)
  checkException(find.threshold(difexp = data[1,5],method = "mutual information"),silent = TRUE)
}

## test for dif.exprs function ##

test_dif.exprs <- function(){
  
  ## Correct cases
  
  t <- c(rep(0,10),rep(1,10))
  
  norm <- read.table(system.file("extdata","expression_example.txt",package = "coexnet"))
  
  # sam
  
  sam <- dif.exprs(eset = norm,treatment = t,fdr = 0.05,DifferentialMethod = "sam")
  
  checkTrue(dim(sam)[1] >= 1)
  checkTrue(is.data.frame(sam))
  checkEqualsNumeric(nrow(sam),24)
  
  # acde
  
  acde <- dif.exprs(eset = norm,treatment = t,fdr = 0.05,DifferentialMethod = "acde")
  
  checkTrue(dim(acde)[1] >= 1)
  checkTrue(is.data.frame(acde))
  checkEqualsNumeric(nrow(acde),14)
  
  ## Errors
  
  checkException(dif.exprs(),silent = TRUE)
  checkException(dif.exprs(eset = norm,treatment = t,fdr = 3,DifferentialMethod = "sam"),silent = TRUE)
  checkException(dif.exprs(eset = norm[1:10,],treatment = t,fdr = 0.05,DifferentialMethod = "sam"),silent = TRUE)
  checkException(dif.exprs(eset = norm,treatment = c(0,1),fdr = 0.05,DifferentialMethod = "sam"),silent = TRUE)
}

## test for get.info function ##

test_get.info <- function(){
  
  test <- get.info(GSE = "GSE8216", GPL = "GPL2025",dir = tempdir())
  dir <- dir()
  
  ## Correct cases
  
  checkEquals(test,"./GPL2025.soft")
  checkTrue(any(dir == "GSE8216"))
  
  ## Errors
  
  checkException(get.info(),silent = TRUE)
  checkException(get.info(GPL = "GPL2025",dir = tempdir()),silent = TRUE)
}

## test gene.symbol function ##

test_gene.symbol <- function(){
  
  gene_table <- gene.symbol(GPL = "GPL2025",d = system.file("extdata",package = "coexnet"))
  gene_na <- na.omit(gene_table)
  final_table <- gene_na[gene_na$ID != "",]
  
  ## Correct cases
  
  checkEquals(final_table$ID[1],"Os03g0669200")
  checkEqualsNumeric(dim(final_table)[1],6)
  checkEquals(rownames(final_table)[1],"118")
  checkTrue(is.factor(final_table$probe))
  
  ## Errors
  
  checkException(gene.smbol(),silent = TRUE)
  checkException(gene.smbol(d = system.file("extdata",package = "coexnet")),silent = TRUE)
  checkException(gene.smbol(GPL = "GPL2025.soft",d = system.file("extdata",package = "coexnet")),silent = TRUE)
  checkException(gene.smbol(GPL = "GPL2025"),silent = TRUE)
}