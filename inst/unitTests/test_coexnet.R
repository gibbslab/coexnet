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
  
  checkEqualsNumeric(find.threshold(difexp = data[1:5,],method = "correlation"),0.01)
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

## test for get.affy function

test_get.affy <- function(){
  
  affy <- get.affy(GSE = "GSE1234",dir = system.file("extdata",package = "coexnet"))
  
  ## Correct cases
  
  checkTrue(is.object(affy))
  checkEqualsNumeric(affy@nrow,500)
  checkEqualsNumeric(affy@ncol,500)
  checkEqualsNumeric(round(affy@assayData$exprs[1],digits = 1),77.3)
  checkEqualsNumeric(length(affy@assayData$exprs),250000)
  checkEquals(affy@protocolData@data$ScanDate,"06/30/04 10:20:15")
  
  ## Errors
  
  checkException(get.affy(),silent = TRUE)
  checkException(get.affy(dir = system.file("extdata",package = "coexnet")),silent = TRUE)
  checkException(get.affy(GSE = "1234",dir = system.file("extdata",package = "coexnet")),silent = TRUE)
}

## test for CCP function

test_CCP <- function(){
  
  data("net1")
  data("net2")
  
  ccp <- CCP(net1,net2)
  
  # Correct cases
  
  checkTrue(is.object(ccp))
  checkTrue(is.list(ccp))
  checkEqualsNumeric(length(ccp),10)
  
  # Errors
  
  checkException(CCP(),silent = TRUE)
  checkException(CCP(c(net1,net2)),silent = TRUE)
}

## test for shared.components function

test_shared.components <- function(){
  
  data("net1")
  data("net2")
  
  share <- shared.components(net1,net2)
  
  # Correct cases
  
  checkTrue(is.vector(share))
  checkTrue(is.character(share))
  checkEqualsNumeric(length(share),4)
  checkEquals(share[1],"P")
  
  # Errors
  
  checkException(shared.components(),silent = TRUE)
  checkException(shared.components(c(net1,net2)),silent = TRUE)
  checkException(shared.components(c("P","M","N")),silent = TRUE)
}

## test for ppi.net function

test_ppi.net <- function(){
  
  ppi <- ppi.net(input = system.file("extdata","ppi.txt",package = "coexnet"))
  
  # Correct cases
  
  checkTrue(is.object(ppi))
  checkTrue(is.list(ppi))
  checkEqualsNumeric(length(ppi),10)
  
  # Errors
  
  checkException(ppi.net(),silent = TRUE)
  checkException(ppi.net(input = c("SNCA","UBC"),species_ID = 0),silent = TRUE)
}