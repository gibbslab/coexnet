## test for corVar function ##

test_cofVar <- function(){
  
  ## Correct cases
  
  treat <- c(0,1,0,1)
  
  a <- seq(0.1,0.5,length.out = 20)
  norm <- data.frame(a,a+0.2,a+0.1,a+0.3)
  
  # complete
  
  result <- cofVar(norm)
  checkEqualsNumeric(round(result$cv[1],7),0.5163978)
  
  # case
  
  result2 <- cofVar(norm,complete = FALSE,treatment = treat,type = "case")
  checkEqualsNumeric(round(result2$cv[1],7),0.2020305)
  
  # control
  
  result3 <- cofVar(norm,complete = FALSE,treatment = treat,type = "control")
  checkEqualsNumeric(round(result3$cv[1],7),0.4714045)
  
  ## Errors
  
  checkException(cofVar(),silent = TRUE)
  checkException(cofVar(norm[1,3]),silent = TRUE)
  checkException(cofVar(norm,complete = FALSE,treatment = treat,type = "tratment"),silent = TRUE)
  checkException(cofVar(norm,complete = FALSE,treatment = c(0,1),type = "control"),silent = TRUE)
}

## test for createNet function ##

test_createNet <- function(){
  
  ## Correct cases
  
  pathfile <- system.file("extdata","expression_example.txt",package = "coexnet")
  data <- read.table(pathfile,stringsAsFactors = FALSE)
  
  # correlation
  
  cor_pearson <- createNet(expData = data,threshold = 0.7,method = "correlation")

  checkEquals(names(cor_pearson["RFC2"][4]),"PTPN21")
  checkEquals(names(cor_pearson["RFC2"][137]),"IL28A")
  checkEqualsNumeric(cor_pearson["RFC2"][137],1)
  
  # mutual information
  
  mut_inf <- createNet(expData = data,threshold = 0.5,method = "mutual information")
  
  checkEquals(names(mut_inf["PXK"][40]),"NEDD1")
  checkEqualsNumeric(mut_inf["PXK"][40],0)
  checkEqualsNumeric(table(mut_inf["RFC2"] == 1)[2],3)
  
  ## Errors
  
  checkException(createNet(),silent = TRUE)
  checkException(createNet(expData = data,threshold = 0.3),silent = TRUE)
  checkException(createNet(expData = data[1,5],method = "mutual information",threshold = 0.5),silent = TRUE)
}

## test for findThreshold ##

test_findThreshold <- function(){
  
  ## Correct cases
  
  pathfile <- system.file("extdata","expression_example.txt",package = "coexnet")
  data <- read.table(pathfile,stringsAsFactors = FALSE)
  
  # correlation
  
  checkEqualsNumeric(findThreshold(expData = data[1:5,],method = "correlation"),0.01)
  checkEqualsNumeric(findThreshold(expData = data[1:100,],method = "mutual information"),0.03)
  
  ## Errors
  
  checkException(findThreshold(),silent = TRUE)
  checkException(findThreshold(expData = data),silent = TRUE)
  checkException(findThreshold(expData = data[1,5],method = "mutual information"),silent = TRUE)
}

## test for difExprs function ##

test_difExprs <- function(){
  
  ## Correct cases
  
  treat <- c(rep(0,10),rep(1,10))
  
  norm <- read.table(system.file("extdata","expression_example.txt",package = "coexnet"))
  
  # sam
  
  sam <- difExprs(expData = norm,treatment = treat,fdr = 0.05,DifferentialMethod = "sam")
  
  checkTrue(dim(sam)[1] >= 1)
  checkTrue(is.data.frame(sam))
  checkEqualsNumeric(nrow(sam),24)
  
  # acde
  
  acde <- difExprs(expData = norm,treatment = treat,fdr = 0.05,DifferentialMethod = "acde")
  
  checkTrue(dim(acde)[1] >= 1)
  checkTrue(is.data.frame(acde))
  checkEqualsNumeric(nrow(acde),14)
  
  ## Errors
  
  checkException(difExprs(),silent = TRUE)
  checkException(difExprs(expData = norm,treatment = treat,fdr = 3,DifferentialMethod = "sam"),silent = TRUE)
  checkException(difExprs(expData = norm[1:10,],treatment = treat,fdr = 0.05,DifferentialMethod = "sam"),silent = TRUE)
  checkException(difExprs(expData = norm,treatment = c(0,1),fdr = 0.05,DifferentialMethod = "sam"),silent = TRUE)
}

## test for getInfo function ##

test_getInfo <- function(){
  
  ## Errors
  
  checkException(getInfo(),silent = TRUE)
  checkException(getInfo(GPL = "GPL2025",directory = tempdir()),silent = TRUE)
}

## test geneSymbol function ##

test_geneSymbol <- function(){
  
  gene_table <- geneSymbol(GPL = "GPL2025",directory = system.file("extdata",package = "coexnet"))
  gene_na <- na.omit(gene_table)
  final_table <- gene_na[gene_na$ID != "",]
  
  ## Correct cases
  
  checkEquals(final_table$ID[1],"Os03g0669200")
  checkEqualsNumeric(dim(final_table)[1],6)
  checkEquals(rownames(final_table)[1],"118")
  checkTrue(is.factor(final_table$probe))
  
  ## Errors
  
  checkException(geneSymbol(),silent = TRUE)
  checkException(geneSymbol(directory = system.file("extdata",package = "coexnet")),silent = TRUE)
  checkException(geneSymbol(GPL = "GPL2025.soft",directory = system.file("extdata",package = "coexnet")),silent = TRUE)
}

## test for getAffy function

test_getAffy <- function(){
  
  affy <- getAffy(GSE = "GSE1234",directory = system.file("extdata",package = "coexnet"))
  
  ## Correct cases
  
  checkTrue(is.object(affy))
  checkEqualsNumeric(affy@nrow,500)
  checkEqualsNumeric(affy@ncol,500)
  checkEqualsNumeric(round(affy@assayData$exprs[1],digits = 1),77.3)
  checkEqualsNumeric(length(affy@assayData$exprs),250000)
  checkEquals(affy@protocolData@data$ScanDate,"06/30/04 10:20:15")
  
  ## Errors
  
  checkException(getAffy(),silent = TRUE)
  checkException(getAffy(directory = system.file("extdata",package = "coexnet")),silent = TRUE)
  checkException(getAffy(GSE = "1234",directory = system.file("extdata",package = "coexnet")),silent = TRUE)
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

## test for sharedComponents function

test_sharedComponents <- function(){
  
  data("net1")
  data("net2")
  
  share <- sharedComponents(net1,net2)
  
  # Correct cases
  
  checkTrue(is.vector(share))
  checkTrue(is.character(share))
  checkEqualsNumeric(length(share),4)
  checkEquals(share[1],"P")
  
  # Errors
  
  checkException(sharedComponents(),silent = TRUE)
  checkException(sharedComponents(c(net1,net2)),silent = TRUE)
  checkException(sharedComponents(c("P","M","N")),silent = TRUE)
}

## test for ppiNet function

test_ppiNet <- function(){
  
  ppi <- ppiNet(file = system.file("extdata","ppi.txt",package = "coexnet"))
  
  # Correct cases
  
  checkTrue(is.object(ppi))
  checkTrue(is.list(ppi))
  checkEqualsNumeric(length(ppi),10)
  
  # Errors
  
  checkException(ppiNet(),silent = TRUE)
  checkException(ppiNet(file = c("SNCA","UBC"),species_ID = 0),silent = TRUE)
}