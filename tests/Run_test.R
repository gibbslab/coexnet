if (require("BioGenerics")) {
  BiocGenerics:::testPackage("coexnet")
} else {
  stop("BiocGenerics package not installed")
} 

