# coexnet: An R package to build CO-EXpression NETworks from Microarray Data

## Overview

**coexnet** run under R version 3.3.

## Install

To install **coexnet** from GitHub you needs *devtools* package:

```
# Install devtools
install.packages("devtools")

# Install coexnet
devtools::install_github("gibbslab/coexnet")
library(coexnet)

```

## Available functions

|||
| :-------- | :------------------ | 
| cof.var | Calculate the coefficient of variation to expression matrix. |
| create.net | Create a co-expression network from expression matrix. |
| dif.exprs | Differential expression analysis using two different methods. |
| expr.mat | Calculate the expression matrix from the raw expression data. |
| find.threshold | Find the threshold value to create a co-expression network. |
| gene.symbol | Create a table relating probesets with genes. |
| get.affy | Charge and create an AffyBatch object |
| get.info | Download raw expression data from GEO DataSet |

## Citation

