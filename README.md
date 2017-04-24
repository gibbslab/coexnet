# coexnet: An R package to build CO-EXpression NETworks from Microarray Data

## Overview

Extracts the gene expression matrix from GEO DataSets as a
    AffyBatch object. Additionally, can make the normalization process using two
    different methods (vsn and rma). The summarization (pass from multi-probe to
    one gene) uses two different criteria (Maximum value and Median of the samples
    expression data) and the process of gene differentially expressed analisys using
    two methods (sam and acde). The construction of the co-expression network can
    be conduced using two different methods, Pearson Correlation Coefficient (PCC)
    or Mutual Information (MI) and choosing a threshold value using a graph theory
    approach.

## Install

To install **coexnet** from GitHub you need *devtools* package:

```
# Install devtools
install.packages("devtools")

# Install coexnet
devtools::install_github("gibbslab/coexnet")
library(coexnet)

```

## Available functions

| Name | Description |
| :-------- | :------------------ | 
| CCP | Obtain the Common Connection patterns for two or more compared networks |
| cof.var | Calculate the coefficient of variation to expression matrix. |
| create.net | Create a co-expression network from expression matrix. |
| dif.exprs | Differential expression analysis using two different methods. |
| expr.mat | Calculate the expression matrix from the raw expression data. |
| find.threshold | Find the threshold value to create a co-expression network. |
| gene.symbol | Create a table relating probesets with genes. |
| get.affy | Charge and create an AffyBatch object |
| get.info | Download raw expression data from GEO DataSet |
| ppi.net | Create a protein-protein interaction network |
| shared.components | Obtain the shared components for two or more compared networks |

## Citation

Juan Henao, Liliana Lopez-Kleine Andres Pinzon-Velasco (2016). coexnet: An R package to build CO-EXpression NETworks from Microarray Data (Version 0.1) [software] Available at https://github.com/gibbslab/coexnet
