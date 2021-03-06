% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CCP.R
\name{CCP}
\alias{CCP}
\title{Finding Common Connection Pattern between different networks.}
\usage{
CCP(..., by = NULL)
}
\arguments{
\item{...}{The networks (igraph objects) to obtain the Common Connection Patterns.}

\item{by}{It can be, a degree value or a range of values (c(min, max)), to defined the Common Connection Pattern. If you pass one value all the nodes above this degree will be used. By default is NULL, it calculates the CCPs using all the network.}
}
\value{
An igraph object with all the Common Connection Pattern in the same network.
}
\description{
From the intersection of two or more networks, it obtains the connected components 
by deleting the solitary nodes in the intersection network. They must be igraph objects.
}
\details{
The Common Connection Pattern (CCP), is a new methodological proposal to identify molecular components linked together and common in several biological networks. The principal assumption behind Common Connection Pattern is that the networks to be compared must have the same molecular information from, i.e., must compare one layer of molecular abstraction at the same time, for example, co-expression layer, protein-protein layer, the gene regulation layer, among others.

For this, the intersection between biological networks is calculated whose result are the sub-networks with diameter greater than zero, being each of them considered as a Common Connection Pattern.
}
\examples{
# Loading data

data("net1")
data("net2")

# Obtaining Common Connection Patterns

ccp <- CCP(net1,net2)
ccp
}
\author{
Juan David Henao <judhenaosa@unal.edu.co>
}
