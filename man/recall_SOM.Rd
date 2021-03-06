% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/somobj_doc.R
\name{recall_SOM}
\alias{recall_SOM}
\title{Recall a trained SOM object}
\usage{
SOMobj$recall_SOM(X)
}
\arguments{
\item{X}{the training data matrix. 
This should be the same matrix that was input to \link{initialize_SOM} 
(checks will be performed on its statistics, if they do not match an error will be returned).}
}
\value{
None, the fields described above are calculated and stored
}
\description{
A SOM recall maps all training data to their representative prototypes (their BMUs). 
Several quantities, as outlined below, result from this mapping.
}
\details{
Several internal fields are set by the recall function: 
\itemize{
\item \code{BMU} matrix (nrow = \code{nX}, ncol = \code{nBMU}) containing the 1st, 2nd, ... BMUs for every training vector in its columns. 
The valued in \code{BMU} are the 1-based indices of the prototypes, in neuron order. 
\item \code{SQE} vector (length = \code{nX}) containing the squared quantization error of each training vector, as quantized by its 1st BMU
\item \code{CADJ} the CADJ matrix (nrow = ncol = \code{nW}) of Cumulative (weighted) Topological Adjacencies of the SOM prototypes. 
See \link{CONN} for details. 
\item \code{RF_size} vector (length = \code{nW}) containing the number of training vectors mapped to each prototype
\item \code{Entropy} the (normalized) entropy of the discrete SOM mapping
\item \code{RF_members} a list (length = \code{nW}) of the training vector indices mapped to each prototype 
\item \code{fences} a data frame containing the information for visualizing the U-matrix fences on the SOM lattice, 
see \link{set_lattice_fences} for details. 
\item \code{RF_label} and \code{RF_label_dist}, the values of \code{X_label} propagated to the Receptive Fields, 
according to the SOM mapping. 
}

Note: the "Receptive Field" of a prototype \eqn{w_j} is the set of all training data for which \eqn{w_j} is the BMU. 

It is usually not necessary to directly call this function as it is invoked inside \link{train_SOM}.
}
