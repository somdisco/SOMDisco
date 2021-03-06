% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/somobj_doc.R
\name{calc_Entropy}
\alias{calc_Entropy}
\title{Calculate the Normalized Entropy of the SOM mapping}
\usage{
SOMobj$calc_Entropy()
}
\value{
None, the field \code{Entropy} is set internally.
}
\description{
The normalized entropy of the SOM quantization is given by 
\deqn{entropy = -sum(F*log(F))/log(nW)}
where \code{F} is a vector of Receptive Field relative frequencies (i.e., \eqn{RF_size / nX}). 

SOMs trained with the CSOM update rule seek to maximize the entropy of the learned mapping, so 
this normalized entropy measure provides a way of comparing different SOMs regardless of their size 
(or the number of training vectors available).

It is usually not necessary to directly call this function as it is invoked inside \link{recall_SOM}.
}
