% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/somobj_doc.R
\name{calc_eta}
\alias{calc_eta}
\title{Calculate the SOM neighborhood function}
\usage{
SOMobj$calc_eta(sigma)
}
\arguments{
\item{sigma}{the maximum neighborhood radius for which the neighborhood functional is applied}
}
\value{
a vector of the eta neighborhood coefficients for a given value of \code{sigma}. 
The vector is ordered such that \code{eta[1]} is the coefficient applied to the prototype update of the BMU, 
\code{eta[2]} is the coefficient applied to prototype updates within a radius=1 of the BMU, 
\code{eta[3]} is the coefficient applied to prototype updates within a radius=2 of the BMU, etc.
}
\description{
Topology preservation in SOM mappings is enforced by the neighborhood function, 
which, during each training step, specifies cooperative updates of a small radius of prototype's whose neurons neighbor the BMU. 
The maximum neighborhood radius at each time is set by the \code{sigma} parameter specified in \link{set_LRAS}.  
The eta coefficients are computed from a logistic decay up to this \code{sigma}, so that during each training step 
the winning prototype is updated with the strongest effect, the prototypes within a radius = 1 on the lattice are updated with 
the next strongest effect, and so on (up to a maximum radius of \code{sigma}).  

It is usually not necessary to directly call this function as it is invoked inside \link{train_SOM}.
}
