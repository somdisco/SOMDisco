% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{summarystat_RF_fwdmap}
\alias{summarystat_RF_fwdmap}
\title{Receptive Field Summary Statistics}
\usage{
summarystat_RF_fwdmap(X_vals, fwdmap, nW, stat = "mean", parallel = TRUE)
}
\arguments{
\item{X_vals}{the (possibly multivariate) data values to be summarized, as a matrix (nrows = nrow(X), ncol = num. dimensions)}

\item{fwdmap}{a vector (length = nrow(X)) of indices (1-based) describing the forward mapping. 
The i-th element of FWDmap should contain the index of the prototype that the i-th data vector is mapped to.}

\item{nW}{the total number of prototypes in the vector quantizer}

\item{stat}{the statistic to compute. Possible values are 
\itemize{
\item count \item sum \item mean \item sd \item q0 (min) \item q25 (first quantile) \item q50 (median) \item q75 (third quantile) \item q100 (max)
}}

\item{parallel}{optional, whether to compute in parallel. Default = TRUE}
}
\value{
a matrix (nrows = nW, ncol = num. data dimensions) containing the requested statistic for each prototype
}
\description{
Compute summary statistics of each Receptive Field of a vector quantizer, given its forward mapping. 
A forward mapping is a vector containing the prototype index to which each data vector is mapped.
}
\details{
The statistics are computed individually for each data dimension. 
The statistics for any prototypes in the vector quantizer whose Receptive Fields are empty (no data mapped to them) are returned as NA. 
The rows of the returned matrix of statistics represent prototype-level summary statistics (row1 contains stats of prototype1, etc.).
}
