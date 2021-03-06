% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{summarystat_RF_revmap}
\alias{summarystat_RF_revmap}
\title{Receptive Field Summary Statistics}
\usage{
summarystat_RF_revmap(X_vals, revmap, stat = "mean", parallel = TRUE)
}
\arguments{
\item{X_vals}{the (possibly multivariate) data values to be summarized, as a matrix (nrows = nrow(X), ncol = num. dimensions)}

\item{revmap}{a list (length = nW) of indices (1-based) describing the reverse mapping. 
The i-th element of REVmap should contain the indices of data vectors that are mapped to the the i-th prototype of the quantizer.}

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
Compute summary statistics of each Receptive Field of a vector quantizer, given its reverse mapping. 
A reverse mapping is a list containing the data indices in the Receptive Field of each prototype.
}
\details{
The statistics are computed individually for each data dimension. 
The statistics for any prototypes in the vector quantizer whose Receptive Fields are empty (no data mapped to them) are returned as NA. 
The rows of the returned matrix of statistics represent prototype-level summary statistics (row1 contains stats of prototype1, etc.).
}
