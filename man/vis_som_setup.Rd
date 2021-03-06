% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis_som.R
\name{vis_som_setup}
\alias{vis_som_setup}
\title{Setup the SOM lattice for visualizations}
\usage{
vis_som_setup(
  SOM,
  mar = 0.1,
  lattice_coords = F,
  coords.cex = 0.75,
  active = T,
  subset = NULL,
  change.par = TRUE
)
}
\arguments{
\item{SOM}{a SOM object}

\item{mar}{Optional, the plot margin around the lattice (will be recycled for all sides). 
Default = 0.1}

\item{lattice_coords}{whether to add the lattice (i,j) coordinates to the left & bottom of the plot}

\item{coords.cex}{size of lattice coordinate text, if requested}

\item{active}{Optional, if the SOM object has been recalled,
restricts plotting only to active neurons (those whose RF_size > 0).
Default = TRUE.}

\item{subset}{Optional, a vector of neuron indices to restrict the plotting to. 
Default = NULL imposes no restriction (plots whole lattice)}

\item{change.par}{whether to allow the vis function to optimally change the \code{par} of the plot. Default = TRUE. 
If \code{par} is allowed to changed internally, it is always reset upon function exit.}
}
\description{
Setup the SOM lattice for visualizations
}
