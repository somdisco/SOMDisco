% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/somobj_doc.R
\name{set_monitoring_freq}
\alias{set_monitoring_freq}
\title{(De-)Active Monitoring of SOM Training}
\usage{
SOMobj$set_monitoring_freq(mtr_freq)
}
\arguments{
\item{mtr_freq}{the incremental training step count at which monitoring snapshots are taken and stored during SOM training. 
If = 0, which is the default set when calling \link{initialize_SOM}, then no monitoring is performed.}
}
\value{
None, field \code{mtr_freq} is set internally.
}
\description{
Various quantities can be computed and stored at regular training step intervals to allow the analyst 
to observe how SOM learning progresses.  This function either activates or de-activates this monitoring.
}
\details{
As monitoring proceeds, the internal field \code{mtr_age} is updated with the time (training) steps 
at which snapshots were taken. 

If activated, monitoring computes and stores the following quantities: 
\itemize{
\item \code{mtr_RMSQE} a matrix (nrow = length(\code{mtr_age}), ncol = \code{nW}+1) containing the 
Root Mean Squared Error of the quantization of the data by its BMU prototype (the Root Mean Quantization Error at the prototype level). 
The last (\code{nW} + 1) column contains the RMSQE of the quantization over all prototypes.
\item \code{mtr_QSI} a vector (length = length(\code{mtr_age})) of the Quantization Stability Index 
computed during each monitoring snapshot. The QSI at any monitoring step is the proportion of data samples 
that have switched BMUs from the previous monitoring snapshot (by convention, QSI at the initial monitoring step = NA). 
QSI can convey when the SOM "settles down", and reaches a stable quantization.  
\item \code{mtr_Entropy} a vector (length = length(\code{mtr_age})) of the normalized Entropies of the SOM mappings 
taken at each monitoring step.  
}

The monitored quantities can be visualized after a training round by calling \link{vis_som_training}.
}
