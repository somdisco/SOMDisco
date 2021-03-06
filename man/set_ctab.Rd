% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/somobj_doc.R
\name{set_ctab}
\alias{set_ctab}
\title{Set the color table}
\usage{
SOMobj$set_ctab(ctab)
}
\arguments{
\item{ctab}{a color table defining the mapping between the unique labels found in \code{X_label} and distinct colors. 
Must be a data frame with columns 'label' and 'color' (where 'color' is in HEX format).}
}
\value{
None, field \code{ctab} is set internally
}
\description{
The color table controls the mapping of the unique labels found in \code{X_label} to distinct colors.
}
\details{
\code{ctab} must contain color mappings for all unique labels found in \code{X_label}. A check will 
be performed internally.
}
