#' A SHGRWalk Synthetic Hyperspectral Image Cube 
#' 
#' The SHGRWalk (SHGR, for short) data suite is a collection of Synthetic Hyperspectral 128 x 128 pixel image cubes 
#' whose individual pixel "spectra" were sampled from a multi-component Gaussian Mixture Model whose component means 
#' were set via a (secondary) Gaussian Random Walk across the synthetic "spectral channels" (the data dimension, \code{d}). 
#' Individual component variances possesses a Toeplitz correlation structure with various noise levels. 
#' After sampling the pixels were labeled according to which mixture component (class) they were sampled from, 
#' and then organized into the 128 x 128 pixel image such that pixels from the same sampling component occupy 
#' contiguous blocks within the overall image.  
#' More information about the SHGRWalk data can be found here \code{need link}.  
#' 
#' @format For demonstration of SOM learning a 100-dimensional SHGR cube with 20 distinct classes (each with correlated noise) 
#' has been included in the \code{SOMDisco} package. These data and their associated metadata are stored in a list 
#' variabled named \code{SHGR} with components: 
#' \describe{
#' \item{X}{data matrix whose 16,384 rows represent the 128 x 128 pixels in the image, 
#' and 100 columns represent the 100 spectral channels at which the synthetic reflectances were measured} 
#' \item{label}{character vector (length = 16,384) containing class labels of each row of \code{X}, 
#' denoted by the letters \code{A - T}}
#' \item{ctab}{color table mapping the unique character labels in the dataset to a pre-defined representative color. 
#' The color table is stored as a data frame with columns \code{label} and \code{color}.} 
#' \item{pxl.coords}{matrix (nrows = 16,384) whose two columns give the \code{(row,col)} pixel indices 
#' of the rows of \code{X} in the image.  This information is needed to map the data matrix \code{X} 
#' back to its data cube format.}
#' \item{identifier}{string detailing the specifics of the example SHGR cube}
#' }
"SHGR"