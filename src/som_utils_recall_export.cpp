#ifndef SOM_UTILS_RECALL_HPP
#include "som_utils_recall.hpp"
#endif

//' Find BMU of data vectors 
//' 
//' @param X a data matrix, vectors in rows 
//' @param W a prototype matrix, vectors in rows
//' @param nBMU optional, the number of BMUs returned. Default = 2. 
//' @param parallel optional, whether to compute in parallel. Default = TRUE
//' @param bias optional, a vector of bias to be applied to the distance calculation when computing BMU. 
//' Default = NULL, meaning no bias is applied. If given, length(bias) must equal nrow(W) (one bias for every prototype).
//' @param Xmin optional, a vector giving the minimum of the range of X used in any linearly scaling before distances are computed. 
//' Default = NULL, meaning no scaling is applied. 
//' @param Xmax optional, a vector giving the maximum of the range of X used in any linearly scaling before distances are computed. 
//' Default = NULL, meaning no scaling is applied. 
//' @param Wmin optional, a vector giving the minimum of the range of W used in any linearly scaling before distances are computed. 
//' Default = NULL, meaning no scaling is applied. 
//' @param Wmax optional, a vector giving the maximum of the range of W used in any linearly scaling before distances are computed. 
//' Default = NULL, meaning no scaling is applied. 
//' 
//' @details If given, bias affects the distance calculation between the rows of X and W by: dist(X[i,], W[j,]) - bias[j]. 
//' 
//' The inputs Xmin, Xmax, Wmin, Wmax control any linearly scaling that is applied prior to distance calculation. 
//' If any of the above are given, they all must be given.  Each can be given as either a vector controlling ranges by dimension (length = ncol(W)), 
//' or a length 1 vector (in which case the single value will be recycled across dimensions).  Specifying the X and W ranges is only useful 
//' if you are trying to match a BMU selection that was performed during SOM training, which utilized an external / internal network range.  
//' 
//' @return a list with components 
//' \itemize{ 
//'   \item BMU a matrix (nrow = nrow(X), ncol = nBMU) whose (i,j) entry lists the index (the row of W, 1-based) of the j-th BMU of the i-th row of X. 
//'   \item SQE a vector (length = nrow(X)) whose i-th element gives the Squared Quantization Error of the i-th row of X
//' }
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List find_BMU(const arma::mat& X, const arma::mat& W, unsigned int nBMU = 2, bool parallel = true, 
                    Rcpp::Nullable<Rcpp::NumericVector> bias = R_NilValue, 
                    Rcpp::Nullable<Rcpp::NumericVector> Xmin = R_NilValue, Rcpp::Nullable<Rcpp::NumericVector> Xmax = R_NilValue,
                    Rcpp::Nullable<Rcpp::NumericVector> Wmin = R_NilValue, Rcpp::Nullable<Rcpp::NumericVector> Wmax = R_NilValue) {
  
  // Initialize worker 
  Rcpp::Rcout << "Initializing BMU worker" << std::endl; 
  SOM_UTILS_RECALL::find_BMU_worker wkr(X, W, nBMU); 
  
  // Check for bias. If given, set it inside the worker. 
  if(bias.isNotNull()) {
    arma::vec tmpbias = Rcpp::as<arma::vec>(bias);
    wkr.set_bias(tmpbias); 
    Rcpp::Rcout << "++ bias set" << std::endl; 
  }
  
  // Check for ranges. If given, set them in the worker. 
  Rcpp::LogicalVector use_ranges(4); 
  use_ranges[0] = Xmin.isNotNull(); use_ranges[1] = Xmax.isNotNull(); use_ranges[2] = Wmin.isNotNull(); use_ranges[3] = Wmax.isNotNull(); 
  if(Rcpp::any(use_ranges==true).is_true()) {
    if(Rcpp::all(use_ranges==true).is_true()) {
      arma::rowvec tmpXmin = Rcpp::as<arma::rowvec>(Xmin); 
      arma::rowvec tmpXmax = Rcpp::as<arma::rowvec>(Xmax); 
      arma::rowvec tmpWmin = Rcpp::as<arma::rowvec>(Wmin); 
      arma::rowvec tmpWmax = Rcpp::as<arma::rowvec>(Wmax); 
      wkr.set_ranges(tmpXmin, tmpXmax, tmpWmin, tmpWmax); 
      Rcpp::Rcout << "++ scaling ranges set" << std::endl;   
    } else {
      Rcpp::Rcout << "++ scaling ranges not set, all must be specified" << std::endl;   
    }
  } 
  
  // Compute and increment BMU indices by 1, since it is being returned to R 
  Rcpp::Rcout << "Finding BMUs ..."; 
  if(parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
  wkr.BMU += 1;   
  Rcpp::Rcout << "done" << std::endl; 

  Rcpp::List out; 
  out["BMU"] = wkr.BMU; 
  out["SQE"] = wkr.SQE; 
  return out; 
}




//' Build a CADJ matrix
//' 
//' @param BMU matrix (nrow = nrow(X), ncol = nBMU) of indices (1-based) of the BMUs of each data vector in X. 
//' @param nW the number of SOM prototypes (CADJ matrices have nrow = ncol = nW)
//' @param parallel optional, whether to compute in parallel. Default = TRUE
//' 
//' @details The input BMU matrix should be the result of calling \code{find_BMU} with \code{nBMU} >= 2.  
//' 
//' @return a CADJ matrix 
//' 
//' @export
// [[Rcpp::export]]
arma::umat build_CADJ(const arma::umat& BMU, unsigned int nW, bool parallel = true) {
  
  if(!(BMU.n_cols >= 2)) Rcpp::stop("Input BMU must have at least two columns");
  
  SOM_UTILS_RECALL::build_CADJ_worker wkr(BMU, nW, false); // last argument is zeroidx=false
  
  if(parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
  
  return wkr.CADJ; 
}

//' Produce a list of data vectors in each prototype's Receptive Field 
//' 
//' @param BMU matrix (nrow = nrow(X), ncol = nBMU) of indices (1-based) of the BMUs of each data vector in X. 
//' @param nW the number of SOM prototypes (CADJ matrices have nrow = ncol = nW)
//' @param parallel optional, whether to compute in parallel. Default = TRUE
//' 
//' @details The input BMU matrix should be the result of calling \code{find_BMU} with \code{nBMU} >= 2.  
//' 
//' @return a list (length = nW), containing the (1-based) indices of the data vectors mapped to each prototype. 
//' 
//' @export
// [[Rcpp::export]]
std::vector<arma::uvec> find_RF_members(const arma::umat& BMU, unsigned int nW, bool parallel = true) {
  
  SOM_UTILS_RECALL::find_RF_members_worker wkr(BMU, nW, false); 
  
  if(parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
  
  return wkr.RF_members; 
}


//' Compute the plurality label of each prototype. 
//' 
//' @description For labeled data mapped to each prototype's Receptive Field, the plurality label
//' is returned. 
//' 
//' @param X_label a vector (length = nrow(X)) of data labels
//' @param RF_members a list containing the (1-based) indices of the data vectors mapped to each prototype, 
//' as returned by \code{find_RF_members}. 
//' @param parallel optional, whether to compute in parallel. Default = TRUE
//' 
//' @return a vector of prototype labels. If a prototype's Receptive Field is empty the returned label is NA
//' 
//' @export
// [[Rcpp::export]]
Rcpp::CharacterVector find_RF_label(const Rcpp::CharacterVector& X_label, const std::vector<arma::uvec>& RF_members, bool parallel = true) {
  
  SOM_UTILS_RECALL::find_RF_label_worker wkr(X_label, RF_members, false); 
  
  if(parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
  
  return wkr.RF_label; 
}


//' Receptive Field Summary Statistics
//' 
//' @description Compute summary statistics of each Receptive Field of a vector quantizer, given its forward mapping. 
//' A forward mapping is a vector containing the prototype index to which each data vector is mapped. 
//' 
//' @param X_vals the (possibly multivariate) data values to be summarized, as a matrix (nrows = nrow(X), ncol = num. dimensions)
//' @param fwdmap a vector (length = nrow(X)) of indices (1-based) describing the forward mapping. 
//' The i-th element of FWDmap should contain the index of the prototype that the i-th data vector is mapped to. 
//' @param nW the total number of prototypes in the vector quantizer 
//' @param stat the statistic to compute. Possible values are 
//' \itemize{
//' \item count \item sum \item mean \item sd \item q0 (min) \item q25 (first quantile) \item q50 (median) \item q75 (third quantile) \item q100 (max)
//' }
//' @param parallel optional, whether to compute in parallel. Default = TRUE
//' 
//' @details The statistics are computed individually for each data dimension. 
//' The statistics for any prototypes in the vector quantizer whose Receptive Fields are empty (no data mapped to them) are returned as NA. 
//' The rows of the returned matrix of statistics represent prototype-level summary statistics (row1 contains stats of prototype1, etc.).
//' 
//' @return a matrix (nrows = nW, ncol = num. data dimensions) containing the requested statistic for each prototype
//' 
//' @export
// [[Rcpp::export]]
arma::mat summarystat_RF_fwdmap(const arma::mat& X_vals, const arma::uvec& fwdmap, unsigned int nW, std::string stat = "mean", bool parallel = true) {
  
  SOM_UTILS_RECALL::groupstat_fwdmap_worker wkr(X_vals, fwdmap, nW, false); // last argumenent is zeroidx = false
  wkr.set_stat(stat);
  
  if(parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
  
  return wkr.stats;
}


//' Receptive Field Summary Statistics
//' 
//' @description Compute summary statistics of each Receptive Field of a vector quantizer, given its reverse mapping. 
//' A reverse mapping is a list containing the data indices in the Receptive Field of each prototype. 
//' 
//' @param X_vals the (possibly multivariate) data values to be summarized, as a matrix (nrows = nrow(X), ncol = num. dimensions)
//' @param revmap a list (length = nW) of indices (1-based) describing the reverse mapping. 
//' The i-th element of REVmap should contain the indices of data vectors that are mapped to the the i-th prototype of the quantizer. 
//' @param stat the statistic to compute. Possible values are 
//' \itemize{
//' \item count \item sum \item mean \item sd \item q0 (min) \item q25 (first quantile) \item q50 (median) \item q75 (third quantile) \item q100 (max)
//' }
//' @param parallel optional, whether to compute in parallel. Default = TRUE
//' 
//' @details The statistics are computed individually for each data dimension. 
//' The statistics for any prototypes in the vector quantizer whose Receptive Fields are empty (no data mapped to them) are returned as NA. 
//' The rows of the returned matrix of statistics represent prototype-level summary statistics (row1 contains stats of prototype1, etc.).
//' 
//' @return a matrix (nrows = nW, ncol = num. data dimensions) containing the requested statistic for each prototype
//' 
//' @export
// [[Rcpp::export]]
arma::mat summarystat_RF_revmap(const arma::mat& X_vals, const std::vector<arma::uvec>& revmap, std::string stat = "mean", bool parallel = true) {
  
  SOM_UTILS_RECALL::groupstat_revmap_worker wkr(X_vals, revmap, false); // last argumenent is zeroidx = false
  wkr.set_stat(stat);
  
  if(parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
  
  return wkr.stats;
}


  
//' Calculate the fences of the SOM lattice 
//' 
//' @description Lattice fences between neighboring neurons on the lattice are the (squared) Euclidean distance between 
//' the neurons' corresponding prototypes. 
//' 
//' @param nu_xy a matrix whose rows contain the (x,y) coordinates (cols 1 and 2, respectively) of each lattice neuron, ordered 
//' from the bottom-left of the lattice. 
//' @param W the prototype matrix (rows ordered the same as nu_xy)
//' @param nu_verts a 3-d array whose slices contain the vertices of each lattice tile (slices ordered the same as nu_xy)
//' @param nu_ADJ the neuron adjacency matrix 
//' @param parallel optional, whether to compute in parallel. Default = TRUE
//' 
//' @return A data frame identifying the fence and its corresponding value between each pair of neighboring lattice neurons, with columns: 
//' \itemize{
//' \item i the index of the 1st neuron comprising the fence 
//' \item j the index of the 2nd neuron comprising the fence
//' \item x0 the lattice x-coord of neuron i 
//' \item y0 the lattice y-coord of neuron i 
//' \item x1 the lattice x-coord of neuron j 
//' \item y1 the lattice y-coord of neuron j
//' \item value the fence value separating neurons i and j 
//' }
//' 
//' @export 
// [[Rcpp::export]]
Rcpp::DataFrame calc_SOM_fences(const arma::mat& nu_xy, const arma::mat& W, const arma::cube& nu_verts, const arma::umat& nu_ADJ, bool parallel = true) {
  SOM_UTILS_RECALL::fence_worker wkr(nu_xy, W, nu_verts, nu_ADJ); 
  if(parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
  return wkr.fence_dataframe();
}
