#ifndef SOM_UTILS_SCALING_HPP
#include "som_utils_scaling.hpp"
#endif

//' Linear scaling of the rows of a data matrix
//' 
//' Each row is mapped linearly from a given range, to a given range
//' 
//' @param X a data matrix with vectors in rows
//' @param from_min a vector defining the min from which X is scaled, by dimension
//' @param from_max a vector defining the max from which X is scaled, by dimension
//' @param to_min a vector defining the min to which X is scaled, by dimension
//' @param to_max a vector defining the max to which X is scaled, by dimension
//' @param parallel boolean, whether to compute in parallel
//' 
//' @return a matrix (same dimensions as X) whose rows are scaled to the requested range
//' 
//' @export
// [[Rcpp::export]]
arma::mat linscale(const arma::mat& X, 
                   const arma::rowvec& from_min, const arma::rowvec& from_max, const arma::rowvec& to_min, const arma::rowvec& to_max, 
                   bool parallel = true) {
  
  arma::mat out = SOM_UTILS_SCALING::linscale(X, from_min, from_max, to_min, to_max, parallel);
  return out; 
}

