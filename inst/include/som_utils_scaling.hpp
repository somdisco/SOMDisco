#ifndef RcppArmadillo_H
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifndef RcppParallel_H
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
#endif

// [[Rcpp::plugins(cpp11)]]


// #ifndef SOM_UTILS_DIST_HPP
// #include "som_utils_dist.hpp"
// #endif

#ifndef SOM_UTILS_SCALING_HPP
#define SOM_UTILS_SCALING_HPP

namespace SOM_UTILS_SCALING {

// Parallel worker to perform linear scaling of the rows of a matrix 
// from a given range, to a given range 
struct linscale_worker : public RcppParallel::Worker {
  
  const arma::mat& X; // data matrix whose rows are to be scaled
  
  // Input/Output ranges, given as a vector
  const arma::rowvec& from_min, from_max, to_min, to_max;
  
  // Internally computed variables; 
  unsigned int N,d; 
  arma::rowvec from_rng, to_rng; 
  
  // output container
  arma::mat Xscaled;
  
  // Constructor, from & to are both vecs
  linscale_worker(const arma::mat& X, 
                  const arma::rowvec& from_min_, const arma::rowvec& from_max_, 
                  const arma::rowvec& to_min_, const arma::rowvec& to_max_) : 
    X(X), from_min(from_min_), from_max(from_max_), to_min(to_min_), to_max(to_max_) {
    
    N = X.n_rows;
    d = X.n_cols; 
    
    if(from_min.n_elem != d) Rcpp::stop("length(from_min) != ncol(X)");
    if(from_max.n_elem != d) Rcpp::stop("length(from_max) != ncol(X)");
    if(to_min.n_elem != d) Rcpp::stop("length(to_min) != ncol(X)");
    if(to_max.n_elem != d) Rcpp::stop("length(to_max) != ncol(X)");
    
    from_rng = from_max - from_min; 
    to_rng = to_max - to_min; 
    
    Xscaled.set_size(N,d);
  }
  
  // Fill up neighbors from a single neuron 
  void scale_row_j(unsigned int j) {
    
    Xscaled.row(j) = (X.row(j) - from_min) / from_rng % to_rng + to_min;
    
  }
  
  // Parallel operator 
  void operator()(std::size_t begin, std::size_t end) {
    for(unsigned int j=begin; j<end; ++j) {
      scale_row_j(j); 
    }
  }
  
  // Parallel invoker 
  void calc_parallel() {
    RcppParallel::parallelFor(0, N, *this);
  }
  
  // Serial invoker 
  void calc_serial() {
    for(unsigned int j=0; j<N; ++j) {
      scale_row_j(j); 
    }
  }
  
};

inline arma::mat linscale(const arma::mat& X, 
                   const arma::rowvec& from_min, const arma::rowvec& from_max, const arma::rowvec& to_min, const arma::rowvec& to_max, 
                   bool parallel) {
  
  unsigned int d = X.n_cols; 
  
  // Check whether each of the min/max ranges are a scalar or a vector. 
  // If a scalar, the value will be recycled across dimensions. If a vector, its length must match the data dimension. 
  arma::rowvec from_min_, from_max_, to_min_, to_max_; 
  if(from_min.n_elem == d) {
    from_min_ = from_min; 
  } else if(from_min.n_elem == 1) {
    from_min_.set_size(d); from_min_.fill(from_min[0]);
  } else {
    Rcpp::stop("length(from_min) != {1,d}");
  }
  
  if(from_max.n_elem == d) {
    from_max_ = from_max; 
  } else if(from_max.n_elem == 1) {
    from_max_.set_size(d); from_max_.fill(from_max[0]);
  } else {
    Rcpp::stop("length(from_max) != {1,d}");
  }
  
  if(to_min.n_elem == d) {
    to_min_ = to_min; 
  } else if(to_min.n_elem == 1) {
    to_min_.set_size(d); to_min_.fill(to_min[0]);
  } else {
    Rcpp::stop("length(to_min) != {1,d}");
  }
  
  if(to_max.n_elem == d) {
    to_max_ = to_max; 
  } else if(to_max.n_elem == 1) {
    to_max_.set_size(d); to_max_.fill(to_max[0]);
  } else {
    Rcpp::stop("length(to_max) != {1,d}");
  }
  
  
  linscale_worker wkr(X, from_min_, from_max_, to_min_, to_max_); 
  
  if(parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
  
  return wkr.Xscaled;
}


} // close namespace 

#endif