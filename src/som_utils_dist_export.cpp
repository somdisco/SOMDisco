#ifndef SOM_UTILS_DIST_HPP
#include "som_utils_dist.hpp"
#endif


//' Compute distances between rows of two matrices
//'
//' @param X1 a data matrix, vectors in rows 
//' @param X2 a data matrix, vectors in rows. ncol(X2) must equal ncol(X1). 
//' @param which_dist optional, the distance to compute. Options are \itemize{\item "L2" for Euclidean, \item "L22" for Squared Euclidean", \item "L1" for L-1 distance, \item "LInf" for L-Inf distance} Default = "L22".
//' @param parallel optional, whether to compute in parallel. Default = TRUE
//' @param bias optional, a vector of bias to be applied to the distance calculation when computing BMU. 
//' Default = NULL, meaning no bias is applied. If given, length(bias) must equal nrow(X2) (one bias for every row of X2).
//' @param X1min optional, a vector giving the minimum of the range of X1 used in any linearly scaling before distances are computed. 
//' Default = NULL, meaning no scaling is applied. 
//' @param X1max optional, a vector giving the maximum of the range of X1 used in any linearly scaling before distances are computed. 
//' Default = NULL, meaning no scaling is applied. 
//' @param X2min optional, a vector giving the minimum of the range of X2 used in any linearly scaling before distances are computed. 
//' Default = NULL, meaning no scaling is applied. 
//' @param X2max optional, a vector giving the maximum of the range of X2 used in any linearly scaling before distances are computed. 
//' Default = NULL, meaning no scaling is applied. 
//' 
//' @details 
//' If given, bias affects the distance calculation between the rows of X1 and X2 by: which_dist(X1[i,], X2[j,]) - bias[j]. 
//' 
//' The inputs X1min, X1max, X2min, X2max control any linearly scaling (from the range of X1, to that of X2) that is applied prior to distance calculation. 
//' If any of the above are given, they all must be given.  Each can be given as either a vector controlling ranges by dimension (length = ncol(X1)), 
//' or a length 1 vector (in which case the single value will be recycled across dimensions).  
//' 
//' Specifying scaling forces distance computation to be performed in the X2 range. This is an optional feature.
//' 
//' @return a matrix whose (i,j) entry contains the requested distance between the i-th row of X1 and and j-th row of X2
//' 
//' @export
// [[Rcpp::export]]
arma::mat distmat(const arma::mat& X1, const arma::mat& X2, std::string which_dist = "L22", bool parallel = true, 
                  Rcpp::Nullable<Rcpp::NumericVector> bias = R_NilValue, 
                  Rcpp::Nullable<Rcpp::NumericVector> X1min = R_NilValue, Rcpp::Nullable<Rcpp::NumericVector> X1max = R_NilValue,
                  Rcpp::Nullable<Rcpp::NumericVector> X2min = R_NilValue, Rcpp::Nullable<Rcpp::NumericVector> X2max = R_NilValue) {
  
  // Initialize worker 
  Rcpp::Rcout << "Initializing distance worker" << std::endl; 
  SOM_UTILS_DIST::distmat_worker wkr(X1, X2, which_dist); 
  
  // Check for bias. If given, set it inside the worker. 
  if(bias.isNotNull()) {
    arma::vec tmpbias = Rcpp::as<arma::vec>(bias);
    wkr.set_bias(tmpbias); 
    Rcpp::Rcout << "++ bias set" << std::endl; 
  }
  
  // Check for ranges. If given, set them in the worker. 
  Rcpp::LogicalVector use_ranges(4); 
  use_ranges[0] = X1min.isNotNull(); use_ranges[1] = X1max.isNotNull(); use_ranges[2] = X2min.isNotNull(); use_ranges[3] = X2max.isNotNull(); 
  if(Rcpp::any(use_ranges==true).is_true()) {
    if(Rcpp::all(use_ranges==true).is_true()) {
      arma::rowvec tmpX1min = Rcpp::as<arma::rowvec>(X1min); 
      arma::rowvec tmpX1max = Rcpp::as<arma::rowvec>(X1max); 
      arma::rowvec tmpX2min = Rcpp::as<arma::rowvec>(X2min); 
      arma::rowvec tmpX2max = Rcpp::as<arma::rowvec>(X2max); 
      wkr.set_ranges(tmpX1min, tmpX1max, tmpX2min, tmpX2max); 
      Rcpp::Rcout << "++ scaling ranges set" << std::endl;   
    } else {
      Rcpp::Rcout << "++ scaling ranges not set, all must be specified" << std::endl;   
    }
  } 
  
  // Compute and increment BMU indices by 1, since it is being returned to R 
  Rcpp::Rcout << "Computing pairwise distances ..."; 
  if(parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
  Rcpp::Rcout << "done" << std::endl; 
  
  return wkr.DIST; 
}


//' Shortest path distances between vertices of a graph
//' 
//' @description This is a wrapper for igraph::distances, see its help for more information. 
//' 
//' @param ADJ the adjacency matrix of graph vertices. Can be weighted. 
//' @param weighted optional, whether to consider the edge weights found in ADJ and construct a weighted graph. Default = FALSE. 
//' @param directed optional, whether to consider the graph as directed. Default = FALSE. 
//' 
//' @return a shortest-path distance matrix between graph vertices (nrow = ncol = num. vertices)
//' 
//' @export
// [[Rcpp::export]]
arma::umat geodesicdist(const arma::umat& ADJ, bool weighted = false, bool directed = false) {
  
  arma::umat SPDIST = SOM_UTILS_DIST::geodesicdist(ADJ, weighted, directed);
  
  return SPDIST; 
  
  // // Calculate the geodesic distances between vertices, given an adjacency matrix 
  // 
  // // Tap into igraph package, setup functions 
  // Rcpp::Environment igpkg = Rcpp::Environment::namespace_env("igraph");
  // Rcpp::Function igfxn_graph_adjacency = igpkg["graph.adjacency"];
  // Rcpp::Function igfxn_distances = igpkg["distances"];
  // 
  // 
  // 
  // // Define the graph with given adjacency
  // Rcpp::List ig; 
  // if(weighted && directed) {
  //   ig = igfxn_graph_adjacency(Rcpp::Named("adjmatrix", ADJ), Rcpp::Named("mode","directed"), Rcpp::Named("diag", false), Rcpp::Named("weighted",true));  
  // } else if(weighted && !directed) {
  //   ig = igfxn_graph_adjacency(Rcpp::Named("adjmatrix", ADJ), Rcpp::Named("mode","undirected"), Rcpp::Named("diag", false), Rcpp::Named("weighted",true));  
  // } else if(!weighted && directed) {
  //   ig = igfxn_graph_adjacency(Rcpp::Named("adjmatrix", ADJ), Rcpp::Named("mode","directed"), Rcpp::Named("diag", false), Rcpp::Named("weighted",R_NilValue));  
  // } else {
  //   ig = igfxn_graph_adjacency(Rcpp::Named("adjmatrix", ADJ), Rcpp::Named("mode","undirected"), Rcpp::Named("diag", false), Rcpp::Named("weighted",R_NilValue));  
  // }
  // 
  // 
  // // Compute shortest paths between all vertices 
  // Rcpp::IntegerMatrix SPDISTRCPP = igfxn_distances(Rcpp::Named("graph", ig)); 
  // 
  // // Convert to umat and return 
  // arma::umat SPDIST = Rcpp::as<arma::umat>(SPDISTRCPP); 
  // 
  // return SPDIST; 
}


