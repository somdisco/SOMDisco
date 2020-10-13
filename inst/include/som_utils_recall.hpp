#ifndef RcppArmadillo_H
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifndef RcppParallel_H
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
#endif

// [[Rcpp::plugins(cpp11)]]

#ifndef SOM_UTILS_DIST_HPP
#include "som_utils_dist.hpp"
#endif


#ifndef SOM_UTILS_RECALL_HPP
#define SOM_UTILS_RECALL_HPP

namespace SOM_UTILS_RECALL {

// *** Parallel worker to find BMU, both with and without bias ***
struct find_BMU_worker : public RcppParallel::Worker {
  
  // Inputs 
  const arma::mat& X;
  const arma::mat& W;
  const unsigned int& nBMU; // num. of BMUs
  
  unsigned int nX, nW, d; 
  arma::rowvec Xmin, Xmax, Xrng, Wmin, Wmax, Wrng; 
  arma::vec bias; 
  
  // output container
  arma::umat BMU;
  arma::vec SQE; // squared quantization error 
  
  // Constructor 
  find_BMU_worker(const arma::mat& X, const arma::mat& W, const unsigned int& nBMU)
    : X(X), W(W), nBMU(nBMU)
  {
    nX = X.n_rows; 
    nW = W.n_rows; 
    d = X.n_cols; 
    if(W.n_cols != d) Rcpp::stop("ncol(X) != ncol(W)");
    
    // Initialize ranges to [0,1]
    Xmin.set_size(d); Xmin.zeros(); 
    Xmax.set_size(d); Xmax.ones(); 
    Xrng = Xmax - Xmin; 
    Wmin.set_size(d); Wmin.zeros(); 
    Wmax.set_size(d); Wmax.ones(); 
    Wrng = Wmax - Wmin; 
    
    // Initialize bias to 0 
    bias.set_size(nW); bias.zeros(); 
      
    BMU.set_size(nX, nBMU);
    SQE.set_size(nX);
  }
  
  // Set bias (optional)
  void set_bias(const arma::vec& bias_) {
    if(bias_.n_elem != nW) Rcpp::stop("length(bias) != nrow(W)");
    this->bias = bias_; 
  }
  
  // Set ranges (optional)
  void set_ranges(const arma::rowvec& Xmin_, const arma::rowvec& Xmax_, const arma::rowvec& Wmin_, const arma::rowvec& Wmax_) {
    
    // Set Xmin 
    if(Xmin_.n_elem == 1) {
      Xmin.fill(Xmin_[0]); 
    } else if(Xmin_.n_elem == d) {
      Xmin = Xmin_; 
    } else {
      Rcpp::stop("length(Xmin) != ncol(X)");
    }

    // Set Xmax 
    if(Xmax_.n_elem == 1) {
      Xmax.fill(Xmax_[0]); 
    } else if(Xmax_.n_elem == d) {
      Xmax = Xmax_; 
    } else {
      Rcpp::stop("length(Xmax) != ncol(X)");
    }
    
    // Set Xrng 
    Xrng = Xmax - Xmin; 
    
    // Set Wmin 
    if(Wmin_.n_elem == 1) {
      Wmin.fill(Wmin_[0]); 
    } else if(Wmin_.n_elem == d) {
      Wmin = Wmin_; 
    } else {
      Rcpp::stop("length(Wmin) != ncol(W)");
    }
    
    // Set Wmax
    if(Wmax_.n_elem == 1) {
      Wmax.fill(Wmax_[0]); 
    } else if(Wmax_.n_elem == d) {
      Wmax = Wmax_; 
    } else {
      Rcpp::stop("length(Wmax) != ncol(W)");
    }
    
    // Set Wrng 
    Wrng = Wmax - Wmin; 
    
    return; 
  }
  
  // Scale a single data vector to the network range 
  arma::rowvec scale_x_to_netrng(unsigned int i) {
    arma::rowvec x = (X.row(i) - Xmin) / Xrng % Wrng + Wmin; 
    return x; 
  }
  
  // Set the BMU information for a single row (i) of X
  void bmu_of_x(unsigned int i) {
    // Container to store the (squared) distance to each codebook vector in W
    arma::vec dist_to_W(nW);
    
    // Scale the data to network range 
    arma::rowvec x = scale_x_to_netrng(i);
    
    for(unsigned int j = 0; j < nW; ++j) {
      dist_to_W(j) = SOM_UTILS_DIST::dist_L22(x, W.row(j), bias(j));
    }
    
    arma::uvec dist_idxsort = arma::sort_index(dist_to_W, "ascend");
    
    for(unsigned int k = 0; k < nBMU; ++k) {
      BMU(i,k) = dist_idxsort(k);
    }
    
    SQE(i) = dist_to_W(dist_idxsort(0));
  }
  
  
  // Parallel operator - find BMU of each row of X in parallel
  void operator()(std::size_t begin, std::size_t end) {
    for(unsigned int i = begin; i < end; i++) {
      bmu_of_x(i);
    }
  }
  
  // Parallel call method
  void calc_parallel() {
    RcppParallel::parallelFor(0, nX, *this);
  }
  
  // Non-parallel call method
  void calc_serial() {
    // Find BMU of each row of X
    for(unsigned int i=0; i<nX; ++i) {
      bmu_of_x(i);
    }
  }
};


// *** Parallel worker to build CADJ from BMU list ***
struct build_CADJ_worker : public RcppParallel::Worker {
  
  // Inputs 
  const arma::umat& BMU; // must have at least 2 columns
  unsigned int nprotos; 
  bool zeroidx; 
  
  // Accumulated results
  arma::umat CADJ; 
  
  // constructors
  build_CADJ_worker(const arma::umat& BMU, unsigned int nprotos_, bool zeroidx_) : 
    BMU(BMU), nprotos(nprotos_), zeroidx(zeroidx_), CADJ(arma::zeros<arma::umat>(nprotos,nprotos)) {
    if(BMU.n_cols < 2) Rcpp::stop("BMU matrix must have at least 2 columns to build CADJ");
  };
  
  build_CADJ_worker(const build_CADJ_worker& dummy, RcppParallel::Split) : 
    BMU(dummy.BMU), nprotos(dummy.nprotos), zeroidx(dummy.zeroidx), CADJ(arma::zeros<arma::umat>(nprotos,nprotos)) {};
  
  
  void build_CADJ_loop(unsigned int i) {
    if(zeroidx) {
      CADJ( BMU(i,0), BMU(i,1) ) += 1;   
    } else {
      CADJ( BMU(i,0)-1, BMU(i,1)-1 ) += 1; 
    }
    
  }
  
  
  // process just the elements of the range I've been asked to
  void operator()(std::size_t begin, std::size_t end) {
    for(unsigned int i = begin; i < end; i++) {
      build_CADJ_loop(i);
    }
  }
  
  
  // join my values with that of another thread
  void join(const build_CADJ_worker& rhs) {
    CADJ += rhs.CADJ;
  }
  
  
  void calc_parallel() {
    RcppParallel::parallelReduce(0, BMU.n_rows, *this);
  }
  
  
  void calc_serial() {
    for(unsigned int i=0; i<BMU.n_rows; ++i) {
      build_CADJ_loop(i);
    }
  }
};


// *** Parallel worker to find the RF members, given a BMU matrix ***
struct find_RF_members_worker : public RcppParallel::Worker {
  
  // Inputs 
  const arma::umat& BMU;
  unsigned int nprotos; 
  bool zeroidx; 
  
  // Output container 
  std::vector<arma::uvec> RF_members;
  
  
  // Constructor with bias 
  find_RF_members_worker(const arma::umat& BMU, unsigned int nprotos_, bool zeroidx_)
    : BMU(BMU)
  {
    nprotos = nprotos_; 
    zeroidx = zeroidx_; 
    RF_members.resize(nprotos); 
  }
  
  
  // Find members of a single RF 
  void members_of_j(unsigned int j) {
    unsigned int lookup = j; 
    if(!zeroidx) lookup += 1; 
    
    RF_members[j] = arma::find(BMU.col(0) == lookup); 
    
    if(!zeroidx) RF_members[j] += 1; 
  }
  
  
  // Parallel operator - find BMU of each row of X in parallel
  void operator()(std::size_t begin, std::size_t end) {
    for(unsigned int j = begin; j < end; j++) {
      members_of_j(j);
    }
  }
  
  // Parallel call method
  void calc_parallel() {
    RcppParallel::parallelFor(0, nprotos, *this);
  }
  
  // Non-parallel call method
  void calc_serial() {
    // Find BMU of each row of X
    for(unsigned int j=0; j<nprotos; ++j) {
      members_of_j(j);
    }
  }
};


// *** Parallel worker to find the plularity label of a RF, given a BMU matrix and data labels 
// struct find_RF_label_worker : public RcppParallel::Worker {
//   
//   // Inputs 
//   const Rcpp::CharacterVector& Xlbl; 
//   const std::vector<arma::uvec>& RF_members; // vector of data indices mapped to each RF 
//   bool zeroidx; // whether the data indices are 1-based or 0-based
//   
//   // Output container 
//   Rcpp::CharacterVector RF_label; 
//   
//   
//   // Constructor with bias 
//   find_RF_label_worker(const Rcpp::CharacterVector& Xlbl, const std::vector<arma::uvec>& RF_members, bool zeroidx_)
//     : Xlbl(Xlbl), RF_members(RF_members), zeroidx(zeroidx_)
//   {
//     RF_label = Rcpp::CharacterVector(RF_members.size()); 
//     RF_label.fill(NA_STRING);
//   }
//   
//   
//   // Find label of a single RF 
//   void label_of_j(unsigned int j) {
//     
//     // If this RF is empty (has no members) leave the string as NA  
//     if(RF_members[j].size()==0) return; 
//     
//     // Otherwise, create a map to hold the label count and fill it up 
//     std::map<std::string,unsigned int> lbl_map;
//     if(zeroidx) {
//       for(unsigned int k=0; k < RF_members[j].size(); ++k) {
//         ++lbl_map[ Rcpp::as<std::string>(Xlbl[ RF_members[j][k] ])  ];
//       } 
//     } else {
//       for(unsigned int k=0; k < RF_members[j].size(); ++k) {
//         ++lbl_map[ Rcpp::as<std::string>(Xlbl[ RF_members[j][k]-1 ])  ];
//       } 
//     }
//     
//     typedef decltype(std::pair<std::string,unsigned int>()) pair_type;
//     
//     auto comp = [](const pair_type &pair1, const pair_type &pair2) -> bool {
//       return pair1.second < pair2.second; };
//     
//     RF_label[j] = std::max_element(lbl_map.cbegin(), lbl_map.cend(), comp)->first;
//     
//     return; 
//   }
//   
//   
//   // Parallel operator - find BMU of each row of X in parallel
//   void operator()(std::size_t begin, std::size_t end) {
//     for(unsigned int j = begin; j < end; j++) {
//       label_of_j(j);
//     }
//   }
//   
//   // Parallel call method
//   void calc_parallel() {
//     RcppParallel::parallelFor(0, RF_members.size(), *this);
//   }
//   
//   // Non-parallel call method
//   void calc_serial() {
//     // Find BMU of each row of X
//     for(unsigned int j=0; j<RF_members.size(); ++j) {
//       label_of_j(j);
//     }
//   }
// };

struct find_RF_label_worker : public RcppParallel::Worker {
  
  // Inputs 
  const Rcpp::CharacterVector& Xlbl; 
  const std::vector<arma::uvec>& RF_members; // vector of data indices mapped to each RF 
  bool zeroidx; // whether the data indices are 1-based or 0-based
  
  // Output container 
  std::vector<std::vector<int>> RF_labeldist_count;
  std::vector<std::vector<std::string>> RF_labeldist_label; 
  Rcpp::CharacterVector RF_label; 

  // Internal vars 
  std::string std_NA; 

  // Constructor with bias 
  find_RF_label_worker(const Rcpp::CharacterVector& Xlbl, const std::vector<arma::uvec>& RF_members, bool zeroidx_)
    : Xlbl(Xlbl), RF_members(RF_members), zeroidx(zeroidx_)
  {
    RF_labeldist_count.resize(RF_members.size()); 
    RF_labeldist_label.resize(RF_members.size()); 

    RF_label = Rcpp::CharacterVector(RF_members.size()); 
    RF_label.fill(NA_STRING);

    std_NA = "NA_STRING";
  }
  
  
  // Find label of a single RF 
  void label_of_j(unsigned int j) {
    
    // If this RF is empty (has no members) leave the string as NA
    if(RF_members[j].n_elem==0) {
      std::vector<int> emptyintvec; 
      std::vector<std::string> emptystrvec;
      RF_labeldist_count[j] = emptyintvec; 
      RF_labeldist_label[j] = emptystrvec; 
      return;
    } 

    // Otherwise, create a map of <label,count>
    std::map<std::string,int> lbl_map;
    for(unsigned int k=0; k<RF_members[j].n_elem; ++k) {
      unsigned int Xlblidx = RF_members[j][k]; 
      if(!zeroidx) Xlblidx -= 1; 
      
      if(Rcpp::CharacterVector::is_na(Xlbl[ Xlblidx ])) {
        ++lbl_map[ std_NA  ];
      } else {
        ++lbl_map[ Rcpp::as<std::string>(Xlbl[ Xlblidx ])  ];  
      }
    }
    
    // Declare a multimap of <count,label>, sorted in decreasing order, 
    // and insert the above map into it 
    std::multimap<int, std::string, std::greater<int>> MM;
    for(auto& it : lbl_map) {
      MM.insert({ it.second, it.first });
    }

    // Initialize vectors to hold the counts & labels from MM, and fill them up
    std::vector<int> counts(lbl_map.size());
    std::vector<std::string> labels(lbl_map.size());

    unsigned int counter = 0;
    //std::multimap<int, std::string>::const_iterator it2;
    //for(it2 = MM.begin(); it2 != MM.end(); ++it2) {
    for(auto& it2 : MM) {
      counts[counter] = it2.first;
      labels[counter] = it2.second;
      counter++;
    }
    
    RF_label[j] = labels[0];
    RF_labeldist_count[j] = counts; 
    RF_labeldist_label[j] = labels; 

    // Rcpp::CharacterVector counts_names(counts.size()); 
    // for(unsigned int k=0; k<counts.size(); ++k) {
    //   if(labels[k] == std_NA) {
    //     counts_names[k] = NA_STRING; 
    //   } else {
    //     counts_names[k] = Rcpp::String(labels[k]);
    //   }
    // }
    //counts.names() = Rcpp::wrap(labels); 
    // 
    // RF_label_dist[j] = counts; 
    // RF_label[j] = labels[0]; 
    
    return; 
    //typedef decltype(std::pair<std::string,unsigned int>()) pair_type;
    
    // auto comp = [](const pair_type &pair1, const pair_type &pair2) -> bool {
    //   return pair1.second < pair2.second; };
    // 
    //RF_label[j] = std::max_element(lbl_map.cbegin(), lbl_map.cend(), comp)->first;
    //std::sort(lbl_map.begin(), lbl_map.end(), comp); 
    //RF_label[j] = lbl_map.begin()->first; 
    
    
    // // Declare vector of pairs 
    // std::vector<std::pair<std::string, unsigned int> > A; 
    // 
    // // Copy key-value pair from map to vector of pairs 
    // for(auto& it : lbl_map) { 
    //   A.push_back(it); 
    // } 
    // 
    // // comparator function 
    // // inline bool cmp(std::pair<std::string, unsigned int>& a, 
    // //          std::pair<std::string, unsigned int>& b) 
    // // { 
    // //   return a.second < b.second; 
    // // } 
    // auto cmp = [](std::pair<std::string, unsigned int>& a, 
    //               std::pair<std::string, unsigned int>& b) 
    // { 
    //   return a.second < b.second; 
    // };
    // 
    // // Sort using comparator function 
    // std::sort(A.rbegin(), A.rend(), cmp); 
    // 
    // // // Print the sorted value 
    // // for (auto& it : A) { 
    // //   
    // //   cout << it.first << ' '
    // //        << it.second << endl; 
    // // } 
    // 
    // RF_label[j] = A.begin()->first; 
    // 
    // return; 
  }
  
  
  // Parallel operator - find BMU of each row of X in parallel
  void operator()(std::size_t begin, std::size_t end) {
    for(unsigned int j = begin; j < end; j++) {
      label_of_j(j);
    }
  }
  
  // Parallel call method
  void calc_parallel() {
    RcppParallel::parallelFor(0, RF_members.size(), *this);
  }
  
  // Non-parallel call method
  void calc_serial() {
    // Find BMU of each row of X
    for(unsigned int j=0; j<RF_members.size(); ++j) {
      label_of_j(j);
    }
  }
  
  // Re-combine the RF_labeldist into vec of Named IntegerVectors
  std::vector<Rcpp::IntegerVector> process_RF_labeldist() {
    std::vector<Rcpp::IntegerVector> out(RF_members.size()); 
    
    Rcpp::IntegerVector emptyvec(0); 
    for(unsigned int j=0; j<RF_members.size(); ++j) {
      if(RF_members[j].n_elem == 0) {out[j] = emptyvec; continue;}
      
      Rcpp::IntegerVector this_dist = Rcpp::wrap(RF_labeldist_count[j]); 
      Rcpp::CharacterVector this_dist_label(RF_labeldist_count[j].size()); 
      for(unsigned int k=0; k<RF_labeldist_count[j].size(); ++k) {
        if(RF_labeldist_label[j][k] == "NA_STRING")
          this_dist_label[k] = NA_STRING; 
        else
          this_dist_label[k] = RF_labeldist_label[j][k];
      }
      
      this_dist.names() = this_dist_label; 
      out[j] = this_dist; 
    }
    
    return out; 
  }
};



// *** Summary of values in each Receptive Field of a Vector Quantizer 
struct groupstat_revmap_worker : public RcppParallel::Worker {
  
  // Inputs
  const arma::mat vals; // values whose summaries will be computed, nX rows x d columns
  std::string stat; // which summary stat to compute 
  const std::vector<arma::uvec>& revmap; // vector of data indices mapped to each RF. Must have length = # prototypes in VQ
  bool zeroidx; // whether the data indices are 1-based or 0-based
  
  // Output container
  arma::mat stats;
  
  // Constructor
  groupstat_revmap_worker(const arma::mat& vals, const std::vector<arma::uvec>& revmap, bool zeroidx_)
    : vals(vals), revmap(revmap)
  {
    
    stat = ""; 
    
    stats.set_size(revmap.size(), vals.n_cols);
    stats.fill(NA_REAL);
    
    zeroidx = zeroidx_;
    
    // Check the revmap
    // Each component of this vector is a vector of indices describing which elements (rows) of vals matrix
    // are mapped to each prototype.  Thus, the number of indices it contains should = nrow(vals).
    unsigned int sum_members = 0, max_idx = 0;
    for(unsigned int i=0; i<revmap.size(); ++i) {
      sum_members += revmap[i].n_elem;
      if(revmap[i].n_elem > 0) {
        unsigned int maxidx_this_group = revmap[i].max();
        if(maxidx_this_group > max_idx) max_idx = maxidx_this_group;
      }
    }
    
    if(sum_members != vals.n_rows) Rcpp::stop("Input vals should have nrows = number of members in revmap");
    if(zeroidx) max_idx += 1;
    if(max_idx > vals.n_rows) Rcpp::stop("Input revmap references data index out of range of input vals");
  }
  
  void set_stat(std::string stat_) {
    // Summary function can only be count, sum, mean, sd, q0, q25, q50, q75, q100
    std::vector<std::string> possible_stats(9);
    possible_stats[0] = "count";
    possible_stats[1] = "sum";
    possible_stats[2] = "mean";
    possible_stats[3] = "sd";
    possible_stats[4] = "q0";
    possible_stats[5] = "q25";
    possible_stats[6] = "q50";
    possible_stats[7] = "q75";
    possible_stats[8] = "q100";
    
    std::vector<std::string>::const_iterator it = std::find(possible_stats.begin(), possible_stats.end(), stat_);
    if(it == possible_stats.end()) {
      Rcpp::stop("Input 'stat' must be one of count, sum, mean, sd, q0, q25, q50, q75, q100");
    } else {
      this->stat = stat_;
    }
    
    stats.fill(NA_REAL);
  }
  
  
  // Find stats of a single group
  void summary_group_j(unsigned int j) {
    
    // If this group is empty (has no members) leave the stat as NA (all stats initialized to NA in constructor)
    if(revmap[j].n_elem==0) return;
    
    // Otherwise strip out the members indices
    arma::uvec membersj = revmap[j];
    if(!zeroidx) membersj = membersj - 1;
    
    // Otherwise, decode which stat is to be calculated, and compute it
    if(stat == "count") {
      stats[j] = double(membersj.n_elem);
    } else if(stat == "sum") {
      stats.row(j) = arma::sum(vals.rows(membersj), 0);
    } else if(stat == "mean") {
      stats.row(j) = arma::mean(vals.rows(membersj), 0);
    } else if(stat == "sd") {
      if(membersj.n_elem == 1) {
        stats.row(j).fill(0);
      } else {
        stats.row(j) = arma::stddev(vals.rows(membersj), 0, 0);
      }
    } else if(stat == "q0") {
      stats.row(j) = arma::min(vals.rows(membersj), 0);
    } else if(stat == "q25") {
      arma::vec probs = {0.25};
      arma::mat quants = arma::quantile(vals.rows(membersj), probs, 0);
      stats.row(j) = quants.row(0);
    } else if(stat == "q50") {
      stats.row(j) = arma::median(vals.rows(membersj), 0);
    } else if(stat == "q75") {
      arma::vec probs = {0.75};
      arma::mat quants = arma::quantile(vals.rows(membersj), probs, 0);
      stats.row(j) = quants.row(0);
    } else {
      stats.row(j) = arma::max(vals.rows(membersj), 0);
    }
    
    return;
  }
  
  
  // Parallel operator - find BMU of each row of X in parallel
  void operator()(std::size_t begin, std::size_t end) {
    for(unsigned int j = begin; j < end; j++) {
      summary_group_j(j);
    }
  }
  
  // Parallel call method
  void calc_parallel() {
    if(stat == "") Rcpp::stop("Must call set_stat before computing stats");
    RcppParallel::parallelFor(0, revmap.size(), *this);
  }
  
  // Non-parallel call method
  void calc_serial() {
    if(stat == "") Rcpp::stop("Must call set_stat before computing stats");
    // Find BMU of each row of X
    for(unsigned int j=0; j<revmap.size(); ++j) {
      summary_group_j(j);
    }
  }
};

struct groupstat_fwdmap_worker : public RcppParallel::Worker {
  
  // Inputs
  const arma::mat vals; // values whose summaries will be computed, nX rows x p columns
  std::string stat; // which summary stat to compute 
  const arma::uvec& fwdmap; // vector of data indices mapped to each RF. Must have length = # prototypes in VQ
  unsigned int ngroups; 
  bool zeroidx; // whether the group indices in fwdmap are 1-based or 0-based
  
  // Output container
  arma::mat stats; // nrows = ngroups, ncols = ncol(vals)
  
  // Constructor
  groupstat_fwdmap_worker(const arma::mat& vals, const arma::uvec& fwdmap, unsigned int ngroups_, bool zeroidx_)
    : vals(vals), fwdmap(fwdmap)
  {
   
    stat = ""; 
    
    ngroups = ngroups_; 
    stats.set_size(ngroups, vals.n_cols);
    stats.fill(NA_REAL);
    
    zeroidx = zeroidx_;
    
    // Check the fwdmap
    // Each component of this vector is an indicator identifying group membership of the rows of vals matrix
    // It should be the same size as vals
    if(fwdmap.n_elem != vals.n_rows) Rcpp::stop("Input vals should have nrows = length(fwdmap)");
    
    // Also, it should not contain any elements > ngroups (if zeroidx=F) or > ngroups-1 (if zeroidx=T)
    unsigned int max_idx = fwdmap.max();
    if(zeroidx) max_idx += 1;
    if(max_idx > ngroups) Rcpp::stop("Input fwdmap references group indices > ngroups");
  }
  
  void set_stat(std::string stat_) {
    // Summary function can only be count, sum, mean, sd, q0, q25, q50, q75, q100
    std::vector<std::string> possible_stats(9);
    possible_stats[0] = "count";
    possible_stats[1] = "sum";
    possible_stats[2] = "mean";
    possible_stats[3] = "sd";
    possible_stats[4] = "q0";
    possible_stats[5] = "q25";
    possible_stats[6] = "q50";
    possible_stats[7] = "q75";
    possible_stats[8] = "q100";
    
    std::vector<std::string>::const_iterator it = std::find(possible_stats.begin(), possible_stats.end(), stat_);
    if(it == possible_stats.end()) {
      Rcpp::stop("Input 'stat' must be one of count, sum, mean, sd, q0, q25, q50, q75, q100");
    } else {
      this->stat = stat_;
    }
    
    stats.fill(NA_REAL);
  }
  
  // Find stats of a single group
  void summary_group_j(unsigned int j) {
    
    // Find elements in group j 
    unsigned int searchj = j; 
    if(!zeroidx) searchj += 1; 
    arma::uvec membersj = arma::find(fwdmap == searchj); 
    
    // If this group is empty (has no members) leave the stat as NA (all stats initialized to NA in constructor)
    if(membersj.n_elem==0) return;
    
    // Otherwise, decode which stat is to be calculated, and compute it
    if(stat == "count") {
      stats[j] = double(membersj.n_elem);
    } else if(stat == "sum") {
      stats.row(j) = arma::sum(vals.rows(membersj), 0);
    } else if(stat == "mean") {
      stats.row(j) = arma::mean(vals.rows(membersj), 0);
    } else if(stat == "sd") {
      if(membersj.n_elem == 1) {
        stats.row(j).fill(0);
      } else {
        stats.row(j) = arma::stddev(vals.rows(membersj), 0, 0);
      }
    } else if(stat == "q0") {
      stats.row(j) = arma::min(vals.rows(membersj), 0);
    } else if(stat == "q25") {
      arma::vec probs = {0.25};
      arma::mat quants = arma::quantile(vals.rows(membersj), probs, 0);
      stats.row(j) = quants.row(0);
    } else if(stat == "q50") {
      stats.row(j) = arma::median(vals.rows(membersj), 0);
    } else if(stat == "q75") {
      arma::vec probs = {0.75};
      arma::mat quants = arma::quantile(vals.rows(membersj), probs, 0);
      stats.row(j) = quants.row(0);
    } else {
      stats.row(j) = arma::max(vals.rows(membersj), 0);
    }
    
    return;
  }
  
  
  // Parallel operator - find BMU of each row of X in parallel
  void operator()(std::size_t begin, std::size_t end) {
    for(unsigned int j = begin; j < end; j++) {
      summary_group_j(j);
    }
  }
  
  // Parallel call method
  void calc_parallel() {
    if(stat == "") Rcpp::stop("Must call set_stat before computing stats");
    RcppParallel::parallelFor(0, ngroups, *this);
  }
  
  // Non-parallel call method
  void calc_serial() {
    if(stat == "") Rcpp::stop("Must call set_stat before computing stats");
    // Find BMU of each row of X
    for(unsigned int j=0; j<ngroups; ++j) {
      summary_group_j(j);
    }
  }
};


// *** Compute the fence summary data frame in parallel 
struct fence_worker : public RcppParallel::Worker {
  
  const arma::mat& nu_xy; 
  const arma::mat& W; 
  const arma::cube& nu_verts; 
  const arma::umat& nu_ADJ;
  
  // Internal variables
  arma::uvec nfences_from_j; 
  unsigned int nfences; 
  
  // output containers
  arma::umat nu_jk;
  arma::mat fence_endpts; 
  arma::vec fence_values; 
  
  // Constructor
  fence_worker(const arma::mat& nu_xy, const arma::mat& W, const arma::cube& nu_verts, const arma::umat& nu_ADJ) : 
    nu_xy(nu_xy), W(W), nu_verts(nu_verts), nu_ADJ(nu_ADJ) {
    
    arma::uvec adjacent = arma::find(arma::trimatu(nu_ADJ, 1));
    nu_jk = arma::ind2sub(arma::size(nu_ADJ), adjacent).t(); 
    
    nfences = nu_jk.n_rows;  
    
    fence_values.set_size(nfences); 
    //fence_endpts.resize(W.n_rows); 
    fence_endpts.set_size(nfences,4); // x0,y0,x1,y1 
  }
  
  // Fill up neighbors from a single neuron 
  void fences_of_tile_j(unsigned int j) {
    
    // Find which fences originate from this j
    // If none, skip this j
    arma::uvec fence_idx = arma::find(nu_jk.col(0) == j);
    if(fence_idx.size() == 0) return;
    
    // Compute midpoints between each consecutive pair of vertices in nu_verts for tile j
    unsigned int nverts = nu_verts.slice(j).n_rows;
    arma::mat vert_mid(nverts, 2);
    for(unsigned int k=0; k<(nverts-1); ++k) {
      vert_mid.row(k) = (nu_verts.slice(j).row(k) + nu_verts.slice(j).row(k+1)) / 2.0;
    }
    vert_mid.row(nverts-1) = (nu_verts.slice(j).row(nverts-1) + nu_verts.slice(j).row(0)) / 2.0;
    
    
    // For each fence from tile j
    for(unsigned int k=0; k<fence_idx.size(); ++k) {
      // Compute distance from W_j to W_kneighbor
      fence_values[fence_idx[k]] = SOM_UTILS_DIST::dist_L22(W.row(j), W.row(nu_jk(fence_idx[k],1)), 0.0);
      
      // Find lattice midpoint between nu_j and nu_kneighbor
      arma::rowvec nu_mid = (nu_xy.row(j) + nu_xy.row(nu_jk(fence_idx[k],1))) / 2.0;
      
      // Determine which vertex pair's midpoint is closes to the neuron midpoint
      arma::mat tmp_mid = vert_mid;
      tmp_mid.each_row() -= nu_mid;
      arma::vec dist_to_nu_mid = arma::sum(arma::square(tmp_mid), 1);
      unsigned int which_vertex_mid = dist_to_nu_mid.index_min();
      
      // Store endpoints
      fence_endpts(fence_idx[k],0) = nu_verts.slice(j)(which_vertex_mid,0);
      fence_endpts(fence_idx[k],1) = nu_verts.slice(j)(which_vertex_mid,1);
      if(which_vertex_mid < nverts-1) {
        fence_endpts(fence_idx[k],2) = nu_verts.slice(j)(which_vertex_mid+1,0);
        fence_endpts(fence_idx[k],3) = nu_verts.slice(j)(which_vertex_mid+1,1);
      } else {
        fence_endpts(fence_idx[k],2) = nu_verts.slice(j)(0,0);
        fence_endpts(fence_idx[k],3) = nu_verts.slice(j)(0,1);
      }
      
    }
    
  }
  
  // Parallel operator 
  void operator()(std::size_t begin, std::size_t end) {
    
    for(unsigned int j=begin; j<end; ++j) {
      fences_of_tile_j(j); 
    }
  }
  
  // Parallel invoker 
  void calc_parallel() {
    RcppParallel::parallelFor(0, W.n_rows, *this);
  }
  
  // Serial invoker 
  void calc_serial() {
    for(unsigned int j=0; j<W.n_rows; ++j) {
      fences_of_tile_j(j); 
    }
  }
  
  Rcpp::DataFrame fence_dataframe() {
    Rcpp::DataFrame out = Rcpp::DataFrame::create(Rcpp::Named("i") = nu_jk.col(0) + 1, 
                                                  Rcpp::Named("j") = nu_jk.col(1) + 1,
                                                  Rcpp::Named("x0") = fence_endpts.col(0), 
                                                  Rcpp::Named("y0") = fence_endpts.col(1), 
                                                  Rcpp::Named("x1") = fence_endpts.col(2), 
                                                  Rcpp::Named("y1") = fence_endpts.col(3), 
                                                  Rcpp::Named("value") = fence_values); 
    return out; 
  }
  
};


} // close namespace 

#endif


/* 
 // Version WITHOUT bias, WITHOUT ranges
 inline void cpp_find_BMU(arma::umat& BMU, arma::vec& SQE, const arma::mat& X, const arma::mat& W, unsigned int nBMU, bool parallel) {
 find_BMU_worker wkr(X, W, nBMU); 
 
 if(parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
 
 BMU = wkr.BMU; 
 SQE = wkr.SQE; 
 return; 
 }
 
 // Version WITH bias, WITHOUT ranges
 inline void cpp_find_BMU(arma::umat& BMU, arma::vec& SQE, const arma::mat& X, const arma::mat& W, unsigned int nBMU, bool parallel, const arma::vec& bias) {
 find_BMU_worker wkr(X, W, nBMU); 
 wkr.set_bias(bias); 
 
 if(parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
 
 BMU = wkr.BMU; 
 SQE = wkr.SQE; 
 return; 
 }
 
 // Version WITHOUT bias, WITH ranges
 inline void cpp_find_BMU(arma::umat& BMU, arma::vec& SQE, const arma::mat& X, const arma::mat& W, unsigned int nBMU, bool parallel, 
 const arma::rowvec& Xmin, const arma::rowvec& Xmax, const arma::rowvec& Wmin, const arma::rowvec& Wmax) {
 find_BMU_worker wkr(X, W, nBMU); 
 wkr.set_ranges(Xmin, Xmax, Wmin, Wmax); 
 
 if(parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
 
 BMU = wkr.BMU; 
 SQE = wkr.SQE; 
 return; 
 }
 
 // Version WITH bias, WITH ranges
 inline void cpp_find_BMU(arma::umat& BMU, arma::vec& SQE, const arma::mat& X, const arma::mat& W, unsigned int nBMU, bool parallel, 
 const arma::vec& bias, const arma::rowvec& Xmin, const arma::rowvec& Xmax, const arma::rowvec& Wmin, const arma::rowvec& Wmax) {
 find_BMU_worker wkr(X, W, nBMU); 
 wkr.set_bias(bias); 
 wkr.set_ranges(Xmin, Xmax, Wmin, Wmax); 
 
 if(parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
 
 BMU = wkr.BMU; 
 SQE = wkr.SQE; 
 return; 
 }
 */
 
 