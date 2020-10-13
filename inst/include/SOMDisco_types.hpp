#ifndef RcppArmadillo_H
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifndef RcppParallel_H
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
#endif

// [[Rcpp::plugins(cpp11)]]

#include <string>
#include <random>
#include <algorithm>

#ifndef SOM_UTILS_DIST_HPP
#include "som_utils_dist.hpp"
#endif

#ifndef SOM_UTILS_LATTICE_HPP
#include "som_utils_lattice.hpp"
#endif

#ifndef SOM_UTILS_SCALING_HPP
#include "som_utils_scaling.hpp"
#endif

#ifndef SOM_UTILS_RECALL_HPP
#include "som_utils_recall.hpp"
#endif


// RcppParallel is defining some variable named FALSE, 
// which causes errors when TRN_MODULE is exported 
#ifdef FALSE
#undef FALSE
#endif


#ifndef SOMDISCO_SOMOBJ_HPP
#define SOMDISCO_SOMOBJ_HPP

typedef std::map< unsigned int, std::tuple<double,double,double,unsigned int> > type_LRAS; 

// Forward declaration for RCPP_MODULES 
class SOM; 
RCPP_EXPOSED_CLASS(SOM);


class SOM{
private:
  
  
public:
  
  // CSOM constructor
  SOM(); 
  
  // Parallel flag 
  bool parallel; 
  void set_parallel(bool parallel_); 
  
  // Variables involving lattice 
  unsigned int som_x,som_y,nW; 
  std::string lattice_type; 
  arma::mat nu_xy; 
  arma::umat nu_ij;
  arma::cube nu_verts; 
  arma::umat nu_ADJ; // lattice adjacency matrix 
  std::vector<std::vector<arma::uvec>> nu_nhblist; // vector of neuron indices with in each distance 
  
  // Methods involving lattice 
  void set_lattice();
  arma::mat tile_interior_point(double theta, double radprop = 1.0); 
  
  // Variables involving training data 
  unsigned int d; // data dimension 
  unsigned int nX; // number of training vectors
  Rcpp::NumericVector X_stats; // 
  Rcpp::CharacterVector X_label; 
  Rcpp::DataFrame ctab; 
  void set_X_label(const Rcpp::CharacterVector& X_label_, Rcpp::DataFrame ctab_); 
  void set_ctab(Rcpp::DataFrame ctab_); 
  void initialize_SOM(const arma::mat& X, unsigned int som_x, unsigned int som_y, std::string lattice_type);
  
  
  // Variables involving network scaling 
  arma::rowvec netrng_ext_min, netrng_ext_max, netrng_ext_rng; 
  double netrng_int_min, netrng_int_max, netrng_int_rng; 
  void set_netrng(arma::rowvec ext_min, arma::rowvec ext_max, double int_min, double int_max); 
  arma::mat map_to_netrng(const arma::mat& x); 
  arma::mat map_from_netrng(const arma::mat& w); 
  
  // Variables involving prototype initialization 
  arma::mat W; 
  //void set_W_runif(unsigned int seed = 0); 
  void set_W_runif(); 
  void set_W(const arma::mat& W_); 
  void set_p_equal(); 
  void set_p(const arma::vec& p_); 
  arma::vec p; // winning frequencies  
  arma::vec bias(); 
  
  
  // Variables involing Learning Rate Annealing 
  type_LRAS LRAS; 
  double alpha, beta, gamma; 
  unsigned int sigma; 
  arma::vec eta; 
  
  
  // Methods involving Learning Rate Annealing
  void set_LRAS(Rcpp::DataFrame LRAS); 
  Rcpp::DataFrame get_LRAS();
  arma::vec calc_eta(unsigned int sigma_); 
  void update_learning_rates(); 
  
  
  // Variables involving training 
  unsigned int age; 
  Rcpp::IntegerVector train_order; 
  bool user_train_order; 
  unsigned int report_freq; 
  unsigned int mtr_freq; 
  std::vector<unsigned int> mtr_age;    // training ages at which monitoring snapshots were taken 
  std::vector<double> mtr_RMSQE;        // Root Mean Squared Quantization Error at monitoring snapshots 
  std::vector<double> mtr_QSI;          // Quantization stability at monitoring snapshots 
  std::vector<double> mtr_Entropy;      // Entropy at monitoring snapshots 
  
  // Methods involving training 
  void set_reporting_freq(unsigned int report_freq); 
  void update_p(unsigned int winner_idx); 
  void update_W(unsigned int winner_idx, const arma::rowvec& x); 
  void train_SOM(unsigned int nsteps, const arma::mat& X);
  void set_monitoring_freq(unsigned int mtr_freq_); 
  void monitor_SOM(const arma::mat& X);
  void set_train_order(const Rcpp::IntegerVector& train_order_); 
  
  
  // Variables involving recall 
  unsigned int nBMU;
  arma::umat BMU, CADJ;
  arma::uvec prevBMU; 
  arma::uvec RF_size; 
  double Entropy; 
  std::vector<arma::uvec> RF_members; 
  arma::vec SQE; 
  Rcpp::CharacterVector RF_label; 
  std::vector<Rcpp::IntegerVector> RF_label_dist; 
  Rcpp::DataFrame fences; 
  
  void set_nBMU(unsigned int k); 
  void recall_SOM(const arma::mat& X); 
  void mtrrecall_SOM(const arma::mat& X); 
  void set_RF_label();
  void set_lattice_fences(); 
  arma::umat CONN(); 
  void calc_Entropy(); 
  void clear_recall(); 
  
  // Control flags 
  bool is_netrng_set;
  bool is_protos_init;
  bool is_winfrq_init; 
  bool is_lras_set; 
  bool is_trained; 
  bool is_recalled; 
  
  // IO 
  void save(std::string rdsfile); 
  void load(std::string rdsfile); 
  void load_list(Rcpp::List SOMlist); 
  Rcpp::List as_list(); 
  
  // Vis 
  Rcpp::List vis_par;
  void set_vis_par(Rcpp::List vis_param_);
  
  Rcpp::CharacterVector vis_tile_bg; 
  void set_vis_tile_bg(Rcpp::CharacterVector tile_bg); 
  
  Rcpp::NumericVector vis_xlim; 
  void set_vis_xlim(Rcpp::NumericVector xlim); 
  
  Rcpp::NumericVector vis_ylim; 
  void set_vis_ylim(Rcpp::NumericVector ylim); 
};


// ***** Constructor 
inline SOM::SOM() {
  // Default initializations
  this->parallel = true; 
  this->som_x = 0; 
  this->som_y = 0; 
  this->nW = 0; 
  this->lattice_type = ""; 
  this->nu_xy.set_size(0,0); 
  this->nu_ij.set_size(0,0);
  this->nu_verts.set_size(0,0,0); 
  this->nu_ADJ.set_size(0,0); 
  this->nu_nhblist.resize(0);
  
  this->d = 0; // data dimension 
  this->nX = 0; // number of training vectors
  this->X_stats = Rcpp::NumericVector(0); // 
  
  this->netrng_ext_min.set_size(0); 
  this->netrng_ext_max.set_size(0);
  this->netrng_ext_rng.set_size(0); 
  this->netrng_int_min = 0;
  this->netrng_int_max = 0;
  this->netrng_int_rng = 0; 
  
  this->W.set_size(0,0); 
  this->p.set_size(0); // winning frequencies  
  
  this->LRAS.clear(); 
  this->alpha = 0.0;
  this->beta = 0.0;
  this->gamma = 0.0; 
  this->sigma = 0; 
  this->eta.set_size(0); 
  
  this->age = 0; 
  this->train_order = Rcpp::IntegerVector(0); 
  this->user_train_order = false; 
  this->report_freq = 10000; 
  this->mtr_freq = 10000; 
  this->mtr_age.resize(0); 
  this->mtr_RMSQE.resize(0);
  this->mtr_QSI.resize(0); // Quantization stability 
  this->mtr_Entropy.resize(0); // Entropy history 
  
  this->nBMU = 0;
  this->BMU.set_size(0);
  this->prevBMU.set_size(0);
  this->CADJ.set_size(0,0);
  this->RF_size.set_size(0); 
  this->Entropy = 0.0; 
  this->RF_members.resize(0); 
  this->SQE.set_size(0); 
  this->RF_label = Rcpp::CharacterVector(0); 
  this->RF_label_dist.resize(0); 
  this->fences = Rcpp::DataFrame::create(); 
  
  // Control flags 
  this->is_netrng_set = false;
  this->is_protos_init = false;
  this->is_winfrq_init = false; 
  this->is_lras_set = false; 
  this->is_trained = false; 
  this->is_recalled = false; 
  

}

inline void SOM::set_parallel(bool parallel_) {
  this->parallel = parallel_; 
}


// ***** Lattice Functions 
inline void SOM::set_lattice() {
  
  if(!(lattice_type=="grid" || lattice_type=="hex")) Rcpp::stop("lattice_type must be 'grid' or 'hex'");
  
  // Calculate & store the neuron grid points 
  Rcpp::Rcout << "++ calculating lattice (x,y) coordinates ... ";
  if(lattice_type=="hex") {
    SOM_UTILS_LATTICE::calc_lattice_hex(this->nu_xy, this->nu_ij, som_x, som_y);
  } else if(lattice_type=="grid") {
    SOM_UTILS_LATTICE::calc_lattice_grid(this->nu_xy, this->nu_ij, som_x, som_y);
  }
  Rcpp::Rcout << "done" << std::endl; 
  this->nW = this->nu_xy.n_rows;
  
  
  // Calculate & store the adjacency matrix of neurons 
  Rcpp::Rcout << "++ calculating neuron lattice adjacencies ... "; 
  this->nu_ADJ = SOM_UTILS_LATTICE::calc_lattice_adjacency(this->nu_xy);
  Rcpp::Rcout << "done" << std::endl;  
  
  // Calculate the geodesic distances of the neurons on the lattice using Floyd-Warshall 
  Rcpp::Rcout << "++ calculating geodesic lattice distances between neurons ... "; 
  arma::umat SPDIST = SOM_UTILS_DIST::geodesicdist(this->nu_ADJ, false, false); 
  Rcpp::Rcout << "done" << std::endl;  
  
  // Calculate & store the neuron distances as a list
  Rcpp::Rcout << "++ assigning geodesic lattice distances to distlist ... ";
  SOM_UTILS_LATTICE::lattice_distlist_worker distlist_worker(SPDIST);
  if(this->parallel) distlist_worker.calc_parallel(); else distlist_worker.calc_serial(); 
  this->nu_nhblist = distlist_worker.out;
  Rcpp::Rcout << "done" << std::endl;
  
  // Calculate & store lattice tile vertices
  Rcpp::Rcout << "++ calculating lattice tile vertices ... ";
  SOM_UTILS_LATTICE::lattice_tile_verts_worker vertsworker(this->nu_xy, this->lattice_type);
  if(this->parallel) vertsworker.calc_parallel(); else vertsworker.calc_serial(); 
  this->nu_verts = vertsworker.verts;
  Rcpp::Rcout << "done" << std::endl;
}

inline arma::mat SOM::tile_interior_point(double theta, double radprop) {
  // Theta should be in degrees 
  
  if(radprop > 0.0) {
    SOM_UTILS_LATTICE::tile_boundary_intersect_worker worker(this->nu_verts, this->nu_xy, theta); 
    if(this->parallel) worker.calc_parallel(); else worker.calc_serial(); 
    
    arma::rowvec ray_direction(2); 
    ray_direction[0] = std::cos(SOM_UTILS_LATTICE::deg2rad(theta)); 
    ray_direction[1] = std::sin(SOM_UTILS_LATTICE::deg2rad(theta)); 
    
    return this->nu_xy + radprop * worker.dist_to_edge * ray_direction; 
  } else {
    return this->nu_xy;
  }
}


// ***** Initialize SOM object and link to training data 
inline void SOM::initialize_SOM(const arma::mat& X, unsigned int som_x, unsigned int som_y, std::string lattice_type) {
  
  // Store general lattice properties   
  Rcpp::Rcout << "Setting lattice quantities" << std::endl; 
  this->som_x = som_x; 
  this->som_y = som_y; 
  this->lattice_type = lattice_type; 
  this->set_lattice(); 
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl; 
  
  // Store data dimension & Xstatistics
  Rcpp::Rcout << "Setting training data" << std::endl; 
  this->d = X.n_cols; 
  this->nX = X.n_rows; 
  Rcpp::Rcout << "++ recording " << this->nX << " observations of dimension " << this->d << std::endl; 
  
  this->X_stats = Rcpp::NumericVector({X.min(), X.max(), X(0), X(X.n_elem-1)});
  this->X_stats.names() = Rcpp::CharacterVector({"min","max","X[1]","X[end]"});
  
  // Initialize the labels to "?"
  Rcpp::CharacterVector tmp_X_label(this->nX); 
  std::fill(tmp_X_label.begin(), tmp_X_label.end(), "?"); 
  this->X_label = tmp_X_label; 
  
  Rcpp::DataFrame tmp_ctab = Rcpp::DataFrame::create(Rcpp::Named("label") = "?", Rcpp::Named("color") = "white"); 
  //Rcpp::Environment Renv = Rcpp::Environment::namespace_env("SOMDisco");
  //Rcpp::Function build_ctab = Renv["build_ctab"];
  //Rcpp::List tmp_ctab = build_ctab(this->X_label); 
  this->ctab = Rcpp::as<Rcpp::DataFrame>(tmp_ctab); 
  Rcpp::Rcout << "++ set all X_label = default '?'" << std::endl; 
  Rcpp::Rcout << "   change via $set_X_label" << std::endl; 
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl; 
  
  // Set default external & internal network ranges
  Rcpp::Rcout << "Setting network ranges" << std::endl; 
  arma::rowvec tmp_ext_min(this->d); tmp_ext_min.fill(this->X_stats["min"]);
  arma::rowvec tmp_ext_max(this->d); tmp_ext_max.fill(this->X_stats["max"]);
  SOM::set_netrng(tmp_ext_min, tmp_ext_max, 0.0, 1.0); 
  
  Rcpp::Rcout << std::fixed << std::setprecision(2) << "++ external = [" << X.min() << ", " << X.max() << "]" << std::endl;
  Rcpp::Rcout << std::fixed << std::setprecision(2) << "++ internal = [" << 0.0 << ", " << 1.0 << "]" << std::endl; 
  Rcpp::Rcout << "++ change defaults via $set_netrng" << std::endl; 
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl; 
  
  // Initialize prototypes to default 
  Rcpp::Rcout << "Initializing prototypes to random uniform" << std::endl; 
  SOM::set_W_runif(); 
  Rcpp::Rcout << "++ to set a particular random seed call set.seed() and then $set_W_runif()" << std::endl; 
  Rcpp::Rcout << "++ to set to specific values call $set_W()" << std::endl; 
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl; 
  
  // Initialize prototype win frequencies 
  Rcpp::Rcout << "Initializing prototype win frequencies to equiprobable" << std::endl; 
  SOM::set_p_equal(); 
  Rcpp::Rcout << "++ to set to specific values call $set_p(values)" << std::endl; 
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl; 
  
  
  // Initialize Learning Rate Schedule 
  Rcpp::Rcout << "Setting default learning rates" << std::endl; 
  double dblnX = double(this->nX); 
  Rcpp::NumericVector LRASt = Rcpp::NumericVector({1*dblnX, 5*dblnX, 10*dblnX, 25*dblnX, 100*dblnX});
  Rcpp::NumericVector LRASalpha = Rcpp::NumericVector({0.50, 0.25, 0.10, 0.05, 0.01});
  Rcpp::NumericVector LRASbeta = Rcpp::NumericVector({0.50/10, 0.25/10, 0.10/10, 0.05/10, 0.01/10});
  Rcpp::NumericVector LRASgamma = Rcpp::NumericVector({0.50/100, 0.25/100, 0.10/100, 0.05/100, 0.01/100});
  Rcpp::NumericVector LRASsigma = Rcpp::NumericVector({3, 2, 1, 1, 1});
  Rcpp::DataFrame tmpLRAS = Rcpp::DataFrame::create(Rcpp::Named("t") = LRASt, Rcpp::Named("alpha") = LRASalpha, Rcpp::Named("beta") = LRASbeta, Rcpp::Named("gamma") = LRASgamma, Rcpp::Named("sigma") = LRASsigma);
  this->set_LRAS(tmpLRAS); 
  Rcpp::Rcout << "++ change via $get_LRAS and $set_LRAS" << std::endl; 
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl; 
  
  this->nBMU = 2;
  this->BMU.set_size(0,2);
  
  // Print defaults message 
  Rcpp::Rcout << "Defaults set" << std::endl; 
  Rcpp::Rcout << "++ parallel = " << this->parallel << ", change via $set_parallel" << std::endl; 
  Rcpp::Rcout << "++ mtr_freq = " << this->mtr_freq << ", change via $set_monitoring_freq" << std::endl; 
  Rcpp::Rcout << "++ nBMU = " << this->nBMU << ", change via $set_nBMU" << std::endl; 
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl; 
  
}

inline void SOM::set_X_label(const Rcpp::CharacterVector& X_label_, Rcpp::DataFrame ctab_) {
  
  // Check that input labels have same dimensions as X 
  if(X_label_.size() != this->nX) Rcpp::stop("Input X_label must have size = nX");
  
  // Store the labels 
  this->X_label = X_label_; 
  
  // Set the color table 
  this->set_ctab(ctab_); 
}

inline void SOM::set_ctab(Rcpp::DataFrame ctab_) {
  
  // Check that input data frame is valid 
  Rcpp::Environment somdisco_env = Rcpp::Environment::namespace_env("SOMDisco");
  Rcpp::Function ctab_check = somdisco_env[".ctab_check"];
  bool ctab_pass = ctab_check(Rcpp::Named("ctab") = ctab_, Rcpp::Named("query_labels") = this->X_label);
  
  if(!ctab_pass) {
    Rcpp::stop("Input ctab invalid. Is it a data frame with columns 'label' and 'color' containing all unique labels in X_label?");
  }
  
  // Store the color table 
  this->ctab = ctab_; 
}


// ***** Scaling Functions 
inline void SOM::set_netrng(arma::rowvec ext_min, arma::rowvec ext_max, double int_min, double int_max) {
  
  // Set external network minimum
  if(ext_min.size()==1) {
    Rcpp::Rcout << "note: recycling ext_min across data dimension" << std::endl; 
    this->netrng_ext_min.set_size(this->d);
    this->netrng_ext_min.fill(ext_min[0]);
  } else if(ext_min.size()==this->d) {
    this->netrng_ext_min = ext_min; 
  } else {
    Rcpp::stop("ext_min must be a vector of length = 1 or d (data dimension)");
  }
  
  // Set external network maximum 
  if(ext_max.size()==1) {
    Rcpp::Rcout << "note: recycling ext_max across data dimension" << std::endl; 
    this->netrng_ext_max.set_size(this->d);
    this->netrng_ext_max.fill(ext_max[0]);
  } else if(ext_max.size()==this->d) {
    this->netrng_ext_max = ext_max; 
  } else {
    Rcpp::stop("ext_max must be a vector of length = 1 or d (data dimension)");
  }
  
  this->netrng_ext_rng = this->netrng_ext_max - this->netrng_ext_min;
  
  // Set internal network minimum & maximum 
  this->netrng_int_min = int_min; 
  this->netrng_int_max = int_max; 
  
  this->netrng_int_rng = this->netrng_int_max - this->netrng_int_min;
  
  this->is_netrng_set = true; 
}

inline arma::mat SOM::map_to_netrng(const arma::mat& x) {
  if(!this->is_netrng_set) Rcpp::stop("Must call set_netrng before calling map_to_netrng");
  
  // Check if rowvec
  if(x.n_rows == 1 && x.n_cols == this->d) {
    return (x.row(0) - this->netrng_ext_min) / this->netrng_ext_rng * this->netrng_int_rng + this->netrng_int_min;  
    // Check if colvec  
  } else if(x.n_rows == this->d && x.n_cols == 1) {
    return (x.col(0).t() - this->netrng_ext_min) / this->netrng_ext_rng * this->netrng_int_rng + this->netrng_int_min;  
    // Check if matrix, assuming data vectors in rows 
  } else if(x.n_rows > 1 && x.n_cols == this->d) {
    // Convert the "to" range to vectors 
    arma::rowvec to_min(1); to_min.fill(this->netrng_int_min);
    arma::rowvec to_max(1); to_max.fill(this->netrng_int_max);
    return SOM_UTILS_SCALING::linscale(x, this->netrng_ext_min, this->netrng_ext_max, to_min, to_max, true); // last is parallel = true 
  } else {
    Rcpp::stop("Input not a matrix with vectors to be scaled in rows");
  }
}

inline arma::mat SOM::map_from_netrng(const arma::mat& w) {
  if(!this->is_netrng_set) Rcpp::stop("Must call set_netrng before calling map_from_netrng");
  //return (w - this->netrng_int_min) / this->netrng_int_rng % this->netrng_ext_rng + this->netrng_ext_min;
  
  // Check if rowvec
  if(w.n_rows == 1 && w.n_cols == this->d) {
    return (w.row(0) - this->netrng_int_min) / this->netrng_int_rng % this->netrng_ext_rng + this->netrng_ext_min; 
    // Check if colvec  
  } else if(w.n_rows == this->d && w.n_cols == 1) {
    return (w.col(0).t() - this->netrng_int_min) / this->netrng_int_rng % this->netrng_ext_rng + this->netrng_ext_min; 
    // Check if matrix, assuming data vectors in rows 
  } else if(w.n_rows > 1 && w.n_cols == this->d) {
    // The "from" range is internal, convert to vectors 
    arma::rowvec from_min(1); from_min.fill(this->netrng_int_min);
    arma::rowvec from_max(1); from_max.fill(this->netrng_int_max);
    return SOM_UTILS_SCALING::linscale(w, from_min, from_max, this->netrng_ext_min, this->netrng_ext_max, true); // last is parallel = true 
  } else {
    Rcpp::stop("Input not a matrix with vectors to be scaled in rows");
  }
}


// ***** Prototype Initialization Functions 
inline void SOM::set_W_runif() {
  //inline void SOM::set_W_runif(unsigned int seed) {
  
  if(!this->is_netrng_set) Rcpp::stop("Must call $set_netrng before calling $set_W_runif");
  
  // Initialize the prototype matrix 
  //this->W.set_size(this->nW, this->d); 
  // if(seed > 0) {
  //   Rcpp::Environment base = Rcpp::Environment::namespace_env("base");
  //   Rcpp::Function setseed = base["set.seed"];
  //   setseed(Rcpp::Named("seed", seed));
  // }
  //this->W.randu();
  
  // Use Rcpp::Runif for reproducibility
  //int size_protos = int(this->nW * this->d);
  //Rcpp::NumericVector rcpp_runif = Rcpp::runif(size_protos, 0, 1); 
  //Rcpp::Rcout << rcpp_runif[0] << " " << rcpp_runif[1] << " " << rcpp_runif[2] << std::endl; 
  //this->W = arma::mat(rcpp_runif.begin(), this->nW, this->d); 
  this->W.set_size(this->nW, this->d); 
  this->W.randu();
  
  
  // // Map the initial protos from [0-1] to [internal_min,internal_max]
  // this->W = this->W * (this->netrng_int_max - this->netrng_int_min) + this->netrng_int_min;
  
  // Map the initial protos from [0-1] to 0.1*[internal_min,internal_max]
  //this->W = this->W * 0.1*(this->netrng_int_max - this->netrng_int_min) + 0.1*this->netrng_int_min;
  this->W = this->W * 0.1*(this->netrng_int_max - this->netrng_int_min) + (this->netrng_int_min + 0.5*this->netrng_int_rng);
  
  this->is_protos_init = true; 
}

inline void SOM::set_W(const arma::mat& W_) {
  // Make sure that the input W_ has correct dimensions
  if(this->W.n_rows != W_.n_rows || this->W.n_cols != W_.n_cols) 
    Rcpp::stop("Input must have nrows = som_x*som_y, ncols = d");
  
  // Make sure range of input W_ is inside network range 
  if(W_.max() > this->netrng_int_max || W_.min() < this->netrng_int_min) 
    Rcpp::warning("Input is outside netrng_int. Likely need to scale to internal network range, then call $set_W()");
  
  this->W = W_; 
  
  this->is_protos_init = true; 
}

inline void SOM::set_p_equal() {
  // Initialize winning proportions
  this->p.set_size(this->nW); 
  this->p.fill(1.0/double(this->nW));
  
  this->is_winfrq_init = true; 
}

inline void SOM::set_p(const arma::vec& p_) {
  if(p_.n_elem != this->nW) Rcpp::stop("Input must have length = nW");
  
  this->p = p_; 
  
  // Make sure they are all non-negative 
  if(arma::any(this->p < 0)) this->p += this->p.min(); 
  // Make sure they sum to 1 
  this->p /= arma::accu(this->p); 
  
  this->is_winfrq_init = true; 
}

inline arma::vec SOM::bias() {
  return this->gamma * (1.0/double(this->nW) - this->p); 
}


// ***** Learning Rate and Annealing Functions
inline void SOM::set_LRAS(Rcpp::DataFrame LRAS) {
  
  // Check that input data frame contains all necessary columns and extract them
  if(!LRAS.containsElementNamed("t")) Rcpp::stop("LRAS must contain column named 't'");
  if(!LRAS.containsElementNamed("alpha")) Rcpp::stop("LRAS must contain column named 'alpha'");
  if(!LRAS.containsElementNamed("beta")) Rcpp::stop("LRAS must contain column named 'beta'");
  if(!LRAS.containsElementNamed("gamma")) Rcpp::stop("LRAS must contain column named 'gamma'");
  if(!LRAS.containsElementNamed("sigma")) Rcpp::stop("LRAS must contain column named 'sigma'");
  
  Rcpp::IntegerVector tmp_t = LRAS["t"];
  Rcpp::NumericVector tmp_alpha = LRAS["alpha"];
  Rcpp::NumericVector tmp_beta = LRAS["beta"];
  Rcpp::NumericVector tmp_gamma = LRAS["gamma"];
  Rcpp::IntegerVector tmp_sigma = LRAS["sigma"];
  if(Rcpp::min(tmp_sigma) < 0) Rcpp::stop("all sigma must be >= 0");
  
  // Clear old LRAS
  this->LRAS.clear(); 
  
  // Initialize lrs map and fill it up
  for(unsigned int i=0; i<LRAS.nrows(); ++i) {
    this->LRAS[tmp_t[i]] = std::make_tuple(tmp_alpha[i], tmp_beta[i], tmp_gamma[i], tmp_sigma[i]);
  }
  
  // // Since the map is now filled & sorted,
  // // add a final entry to contain learning rates to t=infinity
  // type_LRAS::const_iterator it = --this->LRAS.end();
  // this->LRAS[std::numeric_limits<unsigned int>::max()] = it->second;
  
  // Print it out: 
  Rcpp::Rcout << "Storing the annealing schedule as:" << std::endl;
  Rcpp::Rcout << "t" << "\t" << "alpha" << "\t" << "beta" << "\t" << "gamma" << "\t" << "sigma" << std::endl;
  for(type_LRAS::const_iterator it = this->LRAS.begin(); it != this->LRAS.end(); it++) {
    Rcpp::Rcout << it->first << "\t" <<
      std::get<0>( it->second ) << "\t" <<
        std::get<1>( it->second ) << "\t" <<
          std::get<2>( it->second ) << "\t" <<   
            std::get<3>( it->second ) << std::endl;
  }
  
  // Update the internal learning rates 
  this->update_learning_rates();
  
  this->is_lras_set = true; 
}

inline Rcpp::DataFrame SOM::get_LRAS() {
  
  if(!this->is_lras_set) Rcpp::stop("Must call set_LRAS before calling get_LRAS");
  
  unsigned int nLRAS = this->LRAS.size();
  
  arma::uvec tmp_t(nLRAS);
  arma::vec tmp_alpha(nLRAS);
  arma::vec tmp_beta(nLRAS);
  arma::vec tmp_gamma(nLRAS);
  arma::uvec tmp_sigma(nLRAS);
  
  unsigned int counter = 0; 
  
  for(type_LRAS::const_iterator it = this->LRAS.begin(); it != this->LRAS.end(); it++) {
    tmp_t[counter] = it->first; 
    tmp_alpha[counter] = std::get<0>( it->second );
    tmp_beta[counter] = std::get<1>( it->second );
    tmp_gamma[counter] = std::get<2>( it->second );
    tmp_sigma[counter] = std::get<3>( it->second );
    counter++; 
  }
  
  return Rcpp::DataFrame::create(Rcpp::Named("t") = tmp_t, 
                                 Rcpp::Named("alpha") = tmp_alpha, 
                                 Rcpp::Named("beta") = tmp_beta, 
                                 Rcpp::Named("gamma") = tmp_gamma, 
                                 Rcpp::Named("sigma") = tmp_sigma);
}

inline arma::vec SOM::calc_eta(unsigned int sigma_) {
  // // Set the neighborhood update radius based on sigma 
  // arma::vec eta_(sigma_+1); 
  // eta_[0] = 1.0; // sigma is always >= 0
  // if(sigma_ >= 1) {
  //   eta_[1] = 1.0; 
  //   for(unsigned int s = 2; s <= sigma_; ++s) {
  //     eta_[s] = 1.0 - 1.0 / sigma_ * double(s-1);
  //   }
  // }
  
  // // Logistic decay 
  // arma::vec eta_(sigma_+1); 
  // for(unsigned int s=0; s<=sigma_; ++s) {
  //   eta_[s] = 1.0 / (1.0 + std::exp(12.0/double(sigma_+1) * double(s) - 6.0)); 
  // }
  
  
  double half_width = std::ceil( std::min(double(this->som_x), double(this->som_y)) / 2.0 ); 
  double decay_rate = 1.0 - std::pow(0.01, 1/(half_width-1)); // based on y=a(1-b)^m, where m=half width & y=0.01 set b
  arma::vec eta_(sigma_+1); 
  eta_[0] = 1.0; 
  for(unsigned int s=1; s<=sigma_; ++s) {
    eta_[s] = std::pow(1.0 - decay_rate, double(s-1)); 
  }
  
  return eta_; 
}

inline void SOM::update_learning_rates() {
  
  if(!this->is_lras_set) return;
  
  // Extract the current rates from the table based on the SOMs age 
  type_LRAS::const_iterator it = this->LRAS.lower_bound(this->age); 
  if(it == this->LRAS.end()) it--; 
  this->alpha = std::get<0>(it->second); 
  this->beta = std::get<1>(it->second); 
  this->gamma = std::get<2>(it->second); 
  this->sigma = std::get<3>(it->second); 
  
  
  // Set the neighborhood update radius based on sigma 
  this->eta = this->calc_eta(this->sigma); 
  // this->eta.set_size(this->sigma+1); 
  // this->eta[0] = 1.0; // sigma is always >= 0
  // if(this->sigma >= 1) {
  //   this->eta[1] = 1.0; 
  //   for(unsigned int s = 2; s <= this->sigma; ++s) {
  //     this->eta[s] = 1.0 - 1.0 / this->sigma * double(s-1);
  //   }
  // }
  
}



// ***** Monitoring Functions 
inline void SOM::set_monitoring_freq(unsigned int mtr_freq_) {
  this->mtr_freq = mtr_freq_; 
}

inline void SOM::monitor_SOM(const arma::mat& X) {
  
  // Perform a monitoring recall 
  this->mtrrecall_SOM(X); 
  
  // Increment the list of monitoring ages 
  this->mtr_age.push_back(age);
  
  // Compute the mean-squared quant error by proto
  this->mtr_RMSQE.push_back(std::sqrt(arma::mean(this->SQE)));
  
  // Compute & store the Quantization Stability Index
  // If this is the 1st monitoring step, just store NA for this value
  if(this->prevBMU.n_elem == 0) {
    this->mtr_QSI.push_back(NA_REAL);
  } else {
    arma::uvec BMUdiff = arma::find(arma::abs(this->prevBMU - this->BMU.col(0)));
    unsigned int ndiff = BMUdiff.n_elem;
    this->mtr_QSI.push_back(double(ndiff) / double(X.n_rows));
  }
  
  // Increment the Entropy monitoring 
  this->mtr_Entropy.push_back(this->Entropy); 
  
  // Visualize, if this is at least the second monitoring step recorded 
  if(this->mtr_age.size() >= 2) {
    // Obtaining namespace of SOMDisco package
    Rcpp::Environment pkg = Rcpp::Environment::namespace_env("SOMDisco");
    // Picking up vis_som_training function
    Rcpp::Function f = pkg["vis_som_training"];
    f(Rcpp::Named("SOM") = *this); 
  }

}



// ***** Training functions 
inline void SOM::set_reporting_freq(unsigned int report_freq_) {
  this->report_freq = report_freq_; 
}

inline void SOM::update_p(unsigned int winner_idx) {
  
  // All p are multiplied by this term 
  this->p *= (1.0 - this->beta);
  
  // Winners have this term added, according to their win strength (eta)
  // must subtract 1 from nhblist b/c those indices are stored 1-based 
  for(unsigned int k=0; k < this->eta.size(); ++k) {
    
    if(this->nu_nhblist[winner_idx][k].size() == 0) continue; 
    
    this->p.elem( this->nu_nhblist[winner_idx][k]-1 ) += this->beta * this->eta[k]; 
  }
  
  // Renormalize? 
  // this->p /= arma::accu(this->p);
  
}

inline void SOM::update_W(unsigned int winner_idx, const arma::rowvec& x) {
  
  // Update winner proto + its sigma neighbors, according to their win strength (eta)
  // must subtract 1 from nhblist b/c those indices are stored 1-based 
  for(unsigned int k=0; k < this->eta.size(); ++k) {
    
    if(this->nu_nhblist[winner_idx][k].size() == 0) continue; 
    
    double multiplier = this->alpha * this->eta[k]; 
    
    this->W.rows( this->nu_nhblist[winner_idx][k]-1 ) *= (1.0 - multiplier); 
    
    this->W.each_row( this->nu_nhblist[winner_idx][k]-1 ) += multiplier*x; 
  }
  
}

inline void SOM::train_SOM(unsigned int nsteps, const arma::mat& X) {
  
  if(!this->is_netrng_set) Rcpp::stop("Must call set_netrng() before calling train_SOM()");
  if(!this->is_protos_init) Rcpp::stop("Must initialize prototypes before calling train_SOM");
  if(!this->is_lras_set) Rcpp::stop("Must call set_LRAS() before calling train_SOM()");
  
  if(X.n_cols != this->d) Rcpp::stop("ncol(X) != data dimension specified during initialize_SOM()");
  //Rcpp::NumericVector inXstats = data_summary(X);
  Rcpp::NumericVector inXstats = Rcpp::NumericVector({X.min(), X.max(), X(0), X(X.n_elem-1)});
  inXstats.names() = Rcpp::CharacterVector({"min","max","X[1]","X[end]"});
  Rcpp::NumericVector diffXstats = Rcpp::abs(inXstats - this->X_stats);
  if(Rcpp::is_true(Rcpp::any(diffXstats > 1e-10))) Rcpp::stop("Stats of training set X do not match those computed during initialize_SOM()");
  
  // Set the training order
  if(this->user_train_order) {
    nsteps = this->train_order.size(); 
    Rcpp::Rcout << "User-specified training order detected. nsteps reset to " << nsteps << std::endl;
  } else {
    this->train_order = Rcpp::sample(X.n_rows, nsteps, true);  
  }
  this->train_order = this->train_order-1; 
  
  // Clear out any previous recalls
  if(this->BMU.n_elem > 0) this->prevBMU = this->BMU.col(0); 
  this->is_recalled = false; 
  this->clear_recall(); 
  
  // Begin training 
  //Rcpp::Rcout << "Training step (1000s), monitoring every " << this->mtr_freq << " steps:" << std::endl;
  Rcpp::Rcout << "SOM Training:" << std::endl;
  Rcpp::Rcout << "++ Reporting every " << this->report_freq << ", monitoring every " << this->mtr_freq << " steps" << std::endl;
  unsigned int report_counter = 0; // count the number of times reporting has been done to console, add newline after every 8 reports 
  
  for(unsigned int step = 1; step <= nsteps; ++step) {
    
    // Increment the age count and update annealing rates 
    this->age++;
    this->update_learning_rates();
    
    // Perform some checks every report_freq iterations 
    if(this->report_freq > 0 && step % this->report_freq == 0) {
      Rcpp::Rcout << step / this->report_freq << "\t"; 
      Rcpp::checkUserInterrupt();  // Check for user interrupt 
      report_counter++; 
    }
    
    //if(step % 10000 == 0 || step % nsteps==0) Rcpp::Rcout << std::endl;
    if(report_counter == 8) {
      Rcpp::Rcout << std::endl; 
      report_counter = 0; 
    }
    
    // Compute monitoring measures, if monitoring 
    if((this->mtr_freq > 0) && (this->age % this->mtr_freq == 0)) {
      
      // Calc & store monitoring quantities 
      this->is_trained = true; 
      this->monitor_SOM(X);  
    }
    
    // Scale the sample and find its BMU, without bias  
    arma::rowvec x = this->map_to_netrng(X.row(this->train_order[step-1]));
    SOM_UTILS_DIST::distmat_worker wkr(x, this->W, "L22"); 
    if(this->parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
    arma::mat xdist = wkr.DIST; 
    unsigned int winner_idx = xdist.index_min(); 
    
    
    // Update the win proportions 
    this->update_p(winner_idx); 
    
    
    // find new BMU, with bias 
    xdist -= this->bias().t();
    winner_idx = xdist.index_min(); 
    
    
    // Update prototypes 
    this->update_W(winner_idx, x); 
    
  } // close training loop 
  if(report_counter > 0) Rcpp::Rcout << std::endl; 
  Rcpp::Rcout << "End Training (current age = " << this->age << ")" << std::endl; 
  
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl; 
  
  this->is_trained = true; 
  
  this->recall_SOM(X); 
  
  this->user_train_order = false; 
  this->train_order = this->train_order+1; // undo the 0-based indexing above 
}

inline void SOM::set_train_order(const Rcpp::IntegerVector& train_order_) {
  this->train_order = train_order_; 
  this->user_train_order = true; 
}



// ***** Recall functions 
inline void SOM::set_nBMU(unsigned int k) {
  if(k < 2) Rcpp::stop("Cannot set nBMU < 2");
  
  this->nBMU = k; 
}

inline void SOM::calc_Entropy() {
  arma::vec RFp = arma::conv_to<arma::vec>::from(this->RF_size); 
  RFp /= arma::accu(RFp); 
  arma::vec logRFp = arma::log2(RFp); 
  logRFp.elem(arma::find(this->RF_size == 0)).zeros(); 
  this->Entropy = -1.0 * arma::accu(RFp % logRFp) / std::log2(double(this->nW));
}

inline void SOM::recall_SOM(const arma::mat& X) {
  
  bool verbose = true; 
  
  // Checks 
  //if(!this->is_trained) Rcpp::stop("Must call train_SOM before calling recall_SOM()");
  if(X.n_cols != this->d) Rcpp::stop("ncol(X) != data dimension specified during initialize_SOM()");
  
  // Examine input data, make sure it's the same as what was specified during initialization
  Rcpp::NumericVector inXstats = Rcpp::NumericVector({X.min(), X.max(), X(0), X(X.n_elem-1)});
  Rcpp::NumericVector diffXstats = Rcpp::abs(inXstats - this->X_stats);
  if(Rcpp::is_true(Rcpp::any(diffXstats > 1e-10))) Rcpp::stop("Stats of training set X do not match those computed during initialize_SOM()");
  
  if(verbose) {
    Rcpp::Rcout << "SOM Recall:" << std::endl; 
  }
  
  // Find BMU of data 
  if(this->BMU.n_elem > 0) this->prevBMU = this->BMU.col(0); 
  if(verbose) Rcpp::Rcout << "++ finding BMUs of data ... "; 
  SOM_UTILS_RECALL::find_BMU_worker BMUwkr(X, this->W, this->nBMU);
  BMUwkr.set_bias(this->bias()); 
  arma::rowvec net_min(1); net_min.fill(this->netrng_int_min); 
  arma::rowvec net_max(1); net_max.fill(this->netrng_int_max); 
  BMUwkr.set_ranges(this->netrng_ext_min, this->netrng_ext_max, net_min, net_max); 
  if(this->parallel) BMUwkr.calc_parallel(); else BMUwkr.calc_serial(); 
  this->BMU = BMUwkr.BMU + 1; // Store 1-indexed 
  this->SQE = BMUwkr.SQE; // Store quantization error 
  
  if(verbose) Rcpp::Rcout << "done" << std::endl; 
  
  // Build CADJ matrix 
  if(verbose) Rcpp::Rcout << "++ building CADJ matrix ... ";
  SOM_UTILS_RECALL::build_CADJ_worker CADJwkr(BMUwkr.BMU, this->nW, true); // last argument is zeroidx = true
  if(this->parallel) CADJwkr.calc_parallel(); else CADJwkr.calc_serial(); 
  this->CADJ = CADJwkr.CADJ; 
  if(verbose) Rcpp::Rcout << "done" << std::endl; 
  
  // Set the RF_Size = rowsums of CADJ 
  if(verbose) Rcpp::Rcout << "++ setting RF_size ... "; 
  this->RF_size = arma::sum(this->CADJ, 1); 
  if(verbose) Rcpp::Rcout << "done" << std::endl; 
  
  // Calc normalized entropy (log-2)
  if(verbose) Rcpp::Rcout << "++ calculating SOM Entropy ... "; 
  this->calc_Entropy(); 
  if(verbose) Rcpp::Rcout << "done" << std::endl; 
  
  // Set the RF_members 
  if(verbose) Rcpp::Rcout << "++ populating RF_members ... ";
  SOM_UTILS_RECALL::find_RF_members_worker RFmemwkr(this->BMU, this->nW, false); // use the 1-based BMUWs and set zeroidx = F, so member indices are returned as 1-indexed
  if(this->parallel) RFmemwkr.calc_parallel(); else RFmemwkr.calc_serial(); 
  this->RF_members = RFmemwkr.RF_members; 
  if(verbose) Rcpp::Rcout << "done" << std::endl; 
  
  // Set the lattice fences 
  if(verbose) Rcpp::Rcout << "++ setting lattice fences ... ";
  this->set_lattice_fences(); 
  if(verbose) Rcpp::Rcout << "done" << std::endl; 
  
  // Propagate the labels 
  if(verbose) Rcpp::Rcout << "++ propagating X_label to RFs ... "; 
  this->set_RF_label(); 
  if(verbose) Rcpp::Rcout << "done" << std::endl; 
  
  this->is_recalled = true; 
}

inline void SOM::mtrrecall_SOM(const arma::mat& X) {
  
  // Checks 
  if(!this->is_trained) Rcpp::stop("Must call train_SOM before calling mtrrecall_SOM()");
  if(X.n_cols != this->d) Rcpp::stop("ncol(X) != data dimension specified during initialize_SOM()");
  
  // Examine input data, make sure it's the same as what was specified during initialization
  Rcpp::NumericVector inXstats = Rcpp::NumericVector({X.min(), X.max(), X(0), X(X.n_elem-1)});
  Rcpp::NumericVector diffXstats = Rcpp::abs(inXstats - this->X_stats);
  if(Rcpp::is_true(Rcpp::any(diffXstats > 1e-10))) Rcpp::stop("Stats of training set X do not match those computed during initialize_SOM()");
  
  // Find BMU of data 
  if(this->BMU.n_elem > 0) this->prevBMU = this->BMU.col(0); 
  SOM_UTILS_RECALL::find_BMU_worker BMUwkr(X, this->W, this->nBMU);
  BMUwkr.set_bias(this->bias()); 
  arma::rowvec net_min(1); net_min.fill(this->netrng_int_min); 
  arma::rowvec net_max(1); net_max.fill(this->netrng_int_max); 
  BMUwkr.set_ranges(this->netrng_ext_min, this->netrng_ext_max, net_min, net_max); 
  if(this->parallel) BMUwkr.calc_parallel(); else BMUwkr.calc_serial(); 
  this->BMU = BMUwkr.BMU + 1; // Store 1-indexed 
  this->SQE = BMUwkr.SQE; // Store quantization error 
  

  // Build CADJ matrix 
  SOM_UTILS_RECALL::build_CADJ_worker CADJwkr(BMUwkr.BMU, this->nW, true); // last argument is zeroidx = true
  if(this->parallel) CADJwkr.calc_parallel(); else CADJwkr.calc_serial(); 
  this->CADJ = CADJwkr.CADJ; 

  // Set the RF_Size = rowsums of CADJ 
  this->RF_size = arma::sum(this->CADJ, 1); 

  // Calc normalized entropy (log-2)
  this->calc_Entropy(); 

  // Set the RF_members 
  SOM_UTILS_RECALL::find_RF_members_worker RFmemwkr(this->BMU, this->nW, false); // use the 1-based BMUWs and set zeroidx = F, so member indices are returned as 1-indexed
  if(this->parallel) RFmemwkr.calc_parallel(); else RFmemwkr.calc_serial(); 
  this->RF_members = RFmemwkr.RF_members; 
  
  // Set the lattice fences 
  this->set_lattice_fences(); 

}

inline void SOM::set_RF_label() {
  //if(!this->is_recalled) Rcpp::stop("Must call recall_SOM before calling set_RF_label");
  //if(X_label.size() != this->nX) Rcpp::stop("X_label must have length = nX");
  
  SOM_UTILS_RECALL::find_RF_label_worker wkr(this->X_label, this->RF_members, false); 
  if(this->parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
  
  this->RF_label = wkr.RF_label; 
  this->RF_label_dist = wkr.process_RF_labeldist(); 
}

inline void SOM::set_lattice_fences() {
  
  if(!this->is_protos_init) Rcpp::stop("Must initialize prototypes before calling set_lattice_fences");
  
  SOM_UTILS_RECALL::fence_worker fenceworker(this->nu_xy, this->W, this->nu_verts, this->nu_ADJ); 
  if(this->parallel) fenceworker.calc_parallel(); else fenceworker.calc_serial(); 
  this->fences = fenceworker.fence_dataframe(); 
}

inline arma::umat SOM::CONN() {
  return this->CADJ+this->CADJ.t(); 
}

inline void SOM::clear_recall() {
  this->BMU.set_size(0);
  this->CADJ.set_size(0,0);
  this->RF_size.set_size(0); 
  this->Entropy = 0.0; 
  this->RF_members.resize(0); 
  this->SQE.set_size(0); 
  this->RF_label = Rcpp::CharacterVector(0); 
  this->RF_label_dist.resize(0); 
  this->fences = Rcpp::DataFrame::create(); 
}



// ***** Saving & Loading 
inline void SOM::save(std::string rdsfile) {
  
  // Parse save file name 
  std::string rdsend = ".som"; 
  if(rdsend.size() > rdsfile.size() || !std::equal(rdsend.rbegin(), rdsend.rend(), rdsfile.rbegin()))
    Rcpp::stop("Save file string must end with .som");
  
  // Convert SOM object to list 
  Rcpp::List SOMList = this->as_list(); 
  
  // Save list as RDS 
  Rcpp::Environment base = Rcpp::Environment::namespace_env("base");
  Rcpp::Function saveRDS = base["saveRDS"];
  saveRDS(Rcpp::wrap(SOMList), Rcpp::Named("file", rdsfile));
  
  return; 
}

inline void SOM::load(std::string rdsfile) {
  
  // Load the RDS file as a list 
  Rcpp::Environment base = Rcpp::Environment::namespace_env("base");
  Rcpp::Function readRDS = base["readRDS"];
  Rcpp::List SOMList = readRDS(Rcpp::Named("file", rdsfile));
  
  // Populate the SOM object with list fields 
  this->load_list(SOMList); 
  
  return; 
  
}

inline void SOM::load_list(Rcpp::List SOMList) {
  
  
  this->parallel = SOMList["parallel"];
  
  // Variables involving lattice 
  this->som_x = SOMList["som_x"];
  this->som_y = SOMList["som_y"];
  this->lattice_type = Rcpp::as<std::string>(SOMList["lattice_type"]);
  this->set_lattice();
  
  // Variables involving training data 
  this->d = SOMList["d"];
  this->nX = SOMList["nX"];
  this->X_stats = SOMList["X_stats"];
  this->X_label = Rcpp::as<Rcpp::CharacterVector>(SOMList["X_label"]);
  this->ctab = Rcpp::as<Rcpp::DataFrame>(SOMList["ctab"]);
  
  // Variables involving network scaling 
  this->netrng_ext_min = Rcpp::as<arma::rowvec>(SOMList["netrng_ext_min"]);
  this->netrng_ext_max = Rcpp::as<arma::rowvec>(SOMList["netrng_ext_max"]);
  this->netrng_ext_rng = Rcpp::as<arma::rowvec>(SOMList["netrng_ext_rng"]);
  
  this->netrng_int_min = SOMList["netrng_int_min"]; 
  this->netrng_int_max = SOMList["netrng_int_max"]; 
  this->netrng_int_rng = SOMList["netrng_int_rng"];
  
  // Variables involving prototypes & updating
  this->W = Rcpp::as<arma::mat>(SOMList["W"]); 
  this->p = Rcpp::as<arma::vec>(SOMList["p"]);
  
  // Variables involing Learning Rate Annealing 
  Rcpp::DataFrame tmpLRAS = Rcpp::as<Rcpp::DataFrame>(SOMList["LRAS"]);
  this->set_LRAS(tmpLRAS); 
  this->update_learning_rates();
  
  // Variables involving monitoring
  this->mtr_freq = SOMList["mtr_freq"]; 
  this->mtr_age = Rcpp::as<std::vector<unsigned int>>(SOMList["mtr_age"]); 
  this->mtr_RMSQE = Rcpp::as<std::vector<double>>(SOMList["mtr_RMSQE"]); 
  this->mtr_QSI = Rcpp::as<std::vector<double>>(SOMList["mtr_QSI"]); // Quantization stability 
  this->mtr_Entropy = Rcpp::as<std::vector<double>>(SOMList["mtr_Entropy"]); // Entropy history
  
  // Variables involving training 
  this->age = SOMList["age"]; 
  this->train_order = Rcpp::as<Rcpp::IntegerVector>(SOMList["train_order"]);
  
  // Variables involving recall 
  this->nBMU = SOMList["nBMU"];
  this->BMU = Rcpp::as<arma::umat>(SOMList["BMU"]);
  this->CADJ = Rcpp::as<arma::umat>(SOMList["CADJ"]);
  this->RF_size = Rcpp::as<arma::uvec>(SOMList["RF_size"]);
  this->Entropy = SOMList["Entropy"];
  this->RF_members = Rcpp::as<std::vector<arma::uvec>>(SOMList["RF_members"]);
  this->SQE = Rcpp::as<arma::vec>(SOMList["SQE"]); 
  
  this->RF_label = Rcpp::as<Rcpp::CharacterVector>(SOMList["RF_label"]);
  this->RF_label_dist = Rcpp::as<std::vector<Rcpp::IntegerVector>>(SOMList["RF_label_dist"]);
  this->fences = Rcpp::as<Rcpp::DataFrame>(SOMList["fences"]);
  
  
  
  // Control flags 
  this->is_netrng_set = Rcpp::as<bool>(SOMList["is_netrng_set"]);
  this->is_protos_init = Rcpp::as<bool>(SOMList["is_protos_init"]);
  this->is_winfrq_init = Rcpp::as<bool>(SOMList["is_winfrq_init"]);
  this->is_lras_set = Rcpp::as<bool>(SOMList["is_lras_set"]); 
  this->is_trained = Rcpp::as<bool>(SOMList["is_trained"]); 
  this->is_recalled = Rcpp::as<bool>(SOMList["is_recalled"]); 
  
  
}

inline Rcpp::List SOM::as_list() {
  
  Rcpp::List SOMList;
  
  SOMList["parallel"] = this->parallel; 
  SOMList["som_x"] = this->som_x;
  SOMList["som_y"] = this->som_y;
  SOMList["nW"] = this->nW; 
  
  SOMList["lattice_type"] = this->lattice_type; 
  SOMList["nu_xy"] = this->nu_xy; 
  SOMList["nu_ij"] = this->nu_ij;
  SOMList["nu_verts"] = this->nu_verts; 
  SOMList["nu_ADJ"] = this->nu_ADJ; 
  SOMList["nu_nhblist"] = this->nu_nhblist;
  
  SOMList["d"] = this->d; 
  SOMList["nX"] = this->nX; 
  SOMList["X_stats"] = this->X_stats; 
  SOMList["X_label"] = this->X_label; 
  SOMList["ctab"] = this->ctab;
  
  SOMList["netrng_ext_min"] = this->netrng_ext_min; 
  SOMList["netrng_ext_max"] = this->netrng_ext_max; 
  SOMList["netrng_ext_rng"] = this->netrng_ext_rng; 
  SOMList["netrng_int_min"] = this->netrng_int_min; 
  SOMList["netrng_int_max"] = this->netrng_int_max; 
  SOMList["netrng_int_rng"] = this->netrng_int_rng; 
  
  SOMList["W"] = this->W; 
  SOMList["p"] = this->p; 
  
  SOMList["LRAS"] = this->get_LRAS();
  
  SOMList["alpha"] = this->alpha;
  SOMList["beta"] = this->beta; 
  SOMList["gamma"] = this->gamma; 
  SOMList["sigma"] = this->sigma; 
  SOMList["eta"] = this->eta; 
  
  SOMList["mtr_freq"] = this->mtr_freq; 
  SOMList["mtr_age"] = this->mtr_age; 
  SOMList["mtr_RMSQE"] = this->mtr_RMSQE; // Root mean square quantization error 
  SOMList["mtr_QSI"] = this->mtr_QSI; // Quantization stability 
  SOMList["mtr_Entropy"] = this->mtr_Entropy; // Entropy history 
  
  SOMList["age"] = this->age; 
  SOMList["train_order"] = this->train_order; 
  
  SOMList["nBMU"] = this->nBMU;
  SOMList["BMU"] = this->BMU;
  SOMList["CADJ"] = this->CADJ;
  SOMList["RF_size"] = this->RF_size; 
  SOMList["Entropy"] = this->Entropy; 
  SOMList["RF_members"] = this->RF_members; 
  SOMList["SQE"] = this->SQE; 
  
  SOMList["RF_label"] = this->RF_label; 
  SOMList["RF_label_dist"] = this->RF_label_dist; 
  
  SOMList["fences"] = this->fences; 
  
  SOMList["is_netrng_set"] = this->is_netrng_set;
  SOMList["is_protos_init"] = this->is_protos_init; 
  SOMList["is_winfrq_init"] = this->is_winfrq_init;
  SOMList["is_lras_set"] = this->is_lras_set; 
  SOMList["is_trained"] = this->is_trained; 
  SOMList["is_recalled"] = this->is_recalled; 
  
  return SOMList; 
}

// ***** Vis Parameters ***** 
inline void SOM::set_vis_par(Rcpp::List vis_par_) {
  this->vis_par = vis_par_; 
}

inline void SOM::set_vis_tile_bg(Rcpp::CharacterVector tile_bg) {
  if(tile_bg.size() != this->nW) Rcpp::stop("length(tile_bg) != nW");
  this->vis_tile_bg = tile_bg; 
} 

inline void SOM::set_vis_xlim(Rcpp::NumericVector xlim) {
  if(xlim.size() != 2) Rcpp::stop("length(vis_xlim) != 2");
  this->vis_xlim = xlim; 
}

inline void SOM::set_vis_ylim(Rcpp::NumericVector ylim) {
  if(ylim.size() != 2) Rcpp::stop("length(vis_ylim) != 2");
  this->vis_ylim = ylim; 
}






// RCPP_EXPOSED_CLASS(SOM)
//   
// RCPP_MODULE(som_module){
//   using namespace Rcpp; // Added (if not done globally)
//     
//   class_<SOM>("SOM")
//   
//   //.constructor<unsigned int, unsigned int, unsigned int, std::string>()
//   //.constructor<unsigned int, unsigned int, std::string>()
//   .constructor()
//   
//   // Fields involving the lattice 
//   .field_readonly("som_x", &SOM::som_x, "Width of SOM lattice")
//   .field_readonly("som_y", &SOM::som_y, "Height of SOM lattice")
//   .field_readonly("lattice_type", &SOM::lattice_type, "Lattice Type (grid/hex)")
//   .field_readonly("nu_xy", &SOM::nu_xy, "Neuron Lattice x-y Coordinates")
//   .field_readonly("nu_ij", &SOM::nu_ij, "Neuron Lattice i-j Coordinates")
//   .field_readonly("nu_verts", &SOM::nu_verts, "Lattice tile vertices")
//   .field_readonly("nu_ADJ", &SOM::nu_ADJ, "Neuron Lattice Adjacency Matrix")
//   .field_readonly("nu_nhblist", &SOM::nu_nhblist, "List of Neuron Lattice Neighbors by Distance")
//   .method("set_lattice", &SOM::set_lattice, "Set the SOM lattice parameters")
//   .method("tile_interior_point", &SOM::tile_interior_point, "Get the (x,y) coords of an interior tile point in a diretion of angle theta from its center")
//   
//   // Fields involving training data, dimension
//   .field_readonly("d", &SOM::d, "Data Dimension")
//   .field_readonly("nX", &SOM::nX, "Numer of training data vectors")
//   .field_readonly("Xstats", &SOM::Xstats, "Summary statistics identifying training data")
//   .method("initialize_SOM", &SOM::initialize_SOM, "Initialize the SOM and link to the training data matrix")
//   
//   // Fields inolving network scaling 
//   .field_readonly("netrng_ext_min", &SOM::netrng_ext_min, "External (data) min, by dimension")
//   .field_readonly("netrng_ext_max", &SOM::netrng_ext_max, "External (data) max, by dimension")
//   .field_readonly("netrng_ext_rng", &SOM::netrng_ext_rng, "External (data) range, by dimension")
//   .field_readonly("netrng_int_min", &SOM::netrng_int_min, "Internal (prototype) min")
//   .field_readonly("netrng_int_max", &SOM::netrng_int_max, "Internal (prototype) max")
//   .field_readonly("netrng_int_rng", &SOM::netrng_int_rng, "Internal (prototype) range")
//   //.method("set_network_range", &SOM::set_network_range, "Set the network range (external and internal)")
//   .method("set_netrng", &SOM::set_netrng, "Set the network range (external and internal)")
//   .method("map_to_netrng", &SOM::map_to_netrng, "Map a data vector to the network range (external to internal)")
//   .method("map_from_netrng", &SOM::map_from_netrng, "Map a prototype vector from the network range (internal to external)")
//   
//   // Fields involving prototype and their initialization 
//   .field_readonly("nW", &SOM::nW, "Number of Protos/Neurons")
//   .field_readonly("W", &SOM::W, "Prototype matrix")
//   .method("set_W_runif", &SOM::set_W_runif, "Set the prototype matrix to random uniform")
//   .method("set_W", &SOM::set_W, "Set the prototype matrix to specific input")
//   
//   .field_readonly("p", &SOM::p, "vector of winning proportions")
//   .method("set_p_equal", &SOM::set_p_equal, "Set the prototype winning frequencies to 1/nW")
//   .method("set_p", &SOM::set_p, "Set the prototype winning frequencies to specific value")
//   
//   .method("bias", &SOM::bias, "vector of prototype biases")
//   
//   // Variables involing Learning Rate Annealing 
//   .method("set_LRAS", &SOM::set_LRAS, "Set the Learning Rate Annealing Schedule")
//   .method("get_LRAS", &SOM::get_LRAS, "Get the Learning Rate Annealing Schedule")
//   
//   .field_readonly("alpha", &SOM::alpha, "Annealed prototype update multiplier")
//   .field_readonly("beta", &SOM::beta, "Annealed winfreq update multiplier")
//   .field_readonly("gamma", &SOM::gamma, "Annealed bias update multiplier")
//   .field_readonly("sigma", &SOM::sigma, "Annealed neighborhood update radius")
//   .field_readonly("eta", &SOM::eta, "Neighbor update strength")
//   .method("calc_eta", &SOM::calc_eta, "Get the neighborhood decay function, based on a max neighborhood size (sigma)")
//   .method("update_learning_rates", &SOM::update_learning_rates, "Update the current learning rates based on SOM age")
//   
//   // Functions to update the win proportions and prototypes 
//   .method("update_p", &SOM::update_p, "Update the win proportions during training given a winner index")
//   .method("update_W", &SOM::update_W, "Update the prototypes during training given a winner index")
//   
//   // Monitoring 
//   .field_readonly("mtr_freq", &SOM::mtr_freq, "Incremental age at which monitoring measures are calculated") 
//   .field_readonly("mtr_age", &SOM::mtr_age, "A list of the training ages at which monitoring occurred")
//   .field_readonly("mtr_RMSQE", &SOM::mtr_RMSQE, "Root Mean Square Quantization Error, for each proto and overall, at monitoring ages")
//   .field_readonly("mtr_QSI", &SOM::mtr_QSI, "Quantization Stability at monitoring ages")
//   .method("set_monitoring_freq", &SOM::set_monitoring_freq, "Set the incremental age used for monitoring")
//   
//   // Variables involving training 
//   .field_readonly("age", &SOM::age, "The number of training steps performed thus far")
//   .method("train_SOM", &SOM::train_SOM, "Train the SOM")
//   
//   // Variables involving recall 
//   .field_readonly("nBMU", &SOM::nBMU, "The number of BMUs calculated during recall")
//   .method("set_nBMU", &SOM::set_nBMU, "Set the BMU search depth")
//   .field_readonly("BMU", &SOM::BMU, "Matrix of data BMUs")
//   .field_readonly("CADJ", &SOM::CADJ, "CADJ Matrix")
//   .field_readonly("RF_size", &SOM::RF_size, "Size of each prototype's receptive field")
//   .field_readonly("RF_members", &SOM::RF_members, "Data members of each prototype's receptive field")
//   .field_readonly("SQE", &SOM::SQE, "Squared Quantization Error for each data vector")
//   .method("recall_SOM", &SOM::recall_SOM, "Recall the SOM")
//   
//   .field_readonly("RF_label", &SOM::RF_label, "Plurality label of data members of each prototype's receptive field")
//   .method("set_RF_label", &SOM::set_RF_label, "Sets the RF plurality label, given a vector of data labels")
//   
//   .field_readonly("fences", &SOM::fences, "Data frame containing fence information")
//   .method("set_lattice_fences", &SOM::set_lattice_fences, "Sets lattice fence values")
//   
//   .method("CONN", &SOM::CONN, "Returns the CONN matrix")
// 
//   // Control flags   
//   .field_readonly("is_netrng_set", &SOM::is_netrng_set, "Flag whether set_network_range has been called")
//   .field_readonly("is_protos_init", &SOM::is_protos_init, "Flag whether prototypes have been initialized")
//   .field_readonly("is_winfrq_init", &SOM::is_winfrq_init, "Flag whether prototype win frequencies have been initialized")
//   .field_readonly("is_lras_set", &SOM::is_lras_set, "Flag whether set_LRAS has been called")
//   .field_readonly("is_trained", &SOM::is_trained, "Flag whether train_SOM has been called")
//   .field_readonly("is_recalled", &SOM::is_recalled, "Flag whether recall_SOM has been called")
//   
//   .method("save", &SOM::save, "Save a SOM object")
//   .method("load", &SOM::load, "Load a SOM object")
//   ;
// }


#endif

