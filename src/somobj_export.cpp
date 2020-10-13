#ifndef SOMDISCO_SOMOBJ_HPP
#include "SOMDisco_types.hpp"
#endif



//RCPP_EXPOSED_CLASS(SOM)

RCPP_MODULE(som_module){
  using namespace Rcpp; // Added (if not done globally)
  
  class_<SOM>("SOM")
    
    //.constructor<unsigned int, unsigned int, unsigned int, std::string>()
    //.constructor<unsigned int, unsigned int, std::string>()
    .constructor()
    
    .field_readonly("parallel", &SOM::parallel, "Whether SOM training and associated computations are performed in parallel")
    .method("set_parallel", &SOM::set_parallel, "Set the parallel flag for the SOM object")
    
    // Fields involving the lattice 
    .field_readonly("som_x", &SOM::som_x, "Width of SOM lattice")
    .field_readonly("som_y", &SOM::som_y, "Height of SOM lattice")
    .field_readonly("lattice_type", &SOM::lattice_type, "Lattice Type (grid/hex)")
    .field_readonly("nu_xy", &SOM::nu_xy, "Neuron Lattice x-y Coordinates")
    .field_readonly("nu_ij", &SOM::nu_ij, "Neuron Lattice i-j Coordinates")
    .field_readonly("nu_verts", &SOM::nu_verts, "Lattice tile vertices")
    .field_readonly("nu_ADJ", &SOM::nu_ADJ, "Neuron Lattice Adjacency Matrix")
    .field_readonly("nu_nhblist", &SOM::nu_nhblist, "List of Neuron Lattice Neighbors by Distance")
    .method("set_lattice", &SOM::set_lattice, "Set the SOM lattice parameters")
    .method("tile_interior_point", &SOM::tile_interior_point, "Get the (x,y) coords of an interior tile point in a diretion of angle theta from its center")
    
    // Fields involving training data, dimension
    .field_readonly("d", &SOM::d, "Data Dimension")
    .field_readonly("nX", &SOM::nX, "Numer of training data vectors")
    .field_readonly("X_label", &SOM::X_label, "Vector of data labels")
    .field_readonly("ctab", &SOM::ctab, "Color table mapping data labels to colors, used for visualizations")
    .field_readonly("X_stats", &SOM::X_stats, "Summary statistics identifying training data")
    .method("initialize_SOM", &SOM::initialize_SOM, "Initialize the SOM and link to the training data matrix")
    .method("set_X_label", &SOM::set_X_label, "Set the data labels and their corresponding color table")
    .method("set_ctab", &SOM::set_ctab, "Set the color table for the X_labels")
    
    // Fields inolving network scaling 
    .field_readonly("netrng_ext_min", &SOM::netrng_ext_min, "External (data) min, by dimension")
    .field_readonly("netrng_ext_max", &SOM::netrng_ext_max, "External (data) max, by dimension")
    .field_readonly("netrng_ext_rng", &SOM::netrng_ext_rng, "External (data) range, by dimension")
    .field_readonly("netrng_int_min", &SOM::netrng_int_min, "Internal (prototype) min")
    .field_readonly("netrng_int_max", &SOM::netrng_int_max, "Internal (prototype) max")
    .field_readonly("netrng_int_rng", &SOM::netrng_int_rng, "Internal (prototype) range")
    //.method("set_network_range", &SOM::set_network_range, "Set the network range (external and internal)")
    .method("set_netrng", &SOM::set_netrng, "Set the network range (external and internal)")
    .method("map_to_netrng", &SOM::map_to_netrng, "Map a data vector to the network range (external to internal)")
    .method("map_from_netrng", &SOM::map_from_netrng, "Map a prototype vector from the network range (internal to external)")
    
    // Fields involving prototype and their initialization 
    .field_readonly("nW", &SOM::nW, "Number of Protos/Neurons")
    .field_readonly("W", &SOM::W, "Prototype matrix")
    .method("set_W_runif", &SOM::set_W_runif, "Set the prototype matrix to random uniform")
    .method("set_W", &SOM::set_W, "Set the prototype matrix to specific input")
    
    .field_readonly("p", &SOM::p, "vector of winning proportions")
    .method("set_p_equal", &SOM::set_p_equal, "Set the prototype winning frequencies to 1/nW")
    .method("set_p", &SOM::set_p, "Set the prototype winning frequencies to specific value")
    
    .method("bias", &SOM::bias, "vector of prototype biases")
    
    // Variables involing Learning Rate Annealing 
    .method("set_LRAS", &SOM::set_LRAS, "Set the Learning Rate Annealing Schedule")
    .method("get_LRAS", &SOM::get_LRAS, "Get the Learning Rate Annealing Schedule")
    
    .field_readonly("alpha", &SOM::alpha, "Annealed prototype update multiplier")
    .field_readonly("beta", &SOM::beta, "Annealed winfreq update multiplier")
    .field_readonly("gamma", &SOM::gamma, "Annealed bias update multiplier")
    .field_readonly("sigma", &SOM::sigma, "Annealed neighborhood update radius")
    .field_readonly("eta", &SOM::eta, "Neighbor update strength")
    .method("calc_eta", &SOM::calc_eta, "Get the neighborhood decay function, based on a max neighborhood size (sigma)")
    .method("update_learning_rates", &SOM::update_learning_rates, "Update the current learning rates based on SOM age")
    
    // Functions to update the win proportions and prototypes 
    .method("update_p", &SOM::update_p, "Update the win proportions during training given a winner index")
    .method("update_W", &SOM::update_W, "Update the prototypes during training given a winner index")
    
    // Monitoring 
    .field_readonly("mtr_freq", &SOM::mtr_freq, "Incremental age at which monitoring measures are calculated") 
    .field_readonly("mtr_age", &SOM::mtr_age, "A list of the training ages at which monitoring occurred")
    .field_readonly("mtr_RMSQE", &SOM::mtr_RMSQE, "Root Mean Square Quantization Error, for each proto and overall, at monitoring ages")
    .field_readonly("mtr_QSI", &SOM::mtr_QSI, "Quantization Stability at monitoring ages")
    .field_readonly("mtr_Entropy", &SOM::mtr_Entropy, "Entropy monitoring history")
    .method("set_monitoring_freq", &SOM::set_monitoring_freq, "Set the incremental age used for monitoring")
    
    // Variables involving training 
    .field_readonly("age", &SOM::age, "The number of training steps performed thus far")
    .field_readonly("train_order", &SOM::train_order, "The order of data vectors picked during training")
    .field_readonly("user_train_order", &SOM::user_train_order, "Whether a user-specified training order exists")
    .field_readonly("report_freq", &SOM::report_freq, "Frequency of reporting during learning")
    .method("train_SOM", &SOM::train_SOM, "Train the SOM")
    .method("set_train_order", &SOM::set_train_order, "Set the order of data drawn sequentially for SOM training.")
    .method("set_reporting_freq", &SOM::set_reporting_freq, "Set the training reporting frequency")
    
    // Variables involving recall 
    .field_readonly("nBMU", &SOM::nBMU, "The number of BMUs calculated during recall")
    .method("set_nBMU", &SOM::set_nBMU, "Set the BMU search depth")
    .field_readonly("BMU", &SOM::BMU, "Matrix of data BMUs")
    .field_readonly("CADJ", &SOM::CADJ, "CADJ Matrix")
    .field_readonly("RF_size", &SOM::RF_size, "Size of each prototype's receptive field")
    .field_readonly("Entropy", &SOM::Entropy, "Normalized Entropy of the SOM mapping")
    .field_readonly("RF_members", &SOM::RF_members, "Data members of each prototype's receptive field")
    .field_readonly("SQE", &SOM::SQE, "Squared Quantization Error for each data vector")
    .method("recall_SOM", &SOM::recall_SOM, "Recall the SOM")
    
    .field_readonly("RF_label", &SOM::RF_label, "Plurality label of data members of each prototype's receptive field")
    .field_readonly("RF_label_dist", &SOM::RF_label_dist, "Distribution of labels in each prototype's receptive field")
    .method("set_RF_label", &SOM::set_RF_label, "Sets the RF plurality label, given a vector of data labels")
    
    .field_readonly("fences", &SOM::fences, "Data frame containing fence information")
    .method("set_lattice_fences", &SOM::set_lattice_fences, "Sets lattice fence values")
    
    .method("CONN", &SOM::CONN, "Returns the CONN matrix")
    
    // Control flags   
    .field_readonly("is_netrng_set", &SOM::is_netrng_set, "Flag whether set_network_range has been called")
    .field_readonly("is_protos_init", &SOM::is_protos_init, "Flag whether prototypes have been initialized")
    .field_readonly("is_winfrq_init", &SOM::is_winfrq_init, "Flag whether prototype win frequencies have been initialized")
    .field_readonly("is_lras_set", &SOM::is_lras_set, "Flag whether set_LRAS has been called")
    .field_readonly("is_trained", &SOM::is_trained, "Flag whether train_SOM has been called")
    .field_readonly("is_recalled", &SOM::is_recalled, "Flag whether recall_SOM has been called")
    
    // IO 
    .method("save", &SOM::save, "Save a SOM object")
    .method("load", &SOM::load, "Load a SOM object")
    .method("as_list", &SOM::as_list, "Convert a SOM object to an R list")
    .method("load_list", &SOM::load_list, "Populate a SOM object from an R list")
    
    // Vis Parameters
    .field_readonly("vis_par", &SOM::vis_par, "Plot par values")
    .field_readonly("vis_tile_bg", &SOM::vis_tile_bg, "BG color of plotted SOM tiles")
    .field_readonly("vis_xlim", &SOM::vis_xlim, "Xlim of plot window")
    .field_readonly("vis_ylim", &SOM::vis_ylim, "Ylim of plot window")
    
    .method("set_vis_par", &SOM::set_vis_par, "Set the plot par values")
    .method("set_vis_tile_bg", &SOM::set_vis_tile_bg, "Set the BG color of plotted SOM tiles")
    .method("set_vis_xlim", &SOM::set_vis_xlim, "Set the xlim of the plot window")
    .method("set_vis_ylim", &SOM::set_vis_ylim, "Set the ylim of the plot window")

    
    
    ;
}
