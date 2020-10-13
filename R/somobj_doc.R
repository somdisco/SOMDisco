#' @name SOM
#' @title The SOM Object Class
#' 
#' @field age The current age of the SOM (the number of training steps that have been performed thus far).
#' @field alpha The strength of the prototype update, as set in the LRAS, given the current age of the SOM.
#' @field beta The strength of the win frequency update, as set in the LRAS, given the current age of the SOM.
#' @field bias Vector of CSOM biases for each prototype.
#' @field BMU Matrix of BMUs (the best matching prototype indices) for each training data vector. nrow = nX, ncol = nBMU. 
#' @field CADJ The Cumulative Adjacency Matrix of SOM prototypes. 
#' @field ctab The color table which maps unique labels found in \code{X_label} to colors. Must be a data frame with columns 'label' and 'color'. 
#' @field d The training data (and prototype) dimension
#' @field Entropy The normalized entropy of the learned SOM quantization 
#' @field eta Vector of coefficients applied to cooperate lattice neighbor updates during training. 
#' The elements of eta exhibit logistic decay as lattice distance between neurons increases 
#' (and cease to effect any prototype updates for neurons beyond distance \code{sigma} from a BMU). 
#' @field fences A dataframe storing the information required to visualize the U-Matrix fences on the SOM lattice. 
#' @field gamma The strength of the win frequencies on the prototype biases, as set in the LRAS, given the current age of the SOM.
#' @field is_netrng_set Flag indicating whether \code{set_netrng} has been called. 
#' @field is_protos_init Flag indicating whether \code{set_W_runif} or \code{set_W} has been called. 
#' @field is_recalled Flag indicating whether \code{recall_SOM} has been called since last call to \code{train_SOM}. 
#' @field is_trained Flag indicating whether \code{train_SOM} has been called since prototypes have been initialized.
#' @field is_winfrq_init Flag indicating whether \code{set_p_equal} or \code{set_p} has been called. 
#' @field lattice_type String defining the SOM lattice topology, either "hex" or "grid". 
#' @field LRAS The Learning Rate Annealing Schedule, as a data frame with columns \code{t}, \code{alpha}, \code{beta}, \code{gamma}, \code{sigma}. See \link{default_LRAS} for an example. 
#' @field mtr_age A vector of SOM ages at which monitoring snapshots were taken during training. 
#' @field mtr_Entropy A vector of SOM (normailized) entropies at the monitoring snapshots. 
#' @field mtr_freq The incremental number of steps at which monitoring snapshots are taken during training. Set = 0 to disable monitoring. 
#' @field mtr_QSI A vector of the Quantization Stability Indices computed at each monitoring snapshot. 
#' QSI is the proportion of data vectors whose BMUs have changed since the last monitoring snapshot was taken. 
#' @field mtr_RMSQE A matrix of Root Mean Square Quantization Errors at the prototype level (columns) at each monitoring snapshot (rows). 
#' This matrix has \code{nW+1} columns, were the last contains the global RMSQE (average over all prototype-level RMSQEs). 
#' @field nBMU The number of BMUs recorded during \code{recall_SOM}. The first BMU for a training data vector is its closest SOM prototype, the 2nd BMU is the next-closest prototype, etc. 
#' @field netrng_ext_max The max of the training data range, used for network scaling. 
#' @field netrng_ext_min The min of the training data range, used for network scaling. 
#' @field netrng_ext_rng The entire range (max - min) of the external data range. 
#' @field netrng_int_max The max of the internal (network) range, used for network scaling. 
#' @field netrng_int_min The min of the internal (network) range, used for network scaling. 
#' @field netrng_int_rng The entire range (max - min) of the internal range. 
#' @field nu_ADJ Adjacency matrix of lattice neurons, as dictated by the \code{lattice_type}. 
#' @field nu_ij Matrix of \code{(i=row,j=col)} coordinates of neurons on the SOM lattice. 
#' @field nu_nhblist A list of list describing shortest-path neuron distances according to the SOM lattice topology set in \code{lattice_type}. 
#' \code{nu_nhblist[[i]][[j]]} contains all neuron indices that are lattice distance \code{j} from neuron \code{i}. 
#' @field nu_verts Cube whose slices contain the \code{(x,y)} coordinates of each lattice tile. 
#' @field nu_xy Matrix of \code{(x,y)} coordinates of neurons on the SOM lattice, primarily used for visualizations. 
#' @field nW Number of SOM prototypes (and neurons), \code{= som_x * som_y}. 
#' @field nX Number of training data vectors (nrows of the training data matrix \code{X}). 
#' @field p Vector of the win frequencies of each prototype, as defined in the Conscience-SOM algorithm. 
#' @field parallel Flag indicating whether computations should be performed in parallel, using the RcppParallel package. 
#' @field RF_label Character vector of prototype labels, derived from the labels of training data mapped to each prototype. Only valid if the training data possesses labels. 
#' @field RF_members A list of training data indices mapped to each prototype during \code{recall_SOM}. 
#' @field RF_size A vector storing the sizes of each prototype's Receptive Field (the number of training vectors mapped to each prototype during \code{recall_SOM}).
#' @field sigma The maximum lattice distance between a BMU and its lattice neighborhood to which cooperative SOM prototype updates are applied. 
#' @field som_x Width of SOM lattice, in neurons
#' @field som_y Height of SOM lattice, in neurons
#' @field SQE Vector of the Squared Quantization Error resulting from quantizing each training vector by its BMU (length = \code{nW}). 
#' @field train_order Vector of indices of training vectors recording the order in which data were drawn to be presented to the network during training. 
#' @field user_train_order Flag indicating whether a user-specified training order has been set prior to training. 
#' @field W The prototype matrix (nrow = \code{nW}, ncol = \code{d}) storing SOM prototypes in its rows. 
#' The rows are in neuron order, so row1 corresponds to the \code{(i=1,j=1)} lattice neuron, row2 corresponds to the \code{(i=1,j=2)} lattice neuron, and so on. 
#' @field X_label A vector (length = \code{nX}) containing character labels of each training observation. 
#' @field X_stats A vector storing simple statistics (min, max, value of 1st element, value of last element) of the entire training data \code{X} that was specified 
#' during \code{initialize_SOM}. Subsequent methods which require the training data matrix (both \code{train_SOM} and \code{recall_SOM}) check 
#' that their inputs match those computed and stored in \code{X_stats} to ensure the same training data is given to all SOM functions. 
#' @field vis_par Stores the \code{par()} settings for base R graphics required to layer the \code{vis_som_*} functions atop each other. 
#' Usually this field is not needed by the user. 
#' @field vis_tile_bg Stores the current coloring of the lattice tiles, which sets the text annotation color if requested as "auto". 
#' Usually this field is not needed by the user. 
#' @field vis_xlim Stores the \code{xlim} of the SOM plotting window 
#' Usually this field is not needed by the user. 
#' @field vis_ylim Stores the \code{ylim} of the SOM plotting window. 
#' Usually this field is not needed by the user. 
#' 
#' @section Methods:
#' Each class method has its own documentation, accessible via \code{?SOMDisco::<method_name>}. 
#' For completeness, the list is repeated here in entirety.  Additional functionality 
#' for visualizing a trained SOM object is available through the \code{vis_som_*} functions. 
#' See their documentation for more information. 
#' 
#' \describe{
#' \item{\code{calc_eta}}{Calculate the coefficients governing cooperative updates of the lattice neighbors of a BMU.}
#' \item{\code{CONN}}{The symmetric Cumulative Adjacency Matrix of SOM prototypes.}
#' \item{\code{calc_Entropy}}{Compute the normalized entropy of the SOM quantization.}
#' \item{\code{get_LRAS}}{Get the Learning Rate Annealing Schedule.}
#' \item{\code{initialize_SOM}}{Setup a SOM object for training.}
#' \item{\code{new}}{Instantiate a new SOM object.}
#' \item{\code{map_from_netrng}}{Scale the rows of a matrix from the internal network range, to the external network range.}
#' \item{\code{map_to_netrng}}{Scale the rows of a matrix from the external network range, to the internal network range.}
#' \item{\code{recall_SOM}}{Recall a trained SOM object.}
#' \item{\code{set_ctab}}{Set the color table mapping unique data labels to colors.}
#' \item{\code{set_lattice}}{Set the SOM lattice quantities.}
#' \item{\code{set_lattice_fences}}{Compute and store the information needed to visualize U-Matrix fences on the SOM lattice.}
#' \item{\code{set_LRAS}}{Set the Learning Rate Annealing Schedule.}
#' \item{\code{set_monitoring_freq}}{(De-)Activate monitoring of SOM training.}
#' \item{\code{set_nBMU}}{Set the number of BMUs recorded during \code{recall_SOM}.}
#' \item{\code{set_netrng}}{Set the ranges used for linear network scaling.}
#' \item{\code{set_p}}{Set the CSOM win frequencies to user-specified values.}
#' \item{\code{set_p_equal}}{Set the CSOM win frequencies to the equiprobable value 1/\code{nW}.}
#' \item{\code{set_parallel}}{Set the parallel computation flag.}
#' \item{\code{set_RF_label}}{Set the prototype label based on labeled training data.}
#' \item{\code{set_train_order}}{Set the order in which training data are presented to the network in \code{train_SOM}.}
#' \item{\code{set_W}}{Set the SOM prototypes to user-specified values.}
#' \item{\code{set_W_runif}}{Initialize the SOM prototypes to uniform random values.}
#' \item{\code{set_X_label}}{Set the training data labels.}
#' \item{\code{tile_interior_point}}{Get a list of interior points within each SOM lattice tile. Useful mostly for visualization routines.}
#' \item{\code{train_SOM}}{Train the SOM object.}
#' \item{\code{update_learning_rates}}{Update the learning rates based on the current SOM age.}
#' \item{\code{update_p}}{Update the CSOM win frequencies during training.}
#' \item{\code{update_W}}{Update the prototypes during training.}
#' 
#' \item{\code{save}}{Save a trained SOM object to disk.}
#' \item{\code{load}}{Load a saved SOM object from disk into a current R environment.}
#' \item{\code{as_list}}{Convert and return all fields of a SOM object to an R list.}
#' \item{\code{load_list}}{Populate an instance of a SOM object from an R list.}
#' }
NULL 


# ***** Constructor *****

#' @name new
#' @title Create an empty SOM object
#' @return an empty SOM object templated for initialization (via \link{initialize_SOM}) and training (via \link{train_SOM}). 
#' @usage SOM$new()
NULL 


# ***** Control *****

#' @name set_parallel
#' @title Setter function for parallel computation.  
#' @description Sets an internal flag controlling whether computations are performed in parallel. 
#' @param parallel either TRUE or FALSE, as desired
#' @details 
#' Parallel computation is supported via RcppParallel, 
#' see \link[RcppParallel]{setThreadOptions} for details to control threading. 
#' By default, \code{parallel = TRUE} at SOMobj instantiation. 
#' @return None, \code{parallel} field is set internally 
#' @usage SOMobj$set_parallel(parallel)
NULL


# ***** Lattice *****

#' @name set_lattice
#' @title Setter function for the SOM lattice. 
#' @description Sets the neuron lattice coordinates and associated quantities needed for SOM training and visualization. 
#' @details \link{initialize_SOM} must be called prior to calling \code{set_lattice}. 
#' 
#' Internally, the following fields are computed and set: 
#' \itemize{
#' \item \code{nu_xy} the (x,y) coordinates of neurons on the lattice
#' \item \code{nu_ij} the (i=row,j=col) coordinates of neurons on the lattice
#' \item \code{nW} the number of neurons / prototypes
#' \item \code{nu_ADJ} the adjacency matrix of neurons on the lattice (identifying immediate neighbors)
#' \item \code{nu_nhblist} an organized hierarchical list of neurons sorted by their lattice distance from each other. 
#' \code{nu_nhblist[[i]][[j]]} contains a vector of neuron IDs that are distance j from neuron i on the lattice. 
#' \item \code{nu_verts} cube whose slices contain the vertices of each lattice tile. 
#' Slice i contains the vertices of lattice tile i in its rows, in vertex order.  
#' }
#' @return None, lattice related fields are set internally
#' @usage SOMobj$set_lattice()
NULL

#' @name tile_interior_point
#' @title Get a list of interior points in each lattice tile 
#' @description Handy for visualization, this function returns the (x,y) coordinates of a point inside each 
#' lattice tile that is along a ray of angle \code{theta} emanating from the tile center. 
#' The "radius" of a tile is the distance to the tile boundary in direction \code{theta}. 
#' The (x,y) coordinates are of distance = \code{rprop} x radius from the tile center. 
#' @param theta the angle at which the (x,y) coordinates should be placed, relative to the tile center 
#' @param rprop the proportion of the tile radius along angle \code{theta} at which the (x,y) coordinates will be placed. 
#' Should be < 1 to produce a point inside the lattice tile.  
#' @return a matrix (nrow = nW, ncol = 2) containing the (x,y) coordinates of interior points for each tile. 
#' The rows are in neuron-order. 
#' @usage SOMobj$tile_interior_point(theta, rprop)
NULL


# ***** Initialization *****

#' @name initialize_SOM
#' @title Initialize a SOM object
#' @description Sets up a SOM lattice for training, populates default training parameters, and attaches the training data to the specified SOM. 
#' @param X matrix of training data (data vectors in rows)
#' @param som_x width of the SOM (number of neurons)
#' @param som_y height of the SOM (number of neurons)
#' @param lattice_type either "grid" for rectangular lattices or "hex" for hexagonal lattices
#' @details This is a wrapper function to perform all steps necessary to setup a SOM for learning. 
#' Internally it calls \link{set_lattice}, \link{set_netrng}, \link{set_W_runif}, \link{set_p_equal}, and \link{set_LRAS}, all with their default arguments.  
#' It also sets internal variables: 
#' \itemize{
#' \item \code{som_x}, \code{som_y}, \code{nW}, \code{lattice_type}
#' \item \code{d}, \code{nX}, \code{X_stats}
#' \item \code{nBMU}
#' }
#' @return None, data and lattice related fields are set internally
#' @usage SOMobj$initialize_SOM(X, som_x, som_y, lattice_type)
NULL

#' @name set_X_label
#' @title Set the training data labels 
#' @description If the training data possess labels, this method stores them in the \code{X_label} field, which is used 
#' to propagate the \code{RF_label} and \code{RF_label_dist} fields during \code{recall_SOM}.  
#' @param X_label a vector (length = \code{nW}) containing character labels for each training observation. 
#' @param ctab a color table defining the mapping between the unique labels found in \code{X_label} and distinct colors. 
#' Must be a data frame with columns 'label' and 'color' (where 'color' is in HEX format). 
#' @details During \code{initialize_SOM} the data labels are set to a default value of "?". This method allows changing these defaults. 
#' The function \code{build_ctab} can be used to construct a default color table, if needed. 
#' @return None, field \code{X_label} is set internally. 
#' @usage SOMobj$set_X_label(X_label, ctab)
NULL 

#' @name set_ctab
#' @title Set the color table 
#' @description The color table controls the mapping of the unique labels found in \code{X_label} to distinct colors. 
#' @param ctab a color table defining the mapping between the unique labels found in \code{X_label} and distinct colors. 
#' Must be a data frame with columns 'label' and 'color' (where 'color' is in HEX format). 
#' @details \code{ctab} must contain color mappings for all unique labels found in \code{X_label}. A check will 
#' be performed internally. 
#' @return None, field \code{ctab} is set internally
#' @usage SOMobj$set_ctab(ctab)
NULL 


# ***** Scaling *****

#' @name set_netrng
#' @title Set the external/internal network range
#' @description Training (i.e., prototype updating) is performed in an "internal" network range. 
#' When presented to the network for training, each data vector is scaled linearly from the external 
#' data range to the internal data range. This function sets the min/max of both ranges. 
#' @param ext_min min of the external data range (the training data min)
#' @param ext_max max of the external data range (the training data max)
#' @param int_min min of the internal data range (the prototype min)
#' @param int_max max of the internal data range (the prototype max)
#' @details Both \code{ext_min} and \code{ext_max} are stored internally 
#' (as \code{netrng_ext_min} and \code{netrng_ext_max}, respectively) vectors of length \code{d} (data dimension). 
#' Thus, either can be input as either a length \code{d} vector (to specify a different range of linearly scaling across dimension),  
#' or as a single number (which will be recycled across dimension to produce a length \code{d} vector). 
#' \code{int_min} and \code{int_max} (which are stored internally as \code{netrng_int_min} and \code{netrng_int_max}, respectively) 
#' should be given as a single number (not a vector, varying internal ranges across dimension is not currently supported).  
#' 
#' \code{set_netrng} is called by default during \link{initialize_SOM} using the observed min/max of the data as the external range, 
#' and 0/1 as the internal range.  
#' 
#' Scaling can be turned off by setting \code{ext_min = int_min = 0} and \code{ext_max = int_max = 1}. 
#' 
#' @return None, scaling ranges are set internally 
#' @usage SOMobj$set_netrng(ext_min, ext_max, int_min, int_max). 
NULL

#' @name map_to_netrng
#' @title Map a data matrix to internal network range 
#' @description During training, prototype updates and BMU selelction are performed in the internal network range.  
#' This function allows mapping an arbitrary data matrix to this range to, e.g., replicate BMU selection or 
#' be able to compare data to prototypes. 
#' @param X a data matrix (occupying the external network range) to be mapped to the internal network range. 
#' \code{X} can be either a single vector, or a matrix with multiple data vectors in its rows. 
#' It must have the same dimension (number of columns) as the SOM (which is stored internally as \code{d}).
#' @return a matrix (same dimension as input) containing the scaled data vectors in its rows 
#' @usage SOMobj$map_to_netrng(X)
NULL 

#' @name map_from_netrng 
#' @title Map a prototype matrix to external network range 
#' @description Provides reverse functionality of \link{map_to_netrng}, to linearly scale vectors from internal to external network range. 
#' @param W a matrix occupying internal network range to be mapped to external network range. 
#' \code{W} can be either a single vector, or a matrix with multiple vectors in its rows. 
#' It must have the same dimension (number of columns) as the SOM (which is stored internally as \code{d}). 
#' 
#' Note that the SOM prototypes are computed and stored in internal network range; to compare them directly 
#' to the training data they should be mapped to the external range via this function. 
#' 
#' @return a matrix (same dimension as input) containing the scaled data vectors in its rows 
#' @usage SOMobj$map_from_netrng(W)
NULL 


# ***** Prototype Setting ***** 

#' @name set_W_runif
#' @title Random uniform initialization of SOM prototypes
#' @description SOM training typically begins with randomly initialized prototype vectors. 
#' This function performs this random initialization in the range \code{mid} +/- 5% of the internal network range, 
#' where \code{mid} is the midpoint of \code{netrng_int_min} and \code{netrng_int_max}. 
#' 
#' To force random sampling according to a specified seed, call R's \link[base]{set.seed} function just prior to calling \code{set_W_runif}
#' 
#' Note: this function is called automatically by \link{initialize_SOM}.
#' 
#' @return None, prototype matrix \code{W} is initialized internally. 
#' @usage SOMobj$set_W_runif()
NULL 

#' @name set_W
#' @title Setter function for the prototype matrix \code{W}
#' @description Random prototype initialization is performed via \link{set_W_runif}. 
#' This function allows setting \code{W} to user-prescribed values. 
#' If the input values reside outside of the internal network range a warning will be issued 
#' @param W matrix of prototype values, with rows following neuron ordering. 
#' \code{nrow(W)} must equal \code{nW} and \code{ncol(W)} must equal \code{d}, otherwise an error will be returned. 
#' @return None, prototype matrix \code{W} is set internally. 
#' @usage SOMobj$set_W(W)
NULL 


# ***** Winfrq Setting *****

#' @name set_p_equal
#' @title Equi-probable initialization of CSOM win frequencies
#' @description The Conscience-SOM algorithm of DeSieno attempts a Maximum-Entropy SOM mapping 
#' by updating a record of each prototype's win frequency. This roughly mimics the proportion of times during 
#' training that each prototype is selected as some datum's BMU.  
#' This function initializes the win frequencies for all prototypes to 1 / \code{nW}, which are stored internally in field \code{p}.
#' 
#' Note: this function is called automatically by \link{initialize_SOM}.
#' 
#' @return None, win frequency vector \code{p} is stored interally. 
#' @usage SOMobj$set_p_equal()
#' @references 
#' \insertRef{DeSieno1988}{SOMDisco}
NULL 

#' @name set_p 
#' @title User-specified initialization of CSOM win frequencies
#' @description This function allows setting of the internal field \code{p}, which stores the CSOM win frequencies for each prototype. 
#' It is made available mostly for experimentation. 
#' @param p user-specified win frequency vector. Must have length = \code{nW}
#' @return None, win frequency vector \code{p} is stored internally 
#' @usage SOMobj$set_p(p)
NULL 


# ***** Learning Rates ***** 

#' @name set_LRAS
#' @title Set the Learning Rate Annealing Schedule 
#' @description The LRAS is a data frame with columns: 
#' \itemize{
#' \item t \item alpha \item beta \item gamma \item sigma
#' }
#' The rows of the data frame specify the learning parameters \code{alpha}, \code{beta}, \code{gamma} and \code{sigma} 
#' that are in effect up to (but not including) the training step identified in \code{t}.  
#' At each training step alpha controls the strength of the prototype update, beta controls the strength of the win frequency update, 
#' gamma controls the magnitude of the entropy-maximizing bias, and sigma controls the maximum radius of the neighborhood function eta. 
#' 
#' Internally, the columns of the input are extracted and stored in a \code{std::map} for quick lookup. At each training step, effective 
#' values from this map are extracted and stored in fields \code{alpha}, \code{beta}, \code{gamma} and \code{sigma} (the latter, in turn, 
#' also updates the neighborhood function \code{eta}). 
#' 
#' A default LRAS, based on the number of training sampled, can be obtained via \link{default_LRAS}.  
#' 
#' Internally, a final time step \code{t = std::max_element<unsigned int>} is appended to the LRAS, which has the effect of 
#' recycling the last set of learning parameters for indefinite training.  
#' 
#' Kohonen's original SOM algorithm can be achieved by setting both beta and gamma = 0 in the LRAS.  
#' 
#' @param LRAS a data frame, which must have the above structure. 
#' @return None, the LRAS is stored internally 
#' @usage SOMobj$set_LRAS(LRAS)
#' @references 
#' \insertRef{DeSieno1988}{SOMDisco}
NULL 

#' @name get_LRAS
#' @title Get the Learning Rate Annealing Schedule 
#' @description A getter function for the LRAS (described in \link{set_LRAS}) stored in the SOM object. 
#' @return the effective LRAS of the SOM, as a data frame. 
#' @usage SOMobj$get_LRAS() 
NULL 

#' @name calc_eta
#' @title Calculate the SOM neighborhood function 
#' @description Topology preservation in SOM mappings is enforced by the neighborhood function, 
#' which, during each training step, specifies cooperative updates of a small radius of prototype's whose neurons neighbor the BMU. 
#' The maximum neighborhood radius at each time is set by the \code{sigma} parameter specified in \link{set_LRAS}.  
#' The eta coefficients are computed from a logistic decay up to this \code{sigma}, so that during each training step 
#' the winning prototype is updated with the strongest effect, the prototypes within a radius = 1 on the lattice are updated with 
#' the next strongest effect, and so on (up to a maximum radius of \code{sigma}).  
#' 
#' It is usually not necessary to directly call this function as it is invoked inside \link{train_SOM}. 
#' 
#' @param sigma the maximum neighborhood radius for which the neighborhood functional is applied 
#' 
#' @return a vector of the eta neighborhood coefficients for a given value of \code{sigma}. 
#' The vector is ordered such that \code{eta[1]} is the coefficient applied to the prototype update of the BMU, 
#' \code{eta[2]} is the coefficient applied to prototype updates within a radius=1 of the BMU, 
#' \code{eta[3]} is the coefficient applied to prototype updates within a radius=2 of the BMU, etc.  
#' 
#' @usage SOMobj$calc_eta(sigma)
NULL 

#' @name update_learning_rates
#' @title Update the effective learning rates 
#' @description The CSOM learning rates (described in \link{set_LRAS}) should be annealed over time. 
#' This function queries the LRAS given the current training age of the SOM and updates the effective 
#' \code{alpha}, \code{beta}, \code{gamma}, \code{sigma}, and \code{eta} (the last via a call to \link{calc_eta}) parameters accordingly.  
#' 
#' It is usually not necessary to directly call this function as it is invoked inside \link{train_SOM}. 
#' 
#' @return None, the above learning rate fields are set internally. 
#' @usage SOMobj$update_learning_rates()
NULL 


# ***** Training *****

#' @name update_p
#' @title Update the CSOM win frequencies
#' @description The win frequency field \code{p} is updated internally according to DeSieno's CSOM algorithm.
#' 
#' It is usually not necessary to directly call this function as it is invoked inside \link{train_SOM}. 
#' 
#' @param winner_idx the prototype index of the BMU selected during a CSOM training step 
#' 
#' @return None
#' @usage SOMobj$update_p(winner_idx)
NULL 

#' @name update_W
#' @title Update the SOM prototypes
#' @description The prototype matrix \code{W} is updated internally according to DeSieno's CSOM algorithm.
#' 
#' It is usually not necessary to directly call this function as it is invoked inside \link{train_SOM}. 
#' 
#' @param winner_idx the prototype index of the BMU selected during a CSOM training step 
#' 
#' @return None
#' @usage SOMobj$update_p(winner_idx)
NULL 

#' @name train_SOM 
#' @title Training function for a SOM object 
#' @description (C)SOM training is performed and all associated learning products updated internally, 
#' according to DeSieno's algorithm, using the learning rates specified in \link{set_LRAS}.  
#' 
#' @param nsteps the number of training steps to perform
#' @param X the training data matrix. 
#' This should be the same matrix that was input to \link{initialize_SOM} 
#' (checks will be performed on its statistics, if they do not match an error will be returned). 
#'  
#' @details This function calls many of the helper function of the SOM class internally 
#' (e.g., \link{update_p}, \link{update_W}, \link{map_to_netrng}, \link{update_learning_rates}). 
#'  
#' Training is performed online (not batch) according to a sequential random sampling of the rows of \code{X}. 
#' Reproducibility of training results can be achieved by calling R's \link[base]{set.seed} function prior to calling \code{train_SOM}, 
#' which will control the random seed used for sampling \code{X}.  The training order used is during each call to \code{train_SOM} is stored 
#' in the field \code{train_order}.  
#'  
#' Additionally, \link{recall_SOM} is called at the end of training to update the SOM products.  
#'  
#' Parallel computation for BMU selection and prototype updating can be achieved by calling \link{set_parallel} prior to calling \code{train_SOM}. 
#'  
#' Various visualizations exist to examine the results of SOM training. See any of the \code{vis_som_*} functions for more information. 
#' 
#' @usage SOMobj$train_SOM(nsteps, X)
#'  
#' @return None, SOM products are updated internally 
#' @references 
#' \insertRef{DeSieno1988}{SOMDisco}
NULL 

#' @name set_monitoring_freq
#' @title (De-)Active Monitoring of SOM Training 
#' @description Various quantities can be computed and stored at regular training step intervals to allow the analyst 
#' to observe how SOM learning progresses.  This function either activates or de-activates this monitoring. 
#' @param mtr_freq the incremental training step count at which monitoring snapshots are taken and stored during SOM training. 
#' If = 0, which is the default set when calling \link{initialize_SOM}, then no monitoring is performed. 
#' 
#' @details As monitoring proceeds, the internal field \code{mtr_age} is updated with the time (training) steps 
#' at which snapshots were taken. 
#' 
#' If activated, monitoring computes and stores the following quantities: 
#' \itemize{
#' \item \code{mtr_RMSQE} a matrix (nrow = length(\code{mtr_age}), ncol = \code{nW}+1) containing the 
#' Root Mean Squared Error of the quantization of the data by its BMU prototype (the Root Mean Quantization Error at the prototype level). 
#' The last (\code{nW} + 1) column contains the RMSQE of the quantization over all prototypes.
#' \item \code{mtr_QSI} a vector (length = length(\code{mtr_age})) of the Quantization Stability Index 
#' computed during each monitoring snapshot. The QSI at any monitoring step is the proportion of data samples 
#' that have switched BMUs from the previous monitoring snapshot (by convention, QSI at the initial monitoring step = NA). 
#' QSI can convey when the SOM "settles down", and reaches a stable quantization.  
#' \item \code{mtr_Entropy} a vector (length = length(\code{mtr_age})) of the normalized Entropies of the SOM mappings 
#' taken at each monitoring step.  
#' }
#' 
#' The monitored quantities can be visualized after a training round by calling \link{vis_som_training}.  
#' @return None, field \code{mtr_freq} is set internally. 
#' @usage SOMobj$set_monitoring_freq(mtr_freq)
NULL 

#' @name set_train_order
#' @title Set the sample training order 
#' @description Usually during training, sample vectors should be drawn randomly 
#' (which \link{train_SOM} does by default, and stores a record of the drawing order here). 
#' If a specific training order is desired (for, e.g., experimental purposes), 
#' this can be set by calling \code{set_train_order} prior to calling \code{train_SOM}. 
#' This forces data presentation during SOM training to follow the prescribed ordering.  
#' If a specific training order is set, it is only valid during the subsequent call to \code{train_SOM}; 
#' if multiple rounds of training are performed, each with a desired training order, then \code{set_train_order} 
#' must be called again prior to each call to \code{train_SOM}.
#' 
#' @param train_order a vector containing the (1-based) indices (i.e., row indices of training data \code{X}) which 
#' will govern the presentation order during training. 
#' 
#' @details If \code{set_train_order} is called it will override the specification of the \code{nsteps} argument during 
#' a subsequent call to \code{train_SOM}. That is, the entire set of training data specified by input \code{train_order} will be 
#' used during SOM training, regardless of the number of training steps requested during the call to \code{train_SOM}.  
#' 
#' 
#' The internal logical field \code{user_train_order} tracks whether this function has been called. 
#' 
#' @return None, \code{train_order} and \code{user_train_order} are set internally. 
#' @usage SOMobj$set_train_order(train_order)
NULL 


# ***** Recall Functions ***** 

#' @name set_nBMU
#' @title Set the number of BMUs 
#' @description The Best Matching Unit of a sample vector \code{x} is the prototype \eqn{w_j} (and its associated neuron) that 
#' minimizes Euclidean distance \eqn{d(x,w_j)}.  The SOM requires calculating at least the first BMU during training and recall. 
#' Calculation of the \code{CADJ} matrix requires calculating and storing the second BMU as well (the second Best Matching Unit).  
#' This function sets the internal field \code{nBMU}, which specifies the number of BMUs stored in the field \code{BMU} during recall. 
#' Any number > 2 could be useful for further analysis.  Be default, \code{nBMU} is set = 2 during \link{initialize_SOM}. 
#' @param nBMU the number of BMUs requested to compute
#' @return None, \code{nBMU} is set internally.  
#' @usage SOMobj$set_nBMU(nBMU)
NULL 

#' @name recall_SOM 
#' @title Recall a trained SOM object 
#' @description A SOM recall maps all training data to their representative prototypes (their BMUs). 
#' Several quantities, as outlined below, result from this mapping.  
#' @param X the training data matrix. 
#'  This should be the same matrix that was input to \link{initialize_SOM} 
#'  (checks will be performed on its statistics, if they do not match an error will be returned). 
#' @details Several internal fields are set by the recall function: 
#' \itemize{
#' \item \code{BMU} matrix (nrow = \code{nX}, ncol = \code{nBMU}) containing the 1st, 2nd, ... BMUs for every training vector in its columns. 
#' The valued in \code{BMU} are the 1-based indices of the prototypes, in neuron order. 
#' \item \code{SQE} vector (length = \code{nX}) containing the squared quantization error of each training vector, as quantized by its 1st BMU
#' \item \code{CADJ} the CADJ matrix (nrow = ncol = \code{nW}) of Cumulative (weighted) Topological Adjacencies of the SOM prototypes. 
#' See \link{CONN} for details. 
#' \item \code{RF_size} vector (length = \code{nW}) containing the number of training vectors mapped to each prototype
#' \item \code{Entropy} the (normalized) entropy of the discrete SOM mapping
#' \item \code{RF_members} a list (length = \code{nW}) of the training vector indices mapped to each prototype 
#' \item \code{fences} a data frame containing the information for visualizing the U-matrix fences on the SOM lattice, 
#' see \link{set_lattice_fences} for details. 
#' \item \code{RF_label} and \code{RF_label_dist}, the values of \code{X_label} propagated to the Receptive Fields, 
#' according to the SOM mapping. 
#' }
#' 
#' Note: the "Receptive Field" of a prototype \eqn{w_j} is the set of all training data for which \eqn{w_j} is the BMU. 
#' 
#' It is usually not necessary to directly call this function as it is invoked inside \link{train_SOM}. 
#' 
#' @return None, the fields described above are calculated and stored
#' @usage SOMobj$recall_SOM(X)
NULL 

# @name set_RF_label
# @title Propagate training data labels to SOM prototypes 
# @description If the training data possesses labels, the learned SOM mapping induces an associated labeling of the prototypes. 
# These labels are set equal to the plurality winning label of data in each prototype's Receptive Field. 
# 
# For example: Suppose 10 training data vectors are mapped to prototype \eqn{w_j} and  
# 4 of these have label "A", 5 have label "B", while 1 has label "C". The plurality label of \eqn{w_j} = "B" in this case.  
# 
# Tracking and storing the prototype labels provides an easy way to assess the organization of the learned SOM; 
# if the labels actually represent cohesive sub-groupings of the training data then we expect groups of labeled 
# prototypes to cluster near each other on the SOM lattice.  See \link{vis_som_tiles} for more information on visualizing 
# the prototype labels on the lattice.  
# 
# @param Xlabel a character vector (length = \code{nX}) containing a label for each row of the training data matrix \code{X}
# 
# @return None, the field \code{RF_label} is computed and stored internally, which is a vector (length = \code{nW}) containing 
# the plurality label of each prototype. 
# 
# @usage SOMobj$set_RF_label(Xlabel)
# NULL 

#' @name set_lattice_fences
#' @title Compute and store the U-Matrix fences
#' @description The U-Matrix fences, when visualized on the SOM lattice, are an easy way to 
#' determine the similarity between neighboring lattice-neighboring prototypes.  This function 
#' computes a data frame containing the information required to visiualize such fences, which 
#' can be accomplished via \link{vis_som_fences}.  
#' 
#' It is usually not necessary to directly call this function as it is invoked inside \link{recall_SOM}. 
#' 
#' @return None, data frame \code{fences} is set internally 
#' @usage SOMobj$set_lattice_fences()
#' @references 
#' \insertRef{UMatrix}{SOMDisco}
NULL 

#' @name CONN 
#' @title The CONNectivity matrix of SOM prototypes 
#' @description The CONN matrix is calculated from the \code{CADJ} matrix computed during \link{recall_SOM}: 
#' \deqn{CONN = CADJ + t(CADJ)}. 
#' 
#' \code{CADJ} describes the strength of the topological connectivites between SOM prototypes in input space 
#' (higher CADJ values between prototypes \eqn{w_i} and \eqn{w_j} indicate stronger evidence that \code{w_i} and \code{w_j} belong to the same cluster). 
#' The (i,j) values of CADJ are defined: 
#' \deqn{CADJ_{ij} = \#(BMU1(X)=i and BMU2(X)=j)}
#' 
#' By this definition, \code{CADJ} is an asymmetric adjacency matrix. CONN is its symmetric counterpart.  
#' 
#' @return the CONN matrix (nrow = ncol = \code{nW})
#' 
#' @usage SOMobj$CONN()
#' @references 
#' \insertRef{TasdemirMerenyi2009}{SOMDisco}
NULL 

#' @name calc_Entropy
#' @title Calculate the Normalized Entropy of the SOM mapping
#' @description The normalized entropy of the SOM quantization is given by 
#' \deqn{entropy = -sum(F*log(F))/log(nW)}
#' where \code{F} is a vector of Receptive Field relative frequencies (i.e., \eqn{RF_size / nX}). 
#' 
#' SOMs trained with the CSOM update rule seek to maximize the entropy of the learned mapping, so 
#' this normalized entropy measure provides a way of comparing different SOMs regardless of their size 
#' (or the number of training vectors available).
#' 
#' It is usually not necessary to directly call this function as it is invoked inside \link{recall_SOM}.  
#' 
#' @return None, the field \code{Entropy} is set internally. 
#' 
#' @usage SOMobj$calc_Entropy()
NULL 


# ***** Saving / Loading ***** 

#' @name save 
#' @title Save a SOM object 
#' @description All fields in a SOM object can be saved to disk with this function, which allows them to be re-loaded 
#' into a new R environment at a later time (for analysis, or possibly extended training).  
#' @param somfile a string indicating the file path and name in which to save the SOM object.  
#' This must end in extension ".som", otherwise an error is returned. 
#' @details The SOM object is saved to disk as an R list, with each field occupying a corresponding field of the list.  
#' The file is saved in .rds format (can check its details with \link[base]{infoRDS}).  
#' 
#' Saved SOMs can be re-loaded with \link{load}
#' 
#' @return None, the SOM object is saved to disk 
#' @usage SOMobj$save(somfile)
NULL 

#' @name load
#' @title Load an existing SOM object 
#' @description SOM objects previously written to disk via the \link{save} method can be re-loaded into 
#' a new R environment with this function. All fields of the internal C++ class will be populated, and 
#' all methods can be called on the loaded SOM object.  
#' @param somfile a string indicating the file path and name of the saved SOM object.  
#' @details Because the .som file is in .rds format it can, technically, be loaded directly into an R environment 
#' as a list via \link[base]{readRDS}. This can be useful for spot checking the contents of a saved SOM object, but 
#' does not allow use of any of its methods (or visualizations).  The \code{load} methods allows for 
#' proper restoration of a previously saved SOM. 
#' 
#' @return None, the SOM object is loaded 
#' @usage SOMobj$load(somfile)
NULL 

#' @name as_list
#' @title Convert an SOM object to a list 
#' @description This method extracts all fields of an SOM object and places their stored values 
#' into the fields of an R list, with list field names matching the SOM object field names. 
#' @return A list
#' @usage SOMobj$as_list()
NULL 

#' @name load_list 
#' @title Populate a SOM object from a list 
#' @description This method populates all fields of a SOM object from the fields of an R list object. 
#' The list must have field names which exactly match the SOM field names. 
#' @param SOMList a SOM object converted to a list, e.g., with \code{as_list}
#' @return None
#' @usage SOMobj$load_list(SOMList)
NULL 

