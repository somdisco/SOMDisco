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

#ifndef SOM_UTILS_DIST_HPP
#include "som_utils_dist.hpp"
#endif


#ifndef SOM_UTILS_LATTICE_HPP
#define SOM_UTILS_LATTICE_HPP

namespace SOM_UTILS_LATTICE {

// *** Convert between radians, degrees 
inline double rad2deg(double rad) {
  return rad * 180 / M_PI; 
}

inline double deg2rad(double deg) {
  return deg * M_PI / 180; 
}

// *** Convert cartesian coords to polar coords
inline arma::rowvec cart2polar(const arma::rowvec& x) {
  
  arma::rowvec out(2); 
  
  // radius 
  out[0] = std::sqrt(arma::accu(arma::square(x)));
  
  // theta, in degrees 
  out[1] = std::atan2(x[1], x[0]); 
  if(out[1] < 0) out[1] += 2*M_PI; 
  out[1] = rad2deg(out[1]); 
  
  return out; 
}

// *** Convert verts of a lattice tile from cartesian to polar coordinates
inline arma::mat cartverts2polarverts(const arma::mat& cartverts, arma::rowvec nu_xy) {
  
  arma::mat polar(cartverts.n_rows, 2); 
  
  for(unsigned int k=0; k<cartverts.n_rows; ++k) {
    polar.row(k) = cart2polar(cartverts.row(k) - nu_xy); 
  }
  
  return polar; 
}

// *** Calculate the (x,y) and (i,j) coordinates of a hex or grid lattice 
inline void calc_lattice_hex(arma::mat& xy, arma::umat& ij, unsigned int som_x, unsigned int som_y) {
  
  std::vector<double> x,y; 
  std::vector<double> i,j; 
  
  
  for(unsigned int row = 1; row <= som_y; ++row) {
    double this_y = 1 + (row-1) * 2 / std::sqrt(3) * 3/4; // 2/sqrt(3) * 3/4 is vertical distance btw hex centers
    unsigned int min_col = 1, max_col = som_x; 
    double offset = 0; 
    if(row % 2 == 0) {
      //min_col -= 1; // shift every other row's tile to the left 
      max_col += 1;
      offset += 0.5;
    } 
    
    
    for(unsigned int col = min_col; col <= max_col; ++col) {
      x.push_back(double(col)-offset); 
      y.push_back(this_y);
      i.push_back(double(row));
      j.push_back(double(col));
    }
  }
  
  // Assemble (x,y) and (i,j) into matrices 
  xy.resize(x.size(), 2); 
  xy.col(0) = arma::conv_to<arma::vec>::from(x); 
  xy.col(1) = arma::conv_to<arma::vec>::from(y); 
  
  ij.resize(x.size(), 2); 
  ij.col(0) = arma::conv_to<arma::uvec>::from(i); 
  ij.col(1) = arma::conv_to<arma::uvec>::from(j); 
}

inline void calc_lattice_grid(arma::mat& xy, arma::umat& ij, unsigned int som_x, unsigned int som_y) {
  
  std::vector<double> x,y; 
  std::vector<double> i,j; 
  
  for(unsigned int row = 1; row <= som_y; ++row) {
    for(unsigned int col = 1; col <= som_x; ++col) {
      x.push_back(double(col)); 
      y.push_back(double(row));
      i.push_back(double(row));
      j.push_back(double(col));
    }
  }
  
  // Assemble (x,y) and (i,j) into matrices 
  xy.resize(x.size(), 2); 
  xy.col(0) = arma::conv_to<arma::vec>::from(x); 
  xy.col(1) = arma::conv_to<arma::vec>::from(y); 
  
  ij.resize(x.size(), 2); 
  ij.col(0) = arma::conv_to<arma::uvec>::from(i); 
  ij.col(1) = arma::conv_to<arma::uvec>::from(j); 
}

// *** Calculate the nearest-neighbor adjacency matrix of lattice neurons 
inline arma::umat calc_lattice_adjacency(const arma::mat& xy) {
  unsigned int nW = xy.n_rows; 
  
  // Get LInf distance 
  SOM_UTILS_DIST::distmat_worker distwkr(xy, xy, "LInf");
  distwkr.calc_parallel(); 
  
  arma::umat out(nW,nW); 
  out.zeros(); 
  out.elem(arma::find(distwkr.DIST < 1.01)).fill(1); 
  return out; 
}

// *** Build the lattice distance list. 
// This is a list of lists of lattice distances: 
// ex: distlist[i][j] is a vector of neuron indices that are distance "j" from neuron "i"
struct lattice_distlist_worker : public RcppParallel::Worker {
  
  const arma::umat& SPDIST; // Shortest path distance matrix  
  
  // Internal variables 
  unsigned int max_dist;
  unsigned int nverts; 
  
  // output container
  std::vector<std::vector<arma::uvec>> out;
  
  // Constructor
  lattice_distlist_worker(const arma::umat& SPDIST) : SPDIST(SPDIST) {
    nverts = SPDIST.n_rows;
    out.resize(nverts);
    max_dist = SPDIST.max(); 
  }
  
  // Fill up neighbors from a single neuron 
  void distlist_neuron_j(unsigned int j) {
    // NOTE:  The indices stored here are 1-indexed!!! 
    
    // There are max_dist + 1 levels of lattice neighbors
    out[j].resize(max_dist+1);
    
    // The first tier of neighbors has radius = 0 (prototype j is neighbor of itself)
    out[j][0] = j + 1;
    
    // Fill up the rest of the neighborhood tiers 
    for(unsigned int d=1; d<=max_dist; ++d) {
      arma::uvec tmp = arma::find(SPDIST.row(j) == d);
      if(tmp.size() > 0) {
        out[j][d] = tmp + 1;   
      } else {
        out[j][d] = tmp; 
      }
    }
  }
  
  // Parallel operator 
  void operator()(std::size_t begin, std::size_t end) {
    
    for(unsigned int j=begin; j<end; ++j) {
      distlist_neuron_j(j); 
    }
  }
  
  // Parallel invoker 
  void calc_parallel() {
    RcppParallel::parallelFor(0, nverts, *this);
  }
  
  // Serial invoker 
  void calc_serial() {
    for(unsigned int j=0; j<nverts; ++j) {
      distlist_neuron_j(j); 
    }
  }
  
};


// *** Compute the vertices of each lattice tile
// Output is a cube whose i-th slice contains the vertices in counter-clockwise order (starting from the positive x-axis)
struct lattice_tile_verts_worker : public RcppParallel::Worker {
  
  const arma::mat& xy; 
  std::string lattice_type; 
  
  // output container
  arma::cube verts;
  
  // Constructor, with bias 
  lattice_tile_verts_worker(const arma::mat& xy, std::string lattice_type_) : xy(xy) {
    lattice_type = lattice_type_; 
    if(lattice_type == "grid") {
      verts.set_size(4,2,xy.n_rows);
    } else if(lattice_type == "hex") {
      verts.set_size(6,2,xy.n_rows);
    }
  }
  
  // Fill up neighbors from a single neuron 
  void verts_of_tile_j(unsigned int j) {
    
    verts.slice(j).each_row() = xy.row(j); 
    
    if(lattice_type == "grid") {
      
      verts(0,0,j) += 0.5; verts(0,1,j) -= 0.5; 
      verts(1,0,j) += 0.5; verts(1,1,j) += 0.5; 
      verts(2,0,j) -= 0.5; verts(2,1,j) += 0.5; 
      verts(3,0,j) -= 0.5; verts(3,1,j) -= 0.5; 
      
    } else if(lattice_type == "hex") {
      
      for(unsigned int i=0; i<6; ++i) {
        double anglerad = (60 * double(i) - 30) * M_PI / 180;
        verts(i,0,j) += std::cos(anglerad) / std::sqrt(3); 
        verts(i,1,j) += std::sin(anglerad) / std::sqrt(3); 
      }
    }
    
    // *** Once verts are stored, re-order them by increasing angle 
    // By convention, theta = 0 is + x-axis, and theta > 0 moves away in counter-clockwise manner 
    
    // Convert vertices to polar coordinates. 
    // Returns a 2-column matrix (r,theta), with theta in degrees
    arma::mat verts_polar = cartverts2polarverts(verts.slice(j), xy.row(j)); 
    
    // Sort by the theta angles 
    arma::uvec sortidx = arma::sort_index(verts_polar.col(1), "ascend");
    
    verts.slice(j) = verts.slice(j).rows(sortidx); 
    
  }
  
  // Parallel operator 
  void operator()(std::size_t begin, std::size_t end) {
    
    for(unsigned int j=begin; j<end; ++j) {
      verts_of_tile_j(j); 
    }
  }
  
  // Parallel invoker 
  void calc_parallel() {
    RcppParallel::parallelFor(0, xy.n_rows, *this);
  }
  
  // Serial invoker 
  void calc_serial() {
    for(unsigned int j=0; j<xy.n_rows; ++j) {
      verts_of_tile_j(j); 
    }
  }
  
};


// *** Find endpoints of a lattice tile edge in direction "angle" from its center
inline arma::mat get_tile_edge_endpts(const arma::mat& verts, const arma::rowvec& nu_xy, double angle) {
  // NOTE: Angle should be in degrees! 
  // returns a 2x2 matrix with endpts in rows 
  
  // Convert cartesian verts to polar verts
  arma::mat verts_polar = SOM_UTILS_LATTICE::cartverts2polarverts(verts, nu_xy); 
  
  // Initialize output matrix 
  arma::mat out(2,2); 
  
  // Extract the vertex angles of tile j into a vector 
  std::vector<double> vert_angles = arma::conv_to<std::vector<double>>::from(verts_polar.col(1)); 
  unsigned int nverts = verts_polar.n_rows;
  
  // First, check if angle is in the list of angles, within tolerance
  // if yes, then the "edge" spanning input angle is just the vertex point itself. 
  std::vector<double>::const_iterator findit = std::find_if(vert_angles.cbegin(), vert_angles.cend(),
                                                            [angle](double b) { return std::abs(angle - b) < 1e-6; });
  if(findit != vert_angles.cend()) {
    unsigned int idx = std::distance(vert_angles.cbegin(), findit); 
    out.row(0) = verts.row(idx);
    out.row(1) = verts.row(idx);
    return out; 
  }
  
  // Otherwise, get iterator to the lower_bound of angle in the angles vector 
  // lower_bound returns iterator to first element that is GREATER OR EQUAL to value
  // upper_bound returns iterator to first element that is STRICTLY GREATER than value 
  std::vector<double>::const_iterator lowit = std::lower_bound(vert_angles.cbegin(), vert_angles.cend(), angle); 
  //std::vector<double>::const_iterator uppit = std::upper_bound(vert_angles.cbegin(), vert_angles.cend(), angle); 
  
  // If the iterator is at the first position or last position
  // the endpoints are the first and last vertices 
  if(lowit==vert_angles.cbegin() || lowit==vert_angles.cend()) {
    out.row(0) = verts.row(0); 
    out.row(1) = verts.row(nverts-1); 
    return out; 
  }
  
  // Otherwise, we need to find the index of the lower_bound iterator 
  unsigned int idx2 = std::distance(vert_angles.cbegin(), lowit); 
  unsigned int idx1 = idx2 - 1;
  
  out.row(0) = verts.row(idx1);
  out.row(1) = verts.row(idx2);
  
  return out;
  
}

// *** Crossproduct of two length-2 vectors
inline double crossprod_rowvec(const arma::rowvec& x, const arma::rowvec& y) {
  return x[0]*y[1] - x[1]*y[0]; 
}

// *** Find the (x,y) coordinates along a lattice tile boundary where a ray with angle theta intersects the tile boundary
struct tile_boundary_intersect_worker : public RcppParallel::Worker {
  
  const arma::cube& verts; 
  const arma::mat& nu_xy; 
  double theta; // direction of ray emanating from center, in degrees
  
  // Output container 
  arma::vec dist_to_edge; // holds length of line segment connecting tile center to point of intersection with boundary along ray of direction theta
  
  // Constructor, with bias 
  tile_boundary_intersect_worker(const arma::cube& verts, const arma::mat& nu_xy, double theta_) : 
    verts(verts), nu_xy(nu_xy) {
    theta = theta_; 
    dist_to_edge.set_size(verts.n_slices);  
  }
  
  // Re-order verts of a single tile 
  void dist_tile_j(unsigned int j) {
    
    // Get endpoints of the line segment forming the tile edge that intersects a ray in direction of angle 
    arma::mat endpts = get_tile_edge_endpts(verts.slice(j), nu_xy.row(j), theta); 
    
    // Describe this line segment in ray form q + u*s, where u in [0,1]
    arma::rowvec q = endpts.row(0); 
    arma::rowvec s = endpts.row(1) - endpts.row(0); 
    
    // Describe the ray emanating from center, in direction of angle, as center + t*r
    arma::rowvec r(2); 
    r[0] = std::cos(SOM_UTILS_LATTICE::deg2rad(theta)); 
    r[1] = std::sin(SOM_UTILS_LATTICE::deg2rad(theta));
    
    // If angle passes through a vertex, then the "edge" is a single point. 
    // In this case, the distance to the boundary is just the distance from tile center to vertex 
    if(std::sqrt(arma::accu(arma::square(s))) < 1e-6) {
      dist_to_edge[j] = std::sqrt(arma::accu(arma::square(q - nu_xy.row(j))));
      return; 
    }
    
    // Otherwise, the tile edge is a proper line segment. 
    // Compute distance to intersection of this ray and the edge segment 
    dist_to_edge[j] = crossprod_rowvec(q - nu_xy.row(j), s) / crossprod_rowvec(r, s); 
  }
  
  // Parallel operator 
  void operator()(std::size_t begin, std::size_t end) {
    for(unsigned int j=begin; j<end; ++j) {
      dist_tile_j(j); 
    }
  }
  
  // Parallel invoker 
  void calc_parallel() {
    RcppParallel::parallelFor(0, verts.n_slices, *this);
  }
  
  // Serial invoker 
  void calc_serial() {
    for(unsigned int j=0; j<verts.n_slices; ++j) {
      dist_tile_j(j); 
    }
  }
  
};


} // close namespace 

#endif