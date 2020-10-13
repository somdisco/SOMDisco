# SOMDisco

SOMDisco is an R package intended to address several advances in SOM learning and analysis that are currently missing in currently available SOM packages for R ([som](https://cran.r-project.org/web/packages/som/), [kohonen](https://cran.r-project.org/web/packages/kohonen/), and [popsom](https://cran.r-project.org/web/packages/popsom/index.html)).  The main contributions of SOMDisco are: 

+ Integration of DeSieno's Conscience-SOM modifications to Kohonen's original algorithm to attempt maximum entropy SOM learning
+ Fast and efficient C++ implementation of CSOM training (based on [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) and [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html))
+ Optional parallel SOM training and recall (as applicable, via [RcppParallel](https://cran.r-project.org/web/packages/RcppParallel/index.html))
+ Computation of the CADJ matrix, which is a weighted topological adjacency of SOM prototypes helpful in cluster discovery
+ The CONNvis visualization (requires the `TopoRNet` R package) which represents the topological information in CADJ on the SOM lattice
+ Other advanced SOM visualizations (the mU-matrix, propagation of discrete and continuous values on the lattice)

# Installation

```r
devtools::install_github("somdisco/SOMDisco")
```

# Documentation

See the [SOMDisco homepage](https://somdisco.github.io/SOMDiscoR/output/index.html) for more information.
