#' Generate a default Learning Rate Annealing Schedule 
#' 
#' @param nX the number of data vectors used for training
#' 
#' @return the Learning Rate Annealing Schedule as a data frame with columns
#' t (time step), 
#' alpha (prototype update rate), 
#' beta (win frequencies update rate), 
#' gamma (bias update rate), 
#' sigma (geodesic lattice distance to which prototype update is applied)
default_LRAS = function(nX) {
  LRAS = data.frame(t = integer(5), alpha = numeric(5), beta = numeric(5), gamma = numeric(5), sigma = numeric(5))
  
  LRAS$t[1] = 1 * nX
  LRAS$t[2] = 5 * nX
  LRAS$t[3] = 10 * nX
  LRAS$t[4] = 25 * nX
  LRAS$t[5] = 100 * nX
  
  LRAS$alpha[1] = 0.50
  LRAS$alpha[2] = 0.25
  LRAS$alpha[3] = 0.10
  LRAS$alpha[4] = 0.05
  LRAS$alpha[5] = 0.01
  
  LRAS$beta[1] = 0.50 / 10^1
  LRAS$beta[2] = 0.25 / 10^1
  LRAS$beta[3] = 0.10 / 10^1
  LRAS$beta[4] = 0.05 / 10^1
  LRAS$beta[5] = 0.01 / 10^1
  
  LRAS$gamma[1] = 0.50 / 10^2
  LRAS$gamma[2] = 0.25 / 10^2
  LRAS$gamma[3] = 0.10 / 10^2
  LRAS$gamma[4] = 0.05 / 10^2
  LRAS$gamma[5] = 0.01 / 10^2
  
  LRAS$sigma[1] = 3
  LRAS$sigma[2] = 2
  LRAS$sigma[3] = 1
  LRAS$sigma[4] = 1
  LRAS$sigma[5] = 1
  
  return(LRAS)
  
  
}