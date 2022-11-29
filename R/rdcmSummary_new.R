#' @param x data matrix of dimension \eqn{n} x \eqn{p}.
#' @param Omega_t_list \eqn{p} x \eqn{p} x \eqn{T_count} array of \eqn{T_count} number of estimated dynamic precision matrices.
#' @param G_est \eqn{p} x \eqn{p} x \eqn{T_count} array of \eqn{T_count} number of estimated dynamic networks.
#' @param n0 Number of points in each segmentation.

bic_cal <- function(S, Omega_t_list, n0, G_est=NULL) {
  #regular BIC
  p = dim(S)[1]
  T_count = dim(S)[3]
  
  if(is.null(G_est)) G_est = (abs(Omega_t_list) >= 0.001) - array(diag(p), c(p, p, T_count))
  ## 1 x T_count vector of nedges of T_count Omega_t
  nedges = apply(G_est, 3, sum) / 2
  
  bic <- mapply(function(x1, x2, x3, x4) {x1*(sum(diag(x2 %*% x3)) - log(det(x3))) + x4 * log(x1)},
                n0, 
                plyr::alply(S, 3), 
                plyr::alply(Omega_t_list, 3), 
                nedges+p)
  
  return(sum(bic))
}

bic_cal_m <- function(S, Omega_t_list, n0, h, G_est=NULL) {
  #regular BIC
  p = dim(S)[1]
  T_count = dim(S)[3]
  
  if(is.null(G_est)) G_est = (abs(Omega_t_list) >= 0.001) - array(diag(p), c(p, p, T_count))
  ## 1 x T_count vector of nedges of T_count Omega_t
  nedges = apply(G_est, 3, sum) / 2
  
  bic_m <- mapply(function(x1, x2, x3, x4) {x1*(sum(diag(x2 %*% x3)) - log(det(x3))) + x4 * log(x1) / h},
                n0, 
                plyr::alply(S, 3), 
                plyr::alply(Omega_t_list, 3), 
                nedges+p)
  
  return(sum(bic_m))
}

aic_cal <- function(S, Omega_t_list, n0, G_est=NULL) {
  #regular AIC
  p = dim(S)[1]
  T_count = dim(S)[3]
  
  if(is.null(G_est)) G_est = (abs(Omega_t_list) >= 0.001) - array(diag(p), c(p, p, T_count))
  nedges = apply(G_est, 3, sum) / 2
  
  aic <- mapply(function(x1, x2, x3, x4) { 
    x1*sum(diag(x2 %*% x3)) - x1*log(det(x3)) + 2*x4},
    n0, 
    plyr::alply(S, 3), 
    plyr::alply(Omega_t_list, 3), 
    nedges+p)
  
  return(sum(aic))
}

aic_cal_m <- function(S, Omega_t_list, n0, h, G_est=NULL) {
  #regular AIC
  p = dim(S)[1]
  T_count = dim(S)[3]
  
  if(is.null(G_est)) G_est = (abs(Omega_t_list) >= 0.001) - array(diag(p), c(p, p, T_count))
  nedges = apply(G_est, 3, sum) / 2
  
  aic_m <- mapply(function(x1, x2, x3, x4) { 
    x1*sum(diag(x2 %*% x3)) - x1*log(det(x3)) + 2*x4/h},
    n0, 
    plyr::alply(S, 3), 
    plyr::alply(Omega_t_list, 3), 
    nedges+p)
  
  return(sum(aic_m))
}

nll_cal <- function(S, Omega_t_list, n0, G_est=NULL) {
  #negative log-likelihood
  p = dim(S)[1]
  T_count = dim(S)[3]
  
  nll <- mapply(function(x1, x2, x3) {x1*(sum(diag(x2 %*% x3)) - log(det(x3)))},
                n0, 
                plyr::alply(S, 3), 
                plyr::alply(Omega_t_list, 3))
  
  return(sum(nll))
}

bic_penalty_cal <- function(S, Omega_t_list, n0, h, G_est=NULL) {
  #regular BIC
  p = dim(S)[1]
  T_count = dim(S)[3]
  
  if(is.null(G_est)) G_est = (abs(Omega_t_list) >= 0.001) - array(diag(p), c(p, p, T_count))
  ## 1 x T_count vector of nedges of T_count Omega_t
  nedges = apply(G_est, 3, sum) / 2
  
  bic_penalty <- mapply(function(x1, x2, x3, x4) {x4 * log(x1) / h},
                n0, 
                plyr::alply(S, 3), 
                plyr::alply(Omega_t_list, 3), 
                nedges+p)
  
  return(sum(bic_penalty))
}

aic_penalty_cal <- function(S, Omega_t_list, n0, h, G_est=NULL) {
  #regular BIC
  p = dim(S)[1]
  T_count = dim(S)[3]
  
  if(is.null(G_est)) G_est = (abs(Omega_t_list) >= 0.001) - array(diag(p), c(p, p, T_count))
  ## 1 x T_count vector of nedges of T_count Omega_t
  nedges = apply(G_est, 3, sum) / 2
  
  aic_penalty <- mapply(function(x1, x2, x3, x4) {2*x4/h},
                        n0, 
                        plyr::alply(S, 3), 
                        plyr::alply(Omega_t_list, 3), 
                        nedges+p)
  
  return(sum(aic_penalty))
}