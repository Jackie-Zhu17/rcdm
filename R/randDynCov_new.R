#' Random Covariance Dynamic Model
#'
#' This function implements the Random Covariance Dynamic Model (RCDM) for joint estimation of
#' sparse time-varying precision matrices. Optimization is conducted using block
#' coordinate descent.
#' @param x data matrix of dimension \eqn{n} x \eqn{p}.
#' @param h Bandwidth of the kernel smooth function.
#' @param K Number of the knots.
#' @param lambda1 Non-negative scalar. Induces sparsity in dynamic precision matrix.
#' @param lambda2 Non-negative scalar. Induces similarity between dynamic precision matrix and kernel matrices.
#' @param n0 Number of points in each segmentation(bin).
#' @return A list of length 2 containing:
#' \enumerate{
#' \item Group-level precision matrix estimate (Omega0).
#' \item \eqn{p} x \eqn{p} x \eqn{K} array of \eqn{K} subject-level precision matrix estimates (Omegas).
#' }
#' @author
#' Hengcheng Zhu, Lin Zhang
#'
#' @export
#'
#' @references
#' Zhang Lin, Andrew DiLernia, Karina Quevedo, Jazmin Camchong, Kelvin Lim, and Wei Pan.
#' "A Random Covariance Model for Bi-level Graphical Modeling with Application to Resting-state FMRI Data." 2019. https://arxiv.org/pdf/1910.00103.pdf
library(doParallel)
library(abind)
library(plyr)
acomb <- function(...) abind(..., along=3)

randDynCov <- function(x, h, K, max.iter=1000, rho, lambda2, n0, ncore=8, k_list, parallel=FALSE) {

  n = dim(x)[1]
  p = dim(x)[2]

  # bin segmentation
  T1 <- n - n %% n0
  Time_series_array <- sapply(1:(T1/n0), function(t) {x[seq(n0*t-n0+1,n0*t),]}, simplify = 'array')
  t_list <- seq(n0%/%2+1,T1,by=n0)
  T_count <- length(t_list)

  S <- sapply(alply(Time_series_array,3),
              function(Tt) {t(Tt) %*% Tt/n0},
              simplify = 'array')

  # knots location(equally distributed)
  ## ensure knots location doesn't exceed t location
  if (is.null(k_list)) {
    k_list <- seq(min(t_list), max(t_list), by = (max(t_list)-min(t_list))%/%(K-1))
  }
  pos_k <- sapply(1:K,
                function(k) {ifelse(k_list[k]%%n0==0, k_list[k]%/%n0, k_list[k]%/%n0+1)})
  S_k <- S[,,pos_k]

  # criteria for convergence
  delta = 0.0001


  # kernel smooth weight
  # dmat: K x T_count
  # w_h_k_t: K x T_count
  dmat <- sapply(t_list, function(t){abs(t-k_list)})
  w_h_k_t <- apply(dmat, 2, function(dmat_t){dnorm(dmat_t/h)})
  colsum_w = colSums(w_h_k_t)
  rowsum_w = rowSums(w_h_k_t)


  # inital values
  #Omega_k_list = sapply(alply(S_k,3), function(x1) solve(x1+diag(1,p)*0.01), simplify = 'array' )
  Omega_k_list <- array(NA, c(p,p,K))
  for(k in 1:K){
    St_wt <- apply(sapply(1:T_count, function(t){w_h_k_t[k,t]*S[,,t]}, simplify = 'array'),
                   1:2,
                   sum)/(rowsum_w[k])
    Omega_k_list[,,k] <- solve(St_wt + diag(1,p)*0.01)
  }
  Omega_t_list = sapply(alply(S,3), function(x1) solve(x1+diag(1,p)*0.01), simplify = 'array')


  Omega_k_list.old = array(0,c(p,p,K))
  Omega_t_list.old = array(0,c(p,p,T_count))
  count = 0
  state = 'Less than max.iter!'

  if (parallel == TRUE) {
    myCluster <- makeCluster(ncore, # number of cores to use
                             type = "PSOCK") # type of cluster
    registerDoParallel(myCluster)

    # start BCD algorithm
    while(max(abs(Omega_k_list-Omega_k_list.old))>delta | max(abs(Omega_t_list-Omega_t_list.old))>delta) {
      # record current Omega0 & Omegas
      Omega_k_list.old = Omega_k_list
      Omega_t_list.old = Omega_t_list

      # 1st step(parallel) :
      Omega_t_list = foreach(t = 1:T_count, .combine='acomb', .multicombine=TRUE) %dopar% {
        tmp <- apply(sapply(1:K, function(k){w_h_k_t[k,t]*solve(Omega_k_list[,,k])}, simplify = 'array'),
                     1:2,
                     sum)
        st = (S[,,t] + lambda2*tmp)/(1+lambda2)
        Omega_t = glassoFast::glassoFast(st, rho)$wi
        return(Omega_t)
      }

      # 2nd step:
      for(k in 1:K){
        tmp <- apply(sapply(1:T_count, function(t){w_h_k_t[k,t]*Omega_t_list[,,t]}, simplify = 'array'),
                     1:2,
                     sum)
        Omega_k_list[,,k] <- tmp/(rowsum_w[k])
      }

      # record BCD iterations
      count = count+1
      # cat('iteration',count,'done \n')

      if(count>max.iter) {
        #cat('Omegas fail to converge. \n')
        state = 'Exceed max.iter!'
        break
      }
    }
    stopCluster(myCluster)

  } else {
    # start BCD algorithm
    while(max(abs(Omega_k_list-Omega_k_list.old))>delta | max(abs(Omega_t_list-Omega_t_list.old))>delta) {
      # record current Omega0 & Omegas
      Omega_k_list.old = Omega_k_list
      Omega_t_list.old = Omega_t_list

      # 1st step :
      for(t in 1:T_count) {
        tmp <- apply(sapply(1:K, function(k){w_h_k_t[k,t]*solve(Omega_k_list[,,k])}, simplify = 'array'),
                     1:2,
                     sum)
        st = (S[,,t]+lambda2*tmp)/(1+lambda2*colsum_w[t])
        Omega_t_list[,,t] = glassoFast::glassoFast(st, rho)$wi
      }

      # 2nd step:
      for(k in 1:K){
        tmp <- apply(sapply(1:T_count, function(t){w_h_k_t[k,t]*Omega_t_list[,,t]}, simplify = 'array'),
                     1:2,
                     sum)
        Omega_k_list[,,k] <- tmp/(rowsum_w[k])
      }

      # record BCD iterations
      count = count+1
      # cat('iteration',count,'done \n')

      if(count>max.iter) {
        #cat('Omegas fail to converge. \n')
        state = 'Exceed max.iter!'
        break
      }
    }
  }


  res = list(t_list, Omega_t_list, S, state, count)
  names(res)=c("t_list", "Omega_t_list", "S", "state", "count")

  return(res)
}

