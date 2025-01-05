#' @title Aggregation matrix
#'
#' @description This function generates the aggregation matrix considered in the manuscript.
#'
#' @param fine_grid an \eqn{N_D \times 5} data frame; need more details...
#' @param coar_grid an \eqn{M_D \times 2} matrix of coarse-resolution coordinates
#'
#' @return An \eqn{M_D\times N_D} aggregation matrix.
#'   
#' @export
#' 
agg <- function(fine_grid, coar_grid) {
  
  Nd <- nrow(fine_grid)
  Md <- nrow(coar_grid)
  
  aggmat <- array(0, dim = c(Md, Nd))
  
  for (j in 1:Md) {
    j_idx <- which(fine_grid$par_label == j)
    aggmat[j, j_idx] <- 1 / length(j_idx)
  }
  
  aggmat
  
}

#' @title Thinniing of MCMC samples
#'
#' @description This function provides thinning of MCMC samples.
#'
#' @param object an object of class gpbau.
#' @param nburn an integer: number of burn-in samples.
#' @param nthin an integer: thinning degree. 
#'
#' @return An object of class \code{gpbau}
#'   
#' @export
#' 
thinning <- function(object, nburn, nthin) {

    niter <- object$mcmc_settings$niter
    idx <- seq(nburn + 1, niter, by = nthin)
    
    post_sams <- object$post_sams
    post_sams$beta <- t(post_sams$beta)
    post_sams_thin <- lapply(post_sams, function(x) x[idx, ])

    post_sams_thin$beta <- as.matrix(post_sams_thin$beta)

    object$post_sams_thin <- post_sams_thin

    object

}


#' @title Find initial values of GPB parameters
#'
#' @description This function finds initial values of GPB range and variance parameters
#'
#' @param xxb a vector of coarse-resolution data of length \eqn{M_D}.
#' @param mod_param a list of parameters needed for model fitting. Currently the function requires:
#'                  \code{aggmat}: an \eqn{M_D \times N_D} matrix that aggregates fine-resolution data to coarse-resolution data;
#'                  \code{bau_distmat}: an \eqn{N_D \times N_D} matrix about pairwise-locaiton distances.
#' @param phi_vals a vector of gridded values of GPB range parameter \code{phi}.
#' @param n_sigmasq an integer to create a grid of sigmasq.
#'
#' @return The return object is a list that comprises a matrix of log-likelihood values
#'         and a vector of gridded values of sigmasq
#'   
#' @export
#' 
findinit <- function(xxb, mod_param, phi_vals, n_sigmasq) {
  
  data <- list(resp = xxb, covar = as.matrix(rep(1, length(xxb))),
               aggmat = mod_param$aggmat, distmat = mod_param$bau_distmat)
  
  ll_vals <- array(NA, dim = c(length(phi_vals), n_sigmasq))
  sigmasq_vals <- array(NA, dim = c(length(phi_vals), n_sigmasq))
  
  # Calculate GPB's log likelihood over a grid 
  for (k in seq_along(phi_vals)) {
    
    ww <- data$aggmat
    bau_covmat <- exp(-data$distmat / phi_vals[k])
    coar_covmat <- ww %*% bau_covmat %*% t(ww)
    ratio <- diag(coar_covmat)[1]
    sigmasq_lb <- var(xxb) / ratio
    sigmasq_ub <- var(xxb) / ratio * 3
    sigmasq_seq <- seq(sigmasq_lb, sigmasq_ub, length = n_sigmasq)
    sigmasq_vals[k, ] <- sigmasq_seq
    
    for (l in seq_along(sigmasq_seq)) {
      param <- c(phi_vals[k], sigmasq_seq[l], mean(xxb))
      ll_vals[k, l] <- gpbauLoglik(param, data)
      # cat(paste0("done the ", k, "-", l, "th eval for the ", i,"-", j, " iter.\n"))
    }
  }
  
  # Find the largest log likelihood from the grid
  idx <- as.numeric(which(ll_vals == max(ll_vals), arr.ind = TRUE))

  # output
  list(ll_vals = ll_vals, 
       sigmasq_vals = sigmasq_vals,
       phi_init = phi_vals[idx[1]],
       sigmasq_init =  sigmasq_vals[idx[1],idx[2]])
  
}
  
# log loglikhood
gpbauLoglik <- function(param, data) {
  
  phi <- param[1]
  sigmasq <- param[2]
  be <- as.matrix(param[-(1:2)])
  
  yy <- data$resp
  rr <- data$covar
  ww <- data$aggmat
  distmat <- data$distmat
  
  bau_covmat <- sigmasq * exp(-distmat / phi)
  coar_covmat <- ww %*% bau_covmat %*% t(ww) 
  ee <- yy - rr %*% be
  
  logdet <- as.numeric(determinant(coar_covmat)[[1]])
  ll <- -.5 * t(ee) %*% solve(coar_covmat, ee) - .5 * logdet
  
  ll
  
}