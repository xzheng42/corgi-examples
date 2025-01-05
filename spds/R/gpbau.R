#' @title Function for COS-based statistical downscaling
#' 
#' @description This function fits a Gaussian process over basic areal units (GPB) to
#'              downscale coarse-resolution data via Markov Chain Monte Carlo simulations.
#' 
#' @param coarse_resp a vector of coarse-resoultion data of length \eqn{M_D}.
#' @param coarse_covar an \eqn{M_D \times p} matrix of covariates for the downscaling variable 
#'                     of interest. If an intercept is desired, the first column of the matrix should be a vector of ones.
#' @param prior a list of priors. Each element of the list corresponds to a combination of an unknown parameter 
#'               and its prior distribution; for example, "\code{phi_unif}" and "\code{sigmasq_invgamma}".
#' @param starting a list of starting values. Each element of the list corresponds to an unknown parameter.
#' @param tuning a list of tuning parameters. Each element of the list corresponds to an unknown parameter.
#'               The value portion of each element defines the variance of the random walk Metropolis sampler 
#'               proposal distribution for the unknown parameter.
#' @param mod_param a list of parameters needed for model fitting. Currently the function requires:
#'                  \code{aggmat}: an \eqn{M_D \times N_D} matrix that aggregates fine-resolution data to coarse-resolution data;
#'                  \code{bau_distmat}: an \eqn{N_D \times N_D} matrix about pairwise-locaiton distances.
#' @param mcmc_settings a list of Markov chain Monte Carlo (MCMC) simulation parameters. 
#'        \code{n_iter}: number of iterations; 
#'        \code{n_burn}: number of burn-in samples;
#'        \code{n_thin}: thinning degree. If \code{verbose = TRUE};
#'        \code{n_report} defines the interval to report MCMC progress and Metropolis acceptance rates.
#' @param verbose logical; if true, model specification, MCMC progress, and Metropolis acceptance rates
#'                are printed to the screen.
#' 
#' @return
#' An object of class \code{gpbau}. The return object is a list that comprises:
#' 
#' \tabular{ll}{
#' 
#' \code{post_samples} \tab a list of posterior samples. Each element of the list corresponds to an unknown parameter. \cr
#' \tab \cr
#' \code{coarse_data} \tab a list that consists of input data 
#'         (\code{response}, \code{covars}, and \code{coords}). \cr
#' \tab \cr
#' \code{mod_spec}: \tab a list that consists of model specifications 
#'         (\code{mod_param}: parameters needed to fit models; 
#'          \code{priors}: priors for unknown parameters). \cr
#' \tab \cr
#' \code{mcmc_params}: \tab a list that consists of unknown parameters' starting values and Metropolis sampler proposal distribution variances
#'         used in the MCMC simulation (
#'         \code{starting}: starting values;
#'         \code{tuning}: Metropolis sampler proposal distribution variances). \cr
#' \tab \cr
#' \code{mh_data}: \tab a list that consists of each unknown parameter's Metropolis 
#'                      sampler acceptance and rejection indicators.\cr
#' \tab \cr
#' \code{runtime}: \tab running time for MCMC simulation calculated using \code{system.time}. \cr
#' 
#'}
#'
#' @export
#'
gpbau <- function(coarse_resp,
                  coarse_covar,
                  prior,
                  starting = NULL,
                  tuning,
                  mod_param,
                  mcmc_settings,
                  verbose = TRUE) {


    xxb <- coarse_resp
    rrb <- coarse_covar

    #----------------------------------------------------
    # Compute aggregate matrix
    #----------------------------------------------------
    # To be added in...Now the aggregate matrix must be an element of mod_param
    aggmat <- mod_param$aggmat
    Nd <- ncol(aggmat)
    Md <- nrow(aggmat)

    #----------------------------------------------------
    # Check priors
    #----------------------------------------------------
    prior2 <- vector("list")

    prior2$u_phi <- prior$phi_unif[1]
    prior2$v_phi <- prior$phi_unif[2]
    prior2$u_sigmasq <- prior$sigmasq_invgamma[1]
    prior2$v_sigmasq <- prior$sigmasq_invgamma[2]

    #----------------------------------------------------
    # Check tuning parameters
    #----------------------------------------------------
    se_phi <- tuning$se_phi
    if (is.null(se_phi)) stop("error: MH tuning parameter for phi must be specified.")
    if (se_phi <= 0) stop("error: MH tuning parameter for phi must be positive.")    

    #----------------------------------------------------
    # Check starting values
    #----------------------------------------------------
    if (is.null(starting)) {
      
      starting <- list(
        phi = mean(prior$phi_unif), 
        sigmasq = var(xxb), 
        regcoef = matrix(c(mean(xxb)))
        )
    }

    #----------------------------------------------------
    # verbose
    #----------------------------------------------------
    if (verbose) {
        cat("--------------------------------------------------------------\n")
        cat(" << Spatial statistical downscaling using GPB >>\n")
        cat("--------------------------------------------------------------\n")
        cat(paste0("Number of coarse-resolution data: ", Md, "\n"))
        cat(paste0("Number of fine-resolution BAUs: ", Nd, "\n"))
    }
    
    #----------------------------------------------------
    # MCMC
    #----------------------------------------------------
    cat("----------------------------------------\n")
    cat("Running MCMC...\n");
    cat(paste0("Number of iterations: ", mcmc_settings$niter, "\n"))
    cat("Report progress every ", mcmc_settings$nreport, " iterations\n");
    cat("----------------------------------------\n")
    runtime <- system.time(
        mcmc_out <- fitGPBAU(xxb, rrb, prior2, starting, tuning, mod_param, mcmc_settings)
    )
    
    #----------------------------------------------------
    # Output
    #----------------------------------------------------
    out <- list(
        post_sams = mcmc_out$post_sams,
        coarse_data = list(coarse_resp = coarse_resp, 
                           coarse_covar = coarse_covar),
        mod_spec = list(mod_param = mod_param,
                        prior = prior),
        mcmc_params = list(starting = starting, 
                           tuning = tuning),
        mcmc_settings = mcmc_settings,
        mh_data = mcmc_out$mh_data,
        runtime = runtime
    )

    class(out) <- "gpbau"
    
    out

}

#' @title Function for obtaining downscaling samples
#'
#' @description This function produces downscaling samples of the coarse-resolution
#' variable using GPB via posterior predictive distribution.
#'
#' @param object an object of class gpbau.
#' @param fine_covar an \eqn{N_D \times p} matrix of covariates at the fine resolution.
#'                   Each column of the matrix corresponds to the column of the matrix in 
#'                   \code{coarse_covar} when it was used to fit GPB.
#' @param type A quoted keyword that specifies the type of prediction.
#'             Supported keywrods are \code{"independent"} and \code{"joint"}.
#' @param verbose logical; if true, progress of the prediction is printed to the screen.
#' @param nreport If \code{verbose = TRUE}, \code{n_report} defines the interval to report
#'                the progress of prediction.
#' @param ... additional arguments. No additional arguments are supported currently.
#' 
#' @return
#' The return object is an \eqn{N_D \times K} matrix, where \eqn{K} is the number 
#' of posterior predictive samples.
#' 
#' @exportS3Method 
predict.gpbau <- function(object, 
                          fine_covar,
                          type = "independent",
                          verbose = TRUE,
                          nreport = NULL,
                          ...) {

    coarse_resp <- object$coarse_data$coarse_resp

    if (!is.null(object$post_sams_thin)){
        post_sams <- object$post_sams_thin
    } else {
        post_sams <- object$post_sams
        post_sams$beta <- t(post_sams$beta)
    }

    mod_param <- object$mod_spec$mod_param
    niter <- length(post_sams$sigmasq)
    
    if (is.null(nreport)) {
        nreport <- object$mcmc_settings$nreport
        if (is.null(nreport)) {
            nreport <- niter / 10
        }
    }

    if (type == "independent") {
        pred_out <- predGPBAU3(coarse_resp, fine_covar, post_sams, mod_param, nreport)
    } else if (type == "joint") {
        pred_out <- predGPBAU2(coarse_resp, fine_covar, post_sams, mod_param, nreport)
    } else {
        stop("error: the type must be specified: independent or joint.")
    }
    
    pred_out

}

