#include <iostream>
#include <chrono>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <boost/math/special_functions/bessel.hpp>
#include <Rcpp/Benchmark/Timer.h>
#include "utils.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;


// Fit the GPB model
// [[Rcpp::export]]
List fitGPBAU(const arma::colvec& y,
              const arma::mat& X,
              const List& prior,
              const List& starting,
              const List& tuning,
              const List& mod_param,
              const List& mcmc_settings) {
  
  // Prior
  double u_phi = prior["u_phi"];
  double v_phi = prior["v_phi"];
  double u_sigmasq = prior["u_sigmasq"];
  double v_sigmasq = prior["v_sigmasq"];
  
  // Model and MCMC parameters
  arma::mat bau_distmat = mod_param["bau_distmat"];
  arma::mat aggmat = mod_param["aggmat"];
  int niter = mcmc_settings["niter"];
  int nreport = mcmc_settings["nreport"];
  
  // Starting values
  arma::colvec bb = starting["regcoef"];
  double sigmasq = starting["sigmasq"];
  double phi = starting["phi"];
  double logit_phi = logit(phi, u_phi, v_phi);
  
  // Tuning parameters
  double se_phi = tuning["se_phi"];
  
  // Other stuff
  int nn = y.n_rows;
  int pp = X.n_cols;
  
  // Initialization
  double blk_log_det, cur_loglik, cur_logpost, quad_log_det;
  double cand_logit_phi, cand_phi, cand_blk_log_det, cand_quad_log_det, cand_loglik, cand_logpost;
  double cand_aa, aa;
  
  arma::colvec bb_hat, proj_y, ee, cand_bb_hat, cand_proj_y, cand_ee;
  arma::mat bau_corrmat, blk_corrmat, dat_corrmat_inv, quad_corrmat;
  arma::mat cand_bau_corrmat, cand_blk_corrmat, cand_dat_corrmat_inv, cand_quad_corrmat;
  arma::mat bb_corr;
  
  bau_corrmat = expCorr(bau_distmat, phi);
  blk_corrmat = aggmat * bau_corrmat * aggmat.t();
  blk_log_det = arma::log_det_sympd(blk_corrmat);
  quad_corrmat = X.t() * arma::solve(blk_corrmat, X, arma::solve_opts::likely_sympd);
  quad_log_det = arma::log_det_sympd(quad_corrmat);
  
  proj_y = X.t() * arma::solve(blk_corrmat, y, arma::solve_opts::likely_sympd);
  bb_hat = arma::solve(quad_corrmat, proj_y);
  ee = y - X * bb_hat;
  aa = as_scalar(ee.t() * arma::solve(blk_corrmat, ee, arma::solve_opts::likely_sympd));
  
  cur_loglik = -.5 * blk_log_det - .5* quad_log_det;
  cur_loglik += -(.5 * (nn - pp) + u_sigmasq) * log(.5 * aa + v_sigmasq);
  cur_logpost = cur_loglik + log(phi - u_phi) + log(v_phi - phi);
  
  arma::colvec phi_save(niter), sigmasq_save(niter), mh_data(niter);
  arma::mat bb_save(pp, niter);
  
  // MCMC
  int icount = 0;
  int batch_accept = 0;
  int accept = 0;
  double block_time = 0.0;
  double time_diff, ert;
  
  for (int iter = 0; iter < niter; ++iter) {
    
    auto start = std::chrono::steady_clock::now();
    
    // Update phi
    cand_logit_phi= R::rnorm(logit_phi, se_phi);
    cand_phi = logitInv(cand_logit_phi, u_phi, v_phi);
    
    cand_bau_corrmat = expCorr(bau_distmat, cand_phi);
    cand_blk_corrmat = aggmat * cand_bau_corrmat * aggmat.t();
    cand_blk_log_det = arma::log_det_sympd(cand_blk_corrmat);
    cand_quad_corrmat = X.t() * arma::solve(cand_blk_corrmat, X, arma::solve_opts::likely_sympd);
    cand_quad_log_det = arma::log_det_sympd(cand_quad_corrmat);
    
    cand_proj_y = X.t() * arma::solve(cand_blk_corrmat, y, arma::solve_opts::likely_sympd);
    cand_bb_hat = arma::solve(cand_quad_corrmat, cand_proj_y);
    cand_ee = y - X * cand_bb_hat;
    cand_aa = as_scalar(ee.t() * arma::solve(cand_blk_corrmat, ee, arma::solve_opts::likely_sympd));
    
    cand_loglik = -.5 * cand_blk_log_det - .5* cand_quad_log_det;
    cand_loglik += -(.5 * (nn - pp) + u_sigmasq) * log(.5 * cand_aa + v_sigmasq);
    cand_logpost = cand_loglik + log(cand_phi - u_phi) + log(v_phi - cand_phi);
    
    if (log(arma::randu()) <= (cand_logpost - cur_logpost)) {
      
      // Update phi-relevant parameters
      phi = cand_phi;
      logit_phi = cand_logit_phi;
      blk_log_det = cand_blk_log_det;
      quad_log_det = cand_quad_log_det;
      cur_logpost = cand_logpost;
      aa = cand_aa;
      batch_accept++;
      accept++;
      mh_data(iter) = 1;
      
      // Update sigmasq
      sigmasq = 1.0 / R::rgamma(u_sigmasq + .5 * (nn - pp), 1.0 / (v_sigmasq + .5 *aa));
      
      // Update regression coefficients
      bb_corr = arma::inv_sympd(cand_quad_corrmat);
      bb = arma::mvnrnd(cand_bb_hat, sigmasq* bb_corr);
      
    }
    
    auto end = std::chrono::steady_clock::now();
    time_diff = std::chrono::duration<double>(end - start).count();
    block_time += time_diff;
    
    // save
    phi_save(iter) = phi;
    sigmasq_save(iter) = sigmasq;
    bb_save.col(iter) = bb;
    
    icount++;
    if (icount == nreport) {
      Rprintf("Sampled: %i of %i\n", (iter + 1), niter);
      ert = (block_time / nreport) * (niter - iter);
      Rprintf("Estimated remaining time: %3.2f minutes\n", ert / 60.0);
      Rprintf("Interval acceptance rate: %3.2f%%\n", 100.0 * batch_accept / nreport);
      Rprintf("Overall acceptance rate: %3.2f%%\n\n", 100.0 * accept / (iter + 1));
      icount = 0;
      batch_accept = 0;
      block_time = 0.0;
    }
    
  }
  
  List post_sams = List::create(Named("phi") = phi_save,
                                Named("sigmasq") = sigmasq_save,
                                Named("beta") = bb_save);
  
  return List::create(Named("post_sams") = post_sams,
                      Named("mh_data") = mh_data);
}


// BAUGP (joint) prediction
// [[Rcpp::export]]
arma::mat predGPBAU2(const arma::colvec& y,
                     const arma::mat& X,
                     const List& post_sams,
                     const List& mod_param,
                     const int& nreport){
   
   arma::colvec phi_save = post_sams["phi"];
   arma::colvec sigmasq_save = post_sams["sigmasq"];
   arma::mat bb_save = post_sams["beta"];
   
   arma::mat bau_distmat = mod_param["bau_distmat"];
   arma::mat aggmat = mod_param["aggmat"];
   
   double phi, sigmasq;
   int niter = phi_save.n_rows;
   int nn = bau_distmat.n_rows;
   
   arma::mat bau_corrmat, blk_corrmat, bau_blk_corr, cond_var;
   arma::colvec cond_mean, bb, fine_mean, coarse_mean;
   arma::mat latent_save(nn, niter);
   arma::mat blk_corrmat_inv, bau_blk_coef;
   
   int icount = 0;
   double block_time = 0.0;
   double time_diff, ert;
   
   for (int iter = 0; iter < niter; ++iter){
     
     auto start = std::chrono::steady_clock::now();
     
     phi = phi_save(iter);
     sigmasq = sigmasq_save(iter);
     bb = bb_save.row(iter).t();
     
     fine_mean = X * bb;
     coarse_mean = aggmat * fine_mean;
     
     bau_corrmat = sigmasq * expCorr(bau_distmat, phi);
     bau_blk_corr = bau_corrmat * aggmat.t();
     blk_corrmat = aggmat * bau_blk_corr;
     
     // Prediction
     cond_mean = fine_mean + bau_blk_corr * solve(blk_corrmat, (y - coarse_mean), arma::solve_opts::likely_sympd);
     cond_var = bau_corrmat - bau_blk_corr * solve(blk_corrmat, bau_blk_corr.t(), arma::solve_opts::likely_sympd);
     latent_save.col(iter) = arma::mvnrnd(cond_mean, cond_var);
     
     auto end = std::chrono::steady_clock::now();
     time_diff = std::chrono::duration<double>(end - start).count();
     block_time += time_diff;
     
     icount++;
     if (icount == nreport) {
       Rprintf("Sampled: %i of %i\n", (iter + 1), niter);
       ert = (block_time / nreport) * (niter - iter);
       Rprintf("Estimated remaining time: %3.2f minutes\n", ert / 60.0);
       icount = 0;
       block_time = 0.0;
     }
     
   }
   
   return latent_save;
   
 }


// BAUGP (independent) prediction
// [[Rcpp::export]]
arma::mat predGPBAU3(const arma::colvec& y,
                      const arma::mat& X,
                      const List& post_sams,
                      const List& mod_param,
                      const int& nreport){

   arma::colvec phi_save = post_sams["phi"];
   arma::colvec sigmasq_save = post_sams["sigmasq"];
   arma::mat bb_save = post_sams["beta"];

   arma::mat bau_distmat = mod_param["bau_distmat"];
   arma::mat aggmat = mod_param["aggmat"];

   double phi, sigmasq;
   int niter = phi_save.n_rows;
   int nn = bau_distmat.n_rows;

   arma::mat bau_corrmat, blk_corrmat, bau_blk_corr;
   arma::colvec bb, fine_mean, coarse_mean;
   arma::mat latent_save(nn, niter);
   double cond_mean, cond_var;
   arma::mat blk_corrmat_inv;
   arma::rowvec bau_blk_coef;
   arma::rowvec ss;


   int icount = 0;
   double block_time = 0.0;
   double time_diff, ert;

   for (int iter = 0; iter < niter; ++iter){

     auto start = std::chrono::steady_clock::now();

     phi = phi_save(iter);
     sigmasq = sigmasq_save(iter);
     bb = bb_save.row(iter).t();

     fine_mean = X * bb;
     coarse_mean = aggmat * fine_mean;

     bau_corrmat = sigmasq * expCorr(bau_distmat, phi);
     bau_blk_corr = bau_corrmat * aggmat.t();
     blk_corrmat = aggmat * bau_blk_corr;

     blk_corrmat_inv = arma::inv_sympd(blk_corrmat);
     
     for (int i = 0; i < nn; ++i) {
       
       ss = bau_blk_corr.row(i);
       bau_blk_coef = ss * blk_corrmat_inv;
       cond_mean = fine_mean(i) + arma::as_scalar(bau_blk_coef * (y - coarse_mean));
       cond_var = sigmasq - arma::as_scalar(bau_blk_coef * ss.t());
       latent_save(i, iter) = R::rnorm(cond_mean, sqrt(cond_var));

     }

     auto end = std::chrono::steady_clock::now();
     time_diff = std::chrono::duration<double>(end - start).count();
     block_time += time_diff;

     icount++;
     if (icount == nreport) {
       Rprintf("Sampled: %i of %i\n", (iter + 1), niter);
       ert = (block_time / nreport) * (niter - iter);
       Rprintf("Estimated remaining time: %3.2f minutes\n", ert / 60.0);
       icount = 0;
       block_time = 0.0;
     }

   }

   return latent_save;

 }




