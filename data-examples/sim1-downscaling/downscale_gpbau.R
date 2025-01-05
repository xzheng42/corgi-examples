################################################################################
# Statistical downscaling using GPB
################################################################################
rm(list = ls())

inpath <- "output/to/sim1-downscaling/"
outpath <- "output/to/sim1-downscaling/"


library(ggplot2)
library(spds)

#-------------------------------------------------------------------------------
# Load coarse-resolution covariate data
#-------------------------------------------------------------------------------
load(paste0(inpath, "covar_dat.rdata"))

#-------------------------------------------------------------------------------
# Parallel computation settings
#-------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(doFuture)
})

registerDoFuture()

ncores <- parallel::detectCores() - 1
nrange <- length(range_vals)
nsill <- length(sill_vals)

plan(multisession, workers = min(ncores, nrange*nsill))

#-------------------------------------------------------------------------------
# GPB model settings
#-------------------------------------------------------------------------------
Nd <- nrow(fine_grid)
Md <- nrow(coar_grid)
bau_distmat <- as.matrix(dist(fine_grid[, c("lon", "lat")]))

# Aggregate matrix
aggmat <- agg(fine_grid, coar_grid)

# Set up
prior <- list(phi_unif = c(0.1/3, 1),
              sigmasq_invgamma = c(3, 2))
tuning <- list(se_phi = 2.5)
mod_param <- list(bau_distmat = bau_distmat, aggmat = aggmat)
mcmc_settings <- list(niter = 8000, nreport = 1000)

coar_covar <- as.matrix(rep(1, Md))

#-------------------------------------------------------------------------------
# Fit GPB models (Step 2 in the manuscript)
#-------------------------------------------------------------------------------
set.seed(42)

# Grid search to find initial values of range and variance for model fitting 
# Evaluating log likelihood over a grid of range and variance values
phi_vals <- seq(0.1/3, 1, length = 50) 
n_sigma_vals <- 5
findint_runtime <- system.time(
  init_list <- foreach(i = seq_along(range_vals), .options.future = list(packages = "spds", seed = TRUE)) %:%
    foreach(j = seq_along(sill_vals), .options.future = list(packages = "spds", seed = TRUE)) %dofuture% {
      findinit(coar_xx[, i, j], mod_param, phi_vals, n_sigma_vals)
    }
)

phi_init <- array(NA, dim = c(length(range_vals), length(sill_vals)))
sigmasq_init <- array(NA, dim = c(length(range_vals), length(sill_vals)))

for (i in seq_along(range_vals)) {
  for (j in seq_along(sill_vals)) {
    phi_init[i, j] <- init_list[[i]][[j]]$phi_init
    sigmasq_init[i, j] <- init_list[[i]][[j]]$sigmasq_init
  }
}

# Fit GPB
cat("----------------------------------------------------\n")
cat("Fitting GPB models\n")

gpb_runtime <- system.time(
  out_list <- foreach(i = seq_along(range_vals), .options.future = list(packages = "spds", seed = TRUE)) %:%
    foreach(j = seq_along(sill_vals), .options.future = list(packages = "spds", seed = TRUE)) %dofuture%
    {
      gpbau(coarse_resp = coar_xx[,i,j],
            coarse_covar = coar_covar,
            prior = prior,
            starting = list(phi = phi_init[i, j],
                            sigmasq = sigmasq_init[i, j],
                            regcoef = matrix(c(mean(coar_xx[,i,j])))),
            tuning = tuning,
            mod_param = mod_param,
            mcmc_settings = mcmc_settings,
            verbose = FALSE)
    }
)

cat("Done!\n")
# cat(paste0("Running time: ", round(gpb_runtime[3]/60, 3), " minutes.\n"))
cat("----------------------------------------------------\n")


# Thinning the MCMC chain if needed
for (i in seq_along(range_vals)) {
  for (j in seq_along(sill_vals)) {
    out <- out_list[[i]][[j]]
    out <- thinning(out, nburn = 3000, nthin = 10)
    out_list[[i]][[j]] <- out
  }
}


#-------------------------------------------------------------------------------
# Posterior prediction (Step 3 in the manuscript)
#-------------------------------------------------------------------------------
set.seed(42)
fine_covar <- as.matrix(rep(1, Nd))
xx_pred <- vector("list", length = nrange)

runtime <- system.time(
  xx_pred <- foreach(i = seq_along(range_vals), .options.future = list(packages = "spds", seed = TRUE)) %:%
    foreach(j = seq_along(sill_vals), .options.future = list(packages = "spds", seed = TRUE)) %dofuture%
    {
      predict(out_list[[i]][[j]], fine_covar, nreport = 100)
    }
)

save(xx_pred, file = paste0(outpath, "sim_gpbau_pred.RData"))




