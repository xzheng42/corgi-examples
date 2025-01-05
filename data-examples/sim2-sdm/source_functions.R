fitGLMPois <- function(zz, xx, xx_ppd, BB) {
  
  be_est <- vector("list", length = 4)
  nn <- length(zz)
  KK <- ncol(xx_ppd)
  

  #-------------------------------------------------------------------------------
  # GLM-oracle
  #-------------------------------------------------------------------------------
  # Fit a glm
  orig_out <- glm(zz ~ xx, family = poisson)
  orig_be <- orig_out$coefficients

  # Parametric bootstrap
  orig_be_pb_sam <- array(NA, dim = c(length(orig_be), BB))
  orig_yy <- exp(cbind(1, xx) %*% orig_be)
  
  for (b in 1:BB) {
    zz_b <- rpois(nn, orig_yy)
    orig_out_b <- glm(zz_b ~ xx, family = poisson)
    orig_be_pb_sam[, b] <- orig_out_b$coefficients
  }

  be_est[[1]] <- orig_be_pb_sam


  #-------------------------------------------------------------------------------
  # GLM-plugin
  #-------------------------------------------------------------------------------
  # covariate
  plugin_xx <- rowMeans(xx_ppd)

  # Fit a glm
  plugin_out <- glm(zz ~ plugin_xx, family = poisson)
  plugin_be <- plugin_out$coefficients

  # Parametric bootstrap
  plugin_be_pb_sam <- array(NA, dim = c(length(plugin_be), BB))
  plugin_yy <- exp(cbind(1, plugin_xx) %*% plugin_be)

  for (b in 1:BB) {
    zz_b <- rpois(nn, plugin_yy)
    plugin_out_b <- glm(zz_b ~ plugin_xx, family = poisson)
    plugin_be_pb_sam[, b] <- plugin_out_b$coefficients
  }

  be_est[[2]] <- plugin_be_pb_sam


  #-------------------------------------------------------------------------------
  # GLM-ensemble
  #-------------------------------------------------------------------------------
  ppd_be_pb_sam_pool <- vector("list", length = KK)

  for (k in 1:KK) {

    xx_k <- xx_ppd[, k]
    ppd_out_k <- glm(zz ~ xx_k, family = poisson)
    ppd_be_k <- ppd_out_k$coefficients
    yy_k <- exp(cbind(1, xx_k) %*% ppd_be_k)

    # Parametric bootstrap
    ppd_be_pb_sam_k <- array(NA, dim = c(length(ppd_be_k), BB))
    for (b in 1:BB) {
      zz_b <- rpois(nn, yy_k)
      ppd_out_kb <- glm(zz_b ~ xx_k, family = poisson)
      ppd_be_pb_sam_k[, b] <- ppd_out_kb$coefficients
    }

    ppd_be_pb_sam_pool[[k]] <- ppd_be_pb_sam_k

  }

  be_est[[3]] <- ppd_be_pb_sam_pool


  #-------------------------------------------------------------------------------
  # GLM-Berkson
  #-------------------------------------------------------------------------------
  xx_ppd_mu <- rowMeans(xx_ppd)
  ppd_mu_be <- glm(zz ~ xx_ppd_mu, family = poisson)$coefficients

  bks_be_pb_sam_pool <- vector("list", length = KK)
  
  for (k in 1:KK) {

    xx_k <- xx_ppd[, k]
    delta_k <- ppd_mu_be[2] * (xx_k - xx_ppd_mu)
    bks_out_k <- glm(zz ~ xx_ppd_mu, offset = delta_k, family = poisson)
    bks_be_k <- bks_out_k$coefficients
    yy_k <- exp(cbind(1, xx_ppd_mu) %*% bks_be_k + delta_k)

    # Parametric bootstrap
    bks_be_pb_sam_k <- array(NA, dim = c(length(bks_be_k), BB))
    for (b in 1:BB) {
      zz_b <- rpois(nn, yy_k)
      bks_out_kb <- glm(zz_b ~ xx_ppd_mu, offset = delta_k, family = poisson)
      bks_be_pb_sam_k[, b] <- bks_out_kb$coefficients
    }

    bks_be_pb_sam_pool[[k]] <- bks_be_pb_sam_k

  }

  be_est[[4]] <- bks_be_pb_sam_pool
  
  list(all_be = be_est, ppd_mu_be = ppd_mu_be)
  
}



evalMod <- function(be_est, be, zz, yy, xx, xx_ppd, hfunc1, hfunc2) {
  
  all_be <- be_est[[1]]
  ppd_mu_be <- be_est[[2]]
  
  nn <- length(zz)
  KK <- ncol(xx_ppd)
  BB <- ncol(all_be[[1]])
  

  #-------------------------------------------------------------------------------
  # GLM-oracle
  #-------------------------------------------------------------------------------
  orig_be_pb_sam <- all_be[[1]]
  
  # bias
  orig_be_bias <- rowMeans(orig_be_pb_sam)[2] - be[2]
  
  # 95% CI for beta
  orig_be_qq <- quantile(orig_be_pb_sam[2,], probs = c(0.025, 0.975))
  orig_be_ci_ind <- (orig_be_qq[1] <= be[2]) * (orig_be_qq[2] >= be[2])
  
  # rmse and 95% pred CI
  orig_pred_pb_sam <- array(NA, dim = c(nn, BB))
  for (b in 1:BB) {
    yy_b <- hfunc1(xx, orig_be_pb_sam[, b])
    orig_pred_pb_sam[, b] <- yy_b
  }
  # rmse
  orig_pred <- rowMeans(orig_pred_pb_sam)
  orig_pred_rmse <- sqrt(mean((orig_pred - yy)^2)) 
  # 95% CI
  orig_pred_qq <- apply(orig_pred_pb_sam, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
  orig_pred_ci_cover <- mean((orig_pred_qq[1,] <= yy) * (orig_pred_qq[2, ] >= yy))
  

  #-------------------------------------------------------------------------------
  # GLM-plugin
  #-------------------------------------------------------------------------------
  plugin_be_pb_sam <- all_be[[2]]
  
  # bias
  plugin_be_bias <- rowMeans(plugin_be_pb_sam)[2] - be[2]
  # 95$ CI for beta
  plugin_be_qq <- quantile(plugin_be_pb_sam[2, ], probs = c(0.025, 0.975))
  plugin_be_ci_ind <- (plugin_be_qq[1] <= be[2]) * (plugin_be_qq[2] >= be[2])

  # rmse and 95% pred CI
  plugin_xx <- rowMeans(xx_ppd)
  plugin_pred_pb_sam <- array(NA, dim = c(nn, BB))
  for (b in 1:BB) {
    yy_b <- hfunc1(plugin_xx, plugin_be_pb_sam[, b])
    plugin_pred_pb_sam[, b] <- yy_b
  }
  # rmse
  plugin_pred <- rowMeans(plugin_pred_pb_sam)
  plugin_pred_rmse <- sqrt(mean((plugin_pred - yy)^2)) 
  # 95% CI
  plugin_pred_qq <- apply(plugin_pred_pb_sam, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
  plugin_pred_ci_cover <- mean((plugin_pred_qq[1,] <= yy) * (plugin_pred_qq[2, ] >= yy))
  
  
  #-------------------------------------------------------------------------------
  # GLM-ensemble
  #-------------------------------------------------------------------------------
  ppd_be_pb_sam_pool <- all_be[[3]]
  ppd_be_pb_sam_pool <- simplify2array(ppd_be_pb_sam_pool)
  
  # bias
  ppd_be_bias <- mean(ppd_be_pb_sam_pool[2, ,]) - be[2]
  
  # 95% CI for beta
  ppd_be2_qq <- quantile(ppd_be_pb_sam_pool[2, ,], probs = c(0.025, 0.975))
  ppd_be2_ci_ind <- (ppd_be2_qq[1] <= be[2]) * (ppd_be2_qq[2] >= be[2])
  
  # rmse and 95% pred CI
  # rmse
  ppd_pred_pb_sam_pool <- array(NA, dim = c(nn, BB, KK))
  for (k in 1:KK) {
    xx_k <- xx_ppd[, k]
    for (b in 1:BB) {
      yy_kb <- hfunc1(xx_k, ppd_be_pb_sam_pool[, b, k])
      ppd_pred_pb_sam_pool[, b, k] <- yy_kb
    }
  }
  ppd_pred <- rowMeans(ppd_pred_pb_sam_pool)
  ppd_pred_rmse <- sqrt(mean((ppd_pred - yy)^2))
  
  # 95% pred CI
  ppd_pred_qq <- apply(ppd_pred_pb_sam_pool, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
  ppd_pred_ci_cover <- mean((ppd_pred_qq[1,] <= yy) * (ppd_pred_qq[2,] >= yy))
  
  
  #-------------------------------------------------------------------------------
  # GLM-Berkson
  #-------------------------------------------------------------------------------
  bks_be_pb_sam_pool <- all_be[[4]]
  bks_be_pb_sam_pool <- simplify2array(bks_be_pb_sam_pool)
  
  # bias
  bks_be_bias <- mean(bks_be_pb_sam_pool[2, ,]) - be[2]
  
  # 95% CI for beta
  bks_be2_qq <- quantile(bks_be_pb_sam_pool[2, ,], probs = c(0.025, 0.975))
  bks_be2_ci_ind <- (bks_be2_qq[1] <= be[2]) * (bks_be2_qq[2] >= be[2])
  
  # rmse and 95% pred CI
  bks_pred_pb_sam_pool <- array(NA, dim = c(nn, BB, KK))
  xx_ppd_mu <- rowMeans(xx_ppd)
  for (k in 1:KK) {
    xx_k <- xx_ppd[, k]
    delta_k <- ppd_mu_be[2] * (xx_k - xx_ppd_mu)
    for (b in 1:BB) {
      yy_kb <- hfunc2(xx_ppd_mu, bks_be_pb_sam_pool[, b, k], delta_k)
      bks_pred_pb_sam_pool[, b, k] <- yy_kb
    }
  }
  # rmse
  bks_pred <- rowMeans(bks_pred_pb_sam_pool)
  bks_pred_rmse <- sqrt(mean((bks_pred - yy)^2))
  # 95% pred CI
  bks_pred_qq <- apply(bks_pred_pb_sam_pool, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
  bks_pred_ci_cover <- mean((bks_pred_qq[1,] <= yy) * (bks_pred_qq[2,] >= yy))

  
  #-------------------------------------------------------------------------------
  # output
  #-------------------------------------------------------------------------------
  orig <- list(orig_be_bias = orig_be_bias, orig_be_ci_ind = orig_be_ci_ind,
               orig_pred_rmse = orig_pred_rmse, orig_pred_ci_cover = orig_pred_ci_cover)
  
  plugin <- list(plugin_be_bias = plugin_be_bias, plugin_be_ci_ind = plugin_be_ci_ind,
                 plugin_pred_rmse = plugin_pred_rmse, plugin_pred_ci_cover = plugin_pred_ci_cover)
  
  ppd <- list(ppd_be_bias = ppd_be_bias, ppd_be_ci_ind = ppd_be2_ci_ind,
              ppd_pred_rmse = ppd_pred_rmse, ppd_pred_ci_cover = ppd_pred_ci_cover)
  
  bks <- list(bks_be_bias = bks_be_bias, bks_be_ci_ind = bks_be2_ci_ind,
              bks_pred_rmse = bks_pred_rmse, bks_pred_ci_cover = bks_pred_ci_cover)
  
  eval <- list(orig = orig, plugin = plugin, ppd = ppd, bks = bks)
  
  eval
  
  
}


