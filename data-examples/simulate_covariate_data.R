################################################################################
# Simulate covariate data
################################################################################
rm(list = ls())

outpath <- "output/to/sim1-downscaling/"

library(SpatialExtremes)
library(ggplot2)
library(dplyr)

set.seed(42)

#-------------------------------------------------------------------------------
# Generate true fine-resolution covariate fields
#-------------------------------------------------------------------------------
range_vals <- seq(0.1, 0.8, by = 0.1)
sill_vals <- seq(1, 5, by = 1)

x_seq <- seq(0, 1, length = 50)
y_seq <- seq(0, 1, length = 60)
fine_grid <- as.matrix(expand.grid(x_seq, y_seq))
colnames(fine_grid) <- c("lon", "lat")
nn <- nrow(fine_grid)

# Simulate from GPB models
fine_xx <- array(NA, dim = c(nn, length(range_vals), length(sill_vals)))

for (i in seq(range_vals)) {
  
  for (j in seq(sill_vals)) {
    
    xx <- t(rgp(1, fine_grid, cov.mod = "powexp", nugget = 0, smooth = 1,
                sill = sill_vals[j], range = range_vals[i]))
    fine_xx[, i, j] <- xx
    
    dat <- data.frame(fine_grid, xx = xx)
    
  }
}

#-------------------------------------------------------------------------------
# Aggregate to coarse-resolution data
#-------------------------------------------------------------------------------
kk_x <- 5
kk_y <- 6
mm <- kk_x * kk_y

fine_grid <- as.data.frame(fine_grid)
fine_grid$x_cut <- as.numeric(cut(fine_grid[, "lon"], kk_x))
fine_grid$y_cut <- as.numeric(cut(fine_grid[, "lat"], kk_y))
grid_index <- expand.grid(x_cut = 1:kk_x, y_cut = 1:kk_y)
grid_index$par_label <- as.numeric(row.names(grid_index))
fine_grid <- left_join(fine_grid, grid_index)

coar_lon <- aggregate(fine_grid$lon, by = list(fine_grid$par_label), mean)$x
coar_lat <- aggregate(fine_grid$lat, by = list(fine_grid$par_label), mean)$x
coar_grid <- data.frame(lon = coar_lon, lat = coar_lat)

coar_xx <- array(NA, dim = c(mm, length(range_vals), length(sill_vals)))

for (i in seq(range_vals)) {
  
  for (j in seq(sill_vals)) {
    
    xx <- fine_xx[, i, j]
    xxb <- aggregate(xx, by = list(fine_grid$par_label), mean)$x
    coar_xx[, i, j] <- xxb
    
    coar_dat <- data.frame(coar_grid, xxb = xxb)
    
  }
  
}


#-------------------------------------------------------------------------------
# Output covariate data
#-------------------------------------------------------------------------------
save(fine_grid, coar_grid, fine_xx, coar_xx, range_vals, sill_vals,
     file = paste0(outpath, "covar_dat.rdata"))
