##########################################################################################
# Date: Correct Parallel Version for d = 1 to 10 (faithful to your working structure)
##########################################################################################

set.seed(123456789)
library(stats)
library(MASS)
library(doSNOW)
library(foreach)
library(parallel)

# Set folder
folder <- c("/Users/sunjiajing/Desktop/R2025/fda-sequential-20250428/critical-values")
setwd(folder)

# Main parameters
numrep <- 1000   # number of repetitions
T.chan <- 10  # change to 1, 2, 5, 10
m.chan <- 10000 / T.chan
gamma <- 0.15

# Dimensions to compute
d_vec <- 1:10

##########################################################################################
# Parallel settings
##########################################################################################
no.cores <- detectCores()
cl <- makeCluster(no.cores - 3, type = "SOCK")
registerDoSNOW(cl)

iterations <- numrep
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

##########################################################################################
# Start parallel simulation
##########################################################################################

results_list <- foreach(i = 1:numrep, .combine = "rbind", .options.snow = opts) %dopar% {
  
  total.sample.size <- m.chan + m.chan * T.chan
  
  # For each replication, store results for d=1,...,10
  result_row <- numeric(length(d_vec))
  
  for (d_index in seq_along(d_vec)) {
    d <- d_vec[d_index]
    
    # Generate d independent Brownian motions
    W_list <- list()
    for (k in 1:d) {
      dWk <- rnorm(total.sample.size) / sqrt(m.chan)
      Wk <- cumsum(dWk)
      W_list[[k]] <- Wk
    }
    W_mat <- do.call(cbind, W_list)  # total Brownian paths
    
    s <- seq(0 + 1e-6, T.chan, length.out = m.chan * T.chan)
    
    # Compute U.2.s
    U.2.s <- matrix(NA, nrow = d, ncol = length(s))
    for (s.ind in 1:length(s)) {
      U.2.s[, s.ind] <- W_mat[m.chan + s.ind, ] - (1 + s[s.ind]) * W_mat[m.chan, ]
    }
    
    # Compute denominator matrix V
    times <- seq(0 + 1e-6, 1, length.out = m.chan)
    temp_mat <- matrix(NA, nrow = m.chan, ncol = d)
    
    for (k in 1:d) {
      Wk_train <- W_list[[k]][1:m.chan]
      temp_mat[,k] <- Wk_train - times * Wk_train[m.chan]
    }
    
    dt <- times[2] - times[1]
    
    if (d == 1) {
      V <- sum(temp_mat^2) * dt
    } else {
      V <- diag(0, nrow = d)
      for (j in 1:m.chan) {
        V <- V + (temp_mat[j, ] %*% t(temp_mat[j, ])) * dt
      }
    }
    
    # Compute statistic
    if (d == 1) {
      result_row[d_index] <- max((U.2.s^2 / V) / ((1 + s)^2 * (s/(1+s))^(2*gamma)))
    } else {
      result_row[d_index] <- max(
        colSums((MASS::ginv(V) %*% U.2.s) * U.2.s) / ((1 + s)^2 * (s / (1 + s))^(2 * gamma))
      )
    }
  }
  result_row
}

##########################################################################################
# Stop cluster
##########################################################################################
stopCluster(cl)
close(pb)

##########################################################################################
# Save results
##########################################################################################
folder.2 <- paste0(folder, "/output")
if (!dir.exists(folder.2)) dir.create(folder.2)

setwd(folder.2)

# Save each dimension's results separately
for (d_index in seq_along(d_vec)) {
  d <- d_vec[d_index]
  temp.file.name <- paste0("M_chan_parallel_T_", T.chan, "_gamma_", gamma, "_d", d, "_nrep_", numrep, ".csv")
  write.csv(results_list[,d_index], temp.file.name, row.names = FALSE)
  
  # Compute quantiles
  r1_5 <- quantile(results_list[,d_index], 0.95)
  r1_10 <- quantile(results_list[,d_index], 0.90)
  
  cat("\nDimension d =", d, "\n")
  cat("95% quantile:", r1_5, "\n")
  cat("90% quantile:", r1_10, "\n")
}