# Assumes previous environment where 'fdata' has been constructed as an fd object of 102 curves
library(refund)  # for FPCA
library(MASS)    # for ginv
source("Method.R")  

# Parameters
m <- 34
T.chan <-2
total_n <- m + m * T.chan  # 34 + 68 = 102
t_grid <- seq(0, 1, length.out = 301)

# Evaluate functional data
values_all <- t(eval.fd(t_grid, fdata))  # Dimensions: 102 x 301

# Split into training and monitoring sets
values_train <- values_all[1:m, ]
values_monitor <- values_all[(m+1):total_n, ]

# Convert to FPCA format
Ly_train <- split(values_train, row(values_train))
Lt_train <- replicate(length(Ly_train), t_grid, simplify = FALSE)

# FPCA decomposition
fpca_train <- FPCA(
  Ly = Ly_train, Lt = Lt_train,
  optns = list(methodMuCovEst = 'smooth', FVEthreshold = 0.95, methodSelectK = 'FVE')
)

# Get FPCA components
mean_train <- fpca_train$mu
phi_train <- fpca_train$phi
scores_train <- fpca_train$xiEst

# Compute monitoring scores
centered_monitor <- sweep(values_monitor, 2, mean_train)
scores_monitor <- centered_monitor %*% phi_train * (1 / length(t_grid))
scores_all <- rbind(scores_train, scores_monitor)

# Whitening
sample.omega <- var(scores_train)
A_estimate <- ldl(sample.omega)$lower
scores_whitened <- t(ginv(A_estimate) %*% t(scores_all))
scores_whitened <- scores_whitened[, 1:min(ncol(scores_whitened), 8), drop = FALSE]

# Run tests
alpha <- 0.05
gamma <- 0   # Set gamma value for RSMS/CSMS thresholding

rsms <- rsms.statistic.fpca.alt(scores_whitened, m = m, T.chan = T.chan, gamm = gamma, alpha = alpha)
ssms <- ssms.statistic.fpca.alt(scores_all, m = m, T.chan = T.chan, gamma = gamma, alpha = alpha)
csms <- csms.statistic.fpca.alt(scores_all, m = m, T.chan = T.chan, gamm = gamma, alpha = alpha)

# Output results
cat("=== RSMS Test ===\n")
print(rsms)

cat("\n=== SSMS Test ===\n")
print(ssms)

cat("\n=== CSMS Test ===\n")
print(csms)
