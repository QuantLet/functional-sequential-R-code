# Clear environment
rm(list = ls()) 

# Set working directory if in RStudio
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Load libraries
library(fda)
library(ggplot2)
library(tibble)
library(dplyr)

# Load functions
source("genData.R")
source("Method.R")

# Read the data
csv.data <- read.csv("plaindata.csv", header = FALSE)
data <- csv.data[1:38325, 1]

# Functional data setup
evalPoints <- 0:364 / 364
basis <- create.bspline.basis(rangeval = c(0, 1), norder = 2, breaks = evalPoints)
coefs <- matrix(nrow = 365, ncol = 105)

for (i in 1:105) {
  fdata_i <- Data2fd(argvals = evalPoints, y = data[((i-1)*365+1):(i*365)], basisobj = basis, lambda = 0)
  coefs[, i] <- fdata_i$coefs
}
fdata <- fd(coefs, basis)

# Compute mean functions
mean1 <- mean.fd(fdata[1:55])
mean2 <- mean.fd(fdata[56:105])

# Evaluate means
t_grid <- seq(0, 1, length.out = 101)
mean_vals1 <- eval.fd(t_grid, mean1)
mean_vals2 <- eval.fd(t_grid, mean2)

# Create plot data
plot_data <- tibble(
  Time = rep(t_grid, 2),
  Flow = c(mean_vals1, mean_vals2),
  Period = rep(c("1910–1964", "1965–2014"), each = length(t_grid))
)

# Plot
p <- ggplot(plot_data, aes(x = Time, y = Flow, linetype = Period)) +
  geom_line(linewidth = 0.4, color = "black") +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(x = "Time", y = "River Flow") +
  theme_minimal(base_size = 9) +
  theme(
    legend.position = c(0.15, 0.9),
    legend.title = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

ggsave("riverFlowMeans_ggplot.png", p, width = 6, height = 4, units = "in", dpi = 1200)
