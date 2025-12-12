rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/projects/DRPT/code/simulationCpp/stroop")

library(MASS)
library(latex2exp)
library(DRPT)

set.seed(180966) 

## --- load & prep ---
X = c(); Y = c()
data = read.csv("Stro.csv", header = TRUE)
for (i in 1:dim(data)[1]){
  str = strsplit(data[i,], " ")[[1]]
  X = c(X, as.numeric(str[2]))
  Y = c(Y, as.numeric(str[3]))
}
normalize = function(x) (x - min(x)) / (max(x) - min(x))
X  = normalize(X); Y = normalize(Y)

d = 1
gaussian.kernel = function(x, y, lambda = 1){
  return(lambda^(-d) * exp(-sum(((x - y) ^ 2) / (lambda ^ 2))))
}
median.bandwidth = median(as.vector(dist(cbind(X, Y), method = "euclidean", diag = FALSE, upper = FALSE)))
kfun = function(x, y) gaussian.kernel(x, y, lambda = median.bandwidth)

## --- core grid & main analysis ---
xxx = 30
lambdas = seq(from = 0, to = 0.3, length.out = xxx)

pval_main = numeric(length(lambdas))
for (k in seq_along(lambdas)) {
  l = lambdas[k]
  cat("lambda:", l, "\n")
  r = function(x) exp(x / l)
  pval_main[k] = DRPT(X, Y, r = r, kernel = kfun)
}

## --- sensitivity: 10 subsamples of size 100 ---
n_sub = 10
m_sub = 100
n_all = length(X)
idx_all = seq_len(n_all)

pval_subs = matrix(NA_real_, nrow = n_sub, ncol = length(lambdas))

for (b in 1:n_sub) {
  idx_b = sample(idx_all, m_sub, replace = FALSE)
  Xb = X[idx_b]; Yb = Y[idx_b]
  
  median.bandwidth = median(as.vector(dist(cbind(Xb, Yb), method = "euclidean", diag = FALSE, upper = FALSE)))
  kfun = function(x, y) gaussian.kernel(x, y, lambda = median.bandwidth)
  
  for (k in seq_along(lambdas)) {
    l = lambdas[k]
    r = function(x) exp(x / l) # same r for consistency
    pval_subs[b, k] = DRPT(Xb, Yb, r = r, kernel = kfun)
  }
}

## --- save results ---
dir.create("experiments/results", recursive = TRUE, showWarnings = FALSE)
out_main = data.frame(Lambda = lambdas, P_lambda = pval_main)
write.csv(out_main, "experiments/results/stroop_pval.csv", row.names = FALSE)

out_subs = data.frame(Lambda = lambdas, t(pval_subs))
names(out_subs)[-1] = paste0("sub_", seq_len(n_sub))
write.csv(out_subs, "experiments/results/stroop_pval_subsamples.csv", row.names = FALSE)

## --- plotting: main curve + 10 thin overlays ---
dir.create("experiments/pictures", recursive = TRUE, showWarnings = FALSE)

# File paths
file_main = "experiments/results/stroop_pval.csv"
file_subs = "experiments/results/stroop_pval_subsamples.csv"

# Safety checks
if (!file.exists(file_main)) stop("Missing: ", normalizePath(file_main, mustWork = FALSE))
if (!file.exists(file_subs)) stop("Missing: ", normalizePath(file_subs, mustWork = FALSE))

# Load saved results
df_main = read.csv(file_main, check.names = FALSE)
df_subs = read.csv(file_subs, check.names = FALSE)

# Extract series
lambdas_main = df_main$Lambda
pval_main    = df_main$P_lambda

lambdas_subs = df_subs$Lambda
sub_cols     = grep("^sub_", names(df_subs), value = TRUE)
if (length(sub_cols) == 0) stop("No columns matching '^sub_' in subsamples file.")

# If Lambda grids differ, we'll use the union for xlim but plot each on its own grid
same_grid = isTRUE(all.equal(lambdas_main, lambdas_subs))

# -----------------------
# Plot
# -----------------------
png("experiments/pictures/p_val_from_saved.png")
par(mfrow = c(1, 1))

# Axis limits
ylim = c(0, 1)
xlim = range(c(lambdas_main, lambdas_subs), finite = TRUE)

# Main line
plot(lambdas_main, pval_main, type = "b", col = "darkviolet",
     ylim = ylim, xlim = xlim,
     ylab = "p-value", xlab = TeX(r'($\eta$)'),
     pch = 14, cex = 0.8)

abline(h = 0.05, col = "red", lty = 2)

# Subsample overlays (thin, semi-transparent)
thin_col = adjustcolor("black", alpha.f = 0.35)

if (same_grid) {
  # Fast path when grids match
  pval_subs_mat = as.matrix(df_subs[, sub_cols, drop = FALSE])
  # Each column is a subsample curve over the same x-grid
  for (j in seq_len(ncol(pval_subs_mat))) {
    lines(lambdas_main, pval_subs_mat[, j], col = thin_col, lwd = 0.7)
  }
} else {
  # Robust path: draw each subsample on its own grid
  for (nm in sub_cols) {
    lines(lambdas_subs, df_subs[[nm]], col = thin_col, lwd = 0.7)
  }
}

legend("topleft", inset = 0.001, bty = "n",
       legend = c("(E1) - Full sample", "(E1) - Subsamples"),
       col = c("darkviolet", thin_col),
       lty = c(1, 1), pch = c(14, NA), lwd = c(1.5, 0.7))

dev.off()

cat("Saved figure to: experiments/pictures/p_val_from_saved.png\n")