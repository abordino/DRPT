rm(list = ls())
set.seed(1203)

library(MASS)
library(pbapply)
library(ggplot2)
library(future.apply)
library(data.table)

setwd("~/Documents/phd/distr_shift/simulationCpp/diamonds")
source("./DRPTcond2sample.R")

tag = "real_low_dim"
data("diamonds")
data = diamonds

d = 7; alpha = 0.05
gaussian.kernel = function(x, y, lambda = 1){
  return(lambda^(-d) * exp(-sum(((x - y) ^ 2) / (lambda ^ 2))))
}

##### sample data ##############
s = 6
X = as.matrix(data[, c("carat", "depth", "table", "x", "y", "z")], nrow=nrow(data), ncol=s)
colnames(X) = c("V1", "V2", "V3", "V4", "V5", "V6")
Y = data$price

normalize = function(x) {
  (x - min(x)) / (max(x) - min(x))
}

X_norm = apply(X, 2, normalize)
Y_norm = normalize(Y)

sample_data = function(X, Y, n, is_null = TRUE, is_x1 = TRUE) {
  if (is_x1) {
    # uniform sampling for X1
    X_idx = sample(1:nrow(X), nrow(X)%/%2, replace = FALSE)
    X1 = X[X_idx, , drop = FALSE]
    Y_subset = Y[X_idx]
    
    # subsample from X1 to construct x1
    x_idx = sample(1:nrow(X1), n, replace=FALSE)
    x = X1[x_idx,,drop=FALSE]
  } else {
    # biased sampling for X2 based on normal distribution
    feature_to_bias = X[, 1]  
    prob = dnorm(feature_to_bias, 0, 1)
    prob = prob / sum(prob)  
    X_idx = sample(1:nrow(X), nrow(X)%/%2, replace = FALSE, prob = prob)
    X2 = X[X_idx, , drop = FALSE]
    Y_subset = Y[X_idx]  
    
    # subsample from X2 to construct x2
    x_idx = sample(1:nrow(X2), n, replace=FALSE)
    x = X2[x_idx,,drop=FALSE]
  }
  
  if (is_null) {
    # Null hypothesis: uniform sampling from Y values
    y = sample(Y_subset, size=n, replace=FALSE)
  } else {
    # Alternative hypothesis: introduce bias in Y1 and Y2
    if (is_x1) {
      u = dunif(Y_subset, 0, 1)
    } else {
      u = exp(-Y_subset)
    }
    u = u / sum(u)
    y = sample(Y_subset, size = n, prob = u, replace = FALSE)
  }
  
  return(list(x = x, y = y))
}

#####----------------------------------------#####
#####----------------------------------------#####
n_values = c(1600, 2000)
n_sims = 200
estimators = c("KLR")

results_list = list()

print("---------Start of the simulation-----------")
for (n in n_values) {
  for (is_null in c(FALSE, TRUE)) {
    
    h_label = if (is_null) "Null" else "Alternative"
    
    for (est in estimators) {
      result = pbapply::pbsapply(1:n_sims, function(sim) {
        seed = 1203 + sim
        set.seed(seed)
        
        # generate data
        d1 = sample_data(X_norm, Y_norm, n, is_null, TRUE)
        set.seed(seed + n_sims)
        d2 = sample_data(X_norm, Y_norm, n, is_null, FALSE)
        
        # Apply our test
        pvalue = DRPTcond2sample(d1, d2, S = 50, H = 99, kernel = gaussian.kernel,
                                 prop = 0.2, est.method = est)
        # Decision
        rejection = 0
        if (pvalue < alpha){
          rejection = 1
        }
        return(rejection)
      })
      
      mean_result = mean(result)
      results_list[[length(results_list) + 1]] = data.table(
        n = n,
        h_label = h_label,
        estimator = est,
        rejection_rate = mean_result
      )
      
      cat("[DRPT2sample]", "| n:", n, "| Estimator:", est, "|", h_label, "| Rejection Rate:", mean_result, "\n", strrep("-", 80), "\n")
    }
  } 
}

results_dt = rbindlist(results_list)

#####----------------------------------------#####
# Save the results
#####----------------------------------------#####
filename = paste0("./experiments/real_examples/results/KLR2000_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")